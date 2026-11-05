# Internal catalogue lookup, shared by the engine and R-side readers.
.star_catalogue_file <- local({
  source <- destination <- NULL
  function(lib.loc=NULL) {
    archive <- system.file('ephemeris', 'sefstars.txt.xz', package='skyscapeR',
                           lib.loc=lib.loc, mustWork=TRUE)
    if (identical(archive,source) && file.exists(destination)) return(destination)
    # Swiss Ephemeris requires this exact plain-text filename. Recreate it at
    # the same session path if temporary files were removed after setup.
    if (!identical(archive,source) || is.null(destination))
      destination <<- file.path(tempfile("skyscape-stars-"),"sefstars.txt")
    directory <- dirname(destination)
    if (!dir.exists(directory) && !dir.create(directory,mode="0700"))
      stop("Cannot create a temporary stellar catalogue directory.",call.=FALSE)
    input <- xzfile(archive,"rb")
    on.exit(close(input),add=TRUE)
    staging <- tempfile("catalogue-",tmpdir=directory)
    on.exit(unlink(staging),add=TRUE)
    output <- file(staging,"wb")
    tryCatch({
      repeat {
        chunk <- readBin(input,"raw",n=65536L)
        if (!length(chunk)) break
        writeBin(chunk,output)
      }
    },finally=close(output))
    if (!file.rename(staging,destination))
      stop("Cannot prepare the temporary stellar catalogue.",call.=FALSE)
    source <<- archive
    destination
  }
})

# Swiss Ephemeris 2.10 reserves 14 bytes in its 256-byte path buffer.
# On long installations expose the same files through one session-local path.
.short_ephemeris_path <- local({
  sources <- directory <- NULL
  function(paths) {
    if (identical(paths, sources) && dir.exists(directory)) return(directory)
    dest <- tempfile("se-", tmpdir=tempdir())
    if (nchar(dest, type="bytes") > 242L || !dir.create(dest, mode="0700"))
      stop("Cannot create a short ephemeris directory. Set a shorter writable TMPDIR before starting R.", call.=FALSE)
    complete <- FALSE
    on.exit(if (!complete) unlink(dest, recursive=TRUE))
    for (src in paths) {
      files <- list.files(src, pattern="([.]se1$|^(sefstars|seleapsec|seorbel)[.]txt$)")
      files <- files[!file.exists(file.path(dest, files))]
      if (!length(files)) next
      from <- file.path(src, files)
      to <- file.path(dest, files)
      for (subdir in unique(dirname(to)))
        if (!dir.exists(subdir)) dir.create(subdir, recursive=TRUE, showWarnings=FALSE)
      linked <- suppressWarnings(file.link(from, to))
      if (any(!linked)) linked[!linked] <- suppressWarnings(file.symlink(from[!linked], to[!linked]))
      if (any(!linked) && !all(file.copy(from[!linked], to[!linked])))
        stop("Cannot stage ephemeris files in the session temporary directory.", call.=FALSE)
    }
    complete <- TRUE
    sources <<- paths
    directory <<- dest
    dest
  }
})

.configure_ephemeris <- function(lib.loc=NULL) {
  catalogue <- .star_catalogue_file(lib.loc=lib.loc)
  external <- strsplit(Sys.getenv("SE_EPHE_PATH"), .Platform$path.sep, fixed=TRUE)[[1]]
  external <- external[nzchar(external) & dir.exists(external)]
  external <- Filter(function(p) length(list.files(p, pattern="[.]se1$"))>0L, external)
  paths <- unique(c(dirname(catalogue), normalizePath(external),
                    system.file('ephemeris', package='swephRdata'),
                    system.file('ephemeris', package='swephR')))
  paths <- paths[nzchar(paths)]
  path <- paste(paths, collapse=.Platform$path.sep)
  if (nchar(path, type='bytes') > 242L) path <- .short_ephemeris_path(paths)

  old_path <- Sys.getenv('SE_EPHE_PATH', unset=NA_character_)
  on.exit({
    if (is.na(old_path)) Sys.unsetenv('SE_EPHE_PATH') else
      Sys.setenv(SE_EPHE_PATH=old_path)
  }, add=TRUE)
  Sys.unsetenv('SE_EPHE_PATH')
  # set_ephe_path alone retains the fixstar2 catalogue cache in SE 2.10.
  swephR::swe_close()
  swephR::swe_set_ephe_path(path)
  # Swiss files cover 600-year blocks. Planetary reductions also need lunar
  # data; the engine can otherwise report SWIEPH while using a Moshier Moon.
  files <- as.character(unique(unlist(lapply(paths, list.files, pattern="^se(pl|mo)(_|m)[0-9]+[.]se1$"))))
  pairs <- intersect(substring(files[startsWith(files,"sepl")],5L),
                     substring(files[startsWith(files,"semo")],5L))
  years <- as.integer(sub("[.]se1$","",chartr("_m","+-",pairs))) * 100L
  boundary <- function(y) swephR::swe_julday(y,1L,1L,0,as.integer(y>=1600))
  .ephemeris_state$starts <- vapply(sort(years),boundary,numeric(1))
  .ephemeris_state$ends <- vapply(sort(years)+600L,boundary,numeric(1))
  .ephemeris_state$available <- length(pairs)>0L
  invisible(paths)
}

.read_star_catalogue <- function() {
  file <- .star_catalogue_file()
  nfields <- max(utils::count.fields(file, sep=',', quote='', comment.char='#'),
                 na.rm=TRUE)
  if (nfields < 14L) stop('Stellar catalogue must contain at least 14 columns.')
  stars <- utils::read.table(file, header=FALSE, sep=',', quote='',
                            comment.char='#', fill=TRUE, strip.white=TRUE,
                            col.names=paste0('V', seq_len(nfields)),
                            stringsAsFactors=FALSE)
  names(stars)[seq_len(14L)] <- c('name', 'identifier', 'ICRS',
    'RA.1', 'RA.2', 'RA.3', 'Dec.1', 'Dec.2', 'Dec.3',
    'pm.1', 'pm.2', 'radvel', 'plx', 'magV')
  if (nfields >= 16L) {
    names(stars)[15:16] <- c('stellar', 'cat')
  }
  stars
}

# Default sky layers draw individual stars and stellar systems. Use exact
# identifiers for integrated sources: M Vel and NGC 1912 5 62 are stars.
.plot_star_catalogue <- function(catalogue = .read_star_catalogue()) {
  integrated <- c("M 45", "M 31", "NGC 104", "IC 2581", "IC 2395", "M 10",
    "NGC 2281", "RMC 136", "M 19", "M 14", "NGC 6250", "M 5", "M 2",
    "M 3", "M 55", "IC 1848", "M 92", "NGC 362", "NAME Centaurus A",
    "NGC 5986", "M 81", "MCG+02-37-031", "GCentre", "GPole")
  id <- catalogue$identifier
  # These unresolved pairs already have combined photometry and system astrometry.
  # Keep component records available when a custom catalogue omits the system.
  systems <- intersect(c("alf UMa","zet Her"),id)
  components <- unlist(lapply(systems,function(system) paste(system,c("A","B"))),
                       use.names=FALSE)
  catalogue[!(id %in% c(integrated,components) | startsWith(id, "SN ")), , drop=FALSE]
}

.ephemeris_state <- new.env(parent=emptyenv())
.ephemeris_state$notified <- FALSE

.ephemeris_install <- 'install.packages("swephRdata", repos = "https://rstub.r-universe.dev", type = "source")'

.ephemeris_failure <- function(jd, body, detail) {
  date <- swephR::swe_revjul(jd, 1L)
  msg <- sprintf("File-backed Swiss Ephemeris required for %s on %d/%02d/%02d (Gregorian UT). %s\nRun %s and restart R, or set SE_EPHE_PATH to covering Swiss Ephemeris files before loading skyscapeR.",
                 body, date$year, date$month, date$day, detail, .ephemeris_install)
  stop(errorCondition(msg, class="skyscapeR_ephemeris_error"))
}

.file_position <- function(jd, body, flags, label=NULL) {
  block <- findInterval(jd,.ephemeris_state$starts)
  missing <- which(block==0L | jd>=c(-Inf,.ephemeris_state$ends)[block+1L])
  if (length(missing)) {
    if (is.null(label)) label <- swephR::swe_get_planet_name(body)
    .ephemeris_failure(jd[missing[1L]],label,"Planetary and lunar file coverage is incomplete.")
  }
  ans <- swephR::swe_calc_ut(jd, body, bitwOr(as.integer(flags), 2L))
  bad <- which(ans$return < 0L | bitwAnd(ans$return, 7L) != 2L |
               grepl("moshier|not found|damaged", ans$serr, ignore.case=TRUE))
  if (length(bad)) {
    i <- bad[1L]
    if (is.null(label)) label <- swephR::swe_get_planet_name(body)
    .ephemeris_failure(jd[i], label, paste(ans$serr[i], collapse=" "))
  }
  ans
}

# Event routines return a crossing status rather than an ephemeris flag.
.file_crossing <- function(day, body, starname, flags, direction, geo, atm, temp, altitude) {
  .file_position(day, body, 2L)
  ans <- swephR::swe_rise_trans_true_hor(day,body,starname,bitwOr(as.integer(flags),2L),
                                         direction,geo,atm,temp,altitude)
  if (ans$return == 0L) .file_position(ans$tret[1], body, 2L)
  if (any(grepl("moshier",ans$serr,ignore.case=TRUE)))
    .ephemeris_failure(day,swephR::swe_get_planet_name(body),ans$serr)
  ans
}

# Leave ordinary catalogue misses and unreachable targets to the caller, while
# preserving missing-data errors in functions that deliberately try many targets.
.try_target <- function(expr) {
  result <- try(expr, silent=TRUE)
  if (inherits(attr(result,"condition"),"skyscapeR_ephemeris_error"))
    stop(attr(result,"condition"))
  result
}
