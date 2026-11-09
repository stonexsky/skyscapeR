

#' Find celestial targets within declination and time ranges
#' @description Screen solar, lunar and stellar candidates whose declinations overlap a measured
#' orientation range during a proposed period of use.
#' @details The search tests overlap in declination. Solar and lunar limiting declinations are
#' sampled across the years of the interval. Site coordinates enable zenith targets
#' and a horizon enables the spatial equinox. The named equinoxes and solstices use
#' geocentric apparent solar longitudes of 0, 90, 180 and 270 degrees. Dates are UT
#' civil days in the requested calendar. Different dates for the first and last year
#' appear in endpoint order separated by a slash. Intermediate years can differ.
#' The sunrise/set date ranges also use the two endpoint years. Solar and lunar candidates
#' are selected before display rounding. Zenith candidates use the period-specific tropical
#' limit and require a reachable solar passage. For spatial equinoxes the reported range
#' uses rising and setting solutions at the endpoint epochs, not an annual interior search.
#' For each star passing the catalogue magnitude cut, the search computes its
#' declination at noon UT on 1 January Gregorian of every astronomical year,
#' including both endpoints. It tests overlap before rounding and linearly
#' interpolates between annual samples to estimate entry and exit. Separate
#' passages are reported separately, with boundaries rounded to the nearest year.
#' A contact within 0.0000001 degree of a stellar band edge is treated as a closed-boundary
#' match. Shorter variations within a sampling year are unresolved. Magnitude describes
#' catalogue brightness and does not model atmospheric visibility.
#'
#' Silva and Lomsdalen (2025, table 5) give an example of stellar target screening for the Maltese temples.
#'
#' Astronomical positions require Swiss Ephemeris files covering the epoch. We
#' recommend installing \pkg{swephRdata}; see \code{\link{skyscapeR-package}} for setup.
#' @references Silva, F. and Lomsdalen, T. (2025). 'No easy way from the earth to the stars': a new statistical approach to the orientation of the Maltese temples. Archaeological and Anthropological Sciences 17, 96. \doi{10.1007/s12520-025-02208-4}.
#' @seealso \code{\link{star.phases}}, \code{\link{riseset}}, \code{\link{sketchSky}}, \code{\link{spatial.equinox}}, \code{\link{solar.date}}, \code{\link{EFM}}, \code{\link{planet.extremes}}
#' @return A list of candidate solar, lunar and stellar matches. The tables give target names and declination limits in degrees. Stellar matches also carry catalogue identifiers and visual magnitudes. Set dates to TRUE to request stellar date ranges. A category with no matches may be absent.
#' @param decrange Range of declination to consider.
#' @param timerange Earlier argument name for \code{year}, retained for existing scripts.
#' @param year Temporal range to consider. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' @param max.mag Maximum magnitude of stars to consider. Defaults to 2.5
#' @param loc FALSE for geocentric declinations, a vector of latitude, longitude and elevation, or a skyscapeR.horizon to also evaluate the local horizon.
#' @param dates Whether to display time period when identified celestial targets
#' occur within declination range. Default is FALSE.
#' @param calendar Calendar for civil dates. Use \code{"G"} for Gregorian or \code{"J"} for Julian.
#' @section Validation:
#' A08 verified closed-band membership, separate stellar passages, sampled contact, magnitude filtering and comparisons before display rounding. See the \href{https://stonexsky.github.io/skyscapeR/manual/#a08-target-ranges}{A08 check in the technical manual}.
#' @export
#' @examples
#' if (requireNamespace("swephRdata", quietly = TRUE)) {
#'   \dontrun{
#'   findTargets(c(-25,-17.5), c(-2500,-1750))
#'   
#'   # if a location is given then the zenith and anti-zenith sun will also be looked at:
#'   findTargets(c(3,12), c(-2500,-1750), loc=c(8.6, 7.3, 200))
#'   
#'   # if a horizon profile is given then the spatial equinox will also be looked at:
#'   hor <- malta_horizons[[1]]
#'   findTargets(c(-7,2), c(-2500,-1750), loc=hor)
#'   }
#' }
#' @section Candidate selection and visibility:
#' A match is a declination overlap, not evidence that a celestial object was
#' visible, deliberately observed or associated with the archaeological feature.
#' The magnitude filter uses fixed catalogue values; it includes integrated
#' magnitudes for non-stellar entries such as M45. It does not predict naked-eye
#' visibility through extinction or twilight, or resolve the members of a cluster.
#'
#' Choose the dating interval and measurement uncertainty from the archaeological
#' and survey evidence. Then examine candidate stars with \code{\link{star.phases}},
#' \code{\link{riseset}} and \code{\link{sketchSky}}. For a statistically selected
#' local declination band, candidate screening is exploratory; the significance
#' of that band does not identify its cause. Physical plausibility and independent
#' archaeological evidence are separate parts of the interpretation.
#' @section Spatial equinox epochs:
#' When a horizon is supplied, the spatial equinox is evaluated at both endpoints
#' of the requested year range (once for a single year), with the current atmospheric
#' settings. The candidate interval spans the finite sunrise and sunset declinations
#' at those endpoints. This is an endpoint envelope, not an exhaustive search for
#' intermediate changes of terrain crossing over very long periods. Undefined
#' crossings generate warnings; if all are undefined, the candidate is omitted.
#' @section Provenance:
#' Structured results carry a \code{provenance} attribute. Read it with
#' \code{attr(result, "provenance")}; its \code{packages} entry records the loaded
#' skyscapeR and swephR versions, with NA when a namespace version is unavailable.
#' This attribute adds no list entries or table columns. Ordinary subsetting or
#' transformations may discard attributes; save the original result for reference.
#' It does not identify ephemeris files, catalogue revisions or a development
#' commit. Retain those and the analysis script separately for exact reproduction.
findTargets <- function(decrange, year, loc=FALSE, max.mag=2.5, dates=FALSE,
                        calendar=skyscapeR.env$calendar, timerange) {
  if (!missing(timerange)) {
    if (!missing(year)) stop("Supply year only once.",call.=FALSE)
    year <- timerange
  }
  if (is.character(year)) year <- astroYear(year)
  targets <- c()

  years <- seq(min(year), max(year))
  seasons <- vapply(unique(range(years)), .season_dates, character(4), calendar=calendar)
  labels <- apply(seasons,1,function(x) paste(unique(x),collapse=" / "))

  ## solar
  targets$solar <- data.frame(name='Test', min.dec=12, max.dec=12, date1='1 Jan', date2='1 Feb', stringsAsFactors=F)
  targets$solar[1,] <- c('june solstice', range(jS(seq(min(year), max(year),1), loc, verbose=FALSE)), labels[2], NA)
  targets$solar[2,] <- c('december solstice', range(dS(seq(min(year), max(year),1), loc, verbose=FALSE)), labels[4], NA)
  targets$solar[3,] <- c('astronomical equinox', 0, 0, labels[1], labels[3])
  endpoints <- unique(range(years))
  # A solar target at the limit may be reached once, twice, or not at all
  # in the selected civil year. Retain only dates returned for that epoch.
  target.days <- function(dec, y) {
    z <- .try_target(solar.date(dec, y, calendar, verbose=FALSE))
    if (inherits(z,"try-error")) return(character())
    as.character(stats::na.omit(z[-1,1]))
  }
  endpoint.label <- function(x) {
    x <- unique(x[!is.na(x) & nzchar(x)])
    if (length(x)) paste(x,collapse=" / ") else NA_character_
  }
  if (inherits(loc,"skyscapeR.horizon")) {
    spatial <- lapply(endpoints,function(y) unlist(spatial.equinox(loc,year=y)$declination))
    limits <- unlist(spatial,use.names=FALSE)
    if (any(is.finite(limits))) {
      day.ranges <- lapply(seq_along(endpoints),function(i) {
        vals <- spatial[[i]]
        lapply(vals[is.finite(vals)],target.days,y=endpoints[i])
      })
      branches <- do.call(c,day.ranges)
      targets$solar[4,] <- c("spatial equinox",range(limits,finite=TRUE),
        endpoint.label(vapply(branches,function(z)z[1L],character(1))),
        endpoint.label(unlist(lapply(branches,function(z)tail(z,-1L)))))
    }
  }
  if (!is.logical(loc)) {
    lat <- extractLatitude(loc)
    attainable <- years[obliquity(years)>=abs(lat)]
    if (length(attainable)) {
      for (target in c("zenith sun","antizenith sun")) {
        decl <- if (target=="zenith sun") lat else -lat
        first.passage <- function(epochs) {
          for (y in epochs) {
            passages <- target.days(decl,y)
            if (length(passages)) return(passages)
          }
          character()
        }
        first <- first.passage(attainable)
        if (!length(first)) next
        last <- if (length(attainable)==1L) first else first.passage(rev(attainable))
        dates.by.year <- list(first,last)
        targets$solar[if(target=="zenith sun")5L else 6L,] <- c(target,decl,decl,
          endpoint.label(vapply(dates.by.year,function(z)z[1L],character(1))),
          endpoint.label(unlist(lapply(dates.by.year,function(z)tail(z,-1L)))))
      }
    }
  }

  samples <- unique(c(seq(min(decrange),max(decrange),0.01),max(decrange)))
  aux1 <- .try_target(solar.date(samples, min(year), calendar, verbose=FALSE))
  aux2 <- if (length(endpoints)==1L) aux1 else
    .try_target(solar.date(samples, max(year), calendar, verbose=FALSE))

  ttt <- c()
  if (class(aux1)[1] != 'try-error') { ttt <- c(ttt, aux1[-1,,drop=FALSE]) }
  if (class(aux2)[1] != 'try-error') { ttt <- c(ttt, aux2[-1,,drop=FALSE]) }
  ttt <- ttt[!is.na(ttt)]
  if (length(ttt)) {
    # A leap-year reference preserves 29 February and all month ends.
    parts <- strsplit(trimws(ttt)," +")
    day <- as.integer(vapply(parts,`[`,character(1),1L))
    month <- match(vapply(parts,`[`,character(1),2L),month.abb)
    template <- as.Date("2000-01-01") + 0:365
    ordinal <- sort(unique(match(sprintf("%02d-%02d",month,day),format(template,"%m-%d"))))
    groups <- split(ordinal,cumsum(c(TRUE,diff(ordinal)>2)))
    if (length(groups)>1L && 366L %in% ordinal && 1L %in% ordinal) {
      groups[[1]] <- c(groups[[length(groups)]],groups[[1]])
      groups <- head(groups,-1L)
    }
    windows <- vapply(groups,function(g) {
      endpoints <- as.POSIXlt(template[c(g[1],tail(g,1L))],tz="UTC")
      paste(paste(endpoints$mday,month.abb[endpoints$mon+1L]),collapse=" - ")
    },character(1))
    targets$solar[7,] <- c("sunrise/set",range(decrange),windows[1],
      if(length(windows)>1L) paste(windows[-1L],collapse=" / ") else NA_character_)
  }

  targets$solar[,2] <- as.numeric(targets$solar[,2])
  targets$solar[,3] <- as.numeric(targets$solar[,3])

  filter.limits <- function(tab) {
    keep <- !is.na(tab$name) & tab$min.dec <= max(decrange)+1e-10 &
      tab$max.dec >= min(decrange)-1e-10
    tab <- tab[which(keep), ,drop=FALSE]
    if (!nrow(tab)) return(NULL)
    tab$min.dec <- round(tab$min.dec,2)
    tab$max.dec <- round(tab$max.dec,2)
    rownames(tab) <- NULL
    tab
  }
  targets$solar <- filter.limits(targets$solar)
  lunar <- list("southern major lunar extreme"=sMjLX, "southern minor lunar extreme"=smnLX,
                "northern minor lunar extreme"=nmnLX, "northern major lunar extreme"=nMjLX)
  limits <- vapply(lunar,function(f)range(f(years,loc,verbose=FALSE)),numeric(2))
  targets$lunar <- filter.limits(data.frame(name=names(lunar),min.dec=limits[1,],max.dec=limits[2,]))

  ## stellar
  ss <- .read_star_catalogue()
  ss <- ss[which(ss$magV<=max.mag),,drop=FALSE] # magnitude filter

  targets$stellar <- data.frame(name='Test', identifier='Test', magV=1.2, min.dec=12, max.dec=12, daterange=NA, stringsAsFactors=F)[FALSE,,drop=FALSE]
  for (i in seq_len(NROW(ss))) {
    d <- star(ss$identifier[i],years)$coord$Dec
    bounds <- range(d)
    if (bounds[2]<min(decrange)-.stellar_band_tolerance || bounds[1]>max(decrange)+.stellar_band_tolerance) next
    interval <- if (dates) .stellar_dates(years,d,range(decrange)) else NA_character_
    targets$stellar[nrow(targets$stellar)+1L,] <- list(ss$name[i],ss$identifier[i],
      ss$magV[i],round(bounds[1],2),round(bounds[2],2),interval)
  }
  for (i in seq_len(NROW(targets$stellar))) {
    aux <- star.names$Western[which(star.names$identifier == targets$stellar[i,2])]
    if (length(aux)>0) { targets$stellar[i,1] <- aux } else { targets$stellar[i,1] <- '' }
  }

  rownames(targets$stellar) <- NULL
  if (NROW(targets$stellar) ==0) { targets$stellar <- c() }

  ## for "security" reasons
  if (dates==FALSE) targets$stellar$daterange <- NULL

  return(.with_astronomy_provenance(targets))
}
