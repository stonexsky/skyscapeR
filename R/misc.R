#' @noRd
checkYear <- function(year){
  if (sum(year < -13201 | year > 17191)>0) { stop('Year(s) outside of Swiss Ephemeris range (13201 BCE to 17191 CE).')}
}

#' @noRd
extractLatitude <- function(loc){
  if (is(loc,'skyscapeR.horizon')) { latitude <- loc$metadata$georef[1] }
  if (is(loc,'numeric')) { latitude <- loc[1]}
  # check
  if (latitude > 90 | latitude < -90) { stop('Latitude must be between -90 and 90 degrees.')}
  return(latitude)
}


#' @noRd
stars.pval <- function(p.value) pval2stars(p.value)

#' @noRd
dnorm <- function(x, mean, sd) { return(exp(-(x-mean)^2/(2*sd^2)) / sqrt(2*pi*sd^2)) }  ## replaces R function (twice as fast)

#' @noRd
tWS <- function(days, year) {
  WS <- findWS(year)$ind
  aux <- days-WS
  aux[aux<1] <- aux[aux<1]+365
  return(aux)
}


#' @noRd
eq2hor <- function(ra, dec, loc, refraction, atm, temp) {
  if (missing(refraction)) { refraction <- skyscapeR.env$refraction }
  if (missing(atm)) { atm <- skyscapeR.env$atm }
  if (missing(temp)) { temp <- skyscapeR.env$temp }

  jd <- swephR::swe_julday(2000,1,1,12,1)
  xx <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],loc[3]), atm, temp, c(ra,dec))$xaz

  if (refraction) { alt <- xx[3] } else { alt <- xx[2] }
  az <- xx[1]-180
  if (az < 0) { az <- az + 360 }
  if (az > 360) { az <- az - 360 }

  return(list(alt = alt, az = az))
}

#' @noRd
veq2hor <- Vectorize(eq2hor, 'ra', SIMPLIFY=T)

#' @noRd
minmaxdec = function(name, from, to, loc=FALSE) {
  xx <- seq(from, to, 100)
  dd <- c()

  # Stars
  sefstars <- .read_star_catalogue()

  if (sum(star.names$Western == name) | sum(as.character(sefstars$identifier) == name)) {
    return(range(star(name, seq(from,to))$coord$Dec))
  } else {
    # Solar-Lunar targets
    for (i in 1: NROW(xx)) {
      dd[i] <- do.call(name, list(year=xx[i], loc=loc, verbose=F))
    }
  }
  ff <- splinefun(xx,dd)
  xxz <- seq(from,to,0.01)
  dd <- ff(xxz)
  mm <- c(min(dd),max(dd))

  return(mm)
}

.stellar_dates <- function(years, declinations, band) {
  # Intersect each annual line segment with the closed declination band.
  if (length(years)==1L) return(paste(BC.AD(rep(years,2)),collapse=" - "))
  a <- head(declinations,-1L)
  b <- tail(declinations,-1L)
  idx <- which(pmin(a,b)<=band[2] & pmax(a,b)>=band[1])
  if (!length(idx)) return(NA_character_)
  slope <- b[idx]-a[idx]
  first <- last <- years[idx]
  last <- last+1
  moving <- slope!=0
  x <- (band[1]-a[idx][moving])/slope[moving]
  y <- (band[2]-a[idx][moving])/slope[moving]
  first[moving] <- years[idx[moving]]+pmax(0,pmin(x,y))
  last[moving] <- years[idx[moving]]+pmin(1,pmax(x,y))
  gaps <- which(tail(first,-1L)>head(last,-1L))
  starts <- first[c(1L,gaps+1L)]
  ends <- last[c(gaps,length(last))]
  paste(BC.AD(round(starts)),BC.AD(round(ends)),sep=" - ",collapse=" and ")
}

#' Convert degrees, arcminutes and arcseconds to decimal degrees
#' @description Convert a signed field angle in degrees, minutes and seconds to decimal degrees for
#' coordinate and orientation calculations.
#' @details The sign applies to the whole angle. The magnitude is the sum of absolute degrees,
#' minutes divided by 60 and seconds divided by 3600. With numeric components the first
#' nonzero component sets the sign and later components must be nonnegative. A string may
#' contain one to three colon- or space-separated components and an optional leading sign.
#' A string such as "-0:30" preserves the negative sign even when degrees are zero.
#' @seealso \code{\link{reduct.compass}}, \code{\link{reduct.theodolite}}, \code{\link{mag.dec}}, \code{\link{sunAz}}, \code{\link{az2dec}}
#' @param dd Degrees as a finite numeric scalar, or a single character string
#'   containing an angle, for example \code{"-8:30:00"}.
#' @param mm Arcminutes as a finite numeric scalar. Defaults to zero.
#' @param ss Arcseconds as a finite numeric scalar. Defaults to zero.
#' @return A numeric scalar giving the angle in signed decimal degrees.
#' @export
#' @examples
#' ten(24, 52, 16)
#' ten(-8, 30, 0)     # -8.5 degrees
#' ten(0, -30, 0)     # -0.5 degrees
#' ten(0, 0, -30)     # -1/120 degrees
#' ten("-8:30:00")    # -8.5 degrees
#' ten("-0:30:00")    # Preserve the minus sign when degrees are zero
ten <- function(dd, mm=0, ss=0) {
  if (is.character(dd)) {
    if (!missing(mm) || !missing(ss)) {
      stop('Supply only dd when using a character angle.', call.=FALSE)
    }
    if (length(dd) != 1L || is.na(dd) || !nzchar(trimws(dd))) {
      stop('dd must contain a single non-empty angle.', call.=FALSE)
    }
    angle <- trimws(dd)
    direction <- if (startsWith(angle, '-')) -1 else 1
    angle <- sub('^[+-]', '', angle)
    # A sign is allowed only at the beginning of the complete string.
    if (!grepl('^[0-9.]+([[:space:]]*:[[:space:]]*[0-9.]+|[[:space:]]+[0-9.]+){0,2}$', angle)) {
      stop('Use one to three unsigned components after an optional leading sign.', call.=FALSE)
    }
    components <- suppressWarnings(as.numeric(strsplit(angle, '[:[:space:]]+')[[1]]))
    if (any(!is.finite(components))) {
      stop('Angle components must be finite numbers.', call.=FALSE)
    }
    components <- c(components, rep(0, 3L-length(components)))
  } else {
    inputs <- list(dd, mm, ss)
    valid <- vapply(inputs, function(x) is.numeric(x) && length(x)==1L &&
                      is.finite(x), logical(1))
    if (!all(valid)) stop('dd, mm and ss must be finite numeric scalars.', call.=FALSE)
    components <- c(dd, mm, ss)
    nonzero <- which(components != 0)
    direction <- if (length(nonzero)) sign(components[nonzero[1]]) else 1
    if (length(nonzero) > 1L && any(components[nonzero[-1]] < 0)) {
      stop('Only the first non-zero component may be negative.', call.=FALSE)
    }
  }
  direction * sum(abs(components) / c(1, 60, 3600))
}


#' @noRd
lty2dash <- function(lty){
  if (lty==1) { return('solid') }
  if (lty==2) { return('dash') }
  if (lty==3) { return('dot') }
  if (lty==4) { return('dashdot') }
  if (lty==5) { return('longdash') }
  if (lty==6) { return('longdashdot') }
}

#' @noRd
epoch2yr <- function(epoch){
  aux <- as.numeric(sub(" .*", "", epoch))
  sign <- sub(".* ", "", epoch)
  if (sign=='BCE') { aux <- -aux}
  return(aux)
}


#' @noRd
checkbody <- function(body){
  body <- toupper(body)
  SE <- get(data('SE', package='swephR', envir = environment())); nn <- names(SE[3:13]); nn[1] <- 'ECLIPTIC'
  ind <- which(nn==body)+2
  if (length(ind)==0) { stop('Solar System body not recognised. Please check name.')}
  body <- as.numeric(SE[ind])
  return(body)
}



#' @noRd
mag2size <- function(mag) {
  # size <- -log(mag+4.932) + 2.396 # good for upto mag 6
  size <- 4.45-1.57*log(mag+10) # good for upto mag 7
  if (is.na(mag) | mag > 7) { size <- 0 }
  return(size)
}



#' Estimate magnetic declination at a site
#' @description Estimate the correction between magnetic north and true north for the location and date
#' of a field survey.
#' @details Compass reduction depends on when and where the reading was taken. The function
#' evaluates IGRF-14 through oce at sea level and returns east-positive declination.
#' Supply a single site for several dates, a single date for several sites, or
#' one date per site. Results follow input order. Other size combinations raise
#' an error. The calendar date is converted to a decimal year using the actual
#' length of that year at midnight UTC. Supported dates run from 1 January 1900
#' through 1 January 2030 inclusive. Site elevation is ignored by this model call.
#' @return Numeric magnetic declination in degrees. Use it as the correction added to a magnetic bearing in \code{\link{reduct.compass}}.
#' @param loc Latitude and longitude in degrees as a vector, matrix or data frame with one row per site, or a horizon object. An optional third column gives elevation.
#' @param date One or more survey dates as Date objects or YYYY/MM/DD strings.
#' @export
#' @references Kelley, D. and Richards, C. oce: magneticField. \url{https://dankelley.github.io/oce/reference/magneticField.html}.
#' @seealso \code{\link[oce]{magneticField}}, \code{\link{reduct.compass}}
#' @examples
#' # Magnetic Declination for London on April 1st 2016:
#' loc <- c( 51.5074, -0.1278 )
#' mag.dec( loc, "2016/04/01" )
#' mag.dec(c(40, -8), c("2010/01/01", "2026/01/01"))
mag.dec = function(loc, date) {
  if (inherits(loc, "skyscapeR.horizon")) loc <- loc$metadata$georef
  loc <- if (is.null(dim(loc))) matrix(loc,nrow=1L) else as.matrix(loc)
  if (!is.numeric(loc) || !ncol(loc) %in% 2:3 || nrow(loc)==0L ||
      any(!is.finite(loc[,1:2])) || any(abs(loc[,1])>90 | abs(loc[,2])>180))
    stop("Supply latitude and longitude in degrees, one site per row.")
  if (!inherits(date,"Date")) {
    if (!is.character(date) || anyNA(date) || any(!grepl("^[0-9]{4}/[0-9]{1,2}/[0-9]{1,2}$",date)))
      stop("Supply survey dates as Date objects or YYYY/MM/DD strings.")
    date <- as.Date(as.vector(date),format="%Y/%m/%d")
  }
  if (!length(date) || anyNA(date) || any(date<as.Date("1900-01-01") | date>as.Date("2030-01-01")))
    stop("Provide valid survey dates from 1900/01/01 through 2030/01/01 for IGRF-14.")
  n <- max(nrow(loc),length(date))
  if (!nrow(loc) %in% c(1L,n) || !length(date) %in% c(1L,n))
    stop("Supply one site, one date, or one date per site.")
  year <- as.integer(format(date,"%Y"))
  jan <- as.Date(paste0(year,"-01-01"))
  next.jan <- as.Date(paste0(year+1L,"-01-01"))
  epoch <- year + as.numeric(date-jan)/as.numeric(next.jan-jan)
  oce::magneticField(rep(loc[,2],length.out=n),rep(loc[,1],length.out=n),
                     rep(epoch,length.out=n),version=14)$declination
}

#' @noRd
#' @export
subset.skyscapeR.pdf <- function(x, subset, ...) {
  x$metadata$name <- x$metadata$name[subset]
  x$metadata$az <- x$metadata$az[subset]
  x$metadata$unc <- x$metadata$unc[subset]
  x$metadata$horizon <- x$metadata$horizon[subset]
  x$data <- x$data[subset]
  x
}


#' @noRd
sampleList <- function(ll) {
  vapply(ll, function(tab) {
    row <- sample.int(nrow(tab), 1L, prob=tab[[2L]])
    tab[[1L]][row]
  }, numeric(1), USE.NAMES=FALSE)
}

#' Find highest-density intervals for an orientation distribution
#' @description Identify the angular regions containing a chosen proportion of a measured or summed
#' orientation distribution at its highest densities.
#' @details A distribution with several peaks can give several disjoint intervals. The mass
#' parameter applies to the selected regions together; the result describes concentration
#' of the supplied density rather than significance against a null model.
#'
#' Uses the existing discrete density-threshold method: sort sampled
#' densities, find the threshold reaching the requested fraction of their sum,
#' and include all points tied at that threshold. Consequently, the selected
#' mass can exceed the requested value. Densities need not sum to one.
#'
#' Coordinates must be finite and strictly increasing on an equally spaced
#' grid. Missing grid points are permitted when gaps are integer multiples of
#' the smallest observed step; omitted points are treated as zero mass and
#' split intervals. This accommodates cutoff-trimmed azimuth distributions.
#' Other irregular grids are rejected; no unequal-width integration is done.
#' A single point is permitted and returns a zero-width interval.
#'
#' Densities must be finite, nonnegative and have a positive total.
#' If a saved distribution contains negative densities, regenerate it from the inputs with \code{\link{coordtrans}} before estimating intervals.
#' Endpoints are sampled coordinates, not interpolated boundaries. Regions
#' crossing the azimuth wrap are reported separately at the two ends of the
#' grid; they are not joined into a circular interval.
#'
#' Hyndman (1996) defines a highest-density region by a common density threshold.
#' The discrete grid and tie rule above implement that principle for sampled curves.
#' @param x A \emph{skyscapeR.spd} or \emph{skyscapeR.pdf} object,
#' or one data frame/list with numeric components \code{x} and \code{y}.
#' @param mass Probability mass, greater than zero and at most one.
#' Default is 0.954.
#' @return For a PDF object, a list of interval matrices, one per distribution,
#' even for a single distribution. For SPD or x/y input, one matrix.
#' Each matrix has columns \code{start} and \code{end}, with ordered endpoints
#' and one row per contiguous region, sorted by start.
#' @export
#' @references Hyndman, R. J. (1996). Computing and Graphing Highest Density Regions. The American Statistician 50, 120-126. \doi{10.1080/00031305.1996.10474359}.
#' @seealso \code{\link{az.pdf}}, \code{\link{dec.pdf}}, \code{\link{spd}}, \code{\link{plot.skyscapeR.pdf}}
#' @examples
#' xx <- seq(-6, 6, by=0.01)
#' hpdi(data.frame(x=xx, y=stats::dnorm(xx)), mass=0.95)
#'
#' # Two separated peaks produce two ordered intervals:
#' yy <- stats::dnorm(xx, -3, 0.5) + stats::dnorm(xx, 3, 0.5)
#' hpdi(data.frame(x=xx, y=yy), mass=0.8)
#'
#' # PDF objects always return a list:
#' Az <- az.pdf(az=c(90,120), unc=2, verbose=FALSE)
#' hpdi(Az)
hpdi <- function(x, mass=0.954) {
  if (!is.numeric(mass) || length(mass) != 1L || !is.finite(mass) ||
      mass <= 0 || mass > 1) {
    stop("mass must be a finite number greater than zero and at most one.")
  }
  region <- function(grd) {
    if (!is.list(grd) || !is.numeric(grd$x) || !is.numeric(grd$y) ||
        !length(grd$x) || length(grd$x) != length(grd$y) ||
        any(!is.finite(grd$x)) || any(!is.finite(grd$y)) ||
        any(grd$y < 0) || !any(grd$y > 0)) {
      stop("Each distribution must contain finite numeric x/y vectors of equal nonzero length, with nonnegative densities and positive total mass.")
    }
    delta <- diff(grd$x)
    if (any(delta <= 0)) stop("x coordinates must be strictly increasing.")
    if (length(delta)) {
      step <- min(delta)
      multiples <- delta / step
      # Floating-point tolerance relative to the coordinate/grid scale.
      tolerance <- 1e-7 + 100 * .Machine$double.eps * max(1, abs(grd$x)) / step
      if (tolerance >= 0.01 || any(abs(multiples-round(multiples)) > tolerance)) {
        stop("x must lie on an equally spaced grid (integer-multiple gaps are allowed).")
      }
    }
    # Scaling avoids overflow without changing relative discrete mass.
    weights <- grd$y / max(grd$y)
    sorted <- sort(weights, decreasing=TRUE)
    cumulative <- cumsum(sorted)
    threshold <- sorted[which(cumulative >= mass * tail(cumulative, 1L))[1L]]
    indices <- which(weights >= threshold)
    breaks <- diff(indices) > 1
    if (length(indices) > 1L) {
      breaks <- breaks | diff(grd$x[indices]) > step * (1 + tolerance)
    }
    gaps <- which(breaks)
    starts <- indices[c(1L, gaps + 1L)]
    ends <- indices[c(gaps, length(indices))]
    cbind(start=grd$x[starts], end=grd$x[ends])
  }
  if (inherits(x, 'skyscapeR.pdf')) {
    if (!is.list(x$data) || !length(x$data)) stop("PDF object contains no distributions.")
    return(lapply(x$data, region))
  }
  if (inherits(x, 'skyscapeR.spd')) return(region(x$data))
  region(x)
}



#' Format a p-value as a significance symbol
#' @description Annotate a significance-test result with the conventional symbol for its p-value range.
#' @details Use a symbol for compact plot labels while retaining the numeric p-value for
#' interpretation. The cut points are 0.001, 0.01, 0.05 and 0.1, with symbols ***, **, *, +
#' and an open circle in order of increasing p. An equality belongs to the interval ending at that cut
#' point. For legacy strings beginning with a less-than sign and space, the function
#' classifies the printed bound. It does not calculate a test or adjust for multiple
#' comparisons.
#' @seealso \code{\link{print.skyscapeR.sigTest}}, \code{\link{plot.skyscapeR.sigTest}}
#' @examples pval2stars(c(0.0005, 0.008, 0.03, 0.08, 0.2))
#' @return A character vector of significance symbols: ***, **, *, + or an open circle, using cut points 0.001, 0.01, 0.05 and 0.1.
#' @param p.value p-value
#' @export
pval2stars <- function(p.value) {
  if (is(p.value,'character')) { p.value <- as.numeric(substr(p.value, 3,nchar(p.value))) }
  out <- symnum(p.value, corr = FALSE, na = FALSE,
                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                symbols = c("***", "**", "*", "+", "\u25cb"), legend=F)
  return(as.character(out))
}

#' @noRd
period <- function(obj) {
  if (obj=='Mercury') { return(6.989733) }
  if (obj=='Venus') { return(8) }
  if (obj=='Mars') { return(687/365.25*10) }
  if (obj=='Jupiter') { return(4331/365.25) }
  if (obj=='Saturn') { return(10747/365.25) }
}

#' @noRd
constLines <- function(const){

  if (is.data.frame(const)) {
    if (!all(c("name1","HIP") %in% names(const)) || !nrow(const))
      stop("Supply a non-empty constellation table from importConstellations().",call.=FALSE)
    cc <- as.character(const$HIP)
    valid <- vapply(cc,function(x) {
      if (is.na(x)) return(FALSE)
      if (!nzchar(trimws(x))) return(TRUE)
      z <- strsplit(trimws(x),"\\s+")[[1L]]
      n <- suppressWarnings(as.numeric(z))
      length(z)>=2L && length(z) %% 2L==0L && all(grepl("^[0-9]+$",z)) &&
        all(is.finite(n)) && all(n>0) && all(n<=.Machine$integer.max)
    },logical(1))
    if (!all(valid)) stop("Each HIP entry must contain pairs of positive integer identifiers.",call.=FALSE)
  } else {
    if (!is.character(const) || !length(const) || anyNA(const))
      stop("Supply constellation names or an imported constellation table.",call.=FALSE)
    matches <- lapply(const,function(x) which(apply(constellations[,1:5],1L,
      function(row) any(row==x,na.rm=TRUE))))
    unknown <- const[lengths(matches)==0L]
    ind <- unique(unlist(matches))
    if (!length(ind)) stop("Constellation not identified.",call.=FALSE)
    if (length(unknown)) warning("Constellation(s) not identified: ",
      paste(unknown,collapse=", "),". Drawing the recognised figures only.",call.=FALSE)
    selected <- constellations[ind,]
    selected <- selected[!duplicated(selected[,c("name2","HIP")]),]
    cc <- selected$HIP
  }
  cc <- cc[nzchar(trimws(cc))]
  ct <- list()
  for (i in seq_along(cc)){
    ct[[i]] <- unlist(strsplit(trimws(cc[i]), "\\s+"))
    ct[[i]] <- paste('HIP', ct[[i]])
  }
  rm(cc)

  ss <- .plot_star_catalogue()
  # Line endpoints identify HIP systems. Prefer the system entry, then component A.
  cat <- trimws(ss$cat)
  hip <- grepl("^HIP [0-9]+[A-Z]?$",cat)
  cat <- cat[hip]
  identifiers <- ss$identifier[hip]
  order.by.component <- order(grepl("[A-Z]$",cat),sub("^HIP [0-9]+","",cat))
  systems <- sub("[A-Z]$","",cat[order.by.component])
  identifiers <- identifiers[order.by.component]
  lapply(ct,function(endpoints) {
    rows <- match(endpoints,systems)
    for (missing in unique(endpoints[is.na(rows)]))
      warning("Could not draw constellation line ",missing,".",call.=FALSE)
    matrix(identifiers[rows],ncol=2,byrow=TRUE)
  })
}

.sky_message <- function(...) {
  if (isTRUE(getOption("skyscapeR.verbose",TRUE))) message(...)
}

.sky_progress <- function(verbose) {
  interactive() && isTRUE(verbose) && isTRUE(getOption("skyscapeR.verbose",TRUE))
}
