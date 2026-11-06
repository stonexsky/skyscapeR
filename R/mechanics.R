#' Convert azimuth and altitude to declination
#' @description Find the celestial declination corresponding to a terrestrial sightline at a known
#' latitude and horizon altitude.
#' @details The same true bearing can indicate different declinations when the observer or horizon
#' changes. Supply angular altitude directly or interpolate it from one shared horizon or
#' one horizon per observation. Numeric latitude inputs imply zero longitude and elevation;
#' a location matrix can include both. When refraction is enabled, the function adjusts the
#' input altitude using pressure, temperature and elevation, then calls the inverse Swiss
#' Ephemeris horizontal transformation. Bearings are clockwise from north; the routine
#' shifts them to the engine convention and rounds the returned declinations using the
#' input bearing representation.
#' @return Numeric declination in degrees for the supplied azimuth and altitude values.
#' @param az True bearing or vector of bearings in degrees clockwise from north.
#' @param loc One latitude or one per bearing (degrees north), a site matrix with one row per bearing
#' (latitude, longitude and optional elevation), a shared horizon or a list of horizons in measurement order.
#' @param alt Horizon altitude in degrees, shared or one per bearing. Omit when supplying
#' horizons to interpolate altitude at each bearing with \code{\link{hor2alt}}.
#' @param refraction Logical. Apply atmospheric refraction. Defaults to the session setting in \code{\link{skyscapeR.vars}}.
#' @param atm Atmospheric pressure in mbar. Defaults to the session setting in \code{\link{skyscapeR.vars}}.
#' @param temp Air temperature in degrees Celsius. Defaults to the session setting in \code{\link{skyscapeR.vars}}.
#' @import swephR
#' @export
#' @seealso \code{\link[swephR]{swe_azalt_rev}}, \code{\link{hor2alt}}
#' @examples
#' dec <- az2dec(az=92, loc=35, alt=2)
#'
#' dec <- az2dec(az=c(92,89,102), loc=35, alt=c(0,2,4))
#'
#' # flat horizon with 2 degrees of altitude
#' hor <- createHor(az=c(0,360), alt=c(2,2), loc=c(35,-8,25))
#' dec <- az2dec(92, loc=hor)
#'
#' # Can also be used for an array of azimuths:
#' decs <- az2dec(az=c(87,92,110), loc=hor)
az2dec = function(az, loc, alt, refraction=skyscapeR.env$refraction, atm=skyscapeR.env$atm, temp=skyscapeR.env$temp){
  jd <- swephR::swe_julday(2000,1,1,12,1)
  prec <- max(nchar(sub('.*\\.', '', as.character(az))))

  az[which(az < 0)] <- 360 + az[which(az < 0)]
  az[which(az > 360)] <- az[which(az > 360)] - 360

  if (missing(alt)) {
    if (class(loc)[1] == 'skyscapeR.horizon') {
      alt <- hor2alt(loc, az)
    } else if (class(loc)[1] == 'list') {
      alt <- c()
      for (i in 1:NROW(loc)) {
        alt[i] <- hor2alt(loc[[i]], az[i])
      }
    } else { stop('Altitude missing.') } }

  if (length(alt) == 1) { alt <- rep(alt, NROW(az)) }
  if (length(alt) != NROW(az)) { stop('alt parameter needs to be of either length 1 or same length as az') }

  if (class(loc)[1] == 'skyscapeR.horizon') { georefs <- loc$metadata$georef; georefs <- matrix(georefs, ncol=3, nrow=NROW(az), byrow=TRUE) }
  if (class(loc)[1] == 'list') {
    georefs <- matrix(NA, ncol=3, nrow=NROW(loc))
    for (i in 1:NROW(loc)) { georefs[i,] <- loc[[i]]$metadata$georef } }
  if (class(loc)[1] == 'numeric') {
    if (length(loc) != 1 & length(loc) != NROW(az)) { stop('loc parameter needs to be of either length 1 or same length as az') }
    if (length(loc) == 1) { loc <- rep(loc, NROW(az)) }
    # if (length(loc) > 1) {
      georefs <- matrix(NA, ncol=3, nrow=NROW(loc))
      for (i in 1:NROW(loc)) { georefs[i,] <- c(loc[i],0,0) }
    # }
  }
  if (class(loc)[1] == 'matrix') {
    georefs <- matrix(NA, ncol=3, nrow=NROW(loc))
    for (i in 1:NROW(loc)) {
      if (dim(loc)[2]==2) { georefs[i,] <- c(loc[i,1], loc[i,2], 0) }
      if (dim(loc)[2]==3) { georefs[i,] <- loc[i,] }
    }
  }

  dec <- c()
  for (i in 1:NROW(az)) {
    georef <- georefs[i,]
    if (refraction) {
      alt[i] <- alt[i] - (swephR::swe_refrac_extended(alt[i], georef[3], atm, temp, 0, 0)$return - alt[i])
    }
    dec[i] <- round( swephR::swe_azalt_rev(jd, 1, c(georef[2],georef[1],georef[3]), c(az[i]-180, alt[i]))$xout[2], prec)
  }

  return(dec)
}



#' Calculate the obliquity of the ecliptic
#' @description Estimate the tilt of the ecliptic at a proposed epoch to set solar and lunar declination
#' limits.
#' @details The tilt changes over archaeological timescales. For each supplied astronomical year the
#' calculation uses Gregorian 1 January at noon UT, requests the obliquity/nutation record
#' from Swiss Ephemeris and returns its first element, the true obliquity in degrees. An
#' annual value is a reference for limiting declinations rather than a search for the
#' instant of solstice.
#' @return Numeric obliquity of the ecliptic in degrees, one value per astronomical year.
#' @param year Year for which to calculate the obliquity. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' Defaults to present year as given by \emph{Sys.Date}
#' @import swephR
#' @export
#' @references Laskar, J. et al. (2004), A long-term numerical
#' solution for the insolation quantities of the Earth, \emph{Astron.
#'  Astroph.}, 428, 261-285, doi:10.1051/0004-6361:20041335.
#' @seealso \code{\link[swephR]{swe_calc_ut}}, \code{\link{jS}}, \code{\link{dS}}
#' @examples
#' # Obliquity for year 4001 BC:
#' obliquity("4001 BC")
obliquity = function(year = skyscapeR.env$cur.year) {
  if (is.character(year)) year <- astroYear(year)
  aux <- c()
  for (i in 1:length(year)) {
    jd <- swephR::swe_julday(year[i],1,1,12,1)
    aux[i] <- swephR::swe_calc_ut(jd, -1, 0)$xx[1]
  }
  return(aux)
}

#' @noRd
body.position.unvec = function(obj='sun', time, timezone, calendar, dec, loc = NULL, refraction, atm, temp, verbose=TRUE) {

  if (is.null(loc)) {
    if (verbose) { .sky_message('No location given. Forcing output to equatorial coordinates in the geocentric frame of reference.') }
    dec <- 'geo'
    nohoriz <- TRUE
  } else { nohoriz <- FALSE }

  if (inherits(loc, 'skyscapeR.horizon')) { loc <- loc$metadata$georef }

  body <- checkbody(obj)

  coords.eq <- data.frame(RA=NA, Dec=NA)
  if (!nohoriz) { coords.hor <- data.frame(az=NA, alt=NA) }

  if (class(time)[1]=='character') {
    checkYear(year(time))
    jd <- time2jd(time, timezone, calendar)
  } else { jd <- time }
  if (dec == 'geo') {
    coords.eq[1,] <- .file_position(jd, body, 2048)$xx[1:2]

  } else if (dec == 'topo') {
    swephR::swe_set_topo(loc[2],loc[1],loc[3])
    coords.eq[1,] <- .file_position(jd, body, 2048+32*1024+16)$xx[1:2]
  }

  if (!nohoriz) {
    aux <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],loc[3]), atm, temp, as.numeric(coords.eq))$xaz

    if (refraction) {
      coords.hor[1,] <- aux[c(1,3)] ## apparent altitude
      coords.eq[1,] <- swephR::swe_azalt_rev(jd, 1, c(loc[2],loc[1],loc[3]), aux)$xout[1:2] ## apparent declination
    } else { coords.hor[1,] <- aux[1:2] }
    coords.hor[1,1] <- coords.hor[1,1] - 180
    if (coords.hor[1,1] > 360) { coords.hor[1,1] <- coords.hor[1,1] - 360 }
    if (coords.hor[1,1] < 0) { coords.hor[1,1] <- coords.hor[1,1] + 360 }
  }

  out <- c()
  out$equatorial <- coords.eq
  if (!nohoriz) { out$horizontal <- coords.hor }
  return(out)
}

#' @noRd
body.position.vec <- Vectorize(body.position.unvec, 'time', SIMPLIFY = FALSE)


#' Find the equatorial and horizontal position of a celestial body
#' @description Locate a Solar System body in equatorial coordinates at a specified instant, or
#' determine where it appears in the observer's sky.
#' @details An alignment comparison needs an epoch, a reference frame and, for a local sightline, a
#' site. Character times are converted to UT Julian dates using the selected calendar and
#' timezone; numeric times are already UT Julian dates. Swiss Ephemeris supplies right
#' ascension and declination. Without a site the function uses the geocentric frame. With a
#' site it can use a topocentric position and convert to local azimuth and altitude. A
#' vector of times is evaluated one instant at a time and returned in input order.
#'
#' Only the horizon object's location is used; its altitude profile is not used
#' to determine whether the body is visible above the terrain. With a location,
#' horizontal altitude is apparent when refraction is TRUE and geometric when
#' FALSE. Azimuth is measured clockwise from north.
#'
#' The Swiss Ephemeris manual documents planetary coordinates and reference frames.
#'
#' Astronomical positions require Swiss Ephemeris files covering the epoch. We
#' recommend installing \pkg{swephRdata}; see \code{\link{skyscapeR-package}} for setup.
#' @param obj String containing name of the solar system body of interest. Can be
#' any of the planets (inc. Pluto), the Moon, the Sun or the Ecliptic. Defaults to 'sun'.
#' @param time Either a string containing the date and time in the format "YYYY/MM/DD HH:MM:SS"
#'  (see \code{\link{timestring}}), or a numeric containing the julian date (see \code{\link{time2jd}}).
#' @param timezone Timezone for civil timestamps. Use UTC or a fixed offset for dates before 1900. See \code{\link{time2jd}} for named zones and offsets.
#' @param calendar Calendar for civil dates. Use \code{"G"} for Gregorian or \code{"J"} for Julian.
#' @param dec Output declination: \emph{geo} for the geocentric, or \emph{topo} for the topocentric
#' frame of reference. If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param loc Observer location: a numeric vector containing latitude and
#' longitude in degrees and elevation in metres, in that order, or a
#' \emph{skyscapeR.horizon} object whose metadata supplies those coordinates.
#' Required for topocentric coordinates and for horizontal output in either
#' frame. With no location, only geocentric equatorial coordinates are returned.
#' @param refraction Logical. Apply atmospheric refraction. Defaults to the session setting in \code{\link{skyscapeR.vars}}.
#' @param atm Atmospheric pressure in mbar. Defaults to the session setting in \code{\link{skyscapeR.vars}}.
#' @param temp Air temperature in degrees Celsius. Defaults to the session setting in \code{\link{skyscapeR.vars}}.
#' @param verbose Logical. Display progress or informational messages. Default \code{TRUE}.
#' @return A list containing an \code{equatorial} data frame with columns
#' \code{RA} and \code{Dec}, in degrees. When \code{loc} is supplied, the list
#' also contains a \code{horizontal} data frame with columns \code{az} and
#' \code{alt}, in degrees. Each data frame has one row per input time.
#' @import swephR
#' @export
#' @references Astrodienst. Swiss Ephemeris programmer's manual. \url{https://www.astro.com/swisseph/swephprg.htm}.
#' @seealso \code{\link[swephR]{swe_calc_ut}}, \code{\link{timestring}}, \code{\link{time2jd}}
#' @examples
#' if (requireNamespace("swephRdata", quietly = TRUE)) {
#'   # Position of the sun at noon GMT on Christmas day 2018:
#'   body.position(obj='sun', time='2018/12/25 12:00:00', timezone='GMT')
#'   
#'   # Declination of the moon at same time
#'   body.position(obj='moon', time='2018/12/25 12:00:00', timezone='GMT')$equatorial$Dec
#'   
#'   # Geometric position using a horizon object's observer location:
#'   hor <- createHor(az=c(0,360), alt=c(2,2), loc=c(51.5,-0.1,25))
#'   body.position(obj='sun', time='2018/12/25 12:00:00', timezone='UTC',
#'   loc=hor, refraction=FALSE)$horizontal
#' }
body.position <- function(time, obj='sun', loc=NULL, timezone, calendar, dec, refraction, atm, temp, verbose=T) {
  if (missing(timezone)) { timezone <- skyscapeR.env$timezone }
  if (missing(calendar)) { calendar <- skyscapeR.env$calendar }
  if (missing(dec)) { dec <- skyscapeR.env$dec }
  if (missing(refraction)) { refraction <- skyscapeR.env$refraction }
  if (missing(atm)) { atm <- skyscapeR.env$atm }
  if (missing(temp)) { temp <- skyscapeR.env$temp }

  if (length(time)==1) {
    return(body.position.unvec(obj, time, timezone, calendar, dec, loc, refraction, atm, temp, verbose))
  } else {
    aux <- body.position.vec(obj, time, timezone, calendar, dec, loc, refraction, atm, temp, verbose)
    out <- c()
    out$equatorial <- as.data.frame(matrix(unlist(lapply(aux,'[[', 'equatorial')), ncol=2, byrow=TRUE)); names(out$equatorial) <- names(aux[[1]]$equatorial)
    if (!is.null(loc)) { out$horizontal <- as.data.frame(matrix(unlist(lapply(aux,'[[', 'horizontal')), ncol=2, byrow=T)); names(out$horizontal) <- names(aux[[1]]$horizontal) }
    return(out)
  }
}



#' Calculate the illuminated fraction of the Moon
#' @description Determine how much of the lunar disc is illuminated at the time of an observation or
#' proposed lunar event.
#' @details Lunar illumination can be used to screen events for a nearly full Moon. The routine
#' converts calendar input to UT Julian dates when needed and obtains the illuminated
#' fraction from the lunar phenomena calculation in Swiss Ephemeris. Values near zero
#' indicate a nearly new Moon and values near one a nearly full Moon; the fraction alone
#' cannot distinguish waxing from waning.
#'
#' Astronomical positions require Swiss Ephemeris files covering the epoch. We
#' recommend installing \pkg{swephRdata}; see \code{\link{skyscapeR-package}} for setup.
#' @return Numeric illuminated fraction of the lunar disc, one value per input time. At full illumination the fraction is 1.
#' @param time One or more UT Julian dates, or local date-time strings in the form YYYY/MM/DD HH:MM:SS.
#' @param timezone Timezone for civil timestamps. Use UTC or a fixed offset for dates before 1900. See \code{\link{time2jd}} for named zones and offsets.
#' @param calendar Calendar for civil dates. Use \code{"G"} for Gregorian or \code{"J"} for Julian.
#' @import swephR
#' @export
#' @seealso \code{\link[swephR]{swe_pheno_ut}}, \code{\link{EFM}}
#' @examples
#' if (requireNamespace("swephRdata", quietly = TRUE)) {
#'   # Moonphase at noon GMT on Christmas day 2018:
#'   moonphase('2018/12/25 12:00:00', 'GMT')
#' }
moonphase <- function(time, timezone, calendar) {
  if (missing(timezone)) { timezone <- skyscapeR.env$timezone }
  if (missing(calendar)) { calendar <- skyscapeR.env$calendar }

  aux <- c()
  for (i in 1:length(time)) {
    if (class(time)[1]=='character') { jd <- time2jd(time[i], timezone, calendar) } else { jd <- time[i] }
    .file_position(jd, 1L, 2L)
    .file_position(jd, 0L, 2L)
    phase <- swephR::swe_pheno_ut(jd, 1L, 2L)
    if (phase$return < 0L || bitwAnd(phase$return,7L)!=2L)
      .ephemeris_failure(jd,"Moon",phase$serr)
    aux[i] <- phase$attr[2]
  }

  return(aux)
}


#' Find rising and setting times and bearings at a site
#' @description Find when and where a star or Solar System body first rises above or sets below the horizon at
#' an observing site.
#' @details The event is a crossing of the local horizon by the body's centre. A raised ridge can
#' delay the first visible rise and an uneven horizon can give several appearances or
#' disappearances. Choose the horizon model, atmosphere, calendar and time convention to
#' match the observation.
#'
#' Events refer to the centre of the disc. For a Solar System body with a coordinate vector, alt is
#' one constant horizon altitude and the Swiss Ephemeris rise/set solver is used.
#' For stars, the daily track is evaluated against the supplied constant altitude or site skyline. With a horizon object, its profile is linearly interpolated in azimuth and
#' compared with the body's topocentric altitude; do not supply alt as well.
#' Refraction TRUE uses apparent altitude; FALSE uses geometric altitude. The
#' constant-altitude solver likewise disables refraction only when FALSE.
#' For a profile, search the 48 hours after each requested local midnight and
#' return the first upward and first downward crossing separately. They need
#' not enclose one continuous visibility interval or occur on the requested day.
#' Missing crossings return NA fields. A one-minute base grid, quarter-step samples,
#' azimuth turning points and recursively subdivided rapid azimuth motion bracket visits
#' to the profile's corners. Each smooth segment is also searched for an interior altitude
#' extremum, so short emergence/reoccultation pairs need not straddle time samples.
#' A pure tangency without a sign change is not a rise or set. Horizon altitudes at
#' 0 and 360 degrees must agree. Root tolerance is 0.01 seconds; interior optimization
#' uses relative seconds to avoid precision loss from large Julian dates.
#' These settings are not a physical accuracy guarantee: finite precision, unresolved
#' terrain and near-zenith azimuth singularities limit numerical completeness. The dec argument
#' controls returned declination, not the topocentric event search.
#'
#' The Swiss Ephemeris manual documents the fixed-altitude rise/set calculation.
#'
#' Astronomical positions require Swiss Ephemeris files covering the epoch. We
#' recommend installing \pkg{swephRdata}; see \code{\link{skyscapeR-package}} for setup.
#' @param obj Solar System body, Western star name or fixed-star catalogue identifier (including M 45). Defaults to 'sun'.
#' @param date Either a string containing the date in the format "YYYY/MM/DD"
#'  (see \code{\link{timestring}}), or a numeric containing the julian date (see \code{\link{time2jd}}).
#'  Can also be a single year ("YYYY") or a month ("YYYY/MM") to calculate risings and settings for every day in the
#'  year or month, respectively. Not necessary if \emph{jd} is given.
#' @param jd A numeric containing the julian date (see \code{\link{time2jd}}) for which to
#' calculate rising and settings. Only needed if \emph{date} is not given.
#' @param calendar Calendar for civil dates. Use \code{"G"} for Gregorian or \code{"J"} for Julian.
#' @param alt The altitude of the horizon to consider. Defaults to zero degrees.
#' @param timezone Timezone for civil timestamps. Use UTC or a fixed offset for dates before 1900. See \code{\link{time2jd}} for named zones and offsets.
#' @param loc Location, either a \emph{skyscapeR.horizon} or a vector
#' containing the latitude, longitude and elevation of location, in this order.
#' @param dec Output declination: \emph{geo} for the geocentric, or \emph{topo} for the topocentric
#' frame of reference. If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param refraction Logical. Apply atmospheric refraction. Defaults to the session setting in \code{\link{skyscapeR.vars}}.
#' @param atm Atmospheric pressure in mbar. Defaults to the session setting in \code{\link{skyscapeR.vars}}.
#' @param temp Air temperature in degrees Celsius. Defaults to the session setting in \code{\link{skyscapeR.vars}}.
#' @param verbose Logical. Display progress or informational messages. Default \code{TRUE}.
#' @import swephR
#' @return A list containing object, rise and set; rise/set are data frames
#' with date, time, azimuth and declination columns. Missing events have NA fields.
#' @export
#' @references Astrodienst. Swiss Ephemeris programmer's manual. \url{https://www.astro.com/swisseph/swephprg.htm}.
#' @seealso \code{\link[swephR]{swe_calc_ut}}, \code{\link[swephR]{swe_azalt}}, \code{\link{body.position}}, \code{\link{star.phases}}
#' @examples
#' if (requireNamespace("swephRdata", quietly = TRUE)) {
#'   # Rising and setting of the sun on june solstice 2018, from the location of London
#'   riseset(obj='sun', date='2018/06/21', loc=c(51.5, 0.11, 100))
#'   
#'   \dontrun{
#'   # Horizon searches and multi-date calculations can exceed five seconds.
#'   # Rising and setting of the moon on june solstice 2018, using a horizon profile
#'   hor <- malta_horizons[[1]]
#'   riseset(obj='moon', date='2018/06/21', loc=hor)
#'   
#'   # Rising and setting of the sun throughout February 1999, from the location of London
#'   riseset(obj='sun', date='1999/02', loc=c(51.5, 0.11, 100))
#'   
#'   # Rising and setting of the sun throughout 3001 BC, from the location of London
#'   riseset(obj='sun', date=-3000, loc=c(51.5, 0.11, 100))
#'   }
#' }
#' @section Crossing is not first visibility:
#' The result supplies geometry and time. It applies neither atmospheric extinction
#' nor a twilight detection threshold and does not test whether the star can be
#' seen. A centre exactly at the skyline can remain obscured in a plot or too faint
#' to detect. Use \code{\link{star.phases}} to select a candidate season, then
#' inspect moments after the crossing with \code{\link{sketchSky}} and
#' \code{\link{skyStep}}. The extra time needed for clearance depends on track
#' slope, latitude and terrain; a universal minute offset is inappropriate.
#'
#' A numerical root tolerance describes the solver, not observational timing
#' accuracy. Uncertain terrain height, refraction and ancient Earth-rotation
#' models can dominate it. Saved output contains formatted civil times rather
#' than a complete record of the input atmosphere and horizon. Preserve the
#' original call and convert date/time back with matching calendar and timezone.
riseset <- function(loc, date, jd, obj = 'sun', alt=0, timezone, calendar, dec, refraction, atm,
                    temp, verbose=TRUE) {
  if (missing(timezone)) { timezone <- skyscapeR.env$timezone }
  if (missing(calendar)) { calendar <- skyscapeR.env$calendar }
  if (missing(dec)) { dec <- skyscapeR.env$dec }
  if (missing(refraction)) { refraction <- skyscapeR.env$refraction }
  if (missing(atm)) { atm <- skyscapeR.env$atm }
  if (missing(temp)) { temp <- skyscapeR.env$temp }

  out <- c()
  out$object <- obj
  body <- tryCatch(checkbody(obj),error=function(e) NULL)
  stellar <- is.null(body)
  if (stellar) {
    alias <- match(obj,skyscapeR::star.names$Western)
    id <- if (is.na(alias)) obj else skyscapeR::star.names$identifier[alias]
    key <- paste0(",",id)
  }

  hor <- NULL
  if (inherits(loc,'skyscapeR.horizon')) {
    if (!missing(alt)) stop('Supply either a horizon profile or an explicit alt, not both.')
    hor <- loc
    loc <- as.numeric(hor$metadata$georef)
  }
  if (!is.logical(refraction) || length(refraction)!=1L || is.na(refraction))
    stop('refraction must be TRUE or FALSE.')
  if (!is.numeric(loc) || length(loc)!=3L || any(!is.finite(loc)))
    stop('loc must supply finite latitude, longitude and elevation.')
  if (is.null(hor) && (!is.numeric(alt) || length(alt)!=1L || !is.finite(alt)))
    stop('alt must be a finite scalar altitude.')

  if (!missing(date) && is.character(date)) {
    date <- if (grepl("/",date,fixed=TRUE)) .era_date(date) else astroYear(date)
  }
  if (!missing(jd) & !missing(date)) { stop('Both date and JD were found, please provide only one.') }

  if (!missing(jd)) {
    date <- jd2time(jd, timezone, calendar)
    date <- substr(date,1,which(strsplit(date, "")[[1]]==" ")-1)
    date <- paste(date, '00:00:00')
    jd0 <- time2jd(date, timezone, calendar)
  } else
  if (is.character(date) && grepl("^[+-]?[0-9]+/[0-9]{1,2}/[0-9]{1,2}$",date)) {
    date <- paste(date, '00:00:00')
    jd0 <- time2jd(date, timezone, calendar)
  } else
  if (is.numeric(date)) {
    dat <- paste0(date,'/01/01 00:00:00')
    jd0 <- time2jd(dat, timezone, calendar)
    jd <- seq(jd0, jd0+366, 1)
    ind <- which(substr(jd2time(jd, timezone, calendar), 1, which(strsplit(jd2time(jd, timezone, calendar), "")[[1]]=="/")[1]-1)==as.character(date))
    jd0 <- jd[ind]
  } else
  if (is.character(date) && grepl("^[+-]?[0-9]+/[0-9]{1,2}$",date)) {
    dat <- paste0(date,'/01 00:00:00')
    jd0 <- time2jd(dat, timezone, calendar)
    jd <- seq(jd0, jd0+31, 1)
    ind <- which(substr(jd2time(jd, timezone, calendar), 1, which(strsplit(jd2time(jd, timezone, calendar), "")[[1]]=="/")[2]-1)==date)
    jd0 <- jd[ind]
  } else stop("Supply a year, YYYY/MM or YYYY/MM/DD.",call.=FALSE)

  if (stellar) {
    .file_position(jd0[1],0L,2L)
    swephR::swe_set_topo(loc[2],loc[1],loc[3])
    stellar.position <- function(t,topocentric=TRUE) {
      flags <- if (topocentric) 32768L+2048L+2L else 2048L+2L
      eq <- swephR::swe_fixstar2_ut(key,t,flags)
      if (eq$return<0L) stop(eq$serr,call.=FALSE)
      if (bitwAnd(eq$return,7L)!=2L) .ephemeris_failure(t,key,eq$serr)
      pos <- swephR::swe_azalt(t,1L,loc[c(2,1,3)],atm,temp,eq$xx[1:3])$xaz
      height <- pos[if (refraction) 3L else 2L]
      c((pos[1]+180)%%360,height,eq$xx[2])
    }
    profile <- if (is.null(hor)) data.frame(az=c(0,360),alt=rep(alt,2)) else hor$data
    track <- function(t) stellar.position(t)[1:2]
  }
  rises <- data.frame(date=NA, time=NA, azimuth=NA, declination=NA, stringsAsFactors = FALSE); sets <- rises
  if (.sky_progress(verbose) && length(jd0)>1) { pb <- txtProgressBar(max = length(jd0), style=3); on.exit(close(pb),add=TRUE) }

  for (k in 1:length(jd0)) {
    if (stellar) {
      events <- .horizon_crossings(track,jd0[k],profile$az,profile$alt)
    } else if (is.null(hor)) {
      events <- vapply(c(1L,2L),function(direction) {
        result <- .file_crossing(jd0[k],body,'',0,
          direction+256L+ifelse(refraction,0L,512L),
          c(loc[2],loc[1],loc[3]),atm,temp,alt)
        if (result$return == -2L) return(NA_real_)
        if (result$return < 0L) stop(result$serr)
        result$tret[1]
      },numeric(1))
    } else {
      events <- .horizon_events(obj,jd0[k],hor,refraction,atm,temp)
    }
    rows <- lapply(events,function(event) {
      if (!is.finite(event)) return(list(date=NA_character_,time=NA_character_,
                                      azimuth=NA_real_,declination=NA_real_))
      if (stellar) {
        pos <- stellar.position(event)
        stamp <- strsplit(jd2time(event,timezone,calendar)," ",fixed=TRUE)[[1]]
        declination <- if (dec=="geo") stellar.position(event,FALSE)[3] else pos[3]
        return(list(date=stamp[1],time=stamp[2],azimuth=pos[1],declination=declination))
      }
      pos <- body.position(obj=obj,time=event,timezone=timezone,calendar=calendar,dec=dec,loc=loc,refraction,atm,temp,verbose=FALSE)
      stamp <- strsplit(jd2time(event,timezone,calendar),' ',fixed=TRUE)[[1]]
      list(date=stamp[1],time=stamp[2],azimuth=pos$horizontal$az,
           declination=pos$equatorial$Dec)
    })
    rises[k,] <- rows[[1]]
    sets[k,] <- rows[[2]]

    if (.sky_progress(verbose) && length(jd0)>1) { setTxtProgressBar(pb, k) }
  }

  rises$azimuth <- as.numeric(rises$azimuth); rises$declination <- as.numeric(rises$declination)
  sets$azimuth <- as.numeric(sets$azimuth); sets$declination <- as.numeric(sets$declination)

  out$rise = rises
  out$set = sets
  return(out)
}




#' Trace the path of a celestial body at a site
#' @description Trace the daily path of a fixed declination across the sky as seen from a site.
#' @details Latitude and declination determine where a diurnal track can meet the terrain. The
#' function sweeps equatorial right ascension through 360 degrees at the requested angular
#' step and converts each position to horizontal coordinates with the selected refraction
#' and atmosphere settings. It retains altitudes above -20 degrees and sorts the points by
#' azimuth. A horizon object supplies the location; the stored track itself includes no
#' test of terrain obstruction or event times.
#' @seealso \code{\link{plotAzimuth}}, \code{\link{plotBars}}, \code{\link{sketchSky}}, \code{\link{planisphere}}
#' @return A skyscapeR.orbit list with az and alt vectors in degrees, ordered by azimuth, plus the input declination (dec) and location (georef). The sampled track retains altitudes above -20 degrees.
#' @param dec Declination of object.
#' @param loc Location, either a \emph{skyscapeR.horizon} or a vector
#' containing the latitude and longitude of location, in this order.
#' @param res The resolution (in degrees of RA) with which to calculate the path.
#' @param refraction Logical. Apply atmospheric refraction. Defaults to the session setting in \code{\link{skyscapeR.vars}}.
#' @param atm Atmospheric pressure in mbar. Defaults to the session setting in \code{\link{skyscapeR.vars}}.
#' @param temp Air temperature in degrees Celsius. Defaults to the session setting in \code{\link{skyscapeR.vars}}.
#' @import swephR
#' @export
#' @examples
#' # Visible path of sun on june solstice on year 4001 BC from London:
#' sun.dec <- jS(-4000)
#' london.lat <- 51.5074 #N
#' london.lon <- -0.1278 #W
#' loc <- c( london.lat, london.lon, 0 )
#' path <- orbit(sun.dec, loc)
#' plot(path$az, path$alt, ylim=c(0,90), type='l', xlab='AZ', ylab='ALT', col='red', lwd=2)
orbit = function(dec, loc, refraction, atm, temp, res=0.25) {
  if (missing(refraction)) { refraction <- skyscapeR.env$refraction }
  if (missing(atm)) { atm <- skyscapeR.env$atm }
  if (missing(temp)) { temp <- skyscapeR.env$temp }

  if (class(loc)[1]=='skyscapeR.horizon') { loc <- loc$metadata$georef }

  ra <- seq(0, 360, by=res)
  aux <- array(NA, c(NROW(ra),2))

  for (i in 1:NROW(ra)) {
    tmp <- eq2hor(ra[i], dec, loc, refraction, atm, temp)

    aux[i,] <- c(tmp$az,tmp$alt)
  }

  # trim off and sort by aazimuth
  ind <- which (aux[,2] > -20)
  aux <- aux[ind,]

  ind <- sort(aux[,1], index.return=TRUE)$ix
  aux <- aux[ind,]

  # return result
  orbit <- c()
  orbit$az <- aux[,1]
  orbit$alt <- aux[,2]
  orbit$dec <- dec
  orbit$georef <- loc
  class(orbit) <- "skyscapeR.orbit"
  return(orbit)
}



#' Calculate the solar bearing for a field observation
#' @description Establish the true direction of a timed solar sight used to orient an instrument in a
#' field survey.
#' @details A solar reference requires the observer position, timestamp, timezone and the part of
#' the disc that was sighted. Times are parsed using R's civil-time routines and converted
#' to UT Julian dates. Swiss Ephemeris supplies the topocentric solar position and its
#' horizontal coordinates using 1013.25 mbar and 15 degrees Celsius. The returned bearing
#' is measured clockwise from north. Sighting the left or right limb applies a fixed minus
#' or plus 16 arcminutes in azimuth. Optional altitude is the geometric centre altitude.
#'
#' Astronomical positions require Swiss Ephemeris files covering the epoch. We
#' recommend installing \pkg{swephRdata}; see \code{\link{skyscapeR-package}} for setup.
#' @return Solar azimuth in degrees clockwise from north, one value per observation. With alt = TRUE, a data frame with columns az and alt is returned. The altitude is geometric and is expressed in degrees.
#' @param loc Location, either a \emph{skyscapeR.horizon} or a vector
#' containing the latitude, longitude and elevation of location, in this order.
#' @param time String containing the date and time in the following format:
#' "YYYY-MM-DD HH:MM:SS"
#' @param timezone Timezone for civil timestamps. Use UTC or a fixed offset for dates before 1900. See \code{\link{time2jd}} for named zones and offsets.
#' @param limb Measured limb of the sun. Options are \emph{left}, \emph{right}.
#' If missing the center of the sun will be output.
#' @param alt Logical that triggers output of altitude of the sun at exact time.
#' Default is FALSE.
#' @import swephR
#' @export
#' @seealso \code{\link{reduct.theodolite}}
#' @examples
#' if (requireNamespace("swephRdata", quietly = TRUE)) {
#'   sunAz(c(52,-3,100), '2017-10-04 12:32:14', 'Europe/London')
#' }
sunAz = function(loc, time, timezone, limb, alt=FALSE) {
  if (missing(timezone)) { timezone <- skyscapeR.env$timezone }

  if (class(loc)[1]=='skyscapeR.horizon') { loc <- c(loc$metadata$georef, loc$metadata$elevation) }
  if (is.null(dim(loc))) { dim(loc) <- c(1, NROW(loc)) }

  az <- c(); at <- c()
  for (i in 1:NROW(loc)) {
    pb.date <- as.POSIXct(time[i], timezone[i])
    UT <- format(pb.date, tz="UTC",usetz=TRUE)
    UT <- as.POSIXlt(UT, 'UTC')

    swephR::swe_set_topo(loc[i,2], loc[i,1], loc[i,3])
    jd <- swephR::swe_julday(UT$year+1900, UT$mon+1, UT$mday, UT$hour+UT$min/60+UT$sec/3600,1)
    ss <- .file_position(jd, 0, 32*1024+2048)
    xin <- ss$xx[1:2]
    aux <- swephR::swe_azalt(jd, 1, c(loc[i,2],loc[i,1],loc[i,3]), 1013.25, 15, xin)
    az[i] <- aux$xaz[1]+180
    if (az[i] > 360) { az[i] <- az[i]-360 }

    if (!missing(limb)) {
      if (limb=="left") { az[i] <- az[i] - 32/60/2 }
      if (limb=="right") { az[i] <- az[i] + 32/60/2 }
    }
    if (alt) { at[i] <- aux$xaz[2] }
  }
  if (alt) {
    df <- data.frame(az = az, alt = at)
    return(df)
  } else {
    return(az)
  }
}


#' Find dates when the Sun reaches a specified declination
#' @description Estimate the dates on which the Sun could attain a declination indicated by a structure
#' in a selected year.
#' @details The Sun can pass the same declination on either side of a solstice. The search
#' covers the selected civil year from midnight UT, including its first instant and
#' excluding the start of the next year. Daily geocentric apparent solar declinations
#' bracket the annual turning points. The function refines each extremum against the
#' ephemeris and solves separately on either side, then reports the civil day of each
#' crossing. A tangency gives one date. Two crossings can share a calendar day.
#' The calculation uses no terrain or solar limb model.
#'
#' Astronomical positions require Swiss Ephemeris files covering the epoch. We
#' recommend installing \pkg{swephRdata}; see \code{\link{skyscapeR-package}} for setup.
#' @seealso \code{\link{spatial.equinox}}, \code{\link{zenith}}, \code{\link{antizenith}}, \code{\link{findTargets}}, \code{\link{EFM}}, \code{\link{planet.extremes}}
#' @return A character matrix with one column per reachable input, in input order. The rows
#' are Dec, Date1 and Date2, with Date3 added if three passages fall in one civil
#' year. Dates are chronological day/month labels. Unused date cells are missing.
#' The function stops if all requested declinations are unreachable.
#' @param dec One or more declinations in degrees.
#' @param year Year for which to do calculations. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' @param calendar Calendar for sampling and output. G for Gregorian or J for Julian. Defaults to the package calendar, initially Gregorian.
#' @param verbose Logical. Display progress or informational messages. Default \code{TRUE}.
#' @import swephR
#' @export
#' @examples
#' if (requireNamespace("swephRdata", quietly = TRUE)) {
#'   solar.date(-23, 2018)
#'   solar.date(-12, 1200, calendar='G')
#'   solar.date(-12, 1200, calendar='J')
#'   solar.date(14, -2000)
#' }
solar.date <- function(dec, year, calendar, verbose=TRUE){
  if (is.character(year)) year <- astroYear(year)
  checkYear(year)
  if (missing(calendar)) calendar <- skyscapeR.env$calendar
  cal <- .calendar_flag(calendar)
  start <- swephR::swe_julday(year, 1, 1, 0, cal)
  end <- swephR::swe_julday(year+1, 1, 1, 0, cal)
  position <- function(t) .solar_coordinate(start+t, 2048L, 2L)
  days <- 0:(end-start)
  daily <- position(days)
  # Include turning points at the edges of the calendar as well as inside it.
  padded <- c(position(-1), daily, position(end-start+1))
  extrema <- which(diff(sign(diff(padded))) != 0) - 1L
  turns <- vapply(extrema, function(d) {
    optimize(position, c(d-1,d+1), maximum=padded[d+2] > padded[d+1],
             tol=1e-8)[[1L]]
  }, numeric(1))
  knots <- sort(unique(c(days, turns[turns>0 & turns<end-start])))
  values <- position(knots)
  roots <- lapply(dec, function(target) {
    delta <- values-target
    hits <- knots[abs(delta)<=1e-10]
    brackets <- which(head(delta,-1L)*tail(delta,-1L)<0 &
                      abs(head(delta,-1L))>1e-10 & abs(tail(delta,-1L))>1e-10)
    for (i in brackets) {
      hits <- c(hits, uniroot(function(t) position(t)-target, knots[i+0:1],
                              f.lower=delta[i], f.upper=delta[i+1], tol=1e-8)$root)
    }
    sort(hits[hits<end-start])
  })
  keep <- lengths(roots)>0
  if (!any(keep)) stop("Value(s) are outside of solar range.")
  if (verbose && any(!keep)) .sky_message(sum(!keep)," value(s) were outside of solar range and have been excluded.")
  roots <- roots[keep]
  n <- max(2L,lengths(roots))
  out <- matrix(NA_character_, n+1L,length(roots),
                dimnames=list(c("Dec",paste0("Date",seq_len(n))),NULL))
  out[1,] <- as.character(dec[keep])
  for (i in seq_along(roots)) {
    out[1L+seq_along(roots[[i]]),i] <- .solar_day(start+roots[[i]], cal)
  }
  out
}

.solar_coordinate <- function(jd, flags, column) {
  result <- .file_position(jd, 0L, flags)
  if (any(result$return < 0L)) stop(paste(unique(result$serr),collapse=" "))
  matrix(result$xx, nrow=length(jd))[,column]
}

.solar_day <- function(jd, cal) {
  vapply(jd, function(t) {
    date <- swephR::swe_revjul(t, cal)
    paste(date$day,month.abb[date$month])
  }, character(1), USE.NAMES=FALSE)
}

.season_dates <- function(year, calendar) {
  cal <- .calendar_flag(calendar)
  start <- swephR::swe_julday(year,1,1,0,cal)
  end <- swephR::swe_julday(year+1,1,1,0,cal)
  days <- 0:(end-start)
  longitude <- .solar_coordinate(start+days,0L,1L)
  vapply(c(0,90,180,270), function(target) {
    residual <- function(lon) (lon-target+180) %% 360 - 180
    d <- residual(longitude)
    i <- which(head(d,-1L)<=0 & tail(d,-1L)>0)
    if (!length(i)) return(NA_character_)
    times <- vapply(i,function(k) uniroot(function(t)
      residual(.solar_coordinate(start+t,0L,1L)), days[k+0:1],
      f.lower=d[k],f.upper=d[k+1],tol=1e-8)$root,numeric(1))
    paste(.solar_day(start+times,cal),collapse=" / ")
  }, character(1))
}

#' Estimate the parallax correction to declination
#' @description Estimate how the observer's position on Earth shifts a solar or lunar target declination
#' relative to its geocentric value.
#' @details The limiting-declination functions subtract this approximate correction from their
#' geocentric target. The correction in degrees is the supplied parallax multiplied by the
#' cosine of altitude and sine of latitude, with the angular arguments evaluated in
#' degrees. The sign therefore reverses between hemispheres. This scalar approximation is
#' distinct from computing an instantaneous topocentric ephemeris.
#'
#' Ruggles gives an archaeological example of applying mean lunar parallax in a
#' declination calculator.
#' @references Ruggles, C. L. N. Basic declination calculator. \url{https://web.cliveruggles.com/tools/declination-calculator}.
#' @seealso \code{\link{skyObject}}, \code{\link{jS}}, \code{\link{dS}}, \code{\link{eq}}, \code{\link{nmnLX}}, \code{\link{smnLX}}, \code{\link{nMjLX}}, \code{\link{sMjLX}}, \code{\link{obliquity}}
#' @return Numeric declination correction in degrees. The sign depends on latitude.
#' @param parallax Average parallax to correct (e.g. 0.00224 for the Sun, or 0.952 for the Moon)
#' @param loc This can be either the latitude of the location, or a \emph{skyscapeR.horizon} object.
#' @param altitude Altitude of the celestial object. Defaults to 0 degrees.
#' @import swephR
#' @export
#' @examples
#' # Parallax correction for the moon, as seen from latitude 50ºN and at 0º altitude
#' parallax.corr(0.952, 50, 0)
parallax.corr <- function(parallax, loc, altitude=0) {
  latitude <- extractLatitude(loc)
  return(parallax*cos(altitude/180*pi)*sin(latitude/180*pi)) # V Reijs/SE formula
  # return(parallax * cos(altitude/180*pi) * (1-sin(latitude/180*pi)^2/298.3)) # Nautical formula
}
