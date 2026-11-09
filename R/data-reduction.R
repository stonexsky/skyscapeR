#' Reduce theodolite readings using a solar sight
#' @description Orient a theodolite survey to true north using a timed solar sight, then recover the
#' true bearing of a structural axis.
#' @details The observed difference between the instrument readings to the structure and to the Sun
#' is added to the calculated solar azimuth at the observation time. The routine handles a
#' circle-boundary difference, accounts for the selected solar limb through sunAz and
#' rounds the reference bearing according to the input readings. True bearings lie in
#' [0, 360) degrees. A three-coordinate vector supplies one shared site; for multiple
#' sites, supply a matrix with one latitude, longitude and elevation row per reading. An explicit altitude or a
#' supplied horizon allows the corrected bearing to be converted to declination. Correct
#' time, civil timezone, limb and site are all essential to the reduction.
#'
#' Astronomical positions require Swiss Ephemeris files covering the epoch. We
#' recommend installing \pkg{swephRdata}; see \code{\link{skyscapeR-package}} for setup.
#' @return A data frame with one row per observation. True.Azimuth is the corrected azimuth in degrees. The table also records the observed azimuth, solar azimuth, time and site coordinates. Where available it includes Altitude and Declination in degrees and the supplied identifiers.
#' @param loc Location, either a \emph{skyscapeR.horizon} or a vector
#' containing the latitude, longitude and elevation of location, in this order.
#' @param az Array of azimuths. Use \code{\link{ten}} to convert to
#' decimal degrees if necessary.
#' @param date Date of measurements as a string in the format: 'YYYY/MM/DD'
#' @param time Time of sun-sight measurement in the format: 'HH:MM:SS'
#' @param tz Timezone of input either as a known acronym (e.g. "GMT", "CET") or
#' a string with continent followed by country capital (e.g. "Europe/London").
#' @param az.sun Measured azimuth of the sun. Defaults to zero.
#' @param limb Measured limb of the sun. Options are \emph{left}, \emph{right}.
#' If missing the center of the sun will be used for calculations.
#' @param alt Altitude, necessary for automatic declination calculation.
#' If missing and \emph{loc} is a \emph{skyscapeR.horizon} object then the altitude
#' will be automatically read from the horizon profile.
#' @param name Names or labels to identify each measurement.
#' @param ID IDs or codes to identify each measurement.
#' @param HWT.ID HeyWhatsThat IDs relating to a previously generated
#' horizon profile for measurement.
#' @section Validation:
#' W02 recovered both independent circle-difference control bearings and retained each site row in a multi-site solar reduction. See the \href{https://stonexsky.github.io/skyscapeR/manual/#w02-field-to-sky}{W02 check in the technical manual}.
#' @export
#' @seealso \code{\link{sunAz}}, \code{\link{ten}}
#' @references Ruggles, C.L.N. (1999). \emph{Astronomy in Prehistoric Britain and Ireland}.
#' Yale University Press.
#' @examples
#' if (requireNamespace("swephRdata", quietly = TRUE)) {
#'   lat <- ten(35,50,37.8)
#'   lon <- ten(14,34,6.4)
#'   elev <- 100
#'   az <- c( ten(298,24,10), ten(302,20,40))
#'   az.sun <- ten(327,29,50)
#'   date <- "2016/02/20"
#'   time <- "11:07:17"
#'   
#'   data <- reduct.theodolite(c(lat,lon,elev), az, date , time, tz= "Europe/Malta", az.sun)
#'   
#'   # Declination will be automatically calculated if the altitude is also given:
#'   data <- reduct.theodolite(c(lat,lon,elev), az, date, time,
#'   tz="Europe/Malta", az.sun=az.sun, alt=c(2,5))
#'   
#'   # Alternatively, the altitude can be automatically retrieved from a horizon profile:
#'   hor <- malta_horizons[[1]]
#'   data <- reduct.theodolite(hor, az, date, time, tz= "Europe/Malta", az.sun)
#' }
reduct.theodolite = function(loc, az, date, time, tz, az.sun = 0, limb, alt, name, ID, HWT.ID) {
  if (class(loc)[1]=='skyscapeR.horizon') { hor <- loc; loc <- loc$metadata$georef } else { hor <- NULL }

  n <- length(az)
  if (is.null(dim(loc))) loc <- matrix(loc,nrow=1)
  if (!n || ncol(loc)!=3L || !nrow(loc) %in% c(1L,n))
    stop("Supply one latitude, longitude and elevation per site or one site shared by all readings.")
  if (nrow(loc)==1L) loc <- loc[rep(1L,n),,drop=FALSE]
  for (value in list(date,time,tz,az.sun)) {
    if (!length(value) %in% c(1L,n)) stop("Pair the dates, times, timezones and solar readings with the target readings.")
  }
  date <- rep_len(date,n); time <- rep_len(time,n); tz <- rep_len(tz,n)
  if (!missing(alt) && !length(alt) %in% c(1L,n)) stop("Supply one altitude or one per reading.")
  date <- as.Date(date, "%Y/%m/%d")
  time <- paste(date, time)

  diff <- az - az.sun

  prec <- max(nchar(sub('.*\\.', '', as.character(az))))
  az.sun.corr <- round(sunAz(loc, time, tz, limb), prec)
  az.corr <- (az.sun.corr + diff) %% 360

  df <- c()
  if (!missing(ID)) { df$ID <- ID }
  if (!missing(name)) { df$Name <- name }
  df$Latitude=loc[,1]
  df$Longitude=loc[,2]
  if (!missing(HWT.ID)) { df$HWT.ID <- HWT.ID }
  df$Uncorrected.Azimuth=az
  df$Date.Time=time
  df$Sun.Az=az.sun.corr
  df$True.Azimuth=az.corr

  if (!missing(alt)) {
    .sky_message('Altitude values found. Calculating declination...')
    dec <- az2dec(az.corr, loc, alt)
    df$Altitude = alt
    df$Declination <- dec
  } else if (class(hor)[1]=='skyscapeR.horizon') {
    .sky_message('Horizon profile found. Obtaining altitude values and calculating declination...')
    dec <- az2dec(az.corr, hor)
    df$Altitude <- hor2alt(hor, az.corr)
    df$Declination <- dec
  } else { .sky_message('No altitude values or horizon profile found. Declination values were not calculated.') }

  df <- as.data.frame(df)
  return(df)
}


#' Reduce magnetic compass bearings to true bearings
#' @description Recover true bearings from magnetic compass observations and, where horizon altitude is
#' known, estimate the indicated declinations.
#' @details Magnetic north differs from true north at the time and place of survey. Add the supplied
#' magnetic declination to each magnetic bearing, or calculate it from the survey
#' coordinates and date with mag.dec. Corrected bearings are rounded using the input
#' bearing representation. With an explicit altitude or a measured horizon, the corrected
#' bearing also enters the conversion to declination.
#'
#' Supply either a measurement date or a known magnetic declination.
#' If both are given, the supplied magnetic declination is used and the date is
#' retained in the output. Explicit altitude values take precedence over a
#' horizon profile. Otherwise, each measurement's corrected true azimuth is
#' used to retrieve its horizon altitude; the same altitude is reported and
#' used to calculate celestial declination. Supply one site and date for the whole
#' survey or one location row and/or date per bearing, in measurement order.
#' A supplied correction can likewise be shared or paired with each bearing.
#' @param loc Location, either a \emph{skyscapeR.horizon} or a vector
#' containing the latitude, longitude and elevation of location, in this order.
#' @param mag.az Array of magnetic azimuth measurements.
#' @param date Date of measurements as a string in the format: 'YYYY/MM/DD'.
#' Only necessary if \emph{magdec} is not given. When omitted, the output
#' Date column contains missing values (NA).
#' @param magdec Magnetic declination, if known.
#' @param alt Altitude, necessary for automatic declination calculation.
#' If missing and \emph{loc} is a \emph{skyscapeR.horizon} object then the altitude
#' will be automatically read from the horizon profile.
#' @param name Names or labels to identify each measurement.
#' @param ID IDs or codes to identify each measurement.
#' @param HWT.ID HeyWhatsThat IDs relating to a previously generated
#' horizon profile for measurement.
#' @return A data frame with one row per measurement, including location,
#' magnetic azimuth, date, magnetic declination and true azimuth. Altitude and
#' celestial declination are included when altitudes are supplied or a horizon
#' is available. Optional identification columns are retained.
#' @section Validation:
#' W02 checked north-crossing magnetic readings against the supplied corrections before passing true bearings into a skyline analysis. See the \href{https://stonexsky.github.io/skyscapeR/manual/#w02-field-to-sky}{W02 check in the technical manual}.
#' @export
#' @seealso \code{\link{mag.dec}}, \code{\link{az2dec}}, \code{\link{hor2alt}}
#' @examples
#' loc <- c(35,-7, 100)
#' mag.az <- c(89.5, 105, 109.5)
#' data <- reduct.compass(loc, mag.az, "2016/04/02")
#'
#' # Declination will be automatically calculated if the altitude is also given:
#' data <- reduct.compass(loc, mag.az, "2016/04/02", alt=c(1,2,0))
#'
#' # Alternatively, the altitude can be automatically retrieved from a horizon profile:
#' hor <- malta_horizons[[1]]
#' data <- reduct.compass(hor, mag.az, "2016/04/02")
#'
#' # Known magnetic declination: no measurement date is required.
#' reduct.compass(c(35,-7,100), c(89.5,105,109.5), magdec=-2, alt=c(1,2,0))
#'
#' # A varying horizon supplies an altitude for each corrected azimuth.
#' hor <- createHor(c(0,90,180,270,360), c(0,2,6,1,0), loc=c(35,-7,100))
#' reduct.compass(hor, c(89.5,105,109.5), magdec=-2)
reduct.compass = function(loc, mag.az, date, magdec, alt, name, ID, HWT.ID) {
  if (missing(date) && missing(magdec)) {
    stop("Supply either date (to calculate magnetic declination) or magdec.")
  }
  if (missing(date)) { date <- rep(NA_character_, length(mag.az)) }

  if (class(loc)[1]=='skyscapeR.horizon') { hor <- loc; loc <- loc$metadata$georef } else { hor <- NULL }

  if (is.null(dim(loc))) { loc <- matrix(loc, nrow=1) }
  n <- length(mag.az)
  if (!nrow(loc) %in% c(1L,n) || !ncol(loc) %in% 2:3)
    stop("Supply one location or one location per bearing.")
  if (nrow(loc)==1L) loc <- loc[rep(1L,n),,drop=FALSE]
  if (!length(date) %in% c(1L,n)) stop("Supply one date or one date per bearing.")
  date <- rep(date,length.out=n)
  if (!missing(magdec) && !length(magdec) %in% c(1L,n))
    stop("Supply one correction or one correction per bearing.")
  if (!missing(alt)) { if (NROW(alt) < NROW(mag.az)) { alt <- matrix(alt,NROW(mag.az),1, byrow=T) } }

  if (missing(magdec) & !missing(date)) {
    magdec <- mag.dec(loc, date)
  }
  prec <- max(nchar(sub('.*\\.', '', as.character(mag.az))))
  true.az <- round(mag.az + magdec, prec)
  df <- c()
  if (!missing(ID)) { df$ID <- ID }
  if (!missing(name)) { df$Name <- name }
  df$Latitude=loc[,1]
  df$Longitude=loc[,2]

  df$Magnetic.Azimuth=mag.az
  df$Date=date
  df$Mag.Dec=magdec
  df$True.Azimuth=true.az
  if (!missing(HWT.ID)) { df$HWT.ID <- HWT.ID }
  if (!missing(alt)) {
    .sky_message('Altitude values found. Calculating declination...')
    dec <- az2dec(true.az, loc, alt)
    df$Altitude = alt
    df$Declination <- dec
  } else if (class(hor)[1]=='skyscapeR.horizon') {
    .sky_message('Horizon profile found. Obtaining altitude values and calculating declination...')
    alt <- hor2alt(hor, true.az %% 360)
    dec <- az2dec(true.az, hor, alt=alt)
    df$Altitude <- alt
    df$Declination <- dec
  } else { .sky_message('No altitude values or horizon profile found. Declination values were not calculated.') }

  df <- as.data.frame(df)
  return(df)
}
