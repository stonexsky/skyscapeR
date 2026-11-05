################################################
skyscapeR.env <- new.env(parent = emptyenv())

# Timezone and Calendar
skyscapeR.env$timezone <- 'Europe/London'
skyscapeR.env$calendar <- 'Gregorian'

# Atmospheric Refraction
skyscapeR.env$refraction <- TRUE
skyscapeR.env$atm <- 1013.25
skyscapeR.env$temp <- 15
skyscapeR.env$k <- 0.2

# Frame of Reference
skyscapeR.env$dec <- 'topo'

# Current Year (cannot be changed)
skyscapeR.env$cur.year <- as.numeric(format(Sys.Date(), "%Y"))

################################################

#' Inspect or set session defaults for skyscapeR
#' @description Set the observing and time conventions that should be shared across an analysis, or
#' retrieve them for a reproducible methods record.
#' @details The same observation can have different labels or apparent coordinates under different
#' calendar, timezone, reference-frame or atmospheric choices. With no arguments the
#' function returns a named snapshot. Supplied arguments update the package environment for
#' subsequent calls in the current R session. Functions with their own explicit arguments
#' or fixed UT sampling conventions may override these defaults. Record the settings with
#' the survey data.
#' @references Schaefer, B. E. (1986). Atmospheric Extinction Effects on Stellar Alignments. Journal for the History of Astronomy 17, S32-S42. \doi{10.1177/002182868601701003}.
#' @references Schaefer, B. E. (1993). Astronomy and the Limits of Vision. Vistas in Astronomy 36, 311-361. \doi{10.1016/0083-6656(93)90113-X}.
#' @references University of California, Davis. \emph{Physics 157 Observing Guide}, Spring 2022, discussion of extinction coefficients.
#' @seealso \code{\link{skyscapeR-package}}, \code{\link{time2jd}}, \code{\link{body.position}}
#' @examples settings <- skyscapeR.vars()
#' settings$timezone
#' @return A named list of the current session settings, including the current year. Supplying an argument updates that setting for subsequent package calls in this R session.
#' @param timezone Timezone name (for example "Europe/London") or fixed offset in hours east of UTC,
#' or a string such as "UTC+02:00". Initially "Europe/London". Use UTC or a fixed
#' offset for dates before 1900. See \code{\link{time2jd}} for details.
#' @param calendar Calendar for civil dates. Use \code{"G"} for Gregorian or \code{"J"} for Julian.
#' @param refraction  Whether atmospheric refraction is to be taken into account. Default is TRUE.
#' @param atm Air pressure in mbar for the refraction calculation. Default is 1013.25 mbar.
#' @param temp Air temperature in degrees Celsius for the refraction calculation. Default is 15 degrees Celsius.
#' @param dec Output declination: \emph{geo} for the geocentric, or \emph{topo} for the topocentric
#' frame of reference. Defaults to topocentric.
#' @param k Visual extinction in magnitudes per airmass, shared by natural sky sketches and stellar phase calculations. The initial value of 0.2 is a clear-sky baseline, consistent with the approximate V-band coefficient for Kitt Peak in the UC Davis Physics 157 Observing Guide. Extinction varies with site and atmospheric conditions (Schaefer 1986, 1993). Use zero to disable extinction.
#' @export
skyscapeR.vars = function(timezone, calendar, refraction, atm, temp, dec, k) {
  if (!missing(k)) {
    if (!is.numeric(k) || length(k)!=1L || !is.finite(k) || k<0) stop("Supply a non-negative extinction coefficient k.")
    skyscapeR.env$k <- k
  }
  if (!missing(timezone)) { skyscapeR.env$timezone <- timezone }
  if (!missing(calendar)) { skyscapeR.env$calendar <- calendar }
  if (!missing(refraction)) { skyscapeR.env$refraction <- refraction }
  if (!missing(atm)) { skyscapeR.env$atm <- atm }
  if (!missing(temp)) { skyscapeR.env$temp <- temp }
  if (!missing(dec)) { skyscapeR.env$dec <- dec }

  aux <- ls(skyscapeR.env)
  return(mget(aux, skyscapeR.env))
}
globalVariables('star.names')



