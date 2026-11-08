#' @noRd
events <- function(name) {
  return(switch(name,
                'December Solstice' = 'dS',
                'June Solstice' = 'jS',
                'Equinox' = 'eq',
                'Major Lunar Extreme (southern)' = 'sMjLX',
                'Minor Lunar Extreme (southern)' = 'smnLX',
                'Major Lunar Extreme (northern)' = 'nMjLX',
                'Minor Lunar Extreme (northern)' = 'nmnLX',
                name))
}


#' Prepare celestial targets for orientation plots
#' @description Define celestial targets and their declinations so the same hypotheses can be overlaid
#' on several orientation plots.
#' @details An overlay identifies which targets to compare before inspecting a distribution. Names
#' may denote stars, solar or lunar limits, or fixed custom declinations. Group names
#' expand into the individual limits. At one epoch, the function evaluates each target; for
#' two endpoints it uses annual January stellar declinations. Solar and lunar limits
#' retain century-spaced samples and spline interpolation. Styles and epoch travel with the object. Treat range endpoints as
#' model estimates and verify narrow or rapidly changing windows separately.
#' A single object can contain several related targets. Pass a list of objects
#' to plot(hor), plotBars or plotAzimuth to combine groups with independent
#' epochs, colours, line styles and labels. Groups are drawn in list order.
#' Horizon plots use the plotted horizon's location for every track. Compass
#' plots use each group's stored location and a level horizon, drawing both
#' declination boundaries for a year-range group. Declination plots need no
#' location. No common epoch range is imposed and no astronomical positions
#' are recomputed merely by combining objects. The output class remains
#' skyscapeR.object, so previously saved objects remain usable.
#' sky.objects is a silent compatibility alias for skyObject.
#' @seealso \code{\link{plot.skyscapeR.horizon}}, \code{\link{plotBars}},
#' \code{\link{plotAzimuth}}, \code{\link{jS}}, \code{\link{dS}}, \code{\link{eq}}, \code{\link{nmnLX}}, \code{\link{smnLX}}, \code{\link{nMjLX}}, \code{\link{sMjLX}}, \code{\link{obliquity}}, \code{\link{parallax.corr}}
#' @return A skyscapeR.object list. Columns of decs represent targets, in degrees. A single epoch produces one row. For a two-year range the rows are labelled min and max. The other fields record target count, epoch, location and plotting styles.
#' @param names Target names, star names or custom declinations in degrees. Use "solar extremes" for the solstices or "lunar extremes" for the four lunar extremes.
#' @param epoch Earlier argument name for \code{year}, retained for existing scripts.
#' @param year An astronomical year or a pair of year endpoints. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' @param loc This can be either a vector with the latitude and longitude of the
#' location, or a \emph{skyscapeR.horizon} object.
#' @param col The color for plotting, and differentiating these objects.
#' Defaults to ochre for all objects.
#' @param lty Line type (see \code{\link{par}}) used for differentiation.
#' Used for single-year tracks and year-range boundaries.
#' @param lwd Line width (see \code{\link{par}}) used for differentiation.
#' Used for single-year tracks and year-range boundaries.
#' @export
#' @examples
#' # Create an object with solar targets for epoch range 4001-2001 BC:
#' tt <- skyObject('solar extremes', c("4001 BC", "2001 BC"))
#'
#' if (requireNamespace("swephRdata", quietly=TRUE)) {
#' # Create an object with a few stars for same epoch:
#' tt <- skyObject(c('Sirius', 'Betelgeuse', 'Antares'), c("4001 BC", "2001 BC"),
#' col=c('white', 'red', 'orange'))
#' }
#'
#' # Create an object with solstices and a custom declination value:
#' tt <- skyObject(c('December Solstice','June Solstice', -13), c("4001 BC", "2001 BC"))
#' # Mix a single-epoch group with a range, preserving each group's styles:
#' groups <- list(skyObject(12, "2001 BC", col="blue"),
#'                skyObject("solar extremes", c("4001 BC", "2001 BC"), col="orange"))
#' plotBars(c(-24,12,24), unc=1, obj=groups)
#' @section A reference track is not a dated orbit:
#' Each target represents a declination or declination range. A plotted track
#' holds that declination fixed while hour angle varies. It does not follow the
#' changing declination of a moving Solar System body through a particular day.
#' In particular, \code{skyObject("Equinox", year)} is an equinox reference;
#' the Sun at sunrise or sunset on the equinox date need not lie at the same
#' declination as at the equinox instant. A small offset from a solar disc in
#' \code{\link{sketchSky}} can therefore be expected.
#'
#' To compare with a particular instant, obtain its declination with
#' \code{\link{body.position}} using the intended coordinate frame, and pass
#' that numeric value to \code{skyObject}. Even that fixed-declination track is
#' an approximation to a moving body's full daily orbit. Refraction and the
#' plotted site's coordinates are applied when drawing the horizon track.
#' @section Provenance:
#' Structured results carry a \code{provenance} attribute. Read it with
#' \code{attr(result, "provenance")}; its \code{packages} entry records the loaded
#' skyscapeR and swephR versions, with NA when a namespace version is unavailable.
#' This attribute adds no list entries or table columns. Ordinary subsetting or
#' transformations may discard attributes; save the original result for reference.
#' It does not identify ephemeris files, catalogue revisions or a development
#' commit. Retain those and the analysis script separately for exact reproduction.
skyObject = function(names, year, loc=FALSE, col = '#A9542E', lty = 1, lwd = 1, epoch) {
  if (!missing(epoch)) {
    if (!missing(year)) stop("Supply year only once.",call.=FALSE)
    year <- epoch
  }
  if (is.character(year)) year <- astroYear(year)
  if (length(names) > 1 & length(col)==1) { col <- rep(col, length(names)) }
  if (length(names) > 1 & length(lty)==1) { lty <- rep(lty, length(names)) }
  if (length(names) > 1 & length(lwd)==1) { lwd <- rep(lwd, length(names)) }

  i <- which(names=='solar extremes')
  if (length(i)>0) {
    names <- names[-i]
    names <- c('December Solstice','June Solstice',names)
    aux <- col[i]; col <- col[-i]; col <- c(rep(aux,2),col)
    aux <- lty[i]; lty <- lty[-i]; lty <- c(rep(aux,2),lty)
    aux <- lwd[i]; lwd <- lwd[-i]; lwd <- c(rep(aux,2),lwd)
  }

  i <- which(names=='lunar extremes')
  if (length(i)>0) {
    names <- names[-i]
    names <- c('Major Lunar Extreme (southern)', 'Minor Lunar Extreme (southern)', 'Minor Lunar Extreme (northern)', 'Major Lunar Extreme (northern)',names)
    aux <- col[i]; col <- col[-i]; col <- c(rep(aux,4),col)
    aux <- lty[i]; lty <- lty[-i]; lty <- c(rep(aux,4),lty)
    aux <- lwd[i]; lwd <- lwd[-i]; lwd <- c(rep(aux,4),lwd)
  }

  N <- NROW(names)

  if (NROW(col)==1) { col <- rep(col,N) }
  if (NROW(lty)==1) { lty <- rep(lty,N) }
  if (NROW(lwd)==1) { lwd <- rep(lwd,N) }

  tt <- c()
  tt.col <- c()
  tt.lty <- c()
  tt.lwd <- c()

  for (i in 1:N) {
    # stars
    named.limit <- events(as.character(names[i])) %in% c("dS","jS","eq","sMjLX","smnLX","nMjLX","nmnLX")
    numeric.target <- !is.na(suppressWarnings(as.numeric(names[i])))
    if (!numeric.target && !named.limit && !inherits(.try_target(star(names[i],year[1])), "try-error")) {
      aux <- array(NA, c(NROW(year),1))
      for (j in 1:NROW(year)) {
        aux[j,] <- star(names[i], year[j])$coord$Dec
      }
      colnames(aux) <- names[i]
      tt <- cbind(tt, aux)
      tt.col <- c(tt.col, col[i])
      tt.lty <- c(tt.lty, lty[i])
      tt.lwd <- c(tt.lwd, lwd[i])
      next
    } else if (identical(events(as.character(names[i])), "eq")) {
      aux <- matrix(eq(loc = loc, verbose = FALSE), nrow = length(year),
                    ncol = 1L, dimnames = list(NULL, "eq"))
      tt <- cbind(tt, aux)
      tt.col <- c(tt.col, col[i])
      tt.lty <- c(tt.lty, lty[i])
      tt.lwd <- c(tt.lwd, lwd[i])
    } else
      # custom dec
      if (!is.na(suppressWarnings(as.numeric(names[i])))) {
        aux <- array(NA, c(NROW(year),1))
        aux[,1] <- rep(as.numeric(names[i]),NROW(year))
        colnames(aux) <- 'Custom Dec'
        tt <- cbind(tt, aux)
        tt.col <- c(tt.col, col[i])
        tt.lty <- c(tt.lty, lty[i])
        tt.lwd <- c(tt.lwd, lwd[i])
        next
      } else {

        # try calling the functions
        aux <- array(NA, c(NROW(year),1))
        for (j in 1:NROW(year)) {
          aux[j,] <- do.call(events(names[i]), list(year=year[j], loc=loc))
        }
        colnames(aux) <- events(names[i])
        tt <- cbind(tt, aux)
        tt.col <- c(tt.col, col[i])
        tt.lty <- c(tt.lty, lty[i])
        tt.lwd <- c(tt.lwd, lwd[i])
      }

  }
  rownames(tt) <- year

  # check min and max decs
  if (NROW(year)==2) {
    ind <- which(!colnames(tt) %in% c('Custom Dec','eq'))
    for (i in ind) { tt[,i] <- minmaxdec(events(colnames(tt)[i]), from=min(year), to=max(year)) }
    rownames(tt) <- c('min','max')
  }

  # return result
  object <- c()
  object$n <- NCOL(tt)
  object$decs <- tt
  object$epoch <- year
  object$loc <- loc
  object$col <- tt.col
  object$lty <- tt.lty
  object$lwd <- tt.lwd
  class(object) <- "skyscapeR.object"
  return(.with_astronomy_provenance(object))
}


#' @rdname skyObject
#' @export
sky.objects <- skyObject


#' Estimate solar declination at the December solstice
#' @description Estimate the limiting solar declination at the December solstice for comparison with a
#' structural sightline.
#' @details The December solstice is represented by the negative obliquity of the ecliptic for the
#' requested epoch. The result is geocentric when no site is supplied. With a site, the
#' function subtracts the approximate parallax correction using latitude and target
#' altitude. This obliquity-based target does not determine a solstice time or a sunrise or
#' sunset bearing.
#' @return Numeric declination in degrees. With loc = FALSE the result is geocentric. Supplying a location applies the parallax correction for the chosen altitude.
#' @param year Astronomical year. The default is the current year in the package settings. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' @param loc Latitude in degrees north or a \code{skyscapeR.horizon} object. Omit this argument or use \code{FALSE} for geocentric declination.
#' @param parallax Mean horizontal parallax in degrees. The default is 0.00224 degrees for the Sun.
#' @param altitude Altitude of the sun. Defaults to 0 degrees.
#' @param verbose Logical to control output of warnings and messages.
#' Defaults to TRUE.
#' @export
#' @seealso \code{\link{obliquity}}, \code{\link{jS}}, \code{\link{eq}},
#' \code{\link{zenith}}, \code{\link{antizenith}}, \code{\link{spatial.equinox}},
#' \code{\link{parallax.corr}}
#' @examples
#' # December Solstice geocentric declination for year 4001 BC:
#' dS(-4000)
#'
#' # Topocentric declination for same year and latitude of 50 degrees N:
#' dS(-4000, loc=50)
dS = function(year = skyscapeR.env$cur.year, loc = FALSE, parallax = 0.00224, altitude = 0, verbose=TRUE) {
  if (is.character(year)) year <- astroYear(year)
  aux <- -obliquity(year)

  if (is(loc,'logical')) {
    if (verbose) { .sky_message('Returning geocentric declination because no site was supplied.') }
  } else {
    aux <- aux - parallax.corr(parallax, loc, altitude)
  }
  return(aux)
}

#' Estimate solar declination at the June solstice
#' @description Estimate the limiting solar declination at the June solstice for comparison with a
#' structural sightline.
#' @details The June solstice is represented by the positive obliquity of the ecliptic for the
#' requested epoch. The result is geocentric when no site is supplied. With a site, the
#' function subtracts the approximate parallax correction using latitude and target
#' altitude. This obliquity-based target does not determine a solstice time or a sunrise or
#' sunset bearing.
#' @return Numeric declination in degrees. With loc = FALSE the result is geocentric. Supplying a location applies the parallax correction for the chosen altitude.
#' @param year Astronomical year. The default is the current year in the package settings. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' @param loc Latitude in degrees north or a \code{skyscapeR.horizon} object. Omit this argument or use \code{FALSE} for geocentric declination.
#' @param parallax Mean horizontal parallax in degrees. The default is 0.00224 degrees for the Sun.
#' @param altitude Altitude of the sun. Defaults to 0 degrees.
#' @param verbose Logical to control output of warnings and messages.
#' Defaults to TRUE.
#' @export
#' @seealso \code{\link{obliquity}}, \code{\link{dS}}, \code{\link{eq}},
#' \code{\link{zenith}}, \code{\link{antizenith}}, \code{\link{spatial.equinox}},
#' \code{\link{parallax.corr}}
#' @examples
#' # June Solstice geocentric declination for year 4001 BC:
#' jS(-4000)
#'
#' # Topocentric declination for same year and latitude of 50 degres N:
#' jS(-4000, loc=50)
jS = function(year = skyscapeR.env$cur.year, loc = FALSE, parallax = 0.00224, altitude = 0, verbose=TRUE) {
  if (is.character(year)) year <- astroYear(year)
  aux <- obliquity(year)

  if (is(loc,'logical')) {
    if (verbose) { .sky_message('Returning geocentric declination because no site was supplied.') }
  } else {
    aux <- aux - parallax.corr(parallax, loc, altitude)
  }
  return(aux)
}

#' Estimate solar declination at an equinox
#' @description Estimate the solar declination associated with the astronomical equinox, optionally
#' applying a site-dependent parallax correction.
#' @details The geometric reference at an equinox is the celestial equator, with geocentric
#' declination zero. Supplying a location subtracts the package's approximate parallax
#' correction for the latitude and specified altitude. The result is a declination target
#' in degrees; no equinox date, event time or horizon crossing is calculated.
#' @return Numeric declination in degrees. With loc = FALSE the result is geocentric. Supplying a location applies the parallax correction for the chosen altitude.
#' @param loc Latitude in degrees north or a \code{skyscapeR.horizon} object. Omit this argument or use \code{FALSE} for geocentric declination.
#' @param parallax Mean horizontal parallax in degrees. The default is 0.00224 degrees for the Sun.
#' @param altitude Altitude of the sun. Defaults to 0 degrees.
#' @param verbose Logical to control output of warnings and messages.
#' Defaults to TRUE. This affects messages only, not the calculated declination.
#' @export
#' @seealso \code{\link{obliquity}}, \code{\link{jS}}, \code{\link{zenith}}, \code{\link{antizenith}}, \code{\link{spatial.equinox}}, \code{\link{parallax.corr}}, \code{\link{dS}}
#' @examples
#' # Equinoctial geocentric declination:
#' eq()
#'
#' # The same geocentric result without a message:
#' eq(verbose=FALSE)
#'
#' # Topocentric declination at latitude 50 degrees N:
#' eq(loc=50)
eq = function(loc=FALSE, parallax = 0.00224, altitude = 0, verbose=TRUE) {
  aux <- 0

  if (is(loc,'logical')) {
    if (verbose) { .sky_message('Returning geocentric declination because no site was supplied.') }
  } else {
    aux <- aux - parallax.corr(parallax, loc, altitude)
  }
  return(aux)
}

#' Estimate solar declination at the zenith
#' @description Determine the solar declination associated with an overhead passage at the observer's
#' latitude, where such a passage is possible.
#' @details The geometric zenith has declination equal to the site latitude. The function
#' first checks the latitude against its June and December solar limits at the default
#' epoch. Within those limits it subtracts the supplied approximate solar parallax
#' correction from the latitude. It returns NULL outside the permitted band. This supplies
#' a declination criterion rather than the dates of the passages.
#' @return Solar declination in degrees with the requested parallax correction, or NULL when the latitude falls outside the solstitial limits used by the function.
#' @param loc This can be either the latitude of the
#' location, or a \emph{skyscapeR.horizon} object.
#' @param parallax Mean horizontal parallax in degrees. The default is 0.00224 degrees for the Sun.
#' @param altitude Altitude of the sun. Defaults to 0 degrees.
#' @export
#' @seealso \code{\link{jS}}, \code{\link{dS}}, \code{\link{eq}},
#' \code{\link{antizenith}}, \code{\link{spatial.equinox}}, \code{\link{parallax.corr}}
#' @examples
#' # Zenith sun declination for Mexico City:
#' zenith(19.419)
#'
#' # There is no zenith sun phenomenon in London:
#' zenith(51.507)
zenith = function(loc, parallax = 0.00224, altitude = 0) {
  if (is(loc,'skyscapeR.horizon')) {
    lat <- loc$metadata$georef[1]
  } else { lat <- loc[1] }

  if (lat > jS(loc=loc) | lat < dS(loc=loc)) {
    return(NULL)
  } else {
    aux <- lat
    aux <- aux - parallax.corr(parallax, loc, altitude)
    return(aux)
  }
}

#' Estimate solar declination at the anti-zenith
#' @description Determine the solar declination associated with the point directly beneath an observer,
#' where the Sun can attain that declination.
#' @details The direction opposite the zenith corresponds geometrically to minus the site latitude.
#' The function checks the site against the default-epoch solstitial latitude band,
#' then subtracts the approximate solar parallax correction from the negative latitude. It
#' returns NULL outside that band. This is a target declination for comparison with
#' orientations, not an observable solar transit below ground.
#' @return Solar declination in degrees with the requested parallax correction, or NULL when the latitude falls outside the solstitial limits used by the function.
#' @param loc This can be either the latitude of the
#' location, or a \emph{skyscapeR.horizon} object.
#' @param parallax Mean horizontal parallax in degrees. The default is 0.00224 degrees for the Sun.
#' @param altitude Altitude of the sun. Defaults to 0 degrees.
#' @export
#' @seealso \code{\link{jS}}, \code{\link{dS}}, \code{\link{eq}},
#' \code{\link{zenith}}, \code{\link{spatial.equinox}}, \code{\link{parallax.corr}}
#' @examples
#' # Anti-zenith sun declination for Mexico City:
#' antizenith(19.419)
#'
#' # The same location supplied as a horizon object:
#' hor <- createHor(c(0,90,180,270,360), c(0,2,1,3,0),
#'                  loc=c(19.419,-99.133,2240))
#' antizenith(hor)
#'
#' # There is no anti-zenith solar phenomenon in London:
#' antizenith(51.507)
antizenith = function(loc, parallax = 0.00224, altitude = 0) {
  if (is(loc,'skyscapeR.horizon')) {
    lat <- loc$metadata$georef[1]
  } else { lat <- loc[1] }

  if (lat > jS(loc=loc) | lat < dS(loc=loc)) {
    return(NULL)
  } else {
    aux <- -lat
    aux <- aux - parallax.corr(parallax, loc, altitude)
    return(aux)
  }
}

#' Estimate the spatial equinox from a site horizon
#' @description Find the directions halfway between the two solstitial sunrise bearings and between the
#' two solstitial sunset bearings on a particular horizon.
#' @details A spatial midpoint is defined by the local terrain and need not coincide with the
#' astronomical equinox. The function uses the positive and negative obliquity at
#' the requested year as the geocentric solstitial declinations. Each fixed track
#' is followed through one rotation from lower culmination to lower culmination.
#' Horizontal solar parallax is applied to geometric altitude before atmospheric
#' refraction. The terrain is interpolated linearly in azimuth.
#'
#' On the rising half of each track, the first emergence is selected; on the
#' setting half, the final disappearance is selected. The circular midpoint of
#' the two solstitial bearings is then found separately for sunrise and sunset.
#' At each midpoint, refraction and horizontal parallax are inverted to recover
#' the geocentric declination. Atmospheric settings are identical in both stages.
#'
#' Ruggles (1997) distinguishes the spatial midpoint of the solstitial horizon
#' bearings from the midpoint in time between the solstices.
#' @return A list containing azimuth and declination. Each contains rise and set
#' values in degrees for the supplied horizon. Declinations are geocentric.
#' An undefined side returns NA for both values, with a warning; the other side
#' is retained when defined.
#' @param hor A \emph{skyscapeR.horizon} with finite, increasing azimuths covering
#' 0 to 360 degrees and matching endpoint altitudes.
#' @param parallax Mean horizontal parallax in degrees. The default is 0.00224 degrees for the Sun.
#' @param year One astronomical year, defaulting to the package's current-year
#' setting. Historical labels such as \code{"3250 BC"} are accepted.
#' @param refraction Whether to apply atmospheric refraction. Defaults to the
#' package setting, normally TRUE.
#' @param atm Atmospheric pressure in hPa. Defaults to the package setting.
#' Zero requests an estimate from site elevation.
#' @param temp Air temperature in degrees Celsius. Defaults to the package setting.
#' @export
#' @references Ruggles, C. L. N. (1997). Whose Equinox? Journal for the History of Astronomy 28, S45-S96. \doi{10.1177/002182869702802205}.
#' @seealso \code{\link{jS}}, \code{\link{dS}}, \code{\link{eq}},
#' \code{\link{zenith}}, \code{\link{antizenith}}, \code{\link{findTargets}}
#' @examples
#' \dontrun{
#' hor <- createHWT(34.174051531543405, 110.81299818694872, name='Xipo, Lingbao')
#' spatial.equinox(hor, year="3250 BC")
#' }
#' @section Crossing conventions and limitations:
#' Both rising and setting tracks are considered in either hemisphere. The
#' search uses the terrain-aware crossing machinery also used by
#' \code{\link{riseset}}, including terrain vertices and interior extrema.
#' Tangencies without a change of visibility are excluded. Multiple terrain
#' crossings are resolved using the first-emergence/final-disappearance rule;
#' this need not be the locally preferred observational definition of sunrise
#' or sunset. Inspect complex skylines before interpretation.
#'
#' If either solstitial track lacks the required crossing, the corresponding
#' midpoint is undefined. This includes polar non-rising/non-setting cases and
#' complete terrain obstruction. Exactly opposite bearings have no unique
#' circular midpoint and also return NA.
#'
#' These are fixed-declination solar-centre tracks, not a dated solar ephemeris
#' or upper-limb calculation. The mean horizontal parallax is fixed, and the
#' refraction model cannot represent local atmospheric anomalies. Results are
#' limited by horizon resolution and atmospheric assumptions. Solar semidiameter
#' is not included. The result is a spatial midpoint, not the Sun's rising or
#' setting at the astronomical equinox.
#' @section Provenance:
#' Structured results carry a \code{provenance} attribute. Read it with
#' \code{attr(result, "provenance")}; its \code{packages} entry records the loaded
#' skyscapeR and swephR versions, with NA when a namespace version is unavailable.
#' This attribute adds no list entries or table columns. Ordinary subsetting or
#' transformations may discard attributes; save the original result for reference.
#' It does not identify ephemeris files, catalogue revisions or a development
#' commit. Retain those and the analysis script separately for exact reproduction.
spatial.equinox = function(hor, parallax = 0.00224,
                          year = skyscapeR.env$cur.year,
                          refraction = skyscapeR.env$refraction,
                          atm = skyscapeR.env$atm, temp = skyscapeR.env$temp) {
  if (!inherits(hor,'skyscapeR.horizon'))
    stop('Full horizon profile required for spatial equinox calculation.')
  if (is.character(year)) year <- astroYear(year)
  if (!is.numeric(year) || length(year)!=1L || !is.finite(year))
    stop('year must specify one finite astronomical year.')
  if (length(parallax)!=1L || !is.finite(parallax) || parallax<0 || parallax>=90)
    stop('parallax must be a finite angle from 0 to less than 90 degrees.')
  if (!is.logical(refraction) || length(refraction)!=1L || is.na(refraction))
    stop('refraction must be TRUE or FALSE.')
  if (length(atm)!=1L || !is.finite(atm) || atm<0 ||
      length(temp)!=1L || !is.finite(temp) || temp<= -273.15)
    stop('Supply finite atmospheric settings with nonnegative pressure.')
  loc <- as.numeric(hor$metadata$georef)
  if(length(loc)!=3L || any(!is.finite(loc)) || abs(loc[1])>90)
    stop('Horizon georeference must contain latitude, longitude and elevation.')
  horizon <- .horizon_profile(hor$data$az,hor$data$alt)
  pressure <- if(atm==0) 1013.25*(1-.0065*loc[3]/288)^5.255 else atm
  jd <- 2451545
  meridian <- swephR::swe_sidtime(jd)*15+loc[2]
  limits <- c(-1,1)*obliquity(year)
  bearings <- matrix(NA_real_,2,2,dimnames=list(NULL,c('rise','set')))
  for(i in seq_along(limits)) {
    # One rotation from lower culmination to the next, with noon at t=0.5.
    # Horizontal solar parallax is applied before refraction.
    position <- function(t) {
      ra <- meridian+180-360*t
      p <- swephR::swe_azalt(jd,1L,c(loc[2],loc[1],loc[3]),pressure,temp,
                              c(ra,limits[i],1))$xaz
      h <- p[2]*pi/180
      altitude <- atan2(sin(h)-sin(parallax*pi/180),cos(h))*180/pi
      if(refraction) altitude <- swephR::swe_refrac_extended(
        altitude,loc[3],pressure,temp,.0065,0L)$return
      c((p[1]+180)%%360,altitude)
    }
    roots <- .horizon_crossings(position,0,hor$data$az,hor$data$alt,
                               .duration=1,.all=TRUE,horizon=horizon)
    if(length(roots)) {
      clearance <- function(t) {p<-position(t);p[2]-horizon(p[1])}
      direction <- vapply(roots,function(t)
        sign(clearance(t+1e-7)-clearance(t-1e-7)),numeric(1))
      # Ignore touches and retain crossings in the corresponding half-cycle.
      rising <- roots[direction>0 & roots<=.5]
      setting <- roots[direction<0 & roots>=.5]
      if(length(rising)) bearings[i,1] <- position(min(rising))[1]
      if(length(setting)) bearings[i,2] <- position(max(setting))[1]
    }
  }
  midpoint <- function(x) {
    if(anyNA(x)) return(NA_real_)
    delta <- (x[2]-x[1]+180)%%360-180
    if(abs(abs(delta)-180)<1e-8) return(NA_real_)
    (x[1]+delta/2)%%360
  }
  az <- apply(bearings,2,midpoint)
  declination <- vapply(az,function(a) {
    if(is.na(a)) return(NA_real_)
    altitude <- horizon(a)
    if(refraction) altitude <- swephR::swe_refrac_extended(
      altitude,loc[3],pressure,temp,.0065,1L)$return
    # Invert horizontal parallax to recover the geocentric direction.
    altitude <- altitude+asin(sin(parallax*pi/180)*cos(altitude*pi/180))*180/pi
    swephR::swe_azalt_rev(jd,1L,c(loc[2],loc[1],loc[3]),
                        c(a-180,altitude))$xout[2]
  },numeric(1))
  if(anyNA(az)) warning('Spatial equinox undefined for ',
    paste(names(az)[is.na(az)],collapse=' and '),
    ': both solstitial crossings and an unambiguous midpoint are required.',call.=FALSE)
  .with_astronomy_provenance(list(azimuth=as.list(az),declination=as.list(declination)))
}

#' Estimate northern minor lunar standstill declination
#' @description Estimate the northern declination limit of a minor lunar standstill for comparison with
#' a proposed alignment.
#' @details A minor standstill defines a limiting lunar declination at the epoch of interest. The
#' model uses the obliquity minus a fixed 5.145 + 0.145 degrees, then applies the positive
#' sign for the northern limit. With a location it subtracts the parallax correction for
#' latitude and altitude. The default lunar parallax is 0.952 degrees. This approximation
#' defines a target limit; it does not search a time series for an actual lunar event.
#'
#' Sprajc (2016) describes the limits set by obliquity and the mean lunar
#' inclination. This package also includes the fixed 0.145-degree term above.
#' @return Numeric declination in degrees. With loc = FALSE the result is geocentric. Supplying a location applies the parallax correction for the chosen altitude.
#' @param year Astronomical year. The default is the current year in the package settings. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' @param loc Latitude in degrees north or a \code{skyscapeR.horizon} object. Omit this argument or use \code{FALSE} for geocentric declination.
#' @param parallax Mean horizontal parallax in degrees. The default is 0.952 degrees for the Moon.
#' @param altitude Altitude of the Moon. Defaults to 0 degrees.
#' @param verbose Logical to control output of warnings and messages.
#' Defaults to TRUE.
#' @export
#' @references \enc{Šprajc}{Sprajc}, I. (2016). Lunar alignments in Mesoamerican architecture. Anthropological Notebooks 22(3), 61-85.
#' @seealso \code{\link{smnLX}}, \code{\link{nMjLX}}, \code{\link{sMjLX}},
#' \code{\link{parallax.corr}}
#' @examples
#' # Northern minor Lunar Extreme geocentric declination for year 2501 BC:
#' nmnLX(-2500)
#'
#' # Topocentric declination for same year and latitude of 50 degrees N:
#' nmnLX(-2500, loc=50)
nmnLX = function(year = skyscapeR.env$cur.year, loc=FALSE, parallax = 0.952, altitude = 0, verbose = TRUE) {
  if (is.character(year)) year <- astroYear(year)
  aux <- obliquity(year) - (5.145+0.145)

  if (is(loc,'logical')) {
    if (verbose) { .sky_message('Returning geocentric declination because no site was supplied.') }
  } else {
    aux <- aux - parallax.corr(parallax, loc, altitude)
  }
  return(aux)
}


#' Estimate southern minor lunar standstill declination
#' @description Estimate the southern declination limit of a minor lunar standstill for comparison with
#' a proposed alignment.
#' @details A minor standstill defines a limiting lunar declination at the epoch of interest. The
#' model uses the obliquity minus a fixed 5.145 + 0.145 degrees, then applies the negative
#' sign for the southern limit. With a location it subtracts the parallax correction for
#' latitude and altitude. The default lunar parallax is 0.952 degrees. This approximation
#' defines a target limit; it does not search a time series for an actual lunar event.
#'
#' Sprajc (2016) describes the limits set by obliquity and the mean lunar
#' inclination. This package also includes the fixed 0.145-degree term above.
#' @return Numeric declination in degrees. With loc = FALSE the result is geocentric. Supplying a location applies the parallax correction for the chosen altitude.
#' @param year Astronomical year. The default is the current year in the package settings. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' @param loc Latitude in degrees north or a \code{skyscapeR.horizon} object. Omit this argument or use \code{FALSE} for geocentric declination.
#' @param parallax Mean horizontal parallax in degrees. The default is 0.952 degrees for the Moon.
#' @param altitude Altitude of the Moon. Defaults to 0 degrees.
#' @param verbose Logical to control output of warnings and messages.
#' Defaults to TRUE.
#' @export
#' @references \enc{Šprajc}{Sprajc}, I. (2016). Lunar alignments in Mesoamerican architecture. Anthropological Notebooks 22(3), 61-85.
#' @seealso \code{\link{nmnLX}}, \code{\link{nMjLX}}, \code{\link{sMjLX}},
#' \code{\link{parallax.corr}}
#' @examples
#' # Southern minor Lunar Extreme geocentric declination for year 2501 BC:
#' smnLX(-2500)
#'
#' # Topocentric declination for same year and latitude of 50 degrees N:
#' smnLX(-2500, loc=50)
smnLX = function(year = skyscapeR.env$cur.year, loc=FALSE, parallax = 0.952, altitude = 0, verbose = TRUE) {
  if (is.character(year)) year <- astroYear(year)
  aux <- -(obliquity(year) - (5.145+0.145))

  if (is(loc,'logical')) {
    if (verbose) { .sky_message('Returning geocentric declination because no site was supplied.') }
  } else {
    aux <- aux - parallax.corr(parallax, loc, altitude)
  }
  return(aux)
}


#' Estimate northern major lunar standstill declination
#' @description Estimate the northern declination limit of a major lunar standstill for comparison with
#' a proposed alignment.
#' @details A major standstill defines a limiting lunar declination at the epoch of interest. The
#' model uses the obliquity plus a fixed 5.145 + 0.145 degrees, then applies the positive
#' sign for the northern limit. With a location it subtracts the parallax correction for
#' latitude and altitude. The default lunar parallax is 0.952 degrees. This approximation
#' defines a target limit; it does not search a time series for an actual lunar event.
#'
#' Sprajc (2016) describes the limits set by obliquity and the mean lunar
#' inclination. This package also includes the fixed 0.145-degree term above.
#' @return Numeric declination in degrees. With loc = FALSE the result is geocentric. Supplying a location applies the parallax correction for the chosen altitude.
#' @param year Astronomical year. The default is the current year in the package settings. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' @param loc Latitude in degrees north or a \code{skyscapeR.horizon} object. Omit this argument or use \code{FALSE} for geocentric declination.
#' @param parallax Mean horizontal parallax in degrees. The default is 0.952 degrees for the Moon.
#' @param altitude Altitude of the Moon. Defaults to 0 degrees.
#' @param verbose Logical to control output of warnings and messages.
#' Defaults to TRUE.
#' @export
#' @references \enc{Šprajc}{Sprajc}, I. (2016). Lunar alignments in Mesoamerican architecture. Anthropological Notebooks 22(3), 61-85.
#' @seealso \code{\link{nmnLX}}, \code{\link{smnLX}}, \code{\link{sMjLX}},
#' \code{\link{parallax.corr}}
#' @examples
#' # Northern major Lunar Extreme geocentric declination for year 2501 BC:
#' nMjLX(-2500)
#'
#' # Topocentric declination for same year and latitude of 50 degrees N:
#' nMjLX(-2500, loc=50)
nMjLX = function(year = skyscapeR.env$cur.year, loc=FALSE, parallax = 0.952, altitude = 0, verbose = TRUE) {
  if (is.character(year)) year <- astroYear(year)
  aux <- obliquity(year) + (5.145+0.145)

  if (is(loc,'logical')) {
    if (verbose) { .sky_message('Returning geocentric declination because no site was supplied.') }
  } else {
    aux <- aux - parallax.corr(parallax, loc, altitude)
  }
  return(aux)
}


#' Estimate southern major lunar standstill declination
#' @description Estimate the southern declination limit of a major lunar standstill for comparison with
#' a proposed alignment.
#' @details A major standstill defines a limiting lunar declination at the epoch of interest. The
#' model uses the obliquity plus a fixed 5.145 + 0.145 degrees, then applies the negative
#' sign for the southern limit. With a location it subtracts the parallax correction for
#' latitude and altitude. The default lunar parallax is 0.952 degrees. This approximation
#' defines a target limit; it does not search a time series for an actual lunar event.
#'
#' Sprajc (2016) describes the limits set by obliquity and the mean lunar
#' inclination. This package also includes the fixed 0.145-degree term above.
#' @return Numeric declination in degrees. With loc = FALSE the result is geocentric. Supplying a location applies the parallax correction for the chosen altitude.
#' @param year Astronomical year. The default is the current year in the package settings. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' @param loc Latitude in degrees north or a \code{skyscapeR.horizon} object. Omit this argument or use \code{FALSE} for geocentric declination.
#' @param parallax Mean horizontal parallax in degrees. The default is 0.952 degrees for the Moon.
#' @param altitude Altitude of the Moon. Defaults to 0 degrees.
#' @param verbose Logical to control output of warnings and messages.
#' Defaults to TRUE.
#' @export
#' @references \enc{Šprajc}{Sprajc}, I. (2016). Lunar alignments in Mesoamerican architecture. Anthropological Notebooks 22(3), 61-85.
#' @seealso \code{\link{nmnLX}}, \code{\link{smnLX}}, \code{\link{nMjLX}},
#' \code{\link{parallax.corr}}
#' @examples
#' # Southern major Lunar Extreme geocentric declination for year 2501 BC:
#' sMjLX(-2500)
#'
#' # Topocentric declination for same year and latitude of 50 degrees N:
#' sMjLX(-2500, loc=50)
sMjLX = function(year = skyscapeR.env$cur.year, loc=FALSE, parallax = 0.952, altitude = 0, verbose = TRUE) {
  if (is.character(year)) year <- astroYear(year)
  aux <- -(obliquity(year) + (5.145+0.145))

  if (is(loc,'logical')) {
    if (verbose) { .sky_message('Returning geocentric declination because no site was supplied.') }
  } else {
    aux <- aux - parallax.corr(parallax, loc, altitude)
  }
  return(aux)
}

#' Find equinoctial full moonrises and moonsets
#' @description Identify the near-full Moon whose rising or setting position crosses that of the Sun
#' around the spring or autumn equinox.
#' @details Silva and Pimenta (2012) identify the crossover by the first near-full
#' Moon on the seasonal side of the Sun. Spring and autumn refer to the March
#' and September equinoxes. At spring moonrise the lunar azimuth must exceed
#' the solar azimuth, and at spring moonset it must be smaller. Autumn reverses
#' both inequalities. Bearings are topocentric, clockwise from north.
#'
#' The search runs chronologically from 60 days before to 65 days after the
#' astronomical equinox. Moonrise is paired with the first sunrise on the same
#' local civil date. Moonset is paired with the most recent sunset before it,
#' considering the first solar centre crossing on the current and preceding
#' local civil dates. A candidate without such a sunset is skipped. Local-day
#' boundaries account for timezone changes, including daylight saving time.
#' A coordinate vector
#' specifies a level horizon. A horizon object supplies the horizon for both
#' bodies. Illumination is tested at the lunar event itself. A date qualifies
#' as soon as the illumination and strict bearing conditions both hold. A
#' lower illumination threshold can select an earlier lunation. If no date
#' qualifies the function reports the year and event and stops.
#'
#' Dates and clock times use the supplied timezone and calendar. Changing
#' the timezone can change which events share a civil date. Declinations use
#' the package coordinate setting (topocentric by default). For comparison
#' with Figure 3 of Silva and Pimenta, request geocentric declinations with
#' \code{skyscapeR.vars(dec="geo")}. Event selection uses topocentric bearings
#' in either output frame.
#'
#' @section Comparison with the original simulations:
#' With the corrected preceding-sunset pairing, 2,562 of 2,568 moonsets in
#' the supplied 1620--2047 series at latitudes 0, 60 N and 60 S agree within
#' 0.02 degrees in lunar declination. The comparison uses a level horizon,
#' longitude 0, elevation 400 m, UTC, Gregorian dates, geocentric declinations,
#' refraction, pressure 1013.25 mbar, temperature 15 degrees Celsius and the
#' default illumination threshold. Five remaining cases lie just above the
#' 0.99 illumination boundary; one lies about 0.002 degrees beyond the azimuth
#' crossover boundary. The original series selects the following morning in
#' each case. Thresholds have not been tuned to reproduce those selections.
#'
#' The original 40-degree N series retains systematic coordinate differences.
#' Diagnostic comparisons suggest a different coordinate-conversion or timing
#' convention, but the original settings have not been established. Agreement
#' with these simulation tables is not independent observational validation.
#' The paper does not specify all atmospheric settings; the values above describe
#' this comparison, not confirmed settings of the original simulations.
#'
#' Astronomical positions require Swiss Ephemeris files covering the epoch. We
#' recommend installing \pkg{swephRdata}; see \code{\link{skyscapeR-package}} for setup.
#' @references Silva, F. and Pimenta, F. (2012). The Crossover of the Sun and the Moon. Journal for the History of Astronomy 43. \doi{10.1177/002182861204300204}.
#' @seealso \code{\link{CC}}, \code{\link{spatial.equinox}}, \code{\link{zenith}}, \code{\link{antizenith}}, \code{\link{solar.date}}, \code{\link{findTargets}}, \code{\link{planet.extremes}}
#' @return A list with an Event label and Moon and Sun data frames. Each row gives year, date, time, azimuth and declination for the selected rising or setting event. Angles are in degrees.
#' @param season Either 'spring' or 'autumn'. Default is 'spring'.
#' @param rise Logical to choose whether to calculate Equinoctial Full
#' Moon rises or sets. Defaults to \emph{TRUE}.
#' @param year Astronomical year or years. Can be either a single value (the year), Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' two values (range of years), or a vector of years.
#' @param loc This can be either a \emph{skyscapeR.horizon} object, or a vector with the
#' latitude, longitude and elevation of the site, in this order.
#' @param min.phase Minimum illuminated fraction of the Moon, from 0 to 1. The default 0.99 means 99 percent illumination.
#' @param refraction Logical. Apply atmospheric refraction. Defaults to the session setting in \code{\link{skyscapeR.vars}}.
#' @param atm Atmospheric pressure in mbar. Defaults to the session setting in \code{\link{skyscapeR.vars}}.
#' @param temp Air temperature in degrees Celsius. Defaults to the session setting in \code{\link{skyscapeR.vars}}.
#' @param calendar Calendar for civil dates. Use \code{"G"} for Gregorian or \code{"J"} for Julian.
#' @param timezone Timezone for civil timestamps. Use UTC or a fixed offset for dates before 1900. See \code{\link{time2jd}} for named zones and offsets.
#' @export
#' @examples
#' if (requireNamespace("swephRdata", quietly = TRUE)) {
#'   # Spring Full Moon at latitude 35 degrees north in the year 2018
#'   EFM(year=2018, loc=c(35,-8,100))
#'   
#'   # Autumn Full Moons over a 3 year period
#'   \dontrun{
#'   EFM(season='autumn', year=c(2018,2020), loc=c(35,-8,100))
#'   }
#' }
#' @section Scope of the event model:
#' The illumination fraction is a selection criterion, not a naked-eye detection
#' model. The calculation does not include weather, sky brightness or extinction.
#' It uses centre crossings rather than the first appearance of the upper limb.
#' Only the first qualifying event is returned, so changes in illumination
#' threshold, civil-day convention or horizon can change the selected lunation.
#' Remaining differences from the original simulations are described above.
#' Numerical precision does not establish agreement with a published distribution.
#' Preserve all inputs when comparing series.
#' @section Provenance:
#' Structured results carry a \code{provenance} attribute. Read it with
#' \code{attr(result, "provenance")}; its \code{packages} entry records the loaded
#' skyscapeR and swephR versions, with NA when a namespace version is unavailable.
#' This attribute adds no list entries or table columns. Ordinary subsetting or
#' transformations may discard attributes; save the original result for reference.
#' It does not identify ephemeris files, catalogue revisions or a development
#' commit. Retain those and the analysis script separately for exact reproduction.
EFM <- function(year, loc, season='spring', rise=TRUE, min.phase=.99, timezone, calendar,
                refraction, atm, temp) {
  if (is.character(year)) year <- astroYear(year)
  if (missing(timezone)) timezone <- skyscapeR.env$timezone
  if (missing(calendar)) calendar <- skyscapeR.env$calendar
  if (missing(refraction)) refraction <- skyscapeR.env$refraction
  if (missing(atm)) atm <- skyscapeR.env$atm
  if (missing(temp)) temp <- skyscapeR.env$temp
  season <- match.arg(season, c("spring", "autumn"))
  if (!is.logical(rise) || length(rise)!=1L || is.na(rise)) stop("rise must be TRUE or FALSE.")
  if (!is.logical(refraction) || length(refraction)!=1L || is.na(refraction)) stop("refraction must be TRUE or FALSE.")
  if (!is.numeric(min.phase) || length(min.phase)!=1L || !is.finite(min.phase) || min.phase<0 || min.phase>1)
    stop("min.phase must be an illuminated fraction from 0 to 1.")
  if (!is.numeric(year) || !length(year) || any(!is.finite(year)) || any(year!=trunc(year)))
    stop("year must contain finite integer years.")
  checkYear(year)
  if (length(year)==2L) year <- seq(year[1], year[2], by=1)
  geo <- if (inherits(loc,"skyscapeR.horizon")) loc$metadata$georef else loc
  if (!is.numeric(geo) || length(geo)!=3L || any(!is.finite(geo)) || abs(geo[1])>90 || abs(geo[2])>180)
    stop("loc must supply latitude, longitude and elevation.")
  swephR::swe_set_topo(geo[2],geo[1],geo[3])
  terrain <- inherits(loc,"skyscapeR.horizon") && length(unique(loc$data$alt))>1L
  zone <- .timezone_spec(timezone)
  .calendar_flag(calendar)
  spring <- season=="spring"
  target <- if (spring==rise) 1 else -1
  direction <- if (rise) -1 else 1  # .cc_event uses -1 for rising
  event <- function(day, body) .cc_event(day,body,direction,loc,geo,refraction,atm,temp)
  row <- function(y, jd, body, bearing) {
    stamp <- strsplit(jd2time(jd,timezone,calendar)," ",fixed=TRUE)[[1]]
    dec <- body.position(obj=body,time=jd,timezone,calendar,loc=geo,refraction=refraction,
                         atm=atm,temp=temp,verbose=FALSE)$equatorial$Dec
    data.frame(year=y,date=stamp[1],time=stamp[2],azimuth=bearing,declination=dec)
  }
  moons <- suns <- vector("list",length(year))
  for (k in seq_along(year)) {
    # Bracket the astronomical season in Gregorian coordinates irrespective of the output calendar.
    centre <- swephR::swe_julday(year[k],if (spring) 3L else 9L,21L,0,1L)
    longitude <- if (spring) 0 else 180
    equinox <- stats::uniroot(function(t) .cc_wrap(.cc_calc(t,0,2)[1]-longitude),
                              centre+c(-20,20),tol=1e-8)$root
    local.days <- floor(equinox)+(-62:67)-.5
    midnights <- if (!is.null(zone$offset)) local.days-zone$offset/86400 else
      vapply(local.days,.efm_midnight,numeric(1),zone=zone$zone)
    for (i in 2L:(length(midnights)-1L)) {
      day <- midnights[i]; end <- midnights[i+1L]
      if (end<=day) next
      # A civil day spans at most one lunar illumination extremum. Bound its
      # maximum before paying for a terrain crossing, including both endpoints.
      if (terrain && all(moonphase(c(day,end))<min.phase)) {
        peak <- stats::optimize(moonphase,c(day,end),maximum=TRUE,tol=1e-7)
        if (peak$objective+1e-10<min.phase) next
      }
      moon <- event(day,1L)
      if (!is.finite(moon) || moon<day || moon>=end || moon<equinox-60 || moon>equinox+65) next
      if (moonphase(moon)<min.phase) next
      if (rise) {
        sun <- event(day,0L)
        if (!is.finite(sun) || sun<day || sun>=end) next
      } else {
        # Use actual local-day boundaries, including 23/25-hour DST days.
        # Consider today's sunset too, in case it precedes this moonset.
        starts <- c(midnights[i-1L],day)
        stops <- c(day,end)
        sunsets <- vapply(starts,event,numeric(1),body=0L)
        valid <- is.finite(sunsets) & sunsets>=starts & sunsets<stops & sunsets<moon
        if (!any(valid)) next
        sun <- max(sunsets[valid])
      }
      moon.az <- .cc_position(moon,1L,geo,refraction,atm,temp)$az
      sun.az <- .cc_position(sun,0L,geo,refraction,atm,temp)$az
      if (target*(moon.az-sun.az)<=0) next
      moons[[k]] <- row(year[k],moon,"moon",moon.az)
      suns[[k]] <- row(year[k],sun,"sun",sun.az)
      break
    }
    if (is.null(moons[[k]])) stop("No qualifying ",season," full Moon ",
                                  if (rise) "rise" else "set"," for year ",year[k],".")
  }
  .with_astronomy_provenance(list(Event=paste(season,"full moon"),Moon=do.call(rbind,moons),Sun=do.call(rbind,suns)))
}

# First instant of the local date, including dates whose clocks change at midnight.
.efm_midnight <- function(local.day, zone) {
  before <- .named_offset(local.day-1,zone)
  after <- .named_offset(local.day+1,zone)
  if (before==after) return(local.day-before/86400)
  candidates <- local.day-sort(unique(c(before,after)),decreasing=TRUE)/86400
  for (t in candidates)
    if (abs(t+.named_offset(t,zone)/86400-local.day)<.1/86400) return(t)
  # Midnight was skipped. Locate the first valid local clock time on or after it.
  lo <- floor((min(candidates)-2440587.5)*86400)-1
  hi <- ceiling((max(candidates)-2440587.5)*86400)+1
  while (hi-lo>1) {
    mid <- floor((lo+hi)/2)
    t <- 2440587.5+mid/86400
    if (t+.named_offset(t,zone)/86400<local.day) lo <- mid else hi <- mid
  }
  2440587.5+hi/86400
}

#' Extreme declinations for five visible planets
#' @description Estimate the most northerly and southerly declinations reached by a visible planet near
#' an archaeological epoch.
#' @details An epoch alone does not specify the planet's position in its cycle. The search extends
#' three quarters of the package's planet-specific period on either side of the epoch,
#' rounds the limiting years, and samples at daily intervals from noon UT on 1 January of
#' the start year to the corresponding date in the end year. Each sample uses body.position
#' and the requested site if supplied. The returned extrema are the range of these samples;
#' finer extrema can fall between sampled instants.
#'
#' Astronomical positions require Swiss Ephemeris files covering the epoch. We
#' recommend installing \pkg{swephRdata}; see \code{\link{skyscapeR-package}} for setup.
#' @seealso \code{\link{spatial.equinox}}, \code{\link{zenith}}, \code{\link{antizenith}}, \code{\link{solar.date}}, \code{\link{findTargets}}, \code{\link{EFM}}
#' @return A list containing planet, min.dec and max.dec. The two declinations are extrema of daily samples and are expressed in degrees.
#' @param obj Planet name, with this capitalisation: Mercury, Venus, Mars, Jupiter or Saturn.
#' @param epoch Earlier argument name for \code{year}, retained for existing scripts.
#' @param year Year for which to do calculations. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' @param loc This can be either a \emph{skyscapeR.horizon} object, or a vector with the
#' latitude, longitude and elevation of the site, in this order.
#' @export
#' @examples
#' \dontrun{
#' # Searches across multiple planetary cycles; may take more than five seconds.
#' if (requireNamespace("swephRdata", quietly = TRUE)) {
#'   planet.extremes('Jupiter', 2023)
#'   
#'   planet.extremes('Venus', "4001 BC")
#' }
#' }
#' @section Provenance:
#' Structured results carry a \code{provenance} attribute. Read it with
#' \code{attr(result, "provenance")}; its \code{packages} entry records the loaded
#' skyscapeR and swephR versions, with NA when a namespace version is unavailable.
#' This attribute adds no list entries or table columns. Ordinary subsetting or
#' transformations may discard attributes; save the original result for reference.
#' It does not identify ephemeris files, catalogue revisions or a development
#' commit. Retain those and the analysis script separately for exact reproduction.
planet.extremes <- function(obj, year, loc = NULL, epoch) {
  if (!missing(epoch)) {
    if (!missing(year)) stop("Supply year only once.",call.=FALSE)
    year <- epoch
  }
  if (is.character(year)) year <- astroYear(year)
  start <- round(year - period(obj)*3/4,0); end <- round(year + period(obj)*3/4,0)
  if (start-end==0) { end <- start+1 }
  start <- time2jd(paste0(start,'/01/01 12:00:00'), timezone='UTC')
  end <- time2jd(paste0(end,'/01/01 12:00:00'), timezone='UTC')
  tt <- seq(start, end)
  dec <- c()
  for (i in 1:length(tt)) {
    dec[i] <- body.position(obj, time=tt[i], loc=loc, verbose=F)$equatorial$Dec
  }
  out <- c()
  out$planet <- obj
  out$min.dec <- min(dec)
  out$max.dec <- max(dec)
  return(.with_astronomy_provenance(out))
}
