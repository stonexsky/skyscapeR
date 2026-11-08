#' Crescent crossovers of the Sun and Moon
#' @description Identify the first observable seasonal reversal in the relative
#' rising or setting positions of the Sun and a lunar crescent. These events can
#' be compared with archaeological orientations and need not occur close to a solstice.
#' @details Crescent crossovers describe a change in the relative horizon positions
#' of the luminaries across successive lunations. First (waxing) crescents are
#' compared at sunset and moonset; last (waning) crescents at sunrise and moonrise.
#' The method follows Silva and Pimenta (2012), sections 2.1--2.2. Seasons use the
#' paper's northern-hemisphere names: summer means the June-solstice crossover,
#' winter the December-solstice crossover, also for southern-hemisphere locations.
#'
#' For each conjunction, days are searched forwards for the first crescent and
#' backwards for the last crescent. Visibility is evaluated at sunset or sunrise,
#' respectively, using Caldwell and Laney's criterion: lower-limb lunar altitude
#' plus one third of the Sun--Moon angular separation must be at least 11.3 degrees.
#' The lower-limb altitude is topocentric. The arc of light is the geocentric
#' Sun--Moon separation, as explicitly defined by Caldwell and Laney; it excludes
#' atmospheric refraction. The lower limb is calculated from the instantaneous lunar distance
#' and a lunar radius of 1737.4 km. With refraction enabled, the centre altitude is
#' apparent and the angular radius is subtracted from it. This is an empirical
#' visibility prediction, not a guarantee under actual weather or extinction.
#'
#' SIM1 tests only the first qualifying evening or last qualifying morning of
#' each lunation. SIM2 also tests the following evening or preceding morning,
#' respectively. The original visibility day is preferred if it already meets
#' the crossover condition. Dates are searched within eight civil days of
#' conjunction and must remain crescent (elongation below 90 degrees).
#'
#' Each body's own centre-disc horizon event supplies its azimuth. The shortest
#' signed difference Moon minus Sun is negative after a June crossover and positive
#' after a December crossover, for both event types, with bearings clockwise from
#' north. The starting conjunction is two lunations before the last conjunction
#' at or before the solstice (the beginning of the lunation containing it).
#' This can start nearly three synodic months before the solstice; it is not
#' equivalent to starting at the first conjunction after solstice minus two
#' synodic months. From there,
#' successive lunations are searched for the first eligible crescent with the
#' required strict ordering; exact equality does not qualify. The preceding
#' eligible difference is retained for inspection but is not an additional
#' selection condition in the published algorithm. Missing lunations are skipped.
#' Conjunctions are searched until 150 days after the solstice; failure gives NA.
#' Later qualifying lunations are counted in diagnostics but do not replace the
#' first. Observations may fall in the next calendar year.
#'
#' A coordinate vector uses a zero-altitude astronomical horizon. A horizon object
#' uses its horizon for both rise/set events and a separate lower-limb clearance
#' check at the visibility assessment. The visibility formula still uses altitude
#' above the astronomical horizontal, not altitude above the terrain. This horizon
#' extension was not calibrated in the original paper. Events use the existing
#' centre-disc convention of \code{riseset}, including for the Sun.
#'
#' Swiss Ephemeris supplies apparent astronomical positions (light-time and
#' aberration corrections), with atmospheric refraction controlled separately.
#' Returned declinations are geocentric and without atmospheric refraction, as in
#' Figure 4 of the paper. Individual event selections do not yet reproduce every original simulation.
#' The starting-window convention and differences in the 40-degree N series
#' remain under validation. Agreement of a smoothed distribution does not imply
#' agreement of every year's selected crescent. Inspect the returned dates,
#' visibility margins and azimuth differences before interpreting a particular
#' event; the cause of all discrepancies is not established.
#' Years use astronomical numbering (0 is 1 BCE). Two ascending year endpoints
#' expand inclusively; longer vectors are evaluated as supplied. The year labels
#' the Gregorian year of the solstice, so a selected observation may occur in an
#' adjacent calendar year. The calendar argument controls date display and civil-day parsing.
#' No PNG or copy of the input horizon is stored in the output.
#'
#' Astronomical positions require Swiss Ephemeris files covering the epoch. We
#' recommend installing \pkg{swephRdata}; see \code{\link{skyscapeR-package}} for setup.
#' @section Comparison with the original simulations:
#' For SIM1, 5,024 of 5,136 events in the supplied 1620--2047 series at
#' latitudes 0, 60 N and 60 S agree within 0.02 degrees in lunar declination.
#' This comparison uses a level horizon, longitude 0, elevation 400 m, UTC,
#' Gregorian dates, refraction, pressure 1013.25 mbar and temperature 15 degrees
#' Celsius. The original 40-degree N series has additional systematic offsets.
#' Of the 112 discrepancies outside that series, inferred matching dates suggest
#' differences in the lunation search window and first/last crescent selection.
#' Some lie close to the visibility threshold. Alternative window and solar-limb
#' conventions improve some cases but disrupt others that already agree, so no
#' replacement convention has been adopted. The exact original settings remain
#' unresolved; the paper does not specify all atmospheric inputs. Inspect the
#' diagnostics before interpreting individual events. Agreement with simulation
#' tables is not independent observational validation, and these SIM1 counts do
#' not establish the accuracy of SIM2.
#' @param season Either \code{"winter"} (default) or \code{"summer"}.
#' @param crescent Either \code{"first"} (default, waxing/set) or
#'   \code{"last"} (waning/rise).
#' @param year One integer year, two ascending endpoints, or a vector of years. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' @param loc Latitude, longitude and elevation (degrees, degrees, metres), or
#'   a \code{skyscapeR.horizon} object.
#' @param simulation \code{"SIM1"} (default) or \code{"SIM2"}; see Details.
#' @param refraction Logical; use atmospheric refraction. Defaults to package settings.
#' @param atm Atmospheric pressure in mbar. Defaults to package settings.
#' @param temp Air temperature in degrees Celsius. Defaults to package settings.
#' @param timezone Output and observing-day timezone. Defaults to package settings.
#'   Use UTC or a fixed offset before 1900.
#' @param calendar \code{"G"} for Gregorian or \code{"J"} for Julian dates.
#'   Defaults to package settings.
#' @return A list with \code{Event}, paired \code{Moon} and \code{Sun} data
#' frames (year, date, time, azimuth, declination, UT Julian date), and a
#' \code{diagnostics} data frame. Diagnostics include status, solstice and
#' conjunction JDs, visibility-assessment JD, illuminated fraction, lower-limb
#' altitude, elongation, visibility margin, signed azimuth difference, preceding
#' difference, number of qualifying lunations, and day offset (0 for SIM1's day, 1 for SIM2's additional day).
#' All angles are degrees. Missing results retain the requested year,
#' use NA for observations and issue a warning. \code{settings} records conventions.
#' @references Silva, F. and Pimenta, F. (2012). The Crossover of the Sun and the
#' Moon. Journal for the History of Astronomy 43(2), 191--208.
#' \doi{10.1177/002182861204300204}.
#'
#' Caldwell, J. and Laney, C. D. (2001). First visibility of the lunar crescent.
#' African Skies 5, 15--23. The authors' online version defines dalt and arcl:
#' \url{https://astronomycenter.net/pdf/saao_2001.pdf}.
#' @seealso \code{\link{EFM}}, \code{\link{moonphase}},
#'   \code{\link{riseset}}, \code{\link{body.position}}, \code{\link{createHor}}
#' @export
#' @examples
#' if (requireNamespace("swephRdata", quietly = TRUE)) {
#'   CC(year=2026, loc=c(40, 0, 400), timezone="UTC")
#'   \dontrun{
#'   CC(year=c(2000, 2018), loc=c(40, 0, 400), season="summer", crescent="last",
#'   simulation="SIM2", timezone="UTC")
#'   CC(year=-2999, loc=c(40, 0, 400), timezone="UTC")
#'   }
#' }
#' @section Provenance:
#' Structured results carry a \code{provenance} attribute. Read it with
#' \code{attr(result, "provenance")}; its \code{packages} entry records the loaded
#' skyscapeR and swephR versions, with NA when a namespace version is unavailable.
#' This attribute adds no list entries or table columns. Ordinary subsetting or
#' transformations may discard attributes; save the original result for reference.
#' It does not identify ephemeris files, catalogue revisions or a development
#' commit. Retain those and the analysis script separately for exact reproduction.
CC <- function(year, loc, season="winter", crescent="first", simulation="SIM1",
               timezone=skyscapeR.env$timezone, calendar=skyscapeR.env$calendar,
               refraction=skyscapeR.env$refraction, atm=skyscapeR.env$atm, temp=skyscapeR.env$temp) {
  if (is.character(year)) year <- astroYear(year)
  season <- match.arg(season, c("winter", "summer"))
  crescent <- match.arg(crescent, c("first", "last"))
  simulation <- match.arg(simulation, c("SIM1", "SIM2"))
  calendar <- if (.calendar_flag(calendar)==1L) "G" else "J"
  if (!is.numeric(year) || !length(year) || any(!is.finite(year)) ||
      any(year != trunc(year))) stop("year must contain finite integer years.")
  checkYear(year)
  if (length(year)==2L) {
    if (year[2] < year[1]) stop("Year endpoints must be ascending.")
    year <- seq.int(year[1], year[2])
  }
  geo <- if (inherits(loc, "skyscapeR.horizon")) loc$metadata$georef else loc
  if (!is.numeric(geo) || length(geo)!=3L || any(!is.finite(geo)) ||
      abs(geo[1])>=90 || abs(geo[2])>180)
    stop("loc must supply latitude strictly between -90 and 90, longitude in [-180,180], and elevation.")
  if (!is.logical(refraction) || length(refraction)!=1L || is.na(refraction))
    stop("refraction must be TRUE or FALSE.")
  if (!is.numeric(atm) || length(atm)!=1L || !is.finite(atm) || atm<=0 ||
      !is.numeric(temp) || length(temp)!=1L || !is.finite(temp) || temp<= -273.15)
    stop("Supply positive finite pressure and finite temperature above absolute zero.")
  # Validate timezone/calendar using the package's date parser before searching.
  time2jd(paste0(year[1], "/01/01"), timezone, calendar)
  swephR::swe_set_topo(geo[2], geo[1], geo[3])
  direction <- if (crescent=="first") 1 else -1
  target <- if (season=="winter") 1 else -1
  rows <- vector("list", length(year))
  for (k in seq_along(year)) {
    jan <- swephR::swe_julday(year[k], 1, 1, 0, 1L)
    # Locate the actual solstice, not a fixed calendar date.
    centre <- jan + if (season=="summer") 172 else 355
    longitude <- if (season=="summer") 90 else 270
    solstice <- stats::uniroot(function(t) .cc_wrap(.cc_calc(t,0,2)[1]-longitude),
                              c(centre-20, centre+20), tol=1e-8)$root
    conjunctions <- .cc_newmoons(solstice-180, solstice+150)
    start <- tail(conjunctions[conjunctions<=solstice],3)[1]
    candidates <- list(); previous <- NULL
    for (conjunction in conjunctions) {
      base <- .cc_crescent(conjunction, direction, loc, geo, refraction,
                           atm, temp, timezone, calendar)
      if (is.null(base)) { previous <- NULL; next }
      choices <- list(base)
      if (simulation=="SIM2") {
        extra <- .cc_observation(base$day + direction, conjunction, direction,
                                 loc, geo, refraction, atm, temp)
        if (!is.null(extra)) { extra$offset <- 1L; choices[[2]] <- extra }
      }
      for (candidate in choices) {
        if (conjunction>=start && target*candidate$diff > 0) {
          candidate$previous <- if (is.null(previous)) NA_real_ else previous$diff
          candidates[[length(candidates)+1L]] <- candidate
          break
        }
      }
      # The final eligible observation of this lunation establishes the next comparison.
      previous <- choices[[length(choices)]]
    }
    status <- if (!length(candidates)) "not_found" else "ok"
    selected <- if (status=="ok") candidates[[1]] else NULL
    rows[[k]] <- .cc_result(year[k], solstice, selected, status, geo, timezone, calendar)
    rows[[k]]$diagnostics$qualifying.lunations <- length(candidates)
    if (status!="ok") warning(sprintf("CC %s %s %s: %s (%d eligible reversals).",
                           season,crescent,year[k],status,length(candidates)),call.=FALSE)
  }
  .with_astronomy_provenance(list(Event=paste(season,crescent,"crescent crossover"),
       Moon=do.call(rbind,lapply(rows, `[[`, "Moon")),
       Sun=do.call(rbind,lapply(rows, `[[`, "Sun")),
       diagnostics=do.call(rbind,lapply(rows, `[[`, "diagnostics")),
       settings=list(season=season,crescent=crescent,simulation=simulation,
                     refraction=refraction,atm=atm,temp=temp,timezone=timezone,
                     calendar=calendar,declination="geocentric",limb="centre",
                     visibility="Caldwell-Laney at sunset/sunrise",georef=geo)))
}

.cc_wrap <- function(x) (x+180) %% 360 - 180
.cc_calc <- function(jd, body, flags) {
  ans <- .file_position(jd, body, flags)
  if (ans$return<0) stop(ans$serr)
  ans$xx
}
.cc_newmoons <- function(start, end) {
  phase <- function(t) .cc_wrap(.cc_calc(t,1,2)[1]-.cc_calc(t,0,2)[1])
  days <- seq(start,end,by=1)
  values <- vapply(days,phase,numeric(1))
  i <- which(head(values,-1)<=0 & tail(values,-1)>0 & diff(values)<180)
  vapply(i,function(j) stats::uniroot(phase,days[c(j,j+1)],tol=1e-8)$root,numeric(1))
}
.cc_position <- function(jd, body, geo, refraction, atm, temp) {
  eq <- .cc_calc(jd,body,2L+2048L+32768L)
  hor <- swephR::swe_azalt(jd,1,c(geo[2],geo[1],geo[3]),atm,temp,eq[1:3])$xaz
  list(eq=eq,az=(hor[1]+180)%%360,alt=hor[if (refraction) 3 else 2])
}
.cc_event <- function(day, body, direction, loc, geo, refraction, atm, temp) {
  which.event <- if (direction==1) 2L else 1L
  altitude <- 0
  if (inherits(loc,"skyscapeR.horizon")) {
    if (length(unique(loc$data$alt))==1L) altitude <- loc$data$alt[1] else
    return(.horizon_events(if (body==0) "sun" else "moon",day,loc,
                           refraction,atm,temp)[which.event])
  }
  ans <- .file_crossing(day,body,"",2L,
             which.event+256L+if (refraction) 0L else 512L,
             c(geo[2],geo[1],geo[3]),atm,temp,altitude)
  if (ans$return == -2L) return(NA_real_)
  if (ans$return < 0L) stop(ans$serr)
  ans$tret[1]
}
.cc_observation <- function(day, conjunction, direction, loc, geo, refraction, atm, temp) {
  sunjd <- .cc_event(day,0,direction,loc,geo,refraction,atm,temp)
  if (!is.finite(sunjd) || sunjd < day || sunjd >= day+1.01 ||
      direction*(sunjd-conjunction)<=0 || abs(sunjd-conjunction)>9) return(NULL)
  sun <- .cc_position(sunjd,0,geo,refraction,atm,temp)
  moon <- .cc_position(sunjd,1,geo,refraction,atm,temp)
  rad <- pi/180
  sun.geo <- .cc_calc(sunjd,0,2L+2048L)
  moon.geo <- .cc_calc(sunjd,1,2L+2048L)
  cossep <- sin(sun.geo[2]*rad)*sin(moon.geo[2]*rad)+
    cos(sun.geo[2]*rad)*cos(moon.geo[2]*rad)*cos((sun.geo[1]-moon.geo[1])*rad)
  elongation <- acos(max(-1,min(1,cossep)))/rad
  radius <- asin(1737.4/(moon$eq[3]*149597870.7))/rad
  lower <- moon$alt-radius
  margin <- lower+elongation/3-11.3
  horizon <- if (inherits(loc,"skyscapeR.horizon")) stats::approx(loc$data$az,loc$data$alt,xout=moon$az)$y else 0
  if (elongation>=90 || margin<0 || lower<=horizon) return(NULL)
  moonjd <- .cc_event(if (direction==1) sunjd+1e-7 else sunjd-.5,
                      1,direction,loc,geo,refraction,atm,temp)
  if (!is.finite(moonjd) || direction*(moonjd-sunjd)<=0 ||
      abs(moonjd-sunjd)>.5) return(NULL)
  moon.event <- .cc_position(moonjd,1,geo,refraction,atm,temp)
  list(day=day,assessment=sunjd,sunjd=sunjd,moonjd=moonjd,sunaz=sun$az,
       moonaz=moon.event$az,diff=.cc_wrap(moon.event$az-sun$az),
       lower=lower,elongation=elongation,margin=margin,
       illumination=moonphase(sunjd),conjunction=conjunction,offset=0L)
}
.cc_crescent <- function(conjunction, direction, loc, geo, refraction, atm, temp, timezone, calendar) {
  stamp <- strsplit(jd2time(conjunction,timezone,calendar)," ",fixed=TRUE)[[1]][1]
  midnight <- time2jd(stamp,timezone,calendar)
  for (offset in 0:8) {
    ans <- .cc_observation(midnight+direction*offset,conjunction,direction,
                            loc,geo,refraction,atm,temp)
    if (!is.null(ans)) return(ans)
  }
  NULL
}
.cc_result <- function(year, solstice, x, status, geo, timezone, calendar) {
  event <- function(body) {
    jd <- if (is.null(x)) NA_real_ else if (body==0) x$sunjd else x$moonjd
    stamp <- if (is.na(jd)) c(NA_character_,NA_character_) else
      strsplit(jd2time(jd,timezone,calendar)," ",fixed=TRUE)[[1]]
    data.frame(year=year,date=stamp[1],time=stamp[2],
      azimuth=if(is.null(x)) NA_real_ else if(body==0) x$sunaz else x$moonaz,
      declination=if(is.null(x)) NA_real_ else .cc_calc(jd,body,2L+2048L)[2],jd=jd)
  }
  get <- function(name) if (is.null(x)) NA_real_ else x[[name]]
  list(Moon=event(1),Sun=event(0),diagnostics=data.frame(year=year,status=status,
       solstice.jd=solstice,conjunction.jd=get("conjunction"),
       assessment.jd=get("assessment"),illumination=get("illumination"),
       lower.limb.altitude=get("lower"),elongation=get("elongation"),
       visibility.margin=get("margin"),azimuth.difference=get("diff"),
       previous.difference=get("previous"),day.offset=get("offset")))
}
