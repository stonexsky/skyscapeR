# Actual transitions on the closed model-year cycle. Day keeps the existing
# convention: last day before entry, first day after exit from I or V.
.phase_transitions <- function(phase, day=seq_along(phase)-1L) {
  n <- length(phase)
  previous <- c(n,seq_len(n-1L))
  rows <- which(phase != phase[previous])
  out <- data.frame(index=integer(), day=numeric(), event=character())
  enter <- c(I="acronycal setting", V="acronycal rising")
  leave <- c(I="heliacal rising", V="heliacal setting")
  for (i in rows) {
    if (phase[previous[i]] %in% names(leave))
      out <- rbind(out,data.frame(index=i,day=day[i],event=unname(leave[phase[previous[i]]])))
    if (phase[i] %in% names(enter))
      out <- rbind(out,data.frame(index=i,day=day[previous[i]],event=unname(enter[phase[i]])))
  }
  out
}

#' Calculate a fixed star's coordinates for an epoch
#' @description Find a catalogue object's equatorial position at an archaeological epoch and its
#' apparent magnitude for a sky or visibility model.
#' @details A historical comparison requires a stellar coordinate for the epoch under study. The
#' function resolves a Western common name through star.names or accepts a catalogue
#' identifier and checks it against Swiss Ephemeris. It evaluates equatorial coordinates at
#' noon UT on 1 January of each supplied astronomical year in the Gregorian calendar. The
#' magnitude is a catalogue value. These outputs can feed horizon transformations or
#' visibility estimates, which additionally require observer and atmosphere information.
#'
#' The Swiss Ephemeris manual documents fixed-star positions and catalogue magnitudes.
#'
#' Astronomical positions require Swiss Ephemeris files covering the epoch. We
#' recommend installing \pkg{swephRdata}; see \code{\link{skyscapeR-package}} for setup.
#' @return A skyscapeR.star list with name, Bayer designation, apparent magnitude (app.mag) and a coord data frame. Each row contains year, right ascension (RA) and declination (Dec), both angles in degrees.
#' @param string This can be either the traditional name for the star or its Bayer designation.
#' @param year Year for which to calculate the coordinates. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' Defaults to current year.
#' @import swephR
#' @export
#' @references Astrodienst. Swiss Ephemeris programmer's manual. \url{https://www.astro.com/swisseph/swephprg.htm}.
#' @seealso \code{\link[swephR]{swe_fixstar2_ut}}, \code{\link[swephR]{swe_fixstar2_mag}}, \code{\link{star.names}}, \code{\link{star.phases}}
#' @examples
#' # M45 and M 45 identify the Pleiades.
#' if (requireNamespace("swephRdata", quietly=TRUE)) star("M45", -4300)
#'
#' if (requireNamespace("swephRdata", quietly = TRUE)) {
#'   # Retrieve data for Aldebaran:
#'   Aldeb <- star('Aldebaran')
#'   
#'   # Retrieve data for Aldebaran on 3001 BC:
#'   ss <- star('Aldebaran', -3000)
#' }
#' @section Provenance:
#' Structured results carry a \code{provenance} attribute. Read it with
#' \code{attr(result, "provenance")}; its \code{packages} entry records the loaded
#' skyscapeR and swephR versions, with NA when a namespace version is unavailable.
#' This attribute adds no list entries or table columns. Ordinary subsetting or
#' transformations may discard attributes; save the original result for reference.
#' It does not identify ephemeris files, catalogue revisions or a development
#' commit. Retain those and the analysis script separately for exact reproduction.
#' @section Scientific checks:
#' \itemize{
#' \item For six stars compared with independently propagated catalogue directions, the largest separation was 0.0165 degrees, for Sirius in 3500 BCE. \href{https://stonexsky.github.io/skyscapeR/manual/a03-stellar-positions.html}{A03 study}.
#' \item The Sirius target screen recovered two disjoint intervals in the test band and a single contact at the annual maximum in 36 CE. \href{https://stonexsky.github.io/skyscapeR/manual/a08-target-ranges.html}{A08 study}.
#' \item First-morning dates fell within seven days of 22 of 23 sightings and last-evening dates within seven days of 21 of 33. The respective median predictions were two days early and five days late. \href{https://stonexsky.github.io/skyscapeR/manual/heliacal-rising.html#v01-observations}{V01 study}.
#' }
star <- function(string, year=skyscapeR.env$cur.year) {
  if (is.character(year)) year <- astroYear(year)
  star.names <- skyscapeR::star.names
  aux <- star.names$identifier[which(star.names$Western==string)]
  if (length(aux)==1) {
    name <- string
    string <-  aux
  } else {
    key <- function(x) tolower(gsub("[[:space:]]+", "", x))
    i <- match(key(string), key(star.names$identifier))
    name <- if (is.na(i)) "" else star.names$Western[i]
    if (!is.na(i)) string <- star.names$identifier[i]
  }

  string <- paste0(',',string)
  jd <- vapply(year, swephR::swe_julday, numeric(1), month=1L, day=1L, hourd=12, gregflag=1L)
  .file_position(jd, 0L, 2L, label=paste("fixed star", substring(string,2L)))
  info <- swephR::swe_fixstar2_ut(string, jd, 2048+16384+16+2L)
  if (any(info$return < 0L)) stop(paste(unique(info$serr),collapse=" "))
  bad <- which(bitwAnd(info$return,7L)!=2L | grepl("moshier",info$serr,ignore.case=TRUE))
  if (length(bad)) .ephemeris_failure(jd[bad[1]],string,info$serr[bad[1]])
  coords <- matrix(info$xx,nrow=length(year))
  aux <- data.frame(year=year, RA=coords[,1], Dec=coords[,2])

  if (length(name)==0) { name <- ''}

  star <- c()
  star$name <- name
  star$Bayer <- strsplit(info$starname[1],',')[[1]][2]
  star$app.mag <- swephR::swe_fixstar2_mag(paste(string))$mag
  star$coord <- aux
  class(star) <- "skyscapeR.star"

  return(.with_astronomy_provenance(star))
}


#' @noRd
ff <- function(x, loc, atm, temp, pos) { return(swephR::swe_azalt(x, 1, c(loc[2],loc[1],loc[3]), atm, temp, c(pos[1], pos[2]))$xaz) }

#' @noRd
# Under the fixed-coordinate stellar model every model day uses the same track.
# Resolve the terrain windows once, including separate emergence/reoccultation
# pairs. Twilight and visibility thresholds can then vary by solar longitude.
.phase_track <- function(jd0,star,loc,alt.hor,k,limit,res,refraction,atm,temp,horizon) {
  if(is.null(horizon)) horizon <- .horizon_profile(c(0,360),rep(alt.hor,2))
  eq <- as.numeric(star$coord[2:3])
  component <- if(refraction) 3L else 2L
  position <- function(t) {
    p <- ff(t,loc,atm,temp,eq)
    c((p[1]+180)%%360,p[component])
  }
  # Terrain can cross the stellar track only while the altitude is inside
  # the profile's range. Bound those visits with the smooth fixed-star arc.
  limits <- attr(horizon,"range")
  level.roots <- function(level) .horizon_crossings(position,jd0,c(0,360),c(level,level),
                                                   .step=1/24,.duration=1,.all=TRUE)
  if(is.null(limits)) {
    crossings <- .horizon_crossings(position,jd0,.duration=1,.all=TRUE,horizon=horizon)
  } else if(limits[1]==limits[2]) {
    crossings <- level.roots(limits[1])
  } else {
    bounds <- sort(unique(c(jd0,level.roots(limits[1]),level.roots(limits[2]),jd0+1)))
    crossings <- bounds
    for(i in seq_len(length(bounds)-1L)) {
      a <- bounds[i]; b <- bounds[i+1L]
      h <- position((a+b)/2)[2]
      if(h>=limits[1] && h<=limits[2]) {
        roots <- .horizon_crossings(position,a,.duration=b-a,.all=TRUE,horizon=horizon)
        crossings <- c(crossings,roots)
      }
    }
  }
  solar.grid <- sort(unique(c(seq(jd0,jd0+1,by=res),jd0+1)))
  times <- sort(unique(c(seq(jd0,jd0+1,by=1/1440),solar.grid,crossings,jd0+1)))
  height <- function(t) position(t)[2]
  alt <- vapply(times,height,numeric(1))
  slope <- diff(alt)
  turns <- which(head(slope,-1)*tail(slope,-1)<0)+1L
  extrema <- vapply(turns,function(i) {
    a <- times[i-1L]
    opt <- stats::optimize(function(seconds) height(a+seconds/86400),
      c(0,(times[i+1L]-a)*86400),maximum=slope[i-1L]>0,tol=.001)
    a+opt[[1L]]/86400
  },numeric(1))
  if(length(extrema)) {
    times <- sort(unique(c(times,extrema)))
    alt <- vapply(times,height,numeric(1))
  }
  mid <- head(times,-1)+diff(times)/2
  open <- vapply(mid,function(t) {
    p <- position(t)
    p[2]>horizon(p[1])
  },logical(1))
  attenuation <- .stellar_extinction(alt,k)
  list(time=times,solar.grid=solar.grid,alt=alt,open=open,rising=diff(alt)>0,
       setting=diff(alt)<0,extinction=attenuation,
       magnitude.margin=limit-star$app.mag-attenuation)
}

.stellar_extinction <- function(altitude, k) {
  sine <- sin(altitude*pi/180)
  k/(sine + .025*exp(-11*sine))
}


# Clip fractional positions on each segment against one linear inequality.
.phase_clip <- function(bounds,a,b,strict=FALSE) {
  invalid <- !is.finite(a) | !is.finite(b)
  impossible <- if(strict) a<=0 & b<=0 else a<0 & b<0
  bounds$right[invalid | impossible] <- 0
  up <- which(!invalid & a<0 & b>=0)
  down <- which(!invalid & a>=0 & b<0)
  bounds$left[up] <- pmax(bounds$left[up],a[up]/(a[up]-b[up]))
  bounds$right[down] <- pmin(bounds$right[down],a[down]/(a[down]-b[down]))
  bounds
}

# Precompute the eligible intervals and solar interpolation weights once.
.phase_segments <- function(track,alt.rs=Inf) {
  bounds <- list(left=numeric(length(track$open)),right=as.numeric(track$open))
  for(margin in list(track$extinction,track$magnitude.margin))
    bounds <- .phase_clip(bounds,head(margin,-1L),tail(margin,-1L))
  if(is.finite(alt.rs))
    bounds <- .phase_clip(bounds,alt.rs-head(track$alt,-1L),alt.rs-tail(track$alt,-1L))
  i <- which(bounds$right>bounds$left)
  j <- findInterval(track$time[i],track$solar.grid,all.inside=TRUE)
  nodes <- sort(unique(c(j,j+1L)))
  dt <- diff(track$time)[i]
  ds <- diff(track$solar.grid)[j]
  offset <- track$time[i]-track$solar.grid[j]
  delta <- diff(track$alt)[i]
  list(nodes=nodes,j=match(j,nodes),next.j=match(j+1L,nodes),
       w0=(offset+bounds$left[i]*dt)/ds,w1=(offset+bounds$right[i]*dt)/ds,
       alt0=track$alt[i]+bounds$left[i]*delta,alt1=track$alt[i]+bounds$right[i]*delta,
       rising=track$rising[i],setting=track$setting[i])
}

.phase_visible <- function(segments,sun,solar.alt,arcus) {
  base <- sun[segments$j]
  delta <- sun[segments$next.j]-base
  s0 <- base+segments$w0*delta
  s1 <- base+segments$w1*delta
  bounds <- list(left=numeric(length(base)),right=rep(1,length(base)))
  bounds <- .phase_clip(bounds,solar.alt-s0,solar.alt-s1,TRUE)
  bounds <- .phase_clip(bounds,segments$alt0-s0-arcus,segments$alt1-s1-arcus)
  bounds$right>bounds$left
}

calc.phase <- function(year, day, sun.coord, star, loc, arcus_visionis, alt.hor, k, limit, alt.rs, res, refraction, atm, temp, pb, pbi, horizon=NULL, solar.alt=alt.hor-6, track=NULL) {
  jd0 <- time2jd('2000/01/01 00:00:00',timezone='UTC',calendar='G')
  if (is.null(track)) track <- .phase_track(jd0,star,loc,alt.hor,k,limit,res,refraction,atm,temp,horizon)
  horizon.status <- if (!any(track$open)) "always below" else
    if (all(track$open)) "always above" else "crosses"
  # Geometric clearance and visual detectability are separate properties.
  # Return an ordinary invisible sequence for a completely obstructed track.
  if (horizon.status == "always below")
    return(structure(rep("I", length(day)), horizon.status=horizon.status))
  low <- .phase_segments(track,alt.rs)
  any.altitude <- .phase_segments(track)
  component <- if(refraction) 3L else 2L
  phase <- vapply(seq_along(day),function(i) {
    cache <- rep(NA_real_,length(track$solar.grid))
    coord <- sun.coord[day[i],]
    visible <- function(segments) {
      nodes <- segments$nodes
      missing <- nodes[is.na(cache[nodes])]
      cache[missing] <<- vapply(track$solar.grid[missing],function(t)
        ff(t,loc,atm,temp,coord)[component],numeric(1))
      .phase_visible(segments,cache[nodes],solar.alt,arcus_visionis)
    }
    seen <- visible(low)
    rise <- any(seen & low$rising)
    set <- any(seen & low$setting)
    if (!is.null(pb)) setTxtProgressBar(pb,pbi+i)
    if(rise && set) "RS" else if(rise) "R" else if(set) "S" else
      if(any(seen) || any(visible(any.altitude))) "V" else "I"
  },character(1))
  structure(phase, horizon.status=horizon.status)
}

#' Model seasonal changes in stellar visibility
#' @description Estimate seasonal opportunities for naked-eye stellar visibility.
#' @details The default uses Schaefer's directional twilight brightness and Hecht's naked-eye point-source threshold for heliacal rising (first morning visibility) and acronycal setting (last evening visibility). On each local mean-solar date the calculation searches the terrain-clear rising branch before noon and descending branch after noon, with the Sun below the astronomical horizon. A transition marks the first qualifying morning following an unsuccessful morning, or the last qualifying evening preceding an unsuccessful evening. Solar altitude, stellar zenith distance and the difference in azimuth determine the twilight background. The limiting magnitude sets the dark zenith reference. Rozenberg airmass and the extinction coefficient attenuate the star. Stellar clearance and extinction use the selected apparent or geometric altitude. Solar altitude in the twilight equation is geometric.
#' 
#' The passage bands, phase type, heliacal setting and acronycal rising use Purrington's altitude-difference rule on a 365-sample solar-longitude cycle. The solar.alt cutoff and alt.rs passage ceiling apply to this classification. Therefore a Schaefer event marker may lie inside a passage band. Selecting the Purrington method also derives first morning and last evening from the passage transitions. A horizon object must cover all bearings. With coordinates, alt.hor supplies the constant skyline.
#' 
#' The Schaefer fit was developed for twilight near the horizon. Searches at high latitudes or far from the solar azimuth can extend outside its empirical range. The calculation assumes a clear sky and a fixed observer limit. Moonlight, weather and observer adaptation can affect field sightings. Treat an integrated cluster magnitude as a single unresolved target.
#' @return A skyscapeR.starphases object. The data and seasons tables describe the passage classification. In metadata, events records the event, schematic day, solar.longitude, local date.julian, method and numeric model.day. Use date.julian with \code{\link{riseset}} and calendar = "J" to locate a crossing. In Schaefer mode visibility.calendar contains morning and evening flags for each local date in the requested Gregorian astronomical year. local.jd is the Julian day number at 00:00 UT for that date label. The horizon.status field describes clearance from the supplied skyline.
#' @param star Either the star name or a \emph{skyscapeR.star} object.
#' @param year The year of interest. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' @param loc Location, either a \emph{skyscapeR.horizon} or a vector
#' containing the latitude, longitude and elevation of location, in this order.
#' @param alt.hor The altitude of the horizon to consider.
#' Defaults to zero degrees.
#' @param k Visual extinction in magnitudes per airmass (Schaefer 1986, p. S33). The package adopts 0.2 as a clear-sky baseline. Extinction varies with site and atmospheric conditions (Schaefer 1986, 1993). Defaults to the session value from skyscapeR.vars().
#' @param limit The maximum magnitude of a star that can be
#' visible with the naked eye. Defaults to 6.
#' @param alt.rs The maximum altitude of a star's first or last
#' visibility for it to still be considered to be as rising or setting.
#' Defaults to ten degrees.
#' @param res Spacing of solar altitude samples in days, greater than zero and at
#' most one. Default 1/144 day (ten minutes). Terrain crossings are refined separately.
#' @param refraction Logical. Apply atmospheric refraction. Defaults to the session setting in \code{\link{skyscapeR.vars}}.
#' @param atm Atmospheric pressure in mbar. Defaults to the session setting in \code{\link{skyscapeR.vars}}.
#' @param temp Air temperature in degrees Celsius. Defaults to the session setting in \code{\link{skyscapeR.vars}}.
#' @param solar.alt Solar altitude threshold in degrees. NULL uses
#' -6 for a horizon, or alt.hor - 6 for coordinate input.
#' @param method Visibility model for first morning and last evening. Defaults to Schaefer. Select "purrington" for the passage-transition rule.
#' @section Accuracy of first and last visibility dates:
#' Across 56 recorded dates, heliacal rising fell within seven days of 22 of 23 first morning sightings. The median prediction was two days early. Acronycal setting fell within seven days of 21 of 33 last evening sightings. The median prediction was five days late. A week is a practical search window. Larger discrepancies occur, especially at dusk. These results apply to the two tested events. \href{https://stonexsky.github.io/skyscapeR/manual/heliacal-rising.html#v01-observations}{V01 observed dates and scientific checks}.
#' @section Workflow checks:
#' \href{https://stonexsky.github.io/skyscapeR/manual/w02-field-to-sky.html}{W02 field-to-sky replay} checks the Gacrux event at Skorba West. \href{https://stonexsky.github.io/skyscapeR/manual/h07-panorama-coverage.html}{H07 skyline coverage} checks terrain inputs. \href{https://stonexsky.github.io/skyscapeR/manual/a03-stellar-positions.html}{A03 stellar positions} checks catalogue propagation.
#' @export
#' @references Crumey, A. (2014). Human contrast threshold and astronomical
#' visibility. MNRAS 442, 2600-2619. \doi{10.1093/mnras/stu992}.
#' @references Purrington, R. D. (1988). Heliacal Rising and Setting: Quantitative Aspects. Journal for the History of Astronomy 19, S72-S84. \url{https://articles.adsabs.harvard.edu/pdf/1988JHAS...19...72P}.
#' @references Schaefer, B. E. (1986). Atmospheric Extinction Effects on Stellar Alignments. Journal for the History of Astronomy 17, S32-S42. \doi{10.1177/002182868601701003}.
#' @references Schaefer, B. E. (1993). Astronomy and the Limits of Vision. Vistas in Astronomy 36, 311-361. \doi{10.1016/0083-6656(93)90113-X}.
#' @references Schaefer, B. E. (1987). Heliacal Rise Phenomena. Journal for the History of Astronomy 18, S19-S33. \doi{10.1177/002182868701801103}.
#' @references Brady, B. (2015). Star Phases: the Naked-eye Astronomy of the Old Kingdom Pyramid Texts. In F. Silva and N. Campion (eds), Skyscapes: The Role and Importance of the Sky in Archaeology. Oxbow Books, pp. 76-86.
#' @seealso \code{\link{plot.skyscapeR.starphases}}, \code{\link{riseset}}, \code{\link{sketchSky}}, \code{\link{skyStep}}, \code{\link{star}}
#' @examples
#' \dontrun{
#' if (requireNamespace("swephRdata", quietly = TRUE)) {
#'   ss1 <- star.phases('Aldebaran',-4000, c(35,-8,200))
#'   
#'   # One can then look at the star's phase type:
#'   ss1$metadata$type
#'   
#'   # Date range of seasons:
#'   ss1$metadata$seasons
#'   
#'   # Date range of phase-type events:
#'   ss1$metadata$events
#'   
#'   # And plot them:
#'   plot(ss1)
#'   
#'   # You can play with the parameters and see how predictions change:
#'   ss1 <- star.phases('Aldebaran',-4000, c(35,-8,200), alt.hor=2, alt.rs=5)
#'   plot(ss1)
#' }
#' }
#' @section Scientific checks:
#' \itemize{
#' \item A partial photographic horizon returned missing altitudes outside its sector. Full-sky event searches rejected it and accepted the complete profile. \href{https://stonexsky.github.io/skyscapeR/manual/h07-panorama-coverage.html}{H07 study}.
#' \item The zero-declination track crossed a level skyline at 90 and 270 degrees at 50 degrees north. The study also checks irregular skylines and JPL crossings. \href{https://stonexsky.github.io/skyscapeR/manual/a04-stellar-crossings.html}{A04 study}.
#' \item First-morning dates fell within seven days of 22 of 23 sightings and last-evening dates within seven days of 21 of 33. The respective median predictions were two days early and five days late. \href{https://stonexsky.github.io/skyscapeR/manual/heliacal-rising.html#v01-observations}{V01 study}.
#' \item With extinction of 0.2 magnitudes per airmass, transmission was 0.149 at 5 degrees and 0.832 at the zenith under the adopted law. \href{https://stonexsky.github.io/skyscapeR/manual/sky-rendering.html#extinction-check}{V05 study}.
#' \item The field-to-sky replay recovered the Gacrux position between the Skorba West entrance edges, 2.41 degrees above the measured skyline at the illustrated instant. \href{https://stonexsky.github.io/skyscapeR/manual/w02-field-to-sky.html}{W02 study}.
#' }
star.phases <- function(star, year, loc, alt.hor = 0, alt.rs = 10, solar.alt=NULL, k = skyscapeR.vars()$k, limit =
                        6, refraction, atm, temp, res = 1/24/6, method = c("schaefer", "purrington")) {
  method <- match.arg(method)
  if (length(k)!=1L || !is.finite(k) || k<0 || length(limit)!=1L || !is.finite(limit))
    stop("Supply a finite limiting magnitude and nonnegative extinction coefficient.")
  if (is.character(year)) year <- astroYear(year)
  if (missing(refraction)) { refraction <- skyscapeR.env$refraction }
  if (missing(atm)) { atm <- skyscapeR.env$atm }
  if (missing(temp)) { temp <- skyscapeR.env$temp }

  site.name <- if (is(loc,'skyscapeR.horizon')) loc$metadata$name else NULL
  if (is(loc,'skyscapeR.horizon')) {
    lat <- loc$metadata$georef[1]
    lon <- loc$metadata$georef[2]
    elev <- loc$metadata$georef[3]
  } else {
    lat <- loc[1]
    lon <- loc[2]
    elev <- loc[3]
  }
  horizon <- NULL
  if(is(loc,'skyscapeR.horizon')) {
    if(!missing(alt.hor)) stop('Use either a horizon profile or alt.hor, not both.')
    .require_full_skyline(loc)
    az <- loc$data$az; heights <- loc$data$alt
    horizon <- .horizon_profile(az,heights)
  }
  loc <- c(lat,lon,elev)
  if(length(loc)!=3L || any(!is.finite(loc))) stop('Supply latitude, longitude and elevation.')
  if(!is.logical(refraction) || length(refraction)!=1L || is.na(refraction)) stop('refraction must be TRUE or FALSE.')
  if(is.null(solar.alt)) solar.alt <- if(is.null(horizon)) alt.hor-6 else -6
  if(!is.numeric(solar.alt) || length(solar.alt)!=1L || !is.finite(solar.alt)) stop('solar.alt must be finite and scalar.')
  swephR::swe_set_topo(lon, lat, elev)

  # load star data
  if (!is(star, 'skyscapeR.star')) { star <- star(star, year) }
  arcus_visionis <- 2.1*star$app.mag + 10

  # pre-calculations of solar position
  ob <- obliquity(year)
  ll <- seq(0,359.9,0.1)
  Dec <- asin(sin(ob/180*pi)*sin(ll/180*pi))/pi*180
  RA <- atan2(sin(ll/180*pi)*cos(ob/180*pi), cos(ll/180*pi))/pi*180
  RA[RA<0] <- RA[RA<0]+360
  lon <- seq(0,359.99,by=360/365)
  dd <- approx(seq(0,359.9,0.1), Dec, xout=lon)$y
  rr <- approx(seq(0,359.9,0.1), RA, xout=lon)$y
  sun.pos <- cbind(rr,dd)

  if(!is.numeric(res) || length(res)!=1L || !is.finite(res) || res<=0 || res>1)
    stop("res must be a positive fraction of a day no greater than one.")

  pb <- NULL
  if (.sky_progress(TRUE)) {
    pb <- txtProgressBar(max=365,style=3)
    on.exit(close(pb),add=TRUE)
  }
  track <- .phase_track(time2jd("2000/01/01 00:00:00","UTC","G"),star,loc,alt.hor,k,limit,res,refraction,atm,temp,horizon)
  phase <- calc.phase(year,seq_len(365),sun.pos,star,loc,arcus_visionis,
    alt.hor,k,limit,alt.rs,res,refraction,atm,temp,pb,0,horizon=horizon,solar.alt=solar.alt,track=track)
    horizon.status <- attr(phase, "horizon.status")
    phase <- as.vector(phase)
    # find december solstice date and pin calendar to it
    days <- seq(1,365)-which.min(dd)
    days[days<0] <- days[days<0]+365

    ind <- sort(days, index.return=T)$ix
    days <- days[ind]; phase <- phase[ind]

    type <- 'unclassified'
    uu <- unique(phase)
    if (length(uu)==1) {
      if (uu=='V') { type <- if (horizon.status == 'always above') 'circumpolar' else 'visible throughout cycle' }
      if (uu=='I') { type <- 'invisible' }
    } else {
      if ('V' %in% uu) { type <- 'curtailed passage' }
      if ('I' %in% uu) { type <- 'arising and lying hidden' }
      if ('I' %in% uu & 'V' %in% uu) { type <- 'dual phase'}
    }

    # Summarise contiguous seasons on the circular model year. Terrain can
    # produce more than one separated run of the same code.
    starts <- which(phase!=c(tail(phase,1),head(phase,-1)))
    if (!length(starts)) starts <- 1L
    following <- c(starts[-1],starts[1]+length(phase))
    ends <- (following-2L)%%length(phase)+1L
    labels <- c(R="rise only",S="set only",RS="rise and set",
                I="arising and lying hidden",V="curtailed passage")
    seasons <- data.frame(season=unname(labels[phase[starts]]),
                          begin=days[starts],end=days[ends],length=following-starts)
    transitions <- .phase_transitions(phase, days)
    events <- transitions[,c("event","day"),drop=FALSE]
    events <- events[order(events$day),,drop=FALSE]
    rownames(events) <- NULL
    seasons$begin <- dd.to.DD(seasons$begin, char=TRUE, WS=TRUE)
    seasons$end <- dd.to.DD(seasons$end, char=TRUE, WS=TRUE)
    events$solar.longitude <- lon[ind][match(events$day, days)]
    events$date.julian <- .phase_event_dates(events$solar.longitude,year,loc[2])
    events$method <- rep("purrington", nrow(events))
    calendar <- NULL
    if (method == "schaefer") {
      calendar <- .schaefer_calendar(year,star,loc,track,k,limit,atm,temp,refraction)
      n <- nrow(calendar)
      inside <- 2L:(n-1L)
      rise <- inside[calendar$morning[inside] & !calendar$morning[inside-1L]]
      set <- inside[calendar$evening[inside] & !calendar$evening[inside+1L]]
      selected <- c(rise,set)
      labels <- c(rep("heliacal rising",length(rise)),rep("acronycal setting",length(set)))
      events <- events[!events$event %in% c("heliacal rising","acronycal setting"),,drop=FALSE]
      if (length(selected)) {
        local <- calendar$local.jd[selected]
        longitude <- .solar_coordinate(local+.5-loc[2]/360,0L,1L)
        grid.day <- vapply(longitude,function(x) days[which.min(abs((lon[ind]-x+180)%%360-180))],numeric(1))
        dates <- vapply(local,function(j) sub(" .*","",jd2time(j+.5,"UTC","J")),character(1))
        events <- rbind(events,data.frame(event=labels,day=grid.day,solar.longitude=longitude,
                                          date.julian=dates,method="schaefer"))
      }
      calendar <- calendar[inside,,drop=FALSE]
    }
    events <- events[order(events$day),,drop=FALSE]
    events$model.day <- events$day
    events$day <- if(nrow(events)) dd.to.DD(events$day,char=TRUE,WS=TRUE) else character()

    # results
    out <- c()
    out$metadata$star <- star
    out$metadata$year <- year
    out$metadata$loc <- loc
    out$metadata$type <- type
    out$metadata$skyline <- !is.null(horizon)  # Retain the original field and its position.
    out$metadata$solar.alt <- solar.alt
    out$metadata$visibility.search <- "terrain-windows-v1"
    out$metadata$refraction <- refraction
    out$metadata$atm <- atm
    out$metadata$temp <- temp
    out$metadata$method <- method
    out$metadata$visibility.calendar <- calendar
    out$metadata$events <- events
    out$metadata$seasons <- seasons
    out$metadata$k <- k
    out$metadata$horizon <- !is.null(horizon)
    out$metadata$horizon.status <- horizon.status
    out$metadata$limit <- limit
    out$metadata$alt.rs <- alt.rs
    out$metadata$res <- res
    out$metadata$alt.hor <- if (is.null(horizon)) alt.hor else NA_real_
    out$metadata["horizon.profile"] <- list(if (is.null(horizon)) NULL else
      data.frame(az=as.numeric(az), alt=as.numeric(heights)))
    out$metadata$provenance <- .astronomy_provenance()
    out$metadata$site.name <- site.name
    out$data <- data.frame(day = days, phase = phase)
    class(out) <- 'skyscapeR.starphases'
    return(.with_astronomy_provenance(out))
}






# Map seasonal longitude onto the local mean-solar civil day of this epoch.
.phase_event_dates <- function(targets,year,longitude) {
  if (!length(targets)) return(character())
  start <- swephR::swe_julday(year,1L,1L,0,1L)
  end <- swephR::swe_julday(year+1L,1L,1L,0,1L)
  days <- seq(start,end)
  solar <- .solar_coordinate(days,0L,1L)
  wrap <- function(x) (x+180)%%360-180
  vapply(targets,function(target) {
    delta <- wrap(solar-target)
    i <- which(head(delta,-1L)<=0 & tail(delta,-1L)>0)[1]
    if (is.na(i)) return(NA_character_)
    when <- stats::uniroot(function(t) wrap(.solar_coordinate(t,0L,1L)-target),
                           days[i+0:1],tol=1e-7)$root
    sub(" .*","",jd2time(when+longitude/360,"UTC","J"))
  },character(1),USE.NAMES=FALSE)
}
