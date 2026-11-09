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

calc.phase <- function(year, day, sun.coord, star, loc, arcus_visionis, alt.hor, k, limit, alt.rs, res, refraction, atm, temp, pb, pbi, horizon=NULL, solar.alt=alt.hor-6) {
  jd0 <- time2jd('2000/01/01 00:00:00',timezone='UTC',calendar='G')
  track <- .phase_track(jd0,star,loc,alt.hor,k,limit,res,refraction,atm,temp,horizon)
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
#' @description Estimate the seasons in which a star's rising or setting could be visible in twilight at
#' an archaeological site.
#' @details Seasonal visibility depends on solar depression, stellar altitude, apparent magnitude
#' and atmospheric extinction as well as terrain. The model holds the star's position fixed
#' for the chosen year and samples a 365-day solar cycle with uniformly spaced ecliptic
#' longitudes. It requires a star-Sun altitude separation of 2.1 times the catalogue
#' magnitude plus 10 degrees, and compares extinction-adjusted magnitude with the limiting
#' magnitude. A positive-duration overlap of visibility conditions determines the daily
#' code. The stellar altitude trend within a visible interval below alt.rs determines
#' rising or setting. Terrain reappearance of a descending star therefore contributes
#' a setting phase. Disjoint runs are reported as separate seasons on the circular
#' model year, including runs spanning the December-solstice reference.
#'
#' Both constant and varying horizons are evaluated on all 365 model days. A horizon
#' object supplies the coordinates and piecewise-linear horizon. The profile must
#' cover all bearings and close at north. Use alt.hor only with coordinate input.
#' A fixed-coordinate stellar track is shared across model days. The terrain search
#' uses the same corners, azimuth reversals, interior extrema and 0.01-second root
#' tolerance as \code{\link{riseset}}, collecting all crossings during the model day.
#' Visible windows can lie wholly between minute samples. A touch without clearance
#' contributes no visible duration. Source terrain resolution and numerical tolerances
#' limit the shortest resolvable window.
#'
#' Refraction selects apparent or geometric altitudes for both Sun and star.
#' A stellar horizon affects obstruction at the star's azimuth. The solar darkness
#' test uses altitude relative to the astronomical horizon, with default solar.alt
#' of -6 for a profile and alt.hor - 6 for coordinates. Supply the same solar.alt
#' to compare a constant profile and a raised scalar horizon. Extinction, the
#' star-Sun separation and alt.rs use stellar altitude above the astronomical horizon.
#' Solar altitude is interpolated between res samples. Stellar visibility margins
#' are interpolated on a grid containing minute samples, solar samples, terrain
#' crossings and stellar altitude extrema. Decrease res to assess borderline dates.
#' Any positive-duration opportunity in this threshold model qualifies. Actual
#' naked-eye detection over a very short interval depends on observing conditions.
#' The synthetic solar year and solstice calendar mapping
#' are conventions inherited from the phase model.
#'
#' Purrington (1988, p. S74, equation 1 and Fig. 2) fitted the altitude-difference
#' rule to empirical data down to about magnitude 3. It describes rising at dawn and
#' setting at dusk on the same side of the sky as the Sun. Purrington discusses a
#' smaller separation for opposite-horizon events (p. S75). This function applies
#' the same threshold at every relative azimuth. Schaefer (1987, section 4)
#' examines how sky brightness and observing conditions affect heliacal visibility.
#' Brady (2015) gives the archaeological context for the phase terminology.
#'
#' Astronomical positions require Swiss Ephemeris files covering the epoch. We
#' recommend installing \pkg{swephRdata}; see \code{\link{skyscapeR-package}} for setup.
#' @return A skyscapeR.starphases list. The data table gives each model day and its phase code: I for invisible, V for visible without a classified rise or set, R for rising, S for setting and RS for both. Metadata include the star, year, observer coordinates, phase type, season and event tables, refraction setting and solar altitude threshold. For a named horizon, \code{metadata$site.name} retains the site label for comparisons. The logical \code{metadata$horizon} records whether a horizon object was supplied. The separate \code{metadata$horizon.status} reports \code{"always above"}, \code{"always below"} or \code{"crosses"} relative to the supplied skyline. Event rows include solar.longitude in degrees on the synthetic cycle. The date.julian column gives the representative local Julian calendar date at the supplied longitude. Pass this date to \code{\link{riseset}} with \code{calendar = "J"} to time the skyline crossing. For a dawn sketch, advance from that crossing until the star clears the terrain while the sky is still dark enough. The phase date alone does not specify a visible instant.
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
#' @section Validation:
#' Seasonal searches require full skyline coverage. H07 checked rejection of a partial photograph and successful calculation with the full profile. See the \href{https://stonexsky.github.io/skyscapeR/manual/#h07-panorama-coverage}{H07 check in the technical manual}.
#'
#' W02 replayed the Gacrux heliacal-rising date at Skorba West and passed it to the crossing and sky-reconstruction calculations. See the \href{https://stonexsky.github.io/skyscapeR/manual/#w02-field-to-sky}{W02 check in the technical manual}.
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
#' @section Extinction and model scope:
#' The implementation uses Rozenberg's airmass approximation,
#' \eqn{X(h)=1/(\sin h+0.025\exp(-11\sin h))}, multiplying it by \code{k}
#' and adding the result to the fixed catalogue magnitude. Here \code{h} is the
#' stellar altitude selected by \code{refraction}, above the astronomical
#' horizontal. This is distinct from the height above a ridge. The star-Sun
#' altitude-separation rule uses the catalogue magnitude before extinction.
#'
#' The model does not calculate scattered moonlight, artificial light pollution,
#' weather, an observer's adaptation history, or variable-star light curves.
#' The same separation threshold is used at all relative azimuths. It is therefore
#' a screening model for seasonal opportunities, not a probability of detection.
#' Supplying an integrated cluster entry such as M45 treats its total magnitude as
#' one target; cluster size, surface brightness and resolution into individual
#' stars do not enter this calculation. It does not predict when each Pleiad is
#' individually visible.
#'
#' Horizon status is determined from positive-duration clearance intervals along
#' the same stellar track used for the phase calculation, with the supplied
#' skyline and refraction setting. The same calculation applies in both
#' hemispheres. A track with no clearance returns an ordinary all-I result,
#' rather than raising an error. An always-clear track still undergoes the
#' extinction and solar-darkness tests; never setting does not imply visibility.
#' A star that never sets below a level horizontal can cross a raised ridge.
#' Isolated tangencies do not count as intervals of clearance or obstruction.
#' Numerical resolution and the supplied terrain profile limit this status.
#'
#' The separate metadata$horizon.status field is "always above", "always below"
#' or "crosses", relative to the supplied skyline. Phase codes R and S describe
#' the rising or descending part of a visible track below alt.rs; they do not
#' require a physical skyline crossing that day. An always-clear star may thus
#' have seasonal R/S codes. The type "circumpolar" is retained for an all-V
#' result only when the track is always above the skyline; other all-V results
#' are labelled "visible throughout cycle". Use horizon.status, not type, to
#' distinguish clearance from seasonal visibility.
#' @section Dates and a subsequent sky sketch:
#' Events follow adjacent-sample transitions around the closed cycle, including
#' the December-solstice boundary. Entering and leaving I mark acronycal setting
#' and heliacal rising; entering and leaving V mark acronycal rising and heliacal
#' setting. A season continuing across the boundary creates no event there.
#' Constant phase sequences have no transition events.
#'
#' The 365 model samples are equally spaced in solar longitude, not in elapsed
#' civil time. Printed seasonal labels are anchored to the December solstice.
#' \code{metadata$events$date.julian} instead solves for the event solar longitude
#' in the requested year and converts it to a local mean-solar Julian-calendar
#' date using longitude. It is a representative date, not an exact detection
#' instant, and can differ from the synthetic \code{day} label.
#'
#' Pass the representative date to \code{\link{riseset}} with
#' \code{calendar = "J"} and an explicit timezone. With UTC or a fixed offset,
#' inspect neighbouring dates when local mean-solar and civil days differ.
#' The returned rise is a skyline crossing. Follow it by a short time sequence
#' with \code{\link{sketchSky}} and \code{\link{skyStep}} to inspect terrain
#' clearance and the changing brightness of the local sky.
#'
#' A supplied \code{skyscapeR.star} is used as given. Recompute it for the requested
#' year; the phase function does not update an existing object's coordinates.
#' @section How precise is a heliacal-rising date?:
#' In 69 of 70 individual comparisons, the estimated first morning of visibility was within a week of the predicted heliacal-rising date.
#' @section Reproducing a saved result:
#' Metadata retain the resolved star, year, observer coordinates, extinction,
#' atmospheric settings and solar-altitude threshold, together with \code{limit},
#' \code{alt.rs}, \code{res} and \code{alt.hor}. The latter is NA when a terrain
#' profile replaces the constant horizon. \code{metadata$horizon.profile} contains
#' the numerical \code{az}/\code{alt} table used, or NULL for a constant horizon.
#' Landscape images and DEMs are not copied into the result. Existing logical
#' \code{horizon} and \code{skyline} fields are retained.
#'
#' \code{metadata$provenance$packages} records the loaded skyscapeR and swephR
#' versions; an unavailable version is NA, for example when sourcing functions
#' outside a package namespace. Version strings do not uniquely identify
#' development builds. Retain the original script, code revision and ephemeris
#' files as well: the saved metadata do not identify or archive those files.
#' Earlier saved results do not acquire these fields automatically.
#' @references Krisciunas, K. and Schaefer, B. E. (1991). A model of the brightness
#' of moonlight. PASP 103, 1033-1039 (equation 14 for extinction airmass).
#' \doi{10.1086/132921}.
#' @section Provenance:
#' Structured results carry a \code{provenance} attribute. Read it with
#' \code{attr(result, "provenance")}; its \code{packages} entry records the loaded
#' skyscapeR and swephR versions, with NA when a namespace version is unavailable.
#' This attribute adds no list entries or table columns. Ordinary subsetting or
#' transformations may discard attributes; save the original result for reference.
#' It does not identify ephemeris files, catalogue revisions or a development
#' commit. Retain those and the analysis script separately for exact reproduction.
#' The same record remains available in \code{metadata$provenance}.
star.phases <- function(star, year, loc, alt.hor = 0, alt.rs = 10, solar.alt=NULL, k = skyscapeR.vars()$k, limit =
                        6, refraction, atm, temp, res = 1/24/6) {
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
  phase <- calc.phase(year,seq_len(365),sun.pos,star,loc,arcus_visionis,
    alt.hor,k,limit,alt.rs,res,refraction,atm,temp,pb,0,horizon=horizon,solar.alt=solar.alt)
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
