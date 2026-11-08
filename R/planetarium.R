#' Draw a planisphere for a site and epoch
#' @description Map the stars that can enter a site's sky at a chosen epoch, with the celestial equator
#' and ecliptic as optional references.
#' @details The default star layer draws individual Pleiades members and omits integrated clusters and other nonstellar entries. The full catalogue remains available for target searches and phase calculations.
#'
#' The chart compares the accessible part of the celestial sphere with a catalogue sky at
#' the chosen year. Catalogue stars are filtered by magnitude and by the implemented
#' declinations that rise above a level geometric horizon. Northern latitudes use
#' the north celestial pole; southern latitudes use the south celestial pole.
#' Polar distance sets radius, and right ascension sets angle with opposite
#' handedness for the two poles. Magnitude controls symbol size, and hour labels
#' are placed every three hours. The optional ecliptic uses daily solar positions over a
#' year, using geocentric coordinates without atmospheric refraction. The plot is
#' not an instantaneous view
#' through a site's measured horizon. The Milky Way uses the same diffuse-light tracing as sketchSky(), with dark cutouts and two brighter tiers referenced to the whole-band mean. Source: NASA/Goddard SVS and ESA/Gaia/DPAC, Deep Star Maps 2020. It is a fixed sky map rotated to the requested epoch with Swiss Ephemeris precession and nutation. Its width is illustrative. It does not predict visibility or reconstruct the Galaxy’s internal evolution.
#'
#' Astronomical positions require Swiss Ephemeris files covering the epoch. We
#' recommend installing \pkg{swephRdata}; see \code{\link{skyscapeR-package}} for setup.
#' @seealso \code{\link{plotAzimuth}}, \code{\link{plotBars}}, \code{\link{orbit}}, \code{\link{sketchSky}}
#' @param style Analytical display style, either \code{"publication-colour"} (default) or
#' \code{"publication-mono"}. Monochrome uses grayscale and preserves line styles.
#' @param show.legend Show chart pole, latitude, epoch, magnitude limit and the enabled layer key.
#' FALSE centres the bare chart, keeping its data and right-ascension ticks.
#' @param epoch Earlier argument name for \code{year}, retained for existing scripts.
#' @param year The epoch (year) for which the planisphere is calculated. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' @param lat Observer latitude in degrees north, between -90 and 90.
#' @param max.mag Faintest catalogue magnitude drawn as a star symbol (default 6). This does not filter constellation line endpoints.
#' @param show.ecliptic  Whether to display the ecliptic (TRUE or FALSE).
#' @param col.ecliptic Colour of the ecliptic line (default ochre).
#' @param show.milky.way Draw the Milky Way behind the stars (default TRUE). Set to FALSE for a chart without the band.
#' @param show.galactic.equator Draw the Galactic equator as a dash-dot line. Default FALSE.
#' @param show.equator Whether to display the celestial equator (TRUE or FALSE).
#' @param col.equator  Colour of the celestial equator (default muted grey-green).
#' @param star.col Colour of the star symbols (default dark grey-green).
#' @param constellations Optional constellation selection: a bundled set such as
#'   \code{"modern"}, \code{"chinese"} or \code{"chinese_mansions"}, figure names
#'   or identifiers, or a table from \code{\link{importConstellations}}.
#'   The default \code{NULL} draws no constellation lines.
#' @param col.constellations Colour of the constellation lines (default muted teal).
#'
#' @section Constellation patterns:
#' Constellation segments follow the shorter great-circle arc between catalogue
#' endpoints at the requested epoch. They use the same polar projection as the
#' stars and are clipped to the chart's accessible sky. Faint endpoints remain
#' available for line drawing even when no star symbol is shown. Missing catalogue
#' endpoints generate a warning and their segments are omitted. These are line
#' figures; repeated endpoints mark single-star figures with small rings. They
#' are not constellation boundaries or lunar-mansion sectors. See
#' \code{\link{constellations}} for sources, licences and historical limitations.
#'
#' @return Draws the polar chart on the active graphics device and invisibly returns the saved xpd, mar and pty graphics settings.
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace("swephRdata", quietly = TRUE)) {
#'   planisphere(2023, lat=50, max.mag=6, show.milky.way=FALSE, show.galactic.equator=TRUE)
#'   planisphere(2023, lat=-35, max.mag=6)
#'   planisphere(2023, lat=40, max.mag=3, constellations="chinese_mansions")
#' }
#' }
#' @references NASA/Goddard Space Flight Center Scientific Visualization Studio (2020).
#' Deep Star Maps 2020. \url{https://svs.gsfc.nasa.gov/4851/}.
#'
#' Astrodienst. Swiss Ephemeris programmer's manual.
#' \url{https://www.astro.com/swisseph/swephprg.htm}.
#' @section Projection and visibility limits:
#' This is a polar azimuthal equidistant chart. For latitude at or above zero,
#' radius is 90 minus declination and right ascension increases clockwise from
#' the top. For negative latitude, radius is 90 plus declination and right
#' ascension increases anticlockwise. The chart radius is 180 minus absolute
#' latitude; a non-round cutoff is retained exactly. All layers use this same
#' projection. At the equator the north-polar convention is used.
#'
#' Stars are included when their upper culmination is strictly above a level
#' geometric horizon; stars that only touch it are excluded. This is potential
#' annual access to that part of the sky, not visibility at any particular time.
#' There is no date/time of night, local terrain, atmospheric extinction,
#' moonlight or twilight detection model. Magnitudes are fixed catalogue values,
#' and point size is a graphical brightness indicator. For a particular site and
#' instant use \code{\link{sketchSky}} with a georeferenced horizon.
planisphere <- function(year, lat, max.mag=6, show.ecliptic=TRUE, col.ecliptic='#A9542E', show.equator=TRUE, col.equator='#597575', star.col='#23383A', epoch, show.milky.way=TRUE, show.galactic.equator=FALSE, constellations=NULL, col.constellations="#638582", show.legend=TRUE, style="publication-colour") {
  if (!missing(epoch)) {
    if (!missing(year)) stop("Supply year only once.",call.=FALSE)
    year <- epoch
  }
  if (is.character(year)) year <- astroYear(year)

  if (!is.numeric(lat) || length(lat)!=1L || !is.finite(lat) || abs(lat)>90)
    stop("lat must be one finite latitude between -90 and 90 degrees.")
  if (!is.numeric(year) || length(year)!=1L || !is.finite(year))
    stop("year must specify one finite astronomical year.")
  pole <- if (lat < 0) -1 else 1
  outer <- 180-abs(lat)
  limit <- sort(unique(c(seq(0,outer,10),outer)))
  polar.radius <- function(dec) 90-pole*dec
  angle <- function(ra) (pole*ra)%%360
  clipped.radius <- function(dec) {
    r <- polar.radius(dec)
    r[r>outer] <- NA_real_
    r
  }
  ss <- .plot_star_catalogue()
  ind <- which(ss$magV > max.mag | is.na(ss$magV)); if (length(ind)) ss <- ss[-ind,] # magnitude filter

  ra <- c(); dec <- c(); mag <- c(); type <- NA
  for (i in seq_len(NROW(ss))) {
    aux <- star(ss$identifier[i], year=year)
    ra[i] <- aux$coord$RA
    dec[i] <- aux$coord$Dec
    mag[i] <- aux$app.mag
  }

  visible <- abs(dec-lat) < 90
  dec <- dec[visible]
  ra <- ra[visible]
  mag <- mag[visible]

  theme <- .plot_style(style)
  cols <- .plot_colour(c(ecliptic=col.ecliptic,equator=col.equator,stars=star.col,figures=col.constellations),theme)
  names(cols) <- c("ecliptic","equator","stars","figures")
  old <- par(c("xpd","mar","pty")); on.exit(par(old))
  par(mai=rep(.10,4),pty="m",xpd=FALSE)
  plot.new()
  plot.window(outer*if(show.legend) c(-1.17,2.20) else c(-1.19,1.19),
              c(-1.19,1.19)*outer,asp=1,xaxs="i",yaxs="i")
  circle <- function(r,cc,lwd=.7,lty=1) lines(.polar_xy(0:360,r),col=cc,lwd=lwd,lty=lty)
  for (r in seq(0,outer,by=30)) if (r>0) circle(r,theme$grid)
  p <- .polar_xy(seq(0,330,30),outer)
  segments(0,0,p[,1],p[,2],col=theme$grid,lwd=.65)
  xy <- function(ra,dec) .polar_xy(angle(ra),clipped.radius(dec))
  if (show.milky.way || show.galactic.equator) {
    jd <- swephR::swe_julday(year,1L,1L,12,1L)
    frame <- .galactic_frame(jd)
    if (show.milky.way)
      .draw_galactic_layers(t(frame),pole=pole,outer=outer,paper=TRUE)
    if (show.galactic.equator) {
      g <- .galactic_paths(jd,FALSE,TRUE,frame)$galactic.equator
      lines(xy(g[,1],g[,2]),col=theme$muted,lty=4,lwd=1)
    }
  }
  if (show.equator) circle(90,cols["equator"],1.05,2)
  if (show.ecliptic) {
    jd <- swephR::swe_julday(year,1L,1L,0,1L)+0:366
    sun <- body.position(obj="sun",time=jd,loc=NULL,dec="geo",verbose=FALSE)$equatorial
    lines(xy(sun$RA,sun$Dec),col=cols["ecliptic"],lwd=1.2)
  }
  if (!is.null(constellations)) {
    edges <- do.call(rbind,constLines(constellations))
    if (!is.null(edges)) {
      ids <- unique(as.character(edges[!is.na(edges)]))
      pos <- setNames(lapply(ids,function(id) star(id,year=year)$coord),ids)
      for (i in seq_len(nrow(edges))) {
        if (anyNA(edges[i,])) next
        a <- pos[[edges[i,1]]]
        if (edges[i,1]==edges[i,2]) {
          if (abs(a$Dec-lat)<90) points(xy(a$RA,a$Dec),pch=1,cex=.55,col=cols["figures"])
        } else {
          arc <- .planisphere_arc(a,pos[[edges[i,2]]])
          if (nrow(arc)) lines(xy(arc[,1],arc[,2]),col=cols["figures"],lwd=.85)
        }
      }
    }
  }
  if (length(mag)) points(xy(ra,dec),pch=16,col=cols["stars"],
                         cex=vapply(mag,mag2size,numeric(1))*.80)
  circle(outer,theme$muted,.85)
  hours <- seq(0,21,3)
  text(.polar_xy(angle(15*hours),outer*1.09),paste0(hours,"h"),cex=.72,col=theme$muted)
  if (show.legend) {
    left <- outer*1.27
    text(left,outer*.68,if(pole==1) "North pole" else "South pole",
         adj=c(0,.5),font=2,cex=.91,col=theme$ink)
    text(left,outer*.49,sprintf("%.2f\u00b0 %s  |  %s",abs(lat),if(pole==1)"N" else "S",BC.AD(year)),
         adj=c(0,.5),cex=.76,col=theme$muted)
    text(left,outer*.33,paste("Magnitude limit",max.mag),adj=c(0,.5),cex=.76,col=theme$muted)
    on <- c(show.equator,show.ecliptic,show.milky.way,show.galactic.equator,!is.null(constellations))
    if (any(on)) legend(left,outer*.17,
       legend=c("Celestial equator","Ecliptic","Milky Way","Galactic equator","Constellations")[on],
       col=c(cols["equator"],cols["ecliptic"],"grey70",theme$muted,cols["figures"])[on],
       lty=c(2,1,1,4,1)[on],lwd=c(1,1.2,4,1,1)[on],cex=.73,bty="n",x.intersp=.55,
       text.col=theme$ink)
  }
  invisible(old)
}

# Sample the shorter spherical arc, avoiding RA wrap and polar discontinuities.
.planisphere_arc <- function(a, b) {
  vector <- function(p) {
    ra <- p$RA*pi/180; dec <- p$Dec*pi/180
    c(cos(dec)*cos(ra), cos(dec)*sin(ra), sin(dec))
  }
  a <- vector(a); b <- vector(b)
  if (any(!is.finite(c(a,b)))) return(matrix(numeric(),0,2))
  omega <- acos(max(-1,min(1,sum(a*b))))
  if (omega < 1e-12 || abs(pi-omega)<1e-10) return(matrix(numeric(),0,2))
  t <- seq(0,1,length.out=max(3L,ceiling(omega*180/pi/.2)+1L))
  xyz <- outer(sin((1-t)*omega),a) + outer(sin(t*omega),b)
  xyz <- xyz/sqrt(rowSums(xyz^2))
  cbind(atan2(xyz[,2],xyz[,1])*180/pi %% 360, asin(pmax(-1,pmin(1,xyz[,3])))*180/pi)
}
