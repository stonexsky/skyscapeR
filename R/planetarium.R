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
#' are separated by 15 degrees. The optional ecliptic uses daily solar positions over a
#' year, using geocentric coordinates without atmospheric refraction. The plot is
#' not an instantaneous view
#' through a site's measured horizon. The Milky Way is filled with translucent grey behind the stars. Its outline retains the main dark rift and follows the diffuse stellar light in NASA’s Deep Star Maps 2020 (NASA/Goddard Space Flight Center Scientific Visualization Studio, ESA/Gaia/DPAC). It is a fixed sky map rotated to the requested epoch with Swiss Ephemeris precession and nutation. Its width is illustrative. It does not predict visibility or reconstruct the Galaxy’s internal evolution.
#'
#' Astronomical positions require Swiss Ephemeris files covering the epoch. We
#' recommend installing \pkg{swephRdata}; see \code{\link{skyscapeR-package}} for setup.
#' @seealso \code{\link{plotAzimuth}}, \code{\link{plotBars}}, \code{\link{orbit}}, \code{\link{sketchSky}}
#' @param epoch Earlier argument name for \code{year}, retained for existing scripts.
#' @param year The epoch (year) for which the planisphere is calculated. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' @param lat Observer latitude in degrees north, between -90 and 90.
#' @param max.mag Faintest catalogue magnitude drawn as a star symbol (default 6). This does not filter constellation line endpoints.
#' @param show.ecliptic  Whether to display the ecliptic (TRUE or FALSE).
#' @param col.ecliptic Colour of the ecliptic line (default: "orange").
#' @param show.milky.way Draw the Milky Way behind the stars (default TRUE). Set to FALSE for a chart without the band.
#' @param show.galactic.equator Draw the Galactic equator as a dash-dot line. Default FALSE.
#' @param show.equator Whether to display the celestial equator (TRUE or FALSE).
#' @param col.equator  Colour of the celestial equator (default: "red").
#' @param star.col Colour of the star symbols (default black).
#' @param constellations Optional constellation selection: a bundled set such as
#'   \code{"modern"}, \code{"chinese"} or \code{"chinese_mansions"}, figure names
#'   or identifiers, or a table from \code{\link{importConstellations}}.
#'   The default \code{NULL} draws no constellation lines.
#' @param col.constellations Colour of the constellation lines (default grey45).
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
#' @return A polar plot of the visible sky, showing stars, the celestial equator,
#'         and the ecliptic for the specified epoch and location.
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
planisphere <- function(year, lat, max.mag=6, show.ecliptic=TRUE, col.ecliptic='orange', show.equator=TRUE, col.equator='red', star.col='black', epoch, show.milky.way=TRUE, show.galactic.equator=FALSE, constellations=NULL, col.constellations="grey45") {
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

  labels <- c(seq(12,24,1),seq(1,11,1))
  label.pos <- angle(15 * labels)

  DD <- FALSE
  if (show.milky.way || show.galactic.equator) {
    jd <- swephR::swe_julday(year,1L,1L,12,1L)
    paths <- .galactic_paths(jd,show.milky.way,show.galactic.equator)
    # Transform every sky layer to the same pole and RA handedness.
    paths <- lapply(paths, function(p) cbind(angle(p[,1]),pole*p[,2]))
    if (show.milky.way) {
      plotrix::polar.plot(0,0,rp.type="s",point.col=NA,radial.lim=limit,
                          start=-270,clockwise=TRUE,labels=labels,label.pos=label.pos,
                          show.grid.labels=FALSE,main=ifelse(year>0,paste(year,'AD'),paste(abs(year)+1,'BC')))
      .draw_galactic_band(paths$south,paths$north,max(limit))
      DD <- TRUE
    }
  }
  if (show.equator) {
    # draws equator (red path)
    plotrix::polar.plot(rep(90,360), seq(0,359,1), radial.lim=limit,
                        start=-270, clockwise = TRUE, rp.type='p', line.col=col.equator,
                        labels = labels, label.pos = label.pos,
                        show.grid.labels=FALSE, add=DD, main=ifelse(year>0,paste(year,'AD'),paste(abs(year)+1,'BC')))
    DD <- TRUE
  }

  if (show.ecliptic) {
    # draw ecliptic (yellow path)
    jd0 <- swephR::swe_julday(year,1L,1L,0,1L)
    jd <- seq(jd0, jd0+366, 1)

    sun.pos <- body.position(obj='sun', time=jd, loc=NULL, dec='geo', verbose=FALSE)$equatorial

    plotrix::polar.plot(clipped.radius(sun.pos$Dec), angle(sun.pos$RA), rp.type='p', line.col=col.ecliptic,
                        start=-270, clockwise = TRUE, radial.lim=limit,
                        labels = labels, label.pos = label.pos,
                        show.grid.labels=FALSE, main=ifelse(year>0,paste(year,'AD'),paste(abs(year)+1,'BC')),
                        add=DD)
    DD <- TRUE
  }

  if (show.milky.way || show.galactic.equator) {
    for (id in names(paths)) {
      line <- paths[[id]]
      radius <- 90-line[,2]
      radius[radius > max(limit)] <- NA_real_
      plotrix::polar.plot(radius, line[,1], radial.lim=limit, start=-270,
                          clockwise=TRUE, rp.type="l", labels=labels, label.pos=label.pos,
                          show.grid.labels=FALSE, add=DD, line.col="grey45",
                          lty=if (id=="galactic.equator") 4 else 1,
                          main=ifelse(year>0,paste(year,'AD'),paste(abs(year)+1,'BC')))
      DD <- TRUE
    }
  }
  if (!is.null(constellations)) {
    edges <- do.call(rbind, constLines(constellations))
    if (!is.null(edges)) {
      ids <- unique(as.character(edges[!is.na(edges)]))
      positions <- lapply(ids, function(id) star(id, year=year)$coord)
      names(positions) <- ids
      if (!DD) {
        plotrix::polar.plot(0,0,rp.type="s",point.col=NA,radial.lim=limit,
                            start=-270,clockwise=TRUE,labels=labels,label.pos=label.pos,
                            show.grid.labels=FALSE,main=ifelse(year>0,paste(year,'AD'),paste(abs(year)+1,'BC')))
        DD <- TRUE
      }
      for (i in seq_len(nrow(edges))) {
        if (anyNA(edges[i,])) next
        if (edges[i,1] == edges[i,2]) {
          point <- positions[[edges[i,1]]]
          if (abs(point$Dec-lat)<90)
            plotrix::polar.plot(polar.radius(point$Dec), angle(point$RA), rp.type="s",
                                point.col=col.constellations, point.symbols=1, cex=.6,
                                radial.lim=limit, start=-270,clockwise=TRUE,
                                labels=labels,label.pos=label.pos,show.grid.labels=FALSE,add=TRUE)
          next
        }
        path <- .planisphere_arc(positions[[edges[i,1]]], positions[[edges[i,2]]])
        if (!nrow(path)) next
        plotrix::polar.plot(clipped.radius(path[,2]), angle(path[,1]), rp.type="l",
                            line.col=col.constellations, radial.lim=limit,
                            start=-270,clockwise=TRUE,labels=labels,label.pos=label.pos,
                            show.grid.labels=FALSE,add=TRUE)
      }
    }
  }
  # draw stars
  plotrix::polar.plot(polar.radius(dec), angle(ra), rp.type='s', point.symbols=19, cex=sapply(mag, mag2size),
                      start=-270, clockwise = TRUE, radial.lim=limit,
                      labels = labels, label.pos = label.pos,
                      show.grid.labels=FALSE, main=ifelse(year>0,paste(year,'AD'),paste(abs(year)+1,'BC')),
                      add=DD, point.col=star.col)}


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
