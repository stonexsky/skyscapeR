#' Draw a planisphere for a site and epoch
#' @description Map the stars that can enter a site's sky at a chosen epoch, with the celestial equator
#' and ecliptic as optional references.
#' @details The default star layer draws individual Pleiades members and omits integrated clusters and other nonstellar entries. The full catalogue remains available for target searches and phase calculations.
#'
#' The chart compares the accessible part of the celestial sphere with a catalogue sky at
#' the chosen year. Catalogue stars are filtered by magnitude and by the implemented
#' declination cutoff of latitude minus 90 degrees. Right ascension sets polar angle and 90
#' degrees minus declination sets radius. Magnitude controls symbol size, and hour labels
#' are separated by 15 degrees. The optional ecliptic uses daily solar positions over a
#' year. The cutoff is a north-polar construction; the plot is not an instantaneous view
#' through a site's measured horizon. The Milky Way is filled with translucent grey behind the stars. Its outline retains the main dark rift and follows the diffuse stellar light in NASA’s Deep Star Maps 2020 (NASA/Goddard Space Flight Center Scientific Visualization Studio, ESA/Gaia/DPAC). It is a fixed sky map rotated to the requested epoch with Swiss Ephemeris precession and nutation. Its width is illustrative. It does not predict visibility or reconstruct the Galaxy’s internal evolution.
#'
#' Astronomical positions require Swiss Ephemeris files covering the epoch. We
#' recommend installing \pkg{swephRdata}; see \code{\link{skyscapeR-package}} for setup.
#' @seealso \code{\link{plotAzimuth}}, \code{\link{plotBars}}, \code{\link{orbit}}, \code{\link{sketchSky}}
#' @param epoch Earlier argument name for \code{year}, retained for existing scripts.
#' @param year The epoch (year) for which the planisphere is calculated. Historical labels such as \code{"3700 BC"} are accepted by \code{\link{astroYear}}.
#' @param lat Numerical.Latitude of the observer in degrees.
#' @param max.mag  Numerical.Maximum magnitude of stars to be displayed (default: 6).
#' @param show.ecliptic  Whether to display the ecliptic (TRUE or FALSE).
#' @param col.ecliptic Colour of the ecliptic line (default: "orange").
#' @param show.milky.way Draw the Milky Way behind the stars (default TRUE). Set to FALSE for a chart without the band.
#' @param show.galactic.equator Draw the Galactic equator as a dash-dot line. Default FALSE.
#' @param show.equator Whether to display the celestial equator (TRUE or FALSE).
#' @param col.equator  Colour of the celestial equator (default: "red").
#' @param star.col  Colour of the stars (default: "black").
#'
#' @return A polar plot of the visible sky, showing stars, the celestial equator,
#'         and the ecliptic for the specified epoch and location.
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace("swephRdata", quietly = TRUE)) {
#'   planisphere(2023, lat=50, max.mag=6, show.milky.way=FALSE, show.galactic.equator=TRUE)
#' }
#' }
planisphere <- function(year, lat, max.mag=6, show.ecliptic=TRUE, col.ecliptic='orange', show.equator=TRUE, col.equator='red', star.col='black', epoch, show.milky.way=TRUE, show.galactic.equator=FALSE) {
  if (!missing(epoch)) {
    if (!missing(year)) stop("Supply year only once.",call.=FALSE)
    year <- epoch
  }
  if (is.character(year)) year <- astroYear(year)

  ss <- .plot_star_catalogue()
  ind <- which(ss$magV > max.mag | is.na(ss$magV)); if (length(ind)) ss <- ss[-ind,] # magnitude filter

  ra <- c(); dec <- c(); mag <- c(); type <- NA
  for (i in seq_len(NROW(ss))) {
    aux <- star(ss$identifier[i], year=year)
    ra[i] <- aux$coord$RA
    dec[i] <- aux$coord$Dec
    mag[i] <- aux$app.mag
  }

  visible <- dec >= -90+lat
  dec <- dec[visible]
  ra <- ra[visible]
  mag <- mag[visible]

  labels <- c(seq(12,24,1),seq(1,11,1))
  label.pos <- ((15 * labels) %% 360) * pi/180

  DD <- FALSE
  if (show.milky.way || show.galactic.equator) {
    jd <- swephR::swe_julday(year,1L,1L,12,1L)
    paths <- .galactic_paths(jd,show.milky.way,show.galactic.equator)
    limit <- seq(0,180-lat,10)
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
    plotrix::polar.plot(rep(90,360), seq(0,359,1), radial.lim=seq(0,90-(-90+lat),10),
                        start=-270, clockwise = TRUE, rp.type='p', line.col=col.equator,
                        labels = labels, label.pos = label.pos,
                        show.grid.labels=FALSE, add=DD, main=ifelse(year>0,paste(year,'AD'),paste(abs(year)+1,'BC')))
    DD <- TRUE
  }

  if (show.ecliptic) {
    # draw ecliptic (yellow path)
    dat <- paste0(year,'/01/01 00:00:00')
    jd0 <- time2jd(dat, timezone='UTC')
    jd <- seq(jd0, jd0+366, 1)

    sun.pos <- body.position(obj='sun', time=jd, loc=c(lat, 0, 0))$equatorial

    plotrix::polar.plot(90-sun.pos$Dec, sun.pos$RA, rp.type='p', line.col=col.ecliptic,
                        start=-270, clockwise = TRUE, radial.lim=seq(0,90-(-90+lat),10),
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
  # draw stars
  plotrix::polar.plot(90-dec, ra, rp.type='s', point.symbols=19, cex=sapply(mag, mag2size),
                      start=-270, clockwise = TRUE, radial.lim=seq(0,90-(-90+lat),10),
                      labels = labels, label.pos = label.pos,
                      show.grid.labels=FALSE, main=ifelse(year>0,paste(year,'AD'),paste(abs(year)+1,'BC')),
                      add=DD, point.col=star.col)}

