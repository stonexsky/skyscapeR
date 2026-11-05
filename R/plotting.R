# Set the opacity while preserving the RGB channels and any colour names.
.with_alpha <- function(col, alpha) {
  channels <- grDevices::col2rgb(col, alpha=TRUE) / 255
  grDevices::rgb(channels[1,], channels[2,], channels[3,], alpha=alpha,
                 names=colnames(channels))
}

#' Plot orientation measurements on a compass grid
#' @description Inspect the directions of structural axes on a compass and compare them with candidate
#' rise or set bearings.
#' @details Each measurement appears as a radial line of equal length. North is at the top and
#' bearings increase clockwise; repeated bearings can overplot. A data frame supplies the
#' True.Azimuth column. For optional targets the function uses the session atmosphere to compute a diurnal
#' track and seeks a zero-altitude crossing in the eastern and western portions of the
#' track. These reference bearings assume a level horizon and require the corresponding
#' crossings to exist.
#' @seealso \code{\link{plotBars}}, \code{\link{orbit}}, \code{\link{sketchSky}}, \code{\link{planisphere}}
#' @return Draws measured bearings and optional target directions on the active graphics device. North is at the top and azimuth increases clockwise.
#' @param az Azimuths in degrees clockwise from north. Can be omitted if \emph{obj} is given.
#' @param col Single color or color palette to use for plotting measurements.
#' @param lwd Line width for measured bearings. Defaults to 1.5.
#' @param lty Line type to plot measurements. Defaults to 1.
#' @param obj A \emph{skyscapeR.object} from \code{\link{skyObject}}, or
#' a list of these objects. Each group retains its own epochs and styles
#' for displaying the azimuths of celestial objects. Note that this assumes a single
#' location per group and a flat horizon of zero degrees. For a year range,
#' both stored declination boundaries are drawn as reference bearings.
#' @param show.obj.labels Logical to control whether to display celestial objects names.
#' Defaults to TRUE.
#' @export
#' @import utils stats plotrix
#' @examples
#' # Plot some azimuth data:
#' az <- c(120, 100, 93, 97, 88, 115, 112, 67)
#' plotAzimuth(az)
#'
#' # To visualize this data against the common solar and lunar targets:
#' tt <- skyObject('solar extremes', epoch=-2000, loc=c(35,-8), col='red')
#' plotAzimuth(az, obj=tt)
#'
#' # To display only celestial objects
#' plotAzimuth(obj=tt)
plotAzimuth = function(az, obj, show.obj.labels=TRUE, col='blue', lty=1, lwd=1.5) {
  atmosphere <- skyscapeR.vars()[c("refraction","atm","temp")]
  groups <- if (missing(obj)) list() else .sky_groups(obj)
  par(mar=c(1,1,1,1))

  # radial.plot takes radians for both grid positions and the starting angle.
  plotrix::radial.plot(0, 0, radial.lim=c(0,5), start=pi/2, clockwise=TRUE,
                       labels=c('N','NE','E','SE','S','SW','W','NW'),
                       label.pos=(0:7)*pi/4, show.grid.labels=FALSE)

  for (obj in groups) {
    # Add celestial objects
    for (i in seq_len(obj$n)) for (dec in unique(obj$decs[,i])) {
      rise <- c(); set <- c();
      orb <- orbit(dec, obj$loc, refraction=atmosphere$refraction, atm=atmosphere$atm, temp=atmosphere$temp)
      forb <- splinefun(orb$az, orb$alt, ties=mean)
      rise <- uniroot(forb, interval=c(min(orb$az), 180))$root
      set <- uniroot(forb, interval=c(180, max(orb$az)))$root
      tt <- round(c(rise, set),1)

      sl <- TRUE
      for (j in 1:length(tt)) {
        plotrix::polar.plot(5, tt[j], radial.lim=c(0,5), lwd=obj$lwd[i], lty=obj$lty[i], line.col=obj$col[i],
                            start=90, clockwise = TRUE, add=TRUE)

        if (show.obj.labels) {
          radial.plot.labels(5.9, tt[j], radial.lim=c(0,5), units="polar",
                             labels=colnames(obj$decs)[i], start=90, clockwise = TRUE, col=obj$col[i], cex=0.7)
        }
      }
    }
  }

  # Measurements
  if (!missing(az)) {
    if (is(az,'data.frame')) {
      names <- az$Name
      az <- az$True.Azimuth
    } else {
      names <- as.character(1:NROW(az))
    }

    if (length(col)==1) { col <- rep(col, length(az)) }

    plotrix::polar.plot(rep(5,length(az)), az, radial.lim=c(0,5), lwd=lwd, line.col=col,
                        start=90, clockwise=TRUE, add=TRUE)

  }
  mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'), side=3, at=par('usr')[2], cex=0.5, adj=1)
  invisible(list(atmosphere=atmosphere))
}


#' Plot orientation intervals in declination
#' @description Compare central orientations and their angular uncertainties across observations,
#' optionally against celestial target ranges.
#' @details A row represents one measurement. Its horizontal bar spans the supplied value minus and
#' plus its uncertainty, and an optional mark shows the centre. A scalar uncertainty is
#' reused for all rows. Sorting rearranges measurements, labels and widths together.
#' Celestial overlays use the declination axis. The widths reflect the supplied
#' uncertainties and acquire a probability interpretation only from the measurement model.
#'
#' For a two-year epoch range, celestial targets are displayed as bands
#' covering their stored minimum and maximum declinations across the full plot
#' height. Single-year targets are displayed as lines.
#' @return Draws orientation measurements and their uncertainties on the active graphics device.
#' @param val Array of declination or azimuth values.
#' @param unc Single value or array of measurement uncertainty
#' @param names Array of names of measurements in \emph{val}
#' @param unit (Optional). Either 'Declination' or 'Azimuth'. Defaults to 'Declination'.
#' @param col Color to plot measurements in. Defaults to \emph{blue}.
#' @param shade Logical to control whether to shade a polygon of measurements. Defaults
#' to \emph{TRUE}
#' @param mark Logical to control whether to mark the declination value. Defaults to
#' \emph{TRUE}
#' @param sort Logical to control whether to sort the measurements by their declination
#' value. Defaults to \emph{FALSE}
#' @param xrange Array of limits for x-axis.
#' @param yrange Array of limits for y-axis.
#' @param obj A \emph{skyscapeR.object} from \code{\link{skyObject}}, or
#' a list of these objects. Each group retains its own epochs and styles
#' for displaying the declination values of celestial objects.
#' @param show.obj.label Logical to control whether to label the celestial objects in
#' the polar plot. Defaults to \emph{TRUE}.
#' @export
#' @import utils stats graphics
#' @seealso \code{\link{skyObject}}
#' @examples
#' # Plot some declination data:
#' decs <- c(10, 12, -5, 4)
#' plotBars(decs, unc=5)
#'
#' # To visualize this data against the common solar and lunar targets:
#' tt <- skyObject(c('solar extremes','lunar extremes'), epoch=-2000, lty=c(2,3))
#' plotBars(decs, unc=5, obj=tt)
#' @examples
#' # Display a year-range target band:
#' targets <- skyObject('solar extremes', epoch=c("4001 BC", "2001 BC"), col='red')
#' plotBars(c(-24,-23,23,24), unc=0.5, obj=targets)
plotBars <- function(val, unc, names, unit='Declination', obj, show.obj.label=TRUE, xrange, yrange,
                     sort=FALSE, shade=TRUE, mark=TRUE, col='blue') {
  groups <- if (missing(obj)) list() else .sky_groups(obj)
  if (min(val)< -90 | max(val)>90 & unit=='Declination') {
    .sky_message('It appears that val includes azimuth values. Changing unit to Azimuth')
    unit <- 'Azimuth'}
  if (NROW(unc)==1) { unc <- rep(unc, NROW(val)) }
  if (sort) {
    ind <- sort(val, decreasing=TRUE, index.return=TRUE)$ix
    val <- val[ind]
    unc <- unc[ind]
    if (!missing(names)) { names <- names[ind] }
  }

  if (missing(xrange)) xrange <- c(min(val-unc)-5, max(val+unc)+5)
  if (missing(yrange)) yrange <- c(0.5, length(val)+0.5)
  rows <- seq_along(val)
  rows <- rows[rows >= min(yrange) & rows <= max(yrange)]
  labels <- if (missing(names)) rows else names[rows]
  label.width <- max(0, strwidth(as.character(labels),units="inches"), na.rm=TRUE)
  label.size <- min(.85, .40*par("fin")[1]/max(.01,label.width))
  left <- label.width*label.size + .30
  par(mai=c(.80,left,.28,.20), mgp=c(2.2,.65,0), pty="m")
  plot.default(val[1],1,type="n", xlab=paste(unit,'(\u00b0)'), ylab='',
               xlim=xrange, ylim=yrange, axes=FALSE, xaxs="i", yaxs="i")
  # Earlier plots may have enabled drawing in the margins. Restrict each
  # axis explicitly, including minor ticks and zero on a cropped azimuth axis.
  inside <- function(at) at[at >= min(xrange) & at <= max(xrange)]
  major <- pretty(xrange)
  axis(1, at=sort(unique(inside(c(major,0)))), xpd=FALSE)
  step <- abs(diff(major)[1])
  minor <- if (step <= 2) .5 else if (step <= 5) 1 else if (step <= 20) 5 else 10
  axis(1, at=inside(seq(-90,360,by=minor)), labels=FALSE, lwd=.3, xpd=FALSE)
  if (length(rows)) axis(2, at=rows, labels=labels, las=2, cex.axis=label.size,
                          gap.axis=-1, xpd=FALSE)

  box()

  if (shade) { border <- NA } else { border <- col; col <- NA  }

  for (i in 1:NROW(val)) {
    xp <- c(val[i]-unc[i], val[i]+unc[i], val[i]+unc[i], val[i]-unc[i])
    yp <- c(i-0.4,i-0.4,i+0.4,i+0.4)
    polygon(xp, yp, col=.with_alpha(col,0.5), border=border)
    if (mark) { lines(c(val[i],val[i]), c(i-0.4,i+0.4), col=col, lwd=2) }
  }

  .plot_target_bands(groups,unit,show.obj.label)
  mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'), side=3, at=par('usr')[2], cex=0.5, adj=1)
}




#' Plot a site horizon and celestial targets
#' @description Inspect the terrain along a site's horizon and compare it with the apparent paths of
#' candidate celestial targets.
#' @details A track meeting the horizon identifies a possible rise or set direction for that
#' declination. The horizontal axis is true bearing and the vertical axis is angular
#' altitude. The function fills the terrain below the stored profile and can overlay
#' diurnal tracks computed for the observing location. Refraction affects those tracks.
#' Rescaling the display or changing its colours does not change the underlying horizon. Add lines or polygons in azimuth and altitude immediately after plotting. Set \code{par(mfrow=c(2,1))} before plotting to stack two horizons.
#'
#' For a two-year epoch range, draw tracks for the stored minimum and maximum
#' declinations, using the session refraction setting for both tracks. Line
#' styles are taken from the target object; older objects without them use
#' solid lines of width 1.
#' @return Draws the horizon and returns the session atmosphere invisibly in a list.
#' @param x Object of \emph{skyscapeR.horizon} format.
#' @param show.az Logical that controls whether to display azimuth values or cardinal
#' directions. Defaults to TRUE when xrange is supplied.
#' @param xlim Earlier name for xrange.
#' @param xrange Azimuth range for plotting.
#' @param ylim Earlier name for yrange.
#' @param yrange Altitude range for plotting. By default one degree of altitude occupies the same distance as one degree of azimuth on the device.
#' @param obj A \emph{skyscapeR.object} from \code{\link{skyObject}}, or
#' a list of these objects. Each group retains its own epochs and styles
#' for displaying the paths of celestial objects.
#' @param col.ground Color of the ground. Defaults to \emph{darkgreen}.
#' @details A horizon with an attached DEM, PeakFinder or photographic PNG displays the image by default. Add show.skyline = TRUE to inspect its numerical profile, or use.image = FALSE to draw the silhouette alone. With xrange supplied, azimuth labels default to degrees. A partial photograph initially opens on its recorded sector. An image-backed horizon or a supplied sky uses a 2:1 plot box. Copy the PNG with the saved horizon and use relinkHor() after moving it. With automatic altitude limits, a 30-degree azimuth sector shows 15 degrees of altitude.
#' @param sky A colour string or a background from \code{\link{sketchSky}} or \code{\link{stellariumSky}}.
#' @param show.skyline Draw the numerical skyline. By default the line appears only without a panorama PNG.
#' @param fastLoad Use a sampled PNG for quick plots. New photographic imports store a preview. For older horizons the first preview loads the image. The numerical skyline and celestial tracks retain their precision.
#' @param use.image Use an attached external panorama PNG when available. Default TRUE.
#' Missing or changed files warn and fall back to the numerical silhouette.
#' @param col.sky Earlier name for a solid sky colour. Prefer sky.
#' @param ... Additional arguments to be passed to \emph{plot}.
#' @rdname plot.skyscapeR.horizon
#' @export
#' @import utils stats graphics grDevices
#' @seealso \code{\link{downloadHWT}}, \code{\link{skyObject}}
#' @examples
#' # Plot a bundled Maltese horizon:
#' hor <- malta_horizons[[1]]
#' plot(hor)
#'
#' # Add the paths of the solstices and equinoxes sun in the year 1999 BC:
#' tt <- skyObject('solar extremes', epoch=-2000, col='blue')
#' plot(hor, obj=tt)
#' @examples
#' # Compare the bounding tracks across a range of years:
#' hor <- createHor(c(0,360), c(2,2), loc=c(35,-8,25))
#' targets <- skyObject('solar extremes', epoch=c("4001 BC", "2001 BC"), col='red')
#' plot(hor, obj=targets)
plot.skyscapeR.horizon <- function(x, show.az=(!missing(xrange) || !missing(xlim)), xrange, yrange, obj, col.ground='darkgreen', col.sky, use.image=TRUE, ..., show.skyline=NULL, sky=NULL, xlim, ylim, fastLoad=FALSE) {
  if (!missing(xlim)) {
    if (!missing(xrange)) stop("Supply only one of xrange and xlim.")
    xrange <- xlim
  }
  if (!missing(ylim)) {
    if (!missing(yrange)) stop("Supply only one of yrange and ylim.")
    yrange <- ylim
  }
  atmosphere <- skyscapeR.vars()[c("refraction","atm","temp")]
  if (any(c("refraction","atm","temp") %in% names(match.call(expand.dots=FALSE)$...)))
    stop("Set observing conditions with skyscapeR.vars() before plotting.")
  groups <- if (missing(obj)) list() else .sky_groups(obj)
  force(show.az)
  if (missing(xrange)) {
    coverage <- x$metadata$source$coverage
    xrange <- if (identical(x$metadata$source$provider,"Panorama") && diff(coverage)<360-1e-8)
      coverage else c(0,360)
  }
  if (xrange[2]<xrange[1]) xrange[2] <- xrange[2]+360

  old.style <- par("mgp")
  on.exit(par(mgp=old.style), add=TRUE)
  par(mar=c(3,2,1.5,2), mgp=c(2,.65,0))
  if (!missing(col.sky)) {
    if (!is.null(sky)) stop("Supply only one sky background.")
    sky <- col.sky
  }
  if (!is.null(sky) || (use.image && !is.null(x$metadata$image))) {
    pin <- par("pin")
    box <- par("plt")
    fit <- min(pin[1]/2,pin[2])*c(2,1)/pin
    inset <- (1-fit)*c(diff(box[1:2]),diff(box[3:4]))/2
    par(plt=box+c(inset[1],-inset[1],inset[2],-inset[2]))
  }
  if (missing(yrange)) {
    pin <- par("pin")
    span <- min(180,abs(diff(xrange))*pin[2]/pin[1])
    visible <- (x$data$az-xrange[1]) %% 360 <= abs(diff(xrange))
    centre <- if (any(visible)) mean(range(x$data$alt[visible],finite=TRUE)) else mean(range(x$data$alt,finite=TRUE))
    centre <- max(-90+span/2,min(90-span/2,centre))
    yrange <- centre+c(-.5,.5)*span
    # Keep the axis margins when a full-width panorama needs a shorter panel.
    height <- pin[1]*span/abs(diff(xrange))
    if (height < pin[2]) {
      bounds <- par("plt")
      inset <- diff(bounds[3:4]) * (1 - height/pin[2]) / 2
      bounds[3:4] <- bounds[3:4] + c(inset,-inset)
      par(plt=bounds)
    }
  }
  plot(-99999,-99999, xlab = "", ylab = "", yaxs='i', xaxs='i', axes=FALSE, lwd=5, xlim=xrange, ylim=yrange, ...)
  ticks <- function(by) {
    lo <- ceiling(xrange[1]/by); hi <- floor(xrange[2]/by)
    if (lo>hi) numeric() else seq(lo,hi)*by
  }
  scale <- diff(pretty(xrange))[1]
  minor <- if (scale<=1) .1 else if (scale<=2) .5 else if (scale<=5) 1 else 5
  axis(1,at=ticks(minor),lwd=.5,labels=FALSE)
  if (show.az) {
    by <- if (scale>=90) 90 else if (scale>=45) 45 else if (scale>=10) 10 else scale
    at <- ticks(by)
    axis(1,at=at,labels=signif(at %% 360,12),lwd=0,lwd.ticks=1)
  } else {
    at <- ticks(45)
    axis(1,at=at,labels=c("N","NE","E","SE","S","SW","W","NW")[(at/45)%%8+1],lwd=.5)
  }

  # Sky colour
  if (is.character(sky) && length(sky)==1L) {
    box <- par("usr")
    rect(box[1],box[3],box[2],box[4],col=sky,border=NA)
  } else if (inherits(sky,"skyscapeR.sketch")) {
    .draw_sketch_sky(sky)
  } else if (!is.null(sky)) {
    .draw_stellarium_sky(sky)
  }
  draw.track <- function(orb,i,obj) {
    for (turn in seq(floor(xrange[1]/360),ceiling(xrange[2]/360)-1))
      lines(orb$az+360*turn,orb$alt,col=obj$col[i],lty=obj$lty[i],lwd=obj$lwd[i])
  }
  # objects
  for (obj in groups) {
    # Older year-range objects did not store line styles.
    if (is.null(obj$lty)) obj$lty <- rep(1, obj$n)
    if (is.null(obj$lwd)) obj$lwd <- rep(1, obj$n)
    ind <- sort(obj$decs[1,], decreasing=TRUE, index.return=TRUE)$ix
    for (i in ind) {
      if (length(obj$epoch)==1) {
        orb <- orbit(obj$decs[i], x, res=0.5, refraction=atmosphere$refraction, atm=atmosphere$atm, temp=atmosphere$temp)
        draw.track(orb,i,obj)
      } else {
        orb1 <- orbit(obj$decs[1,i], x, res=0.5, refraction=atmosphere$refraction, atm=atmosphere$atm, temp=atmosphere$temp)
        orb2 <- orbit(obj$decs[2,i], x, res=0.5, refraction=atmosphere$refraction, atm=atmosphere$atm, temp=atmosphere$temp)
        draw.track(orb1,i,obj)
        draw.track(orb2,i,obj)
      }
    }
  }

  .draw_horizon(x,col.ground,use.image,show.skyline,fastLoad)

  box()
  mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'), side=3, at=par('usr')[2], cex=0.5, adj=1)
  invisible(list(atmosphere=atmosphere))
}


#' Plot orientation probability distributions
#' @description Compare the uncertainty distributions of individual orientations or examine how a single
#' sightline was transformed to declination.
#' @details Individual curves reveal whether measurements overlap or have disjoint high-density
#' regions. In the stacked view, each selected density is scaled by its own peak for
#' display and its highest-density regions are shaded using hpdi. Heights across different
#' rows therefore do not measure relative probability mass. Use the linked panels for one
#' transformed observation to inspect how its source azimuth distribution and horizon
#' produce the declination distribution. Its annotation gives the measurement name
#' and the highest-density declination interval or intervals containing the requested
#' probability mass, labelled as the confidence range.
#'
#' With neither source input supplied, plot the distributions alone. With both
#' \code{az} and \code{hor}, show the information, azimuth, transformation and
#' declination panels. A single \code{index} is required for combined panels.
#' Supplying only one source input is an error, as is supplying source inputs
#' for an azimuth distribution. Source inputs must be named explicitly.
#'
#' Supply the original sources in measurement order. Counts and structure are checked,
#' but the function cannot verify their identity. The transformation curve uses the
#' refraction, pressure and temperature recorded with the result.
#' @return Draws the selected probability distributions and requested intervals on the active graphics device.
#' @param x A \emph{skyscapeR.pdf} object created with either \code{\link{az.pdf}} or \code{\link{coordtrans}}
#' @param index A value to indicate which distribution to plot, when only one is desired.
#' @param hdr High density region to highlight. Defaults to 0.954
#' @param show.az Deprecated and ignored, with a warning when supplied. Provide both
#' \code{az} and \code{hor} for the combined panels, or omit both for distributions alone.
#' @param xlim Earlier name for xrange.
#' @param xrange Range of x-axis to plot.
#' @param col Color of high density region to highlight. Defaults to blue with 0.5 alpha.
#' @param ... Additional arguments to be passed to \emph{plot}.
#' @param az Original azimuth \emph{skyscapeR.pdf}, in the same order and with
#' the same number of distributions as \code{x}. Supply together with \code{hor}
#' to show combined panels for a single declination distribution.
#' @param hor Original horizon object, or a list of horizons in the same order
#' and with the same length as the distributions in \code{x}. Supply with \code{az}.
#' @seealso \code{\link{az.pdf}}, \code{\link{dec.pdf}}, \code{\link{coordtrans}}, \code{\link{hpdi}}
#' @examples
#' Az <- az.pdf(az=c(87,93), unc=3, verbose=FALSE)
#' hor <- createHor(c(0,360), c(2,2), loc=c(35,-8,25))
#' Dec <- coordtrans(Az, hor, verbose=FALSE)
#' plot(Dec)
#' plot(Dec, index=1, az=Az, hor=hor)
#' @rdname plot.skyscapeR.pdf
#' @export
plot.skyscapeR.pdf <- function(x, index, hdr=0.954, show.az, xrange, col='#0000FF80', ..., az, hor, xlim) {
  if (!missing(xlim)) {
    if (!missing(xrange)) stop("Supply only one of xrange and xlim.")
    xrange <- xlim
  }
  if (!missing(show.az)) {
    warning("show.az is deprecated and ignored. Supply both az and hor for the combined plot, or omit both.", call.=FALSE)
  }
  if (!inherits(x, 'skyscapeR.pdf') || !x$metadata$coord %in% c('az', 'dec')) {
    stop("x must be an azimuth or declination skyscapeR.pdf object.")
  }
  n <- length(x$data)
  if (missing(index)) { index <- seq_len(n) }
  if (!is.numeric(index) || !length(index) || any(!is.finite(index)) ||
      any(index != floor(index)) || any(index < 1 | index > n)) {
    stop("index must contain valid positive integer distribution indices.")
  }
  if (xor(missing(az), missing(hor))) {
    stop("Please supply both az and hor for the combined plot, or omit both.")
  }
  combined <- !missing(az) && !missing(hor)
  if (combined) {
    if (x$metadata$coord != 'dec') stop("az and hor are only used with declination distributions.")
    if (length(index) != 1L) stop("Select a single index for the combined plot.")
    if (!inherits(az, 'skyscapeR.pdf') || !identical(az$metadata$coord, 'az') ||
        length(az$data) != n) {
      stop("az must be the original azimuth PDF with the same number and order of distributions as x.")
    }
    if (inherits(hor, 'skyscapeR.horizon')) {
      hh <- hor
    } else if (is.list(hor) && length(hor) == n &&
               all(vapply(hor, inherits, logical(1), 'skyscapeR.horizon'))) {
      hh <- hor[[index]]
    } else {
      stop("hor must be a horizon object or a list of horizons matching the distributions in x.")
    }
    a <- az$data[[index]]
    if (!is.data.frame(a) || !all(c('x','y') %in% names(a)) ||
        !is.numeric(a$x) || !is.numeric(a$y) || nrow(a) < 2L ||
        any(!is.finite(a$x)) || any(!is.finite(a$y)) || max(a$y) <= 0) {
      stop("The selected azimuth distribution must contain finite numeric x and y values.")
    }
    h <- hh$data
    if (!is.data.frame(h) || !all(c('az','alt','alt.unc') %in% names(h)) ||
        !all(vapply(h[c('az','alt','alt.unc')], is.numeric, logical(1))) ||
        nrow(h) < 2L || any(!is.finite(as.matrix(h[c('az','alt','alt.unc')]))) ||
        !is.numeric(hh$metadata$georef) || length(hh$metadata$georef) != 3L ||
        any(!is.finite(hh$metadata$georef)) ||
        (is.null(hh$metadata$ID) && (min(a$x) < min(h$az) || max(a$x) > max(h$az)))) {
      stop("The selected horizon must have valid coordinates and finite profile data covering the azimuth distribution.")
    }
    settings <- x$metadata$param
    if (!is.logical(settings$refraction) || length(settings$refraction) != 1L ||
        is.na(settings$refraction) ||
        !is.numeric(settings$atm) || length(settings$atm) != 1L || !is.finite(settings$atm) ||
        !is.numeric(settings$temp) || length(settings$temp) != 1L || !is.finite(settings$temp)) {
      stop("x is missing valid stored refraction, pressure or temperature settings.")
    }
  }

  if (combined) {
    # Coordinate Transformation Panels
    curve <- .unwrap_north_pdf(az$data[[index]])
    xx <- curve$x
    dens <- curve$y
    xl <- range(x$data[[index]]$x)
    yl <- range(xx)

    # Advance the caller's layout once. Draw the component plots inside this figure.
    par(mar=rep(0,4))
    plot.new()
    figure <- par("fin")
    line <- par("csi") * par("mex")
    panel <- function(left, bottom, right, top, margins) {
      insets <- c(bottom*figure[2], left*figure[1],
                  (1-top)*figure[2], (1-right)*figure[1])
      margin.scale <- min(1, (top-bottom)*figure[2]/(line*(sum(margins[c(1,3)])+1)),
                          (right-left)*figure[1]/(line*(sum(margins[c(2,4)])+1)))
      par(mai=insets + margin.scale*line*margins, new=TRUE)
    }
    hd <- hpdi(x$data[[index]],hdr)
    intervals <- sprintf("%.2f\u00b0 to %.2f\u00b0",hd[,1],hd[,2])
    heading <- paste(strwrap(x$metadata$name[index],width=32),collapse="\n")
    caption <- paste0(paste(intervals,collapse="\n")," (",format(100*hdr),"% confidence range)")

    # Azimuth Distribution
    panel(0, .25, .25, 1, c(0,4,1,0))
    plot(dens, xx, type='l', ylab='Azimuth (\u00b0)', xlab='', lwd=2, col='black', xlim=c(0, 1.2*max(dens)), ylim=yl, xaxs='i', yaxs='i', axes=FALSE, ...); box()
    polygon(c(dens,rep(0,length(dens))), c(xx, rev(xx)), col='grey', border=NA)
    axis(2, at=pretty(seq(par('usr')[3],par('usr')[4])))
    axis(2, at=0, labels = 0)
    scale <- mean(diff(pretty(seq(par('usr')[3],par('usr')[4]))))
    if (scale <= 1) { axis(2, at=seq(-45,360+45,0.1), lwd=0.2, labels=FALSE) }
    if (scale <= 2 & scale > 1) { axis(2, at=seq(-45,360+45,0.5), lwd=0.2, labels=FALSE) }
    if (scale <= 5 & scale > 2) { axis(2, at=seq(-45,360+45,1), lwd=0.5, labels=FALSE) }
    if (scale <= 20 & scale > 5) { axis(2, at=seq(-45,360+45,5), lwd=0.5, labels=FALSE) }
    if (scale > 10) { axis(2, at=seq(-45,360+45,10), lwd=0.5, labels=FALSE) }

    # Transformation Curve
    panel(.25, .25, 1, 1, c(0,0,1,1))
    plot(-999,999, axes=FALSE, xlim=xl, ylim=yl, xaxs='i', yaxs='i', ...); box()

    ## process horizon profile
    profile <- .coord_horizon(hh)
    query <- .profile_bearings(xx, profile$az)
    alt <- approx(profile$az, profile$alt, xout=query)$y
    alt.unc <- approx(profile$az, profile$alt.unc, xout=query)$y
    refraction <- x$metadata$param$refraction
    atm <- x$metadata$param$atm
    temp <- x$metadata$param$temp

    dec0 <- az2dec(xx, hh, alt, refraction=refraction, atm=atm, temp=temp)
    dec1 <- az2dec(xx, hh, alt+alt.unc, refraction=refraction, atm=atm, temp=temp)
    dec2 <- az2dec(xx, hh, alt-alt.unc, refraction=refraction, atm=atm, temp=temp)
    dec3 <- az2dec(xx, hh, alt+2*alt.unc, refraction=refraction, atm=atm, temp=temp)
    dec4 <- az2dec(xx, hh, alt-2*alt.unc, refraction=refraction, atm=atm, temp=temp)

    yp <- c(xx, rev(xx))
    xp <- c(dec3, rev(dec4))
    polygon(xp, yp, border=NA, col=.with_alpha('blue', 0.25))
    yp <- c(xx, rev(xx))
    xp <- c(dec1, rev(dec2))
    polygon(xp, yp, border=NA, col=.with_alpha('blue', 0.5))
    lines(dec0, xx, lwd=0.5, col='black')


    mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'), side=3, at=par('usr')[2], line=.1, cex=0.5, adj=1)
    # Declination Distribution
    panel(.25, 0, 1, .25, c(4,0,0,1))
    xx <- x$data[[index]]
    plot.new()
    size <- min(.8,.9*par("pin")[1]/max(strwidth(heading,font=2,units="inches"),strwidth(caption,units="inches")),
                .45*par("pin")[2]/(strheight(heading,font=2,units="inches")+strheight(caption,units="inches")+.08))
    w <- max(strwidth(heading,font=2,units="inches",cex=size),strwidth(caption,units="inches",cex=size))/par("pin")[1]*diff(xl)
    heading.height <- (strheight(heading,font=2,units="inches",cex=size)+.08*size)/par("pin")[2]
    h <- heading.height+strheight(caption,units="inches",cex=size)/par("pin")[2]
    pad <- .02*diff(xl)
    left.end <- min(xl[2],xl[1]+w+2*pad)
    right.start <- max(xl[1],xl[2]-w-2*pad)
    peak <- function(lo,hi) max(xx$y[xx$x>=lo & xx$x<=hi],approx(xx$x,xx$y,c(lo,hi),rule=2)$y)
    peaks <- c(peak(xl[1],left.end),peak(right.start,xl[2]))
    side <- which.min(peaks)
    ymax <- max(1.08*max(xx$y),peaks[side]/(1-h-.08))
    plot.window(xlim=xl,ylim=c(0,ymax),xaxs="i",yaxs="i")
    title(xlab="Declination (\u00b0)")
    box()
    label.x <- if (side==1L) xl[1]+pad else xl[2]-pad
    text(label.x,.97*ymax,heading,font=2,adj=c(if (side==1L) 0 else 1,1),cex=size)
    text(label.x,(.97-heading.height)*ymax,caption,
         adj=c(if (side==1L) 0 else 1,1),cex=size)


    for (j in 1:NROW(hd)) {
      ind <- which(x$data[[index]]$x >= min(hd[j,]) & x$data[[index]]$x <= max(hd[j,]))
      xvals <- x$data[[index]]$x[ind]; xvals <- c(xvals, rev(xvals))
      yvals <- x$data[[index]]$y[ind]; yvals <- c(yvals, rep(0, length(yvals)))
      polygon(xvals, yvals, col=col, border=NA)

    }
    lines(x$data[[index]]$x, x$data[[index]]$y)

    axis(1, at=pretty(seq(par('usr')[1],par('usr')[2])))
    axis(1, at=0, labels = 0)
    scale <- mean(diff(pretty(seq(par('usr')[1],par('usr')[2]))))
    if (scale <= 1) { axis(1, at=seq(-90,90,0.1), lwd=0.2, labels=FALSE) }
    if (scale <= 2 & scale > 1) { axis(1, at=seq(-90,90,0.5), lwd=0.2, labels=FALSE) }
    if (scale <= 5 & scale > 2) { axis(1, at=seq(-90,90,1), lwd=0.5, labels=FALSE) }
    if (scale <= 20 & scale > 5) { axis(1, at=seq(-90,90,5), lwd=0.5, labels=FALSE) }
    if (scale > 10) { axis(1, at=seq(-90,90,10), lwd=0.5, labels=FALSE) }

    box()


  } else {
    ## Single or Multiple Distributions
    if (length(col) == 1 & length(index)>1) { col <- rep(col, length(index))}
    if (x$metadata$coord == 'dec') { lab <- 'Declination (\u00b0)'; xl <- c(-90,90) }
    if (x$metadata$coord == 'az') { lab <- 'Azimuth (\u00b0)'; xl <- c(0,360)}
    par(mar=c(4,1,1,1))
    if(!missing(xrange)) { xl <- xrange } # else { xl <- range(hpdi(x, 1)) }
    plot(-9999, -9999, type='l', main='', xlab=lab, ylab='', lwd=2, col='black', xlim=xl, ylim=c(0,length(index)), xaxs='i', yaxs='i', axes=F, ...); box()
    axis(1, at=pretty(seq(par('usr')[1],par('usr')[2])))
    axis(1, at=0, labels = 0)
    scale <- mean(diff(pretty(seq(par('usr')[1],par('usr')[2]))))
    if (scale <= 1) { axis(1, at=seq(-90,360+50,0.1), lwd=0.2, labels=FALSE) }
    if (scale <= 2 & scale > 1) { axis(1, at=seq(-90,360+50,0.5), lwd=0.2, labels=FALSE) }
    if (scale <= 5 & scale > 2) { axis(1, at=seq(-90,360+50,1), lwd=0.5, labels=FALSE) }
    if (scale <= 20 & scale > 5) { axis(1, at=seq(-90,360+50,5), lwd=0.5, labels=FALSE) }
    if (scale > 10) { axis(1, at=seq(-90,360+50,10), lwd=0.5, labels=FALSE) }

    for (i in 1:length(index)) {
      xx <- x$data[[index[i]]]
      hd <- hpdi(xx, hdr)
      for (j in 1:NROW(hd)) {
        ind <- which(xx$x >= min(hd[j,]) & xx$x <= max(hd[j,]))
        xvals <- xx$x[ind]; xvals <- c(xvals, rev(xvals))
        yvals <- xx$y[ind]; yvals <- yvals/max(yvals)*.8
        yvals <- c(yvals, rep(0, length(yvals)))
        polygon(xvals, yvals + (i-1) + 0.1, col=col[i], border=NA)
      }
      ind <- which(diff(xx$x) > 5*mean(diff(xx$x)))
      if (length(ind)>0) { xx[ind,2] <- NA }
      lines(xx$x, xx$y/max(xx$y, na.rm=TRUE)*.8 + (i-1) + 0.1, lwd=1.2, col='black')
      # text(par('usr')[1] - 20, par('usr')[3]+i-0.5, labels=x$metadata$name[index[i]], xpd=T)
      label <- x$metadata$name[index[i]]
      if (length(label) && !is.na(label) && nzchar(label)) {
        limits <- par("usr")[1:2]
        pad <- .01*diff(limits)
        size <- min(1,(diff(limits)-2*pad)/max(strwidth(label),pad))
        width <- strwidth(label,cex=size)
        left <- if (limits[2]-max(hd) >= width+pad) max(hd)+pad else min(hd)-pad-width
        text(max(limits[1]+pad,min(left,limits[2]-pad-width)),i-.5,label,adj=0,cex=size)
      }
      abline(h=i, col='grey', lty=3)
      mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'), side=3, at=par('usr')[2], cex=0.5, adj=1)
    }

  }
}


#' Plot a summed orientation density
#' @description Visualise concentrations of orientation density across a group of structures.
#' @details A peak can arise where several uncertain sightlines coincide. The horizontal axis uses
#' the summed curve's stored coordinate and the vertical axis its density. Optional shading
#' fills the area beneath the curve; changing colours or axis limits does not renormalise
#' it. Whether its integral is one depends on the normalisation used when creating the
#' object. Assess departures from a specified null using a significance test.
#' @seealso \code{\link{spd}}, \code{\link{plot.skyscapeR.pdf}}
#' @examples az <- az.pdf(az = c(85, 92), unc = 2, verbose = FALSE)
#' plot(spd(az))
#' @return Draws the summed density against its angular coordinate on the active graphics device.
#' @param x A \emph{skyscapeR.spd} object created with \code{\link{spd}}
#' @param xlim Earlier name for xrange.
#' @param xrange Range of x-axis to plot.
#' @param ylim Earlier name for yrange.
#' @param yrange Range of y-axis to plot.
#' @param title Title to add to the plot.
#' @param col Color of summed probability density. Defaults to blue.
#' @param shading Logical to control whether to color the entire distribution. Defaults to TRUE
#' @param ... Additional arguments to be passed to \emph{plot}.
#' @param obj A skyObject or list of target groups for a density or significance plot. Declination plots show single epochs as lines and epoch ranges as bands. Azimuth plots use each group's site latitude to mark rising and setting boundary bearings at a geometric level horizon.
#' @rdname plot.skyscapeR.spd
#' @export
plot.skyscapeR.spd <- function(x, xrange, yrange, title=NULL, col='blue', shading=TRUE, ..., obj=NULL, xlim, ylim) {
  if (!missing(xlim)) {
    if (!missing(xrange)) stop("Supply only one of xrange and xlim.")
    xrange <- xlim
  }
  if (!missing(ylim)) {
    if (!missing(yrange)) stop("Supply only one of yrange and ylim.")
    yrange <- ylim
  }
  par(mar=c(5, 4, 1, 1) + 0.1)
  if (missing(xrange)) { xrange <- sort(x$data$x[c(min(which(x$data$y >= 1e-12)), max(which(x$data$y >= 1e-12)))]) }
  if (missing(yrange)) { yrange <- c(0, max(x$data$y)) }

  if (x$metadata$coord == 'dec') { label <- 'Declination (\u00b0)' }
  if (x$metadata$coord == 'az') {
    label <- 'Azimuth (\u00b0)'
    xrange[1] <- max(xrange[1], 0)
    xrange[2] <- min(xrange[2], 360)
  }

  plot(x$data$x, x$data$y, type='l', main=title, xlab=label, ylab='Density', lwd=2, col=col, xlim=xrange, ylim=yrange, xaxs='i', yaxs='i', axes=FALSE, ...); box()
  if (shading) {
    xp <- c(x$data$x, rev(x$data$x))
    yp <- c(x$data$y, rep(0, length(x$data$x)))
    polygon(xp, yp, col=.with_alpha(col,0.5), border=NA)
  }
  if (x$metadata$coord == 'dec') {
    axis(1, at=pretty(seq(par('usr')[1],par('usr')[2])))
    axis(1, at=0, labels = 0)
  } else if (x$metadata$coord == 'az') {
    axis(1, at=c(0, 45, 90, 135, 180, 225, 270, 315, 360))
  }
  scale <- mean(diff(pretty(seq(par('usr')[1],par('usr')[2]))))
  if (scale <= 1) { axis(1, at=seq(-90,360,0.1), lwd=0.2, labels=FALSE) }
  if (scale <= 2 & scale > 1) { axis(1, at=seq(-90,360,0.5), lwd=0.2, labels=FALSE) }
  if (scale <= 5 & scale > 2) { axis(1, at=seq(-90,360,1), lwd=0.5, labels=FALSE) }
  if (scale <= 20 & scale > 5) { axis(1, at=seq(-90,360,5), lwd=0.5, labels=FALSE) }
  if (scale > 10) { axis(1, at=seq(-90,360,10), lwd=0.5, labels=FALSE) }
  .plot_density_targets(obj,x$metadata$coord)
  axis(2)
  mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'), side=3, at=par('usr')[2], cex=0.5, adj=1)

}


.plot_density_targets <- function(obj,coord) {
  groups <- .sky_groups(obj)
  if (coord == "dec") {
    .plot_target_bands(groups,"Declination",TRUE)
  } else for (group in groups) {
    site <- group$loc
    if (inherits(site,"skyscapeR.horizon")) site <- as.numeric(site$metadata$georef)
    if (!is.numeric(site) || length(site)<2L) stop("Supply a site in each skyObject for an azimuth plot.", call.=FALSE)
    for (i in seq_len(group$n)) {
      cosine <- sin(group$decs[,i]*pi/180)/cos(site[1]*pi/180)
      bearings <- acos(cosine[abs(cosine)<=1])*180/pi
      for (v in unique(c(bearings,360-bearings))) {
        abline(v=v,col=group$col[i],lty=group$lty[i],lwd=group$lwd[i])
      }
    }
  }
}

#' Plot an orientation significance test
#' @description Locate the orientations at which the observed summed density departs from the simulated
#' comparison envelope. New tests use an envelope pooled from the observation and simulations.
#' @details The observed curve and envelope use the azimuth or declination grid stored in the test
#' result. The two-tailed display shades between its lower and upper pointwise bounds; the
#' one-tailed display shades up to the upper bound from zero. Optional annotations show the
#' stored global p-value and the local intervals, coloured by deficit or excess and
#' labelled by p-value symbols. The plot reads the fitted result without repeating the
#' simulations. New local values use maximum-cluster adjustment under the complete null;
#' cluster endpoints remain descriptive, not confidence limits. Legacy fits retain their
#' original unadjusted values. Interpret local regions alongside the whole-curve test.
#' @return Draws the empirical density and pointwise comparison envelope. Optional annotations display the global p-value and local regions.
#' @param x A \emph{skyscapeR.sigTest} object created with \code{\link{randomTest}}
#' @param xlim Earlier name for xrange.
#' @param xrange Range of x-axis to plot.
#' @param title Title to add to the plot.
#' @param col Color of summed probability density. Defaults to blue.
#' @param show.pval Logical to control whether to print the global p-value. Default is TRUE.
#' @param show.local Show local regions and their significance symbols. An open circle marks p-values above 0.1. Default is FALSE
#' @param pal Two colours used to mark local regions. The default uses the first and fifth colours of the PRGn palette.
#' @param ... Additional arguments to be passed to \emph{plot}.
#' @inheritParams plot.skyscapeR.spd
#' @rdname plot.skyscapeR.sigTest
#' @export
#' @seealso \code{\link{randomTest}}, \code{\link{modelTest}}, \code{\link{print.skyscapeR.sigTest}}
#' @examples
#' az <- az.pdf(az = c(85, 95, 112), unc = 3, verbose = FALSE)
#' set.seed(42)
#' fit <- randomTest(az, nsims = 10, ncores = 1, .res = 2, verbose = FALSE)
#' plot(fit)
plot.skyscapeR.sigTest <- function(x, xrange, title=NULL, col='blue', show.pval=TRUE, show.local=FALSE, pal=c('#7B3294','#008837'), ..., obj=NULL, xlim) {
  if (!missing(xlim)) {
    if (!missing(xrange)) stop("Supply only one of xrange and xlim.")
    xrange <- xlim
  }
  par(mar=c(5, 4, 1, 1) + 0.1)
  # empirical spd
  spd <- x$data$empirical
  if (missing(xrange)) { xrange <- sort(spd$x[c(min(which(spd$y >= 1e-12)), max(which(spd$y >= 1e-12)))]) }
  if (show.local) { ylim <- c(-max(spd$y)*.05,max(spd$y)) } else { ylim <- c(0,max(spd$y))}
  if (x$metadata$coord == 'dec') { label <- 'Declination (\u00b0)' }
  if (x$metadata$coord == 'az') {
    label <- 'Azimuth (\u00b0)'
    xrange[1] <- max(xrange[1], 0)
    xrange[2] <- min(xrange[2], 360)
  }
  plot(spd$x, spd$y, type='l', main=title, xlab=label, ylab='Density', lwd=2, col=col, xlim=xrange, ylim=ylim, xaxs='i', yaxs='i', axes=FALSE, ...); box()
  xp <- c(spd$x, rev(spd$x))
  yp <- c(spd$y, rep(0, length(spd$x)))
  polygon(xp, yp, col=.with_alpha(col,0.5), border=NA)
  if (x$metadata$coord == 'dec') {
    axis(1, at=pretty(seq(par('usr')[1],par('usr')[2])))
    axis(1, at=0, labels = 0)
  } else if (x$metadata$coord == 'az') {
    axis(1, at=c(0, 45, 90, 135, 180, 225, 270, 315, 360))
  }
  scale <- mean(diff(pretty(seq(par('usr')[1],par('usr')[2]))))
  if (scale <= 1) { axis(1, at=seq(-90,360,0.1), lwd=0.2, labels=FALSE) }
  if (scale <= 2 & scale > 1) { axis(1, at=seq(-90,360,0.5), lwd=0.2, labels=FALSE) }
  if (scale <= 5 & scale > 2) { axis(1, at=seq(-90,360,1), lwd=0.5, labels=FALSE) }
  if (scale <= 20 & scale > 5) { axis(1, at=seq(-90,360,5), lwd=0.5, labels=FALSE) }
  if (scale > 10) { axis(1, at=seq(-90,360,10), lwd=0.5, labels=FALSE) }
  axis(2)

  # sigtest
  xp <- c(spd$x, rev(spd$x))
  if (x$metadata$tails==1) {
    yp <- c(x$data$null.hyp$CE.upper, rep(0, length(spd$x)))
  } else {
    yp <- c(x$data$null.hyp$CE.upper, rev(x$data$null.hyp$CE.lower))
  }
  polygon(xp, yp, col=.with_alpha('grey', 0.5), border=NA)

  .plot_density_targets(obj,x$metadata$coord)

  # global p-value
  if (show.pval) {
    if (x$metadata$global.pval == 0 ) {
      pval <- paste0("p-value < ", round(1/(x$metadata$nsims+1),4))
    } else {
      pval <- paste0('p-value = ', signif(x$metadata$global.pval,5))
    }
    bounds <- par("usr")
    text(bounds[2]-graphics::xinch(.07), bounds[4]-graphics::yinch(.07),
         pval, adj=c(1,1), cex=1.1, font=2)
  }

  # regions of significance
  if (show.local) {
    abline(0,0, lwd=1)
    aux <- as.matrix(x$metadata$local.pval[,2:4])
    if (NROW(aux)>0) {
      for (i in 1:NROW(aux)) {
        if (x$metadata$local.pval[i,1] == '+') { col <- pal[2] } else { col <- pal[1] }
        if (aux[i,1] < xrange[1] & aux[i,2] > xrange[1]) { aux[i,1] <- xrange[1] }
        if (aux[i,2] > xrange[2] & aux[i,1] < xrange[2]) { aux[i,2] <- xrange[2] }
        xp <- c(aux[i,1], aux[i,1], aux[i,2], aux[i,2])
        yp <- c(ylim[1], 0, 0, ylim[1])
        polygon(xp, yp, col=.with_alpha(col, .4), border=NA)
        symbol <- pval2stars(aux[i,3])
        centre <- mean(aux[i,1:2])
        if (identical(symbol,"\u25cb")) {
          graphics::points(centre,ylim[1]/2,pch=1,cex=.75)
        } else {
          text(centre,ylim[1]/2,labels=symbol,cex=.9,srt=90)
        }
      }
    }
  }

  mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'), side=3, at=par('usr')[2], cex=0.5, adj=1)

}


#' Plot seasonal stellar visibility
#' @description Compare the seasonal visibility of one or more stars at the modelled site and epoch. Each row shows the stellar name, designation and epoch. To compare one star across epochs, compute a phase object for each year and pass the second object to plot as its second argument.
#' @details Read each colour as a daily model classification: visible rising and setting, rising
#' only, setting only, visible without a qualifying near-horizon passage, or invisible. The
#' function plots the stored phase sequence on a 365-day scale and can label seasons or
#' display a legend. Supplying a second object or a list compares sequences in rows;
#' optional sorting uses the stored apparent magnitude or seasonal onset. Changing the
#' display does not recalculate the visibility thresholds. Use \code{plot(first, list(second, third))} to compare results.
#' @return Draws a seasonal visibility diagram on the active graphics device, using the phase classifications in the input object.
#' @param x A star phase result.
#' @param y A second object of \emph{skyscapeR.starphases} format or
#' a list of \emph{skyscapeR.starphases} objects for plotting (see Examples).
#' @param show.grid Logical to control whether to show a grid of months.
#' Default is FALSE
#' @param show.legend Show a legend for seasons present among the plotted
#' stars (see example below). Default is TRUE.
#' @param pal Colour palette for displaying phases/seasons.
#' @param ... Additional arguments to be passed to \emph{plot}.
#' @param show.labels Logical. Mark annual events on a single-star strip. Defaults to TRUE.
#' @param sort.by For multiple stars, sort by apparent magnitude ("mag"), rising ("rise") or setting ("set").
#' @rdname plot.skyscapeR.starphases
#' @export
#' @seealso \code{\link{star.phases}}
#' @examples
#' if (requireNamespace("swephRdata", quietly = TRUE)) {
#'   # Plot the seasonality of Aldebaran for 4001 BCE:
#'   s1 <- star.phases('Aldebaran',-4000, c(35,-8, 200))
#'   plot(s1)
#'   # Plot seasonality of several stars
#'   s2 <- star.phases('Alnilam',-4000, c(35,-8, 200))
#'   s3 <- star.phases('Elnath',-4000, c(35,-8, 200))
#'   plot(s1, s2)
#'   plot(s1, list(s2,s3))
#'   plot(s1, show.legend=FALSE, show.grid=TRUE)
#' }
plot.skyscapeR.starphases = function(x, y, show.labels=TRUE, show.legend=TRUE, show.grid=FALSE, sort.by, pal, ...) {
  if (missing(pal)) {
    # First four ColorBrewer Accent colours, followed by the invisibility grey.
    pal <- c('#7FC97F','#BEAED4','#FDC086','#FFFF99','#e0e0e0')
  }
  code <- c('RS', 'R', 'S', 'V', 'I')
  seasons <- c("Rise and Set","Rise Only","Set Only","Curtailed Passage","Arising and Lying Hidden")
  events <- c('acronycal rising','heliacal setting','acronycal setting','heliacal rising')
  if (inherits(x, "skyscapeR.starphases")) x <- list(x)
  if (!missing(y)) x <- c(x, if (inherits(y,"skyscapeR.starphases")) list(y) else y)
  if (!is.list(x) || !length(x) || !all(vapply(x,inherits,logical(1),"skyscapeR.starphases")))
    stop("Supply one or more star.phases() results.")
  if (!missing(sort.by)) {
    sort.by <- match.arg(sort.by, c("mag","rise","set"))
    key <- vapply(x, function(z) {
      if (sort.by == "mag") return(z$metadata$star$app.mag)
      wanted <- if (sort.by == "rise") "R" else "S"
      idx <- which(z$data$phase == wanted)
      if (length(idx)) z$data$day[idx[1]] else NA_real_
    }, numeric(1))
    x <- x[order(key, na.last=TRUE)]
  }
  n <- length(x)
  labels <- vapply(x, function(z) {
    label <- z$metadata$star$name
    if (nzchar(label)) label else z$metadata$star$Bayer
  }, character(1))
  if (!is.null(names(x))) {
    use <- nzchar(names(x))
    labels[use] <- names(x)[use]
  }
  bayer <- vapply(x,function(z) z$metadata$star$Bayer,character(1))
  bayer[bayer==labels] <- ""
  epoch <- rev(vapply(x, function(z) {
    year <- z$metadata$year
    if (length(year) == 1L) as.character(BC.AD(year)) else ""
  }, character(1)))
  labels <- rev(labels)
  bayer <- rev(bayer)
  x <- rev(x)
  old <- par(c("mai","mgp","xpd"))
  on.exit(par(old),add=TRUE)
  width <- par("fin")[1]
  narrow <- width < 6.2
  label.size <- .85
  # Wrap at word boundaries before shrinking to fit an unusually long word.
  wrap.name <- function(label) {
    words <- strsplit(label," ",fixed=TRUE)[[1]]
    rows <- ""
    for (word in words) {
      last <- length(rows)
      trial <- trimws(paste(rows[last],word))
      if (nzchar(rows[last]) && strwidth(trial,units="inches",font=2,cex=.85) > width*.29)
        rows <- c(rows,word)
      else rows[last] <- trial
    }
    paste(rows,collapse="\n")
  }
  labels <- vapply(labels,wrap.name,character(1))
  name.width <- max(strwidth(labels,units="inches",font=2),strwidth(c(bayer,epoch),units="inches"))
  label.size <- min(label.size,width*.30/max(name.width,.01))
  present <- code %in% unlist(lapply(x,function(z) as.character(z$data$phase)),use.names=FALSE)
  legend.labels <- seasons[present]
  legend.fill <- pal[which(present)]
  show.legend <- show.legend && length(legend.labels)>0L
  legend.width <- max(.01,strwidth(legend.labels,units="inches"))
  legend.size <- min(.75,width*.27/max(legend.width,.01))
  right <- if (show.legend && !narrow) legend.width*legend.size+.35 else .40
  below <- narrow && show.legend
  if (below) legend.size <- min(.75,(width/2-.35)/legend.width)
  par(mai=c(if (below) 1.05 else .62,name.width*label.size+.15,.32,right),
      mgp=c(2,.85,0),xpd=FALSE)
  height <- par("pin")[2]
  pad <- max(0,(height/min(.95,height/n)-n)/2)
  plot(NA_real_,NA_real_,xlim=c(.5,365.5),ylim=c(-pad,n+pad),main="",xlab="",ylab="",
       axes=FALSE,xaxs="i",yaxs="i",...)
  aux <- seq(0,365,365/12); aux[1] <- 1
  axis(1, pos=0, at=aux[c(1,4,7,10,13)], cex.axis=min(.8,par('pin')[1]/4), labels=c('December\nSolstice','March\nEquinox','June\nSolstice','September\nEquinox','December\nSolstice'), lwd.ticks=2)
  axis(1, pos=0, at=aux, labels=NA)
  axis(3, pos=n, at=aux[c(1,4,7,10,13)], labels=NA, lwd.ticks = 2)
  axis(3, pos=n, at=aux, labels=NA)

  for (k in 1:n){
    for (i in 1:length(code)) {
      ind <- which(x[[k]]$data$phase == code[i])
      if (!length(ind)) next
      ind.bb <- split(ind,cumsum(c(1,abs(diff(ind))>1)))

      for (j in 1:NROW(ind.bb)) {
        ind.i <- ind.bb[[j]]
        if (length(ind.i)>0) {
          x.poly <- c(ind.i[1]-.5,tail(ind.i,1)+.5,tail(ind.i,1)+.5,ind.i[1]-.5)
          y.poly <- c(k-1,k-1,k,k)
          polygon(x.poly, y.poly, col=pal[i], border=NA)
        }
      }
    }
  }
  row.height <- min(.95,height/n)
  lines.per.name <- lengths(strsplit(labels,"\n",fixed=TRUE))
  label.size <- min(label.size,.85*row.height/(max(lines.per.name+nzchar(bayer)+nzchar(epoch))*par("csi")))
  mid <- seq_len(n)-.5
  line.height <- graphics::yinch(par("csi")*label.size)
  name.x <- .5-graphics::xinch(.06)
  for (i in seq_len(n)) {
    parts <- strsplit(labels[i],"\n",fixed=TRUE)[[1]]
    yy <- mid[i]+((length(parts)+nzchar(bayer[i])+nzchar(epoch[i])-1)/2-seq_along(parts)+1)*line.height
    text(name.x,yy,parts,adj=c(1,.5),font=2,cex=label.size,xpd=NA)
    if (nzchar(bayer[i]))
      text(name.x,tail(yy,1)-line.height,bayer[i],adj=c(1,.5),cex=label.size,xpd=NA)
    if (nzchar(epoch[i]))
      text(name.x,tail(yy,1)-(1+nzchar(bayer[i]))*line.height,epoch[i],
           adj=c(1,.5),cex=label.size,xpd=NA)
  }

  for (i in 1:n) {
    lines(c(1,365),c(i,i))
  }
  lines(c(1,1),c(0,n))
  lines(c(365,365),c(0,n))
  lines(c(1,365),c(0,0))
  lines(c(1,365),c(n,n))

  if (show.grid) {
    for (i in 1:length(aux)) {
      lines(c(aux[i],aux[i]),c(0,n), col='black', lwd=0.3, lty=3)
    }
  }

  if (show.legend) {
    if (below) {
      legend(.5-graphics::xinch(par("mai")[2]-.05),-graphics::yinch(.48),
             legend=legend.labels,fill=legend.fill,cex=legend.size,ncol=min(2,length(legend.labels)),xpd=NA,bty="n",
             x.intersp=.55,y.intersp=.9)
    } else {
      legend(365.5,n/2,yjust=.5,xpd=NA,legend=legend.labels,fill=legend.fill,cex=legend.size,
             bty="n",x.intersp=.55,y.intersp=.9)
    }
  }
  if (n!=1L) return(invisible(NULL))
  x <- x[[1]]
  event.size <- min(.8,min(.95,height)/max(strwidth(events,units="inches"))*.92)
  event.gap <- graphics::xinch(.04)
  if (show.labels & (x$metadata$type == 'curtailed passage' | x$metadata$type == 'arising and lying hidden' | x$metadata$type == 'dual phase')) {
    ind <- which(x$data$phase == code[4])
    if (length(ind) > 0) {
      lines(rep(min(ind)-1,2), c(0,1), lty=2, lwd=1.5, col='black')
      text(min(ind)-1-event.gap,0.5, events[1], cex=event.size, font=1, srt=90 ,adj=c(.5,0))

      lines(rep(max(ind)+1,2), c(0,1), lty=2, lwd=1.5, col='black')
      text(max(ind)+1+event.gap,0.5, events[2], cex=event.size, font=1, srt=90 ,adj=c(.5,1))
    }

    ind <- which(x$data$phase == code[5])
    if (length(ind) > 0) {
      lines(rep(min(ind)-1,2), c(0,1), lty=2, lwd=1.5, col='black')
      text(min(ind)-1-event.gap,0.5, events[3], cex=event.size, font=1, srt=90 ,adj=c(.5,0))

      lines(rep(max(ind)+1,2), c(0,1), lty=2, lwd=1.5, col='black')
      text(max(ind)+1+event.gap,0.5, events[4], cex=event.size, font=1, srt=90 ,adj=c(.5,1))
    }
  }
}


#' Print global and local significance results
#' @description Report how strongly the observed orientation distribution departs from the fitted null
#' model, overall and within local angular intervals.
#' @details Start with the stored global Monte Carlo result, then examine the direction, limits and
#' local p-value of each reported interval. The summary gives the tail setting, confidence
#' level and number of simulations. It counts unequal observations at constant null density
#' separately and reports when there are no local regions. The method formats a saved
#' result on the console without recomputing it. Read interval coordinates in the units and
#' coordinate system of the fitted test.
#' @return Prints the test summary and global and local results to the console.
#' @param x Object of \emph{skyscapeR.sigTest} format.
#' @param ... Additional arguments to be passed to \emph{print}.
#' @export
#' @seealso \code{\link{randomTest}}, \code{\link{modelTest}}, \code{\link{plot.skyscapeR.sigTest}}
#' @examples
#' az <- az.pdf(az = c(85, 95, 112), unc = 3, verbose = FALSE)
#' set.seed(42)
#' fit <- randomTest(az, nsims = 10, ncores = 1, .res = 2, verbose = FALSE)
#' print(fit)
print.skyscapeR.sigTest <- function (x, ...) {
  cat("\n*** Results of Significance Test ***\n\n")
  cat(paste0(x$metadata$tails,'-tailed test at ',x$metadata$conf*100,'% confidence, based on ', x$metadata$nsims, ' simulations.\n'))

  if (x$metadata$global.pval == 0) {
    p.value <- paste0("< ", round(1/(x$metadata$nsims+1),3))
  } else { p.value <- x$metadata$global.pval }
  cat(paste0('p-value: ', format(p.value,digits=5,trim=TRUE), ' (',pval2stars(p.value),')\n'))
  if (NROW(x$data$constant.departures) > 0L) {
    cat("Grid points differing from a constant simulated density:", nrow(x$data$constant.departures),
        if (identical(x$metadata$inference,"pooled-exchangeable-v1"))
          "(included in pooled scoring; see $data$constant.departures).\n" else
          "(excluded from legacy scoring; see $data$constant.departures).\n")
  }
  if (identical(x$metadata$local.inference,'maximum-cluster-v1'))
    cat('local cluster p-values (maximum-cluster adjusted; boundaries descriptive):\n') else
    cat('local p-values (legacy exploratory; unadjusted interval search):\n')
  if (NROW(x$metadata$local.pval) == 0L) {
    cat('      No local regions detected.\n')
  }
  coord <- if (identical(x$metadata$coord, "az")) "az" else "dec"
  for (i in seq_len(NROW(x$metadata$local.pval))) {
    if (x$metadata$local.pval[i,]$p.value == 0) {
      p.value <- paste0("< ", round(1/(x$metadata$nsims+1),3))
    } else { p.value <- x$metadata$local.pval[i,]$p.value }
    cat(paste0('      ',x$metadata$local.pval[i,]$type,'  ',coord,' range [',round(x$metadata$local.pval[i,]$start,2), ', ',round(x$metadata$local.pval[i,]$end,2),'] :: p-value: ', format(p.value,digits=5,trim=TRUE), ' (',pval2stars(p.value),')\n'))
  }
}

.plot_target_bands <- function(groups,unit,show.obj.label) {
  # objects
  for (obj in groups) {
    if (unit == 'Azimuth') { warning('Celestial overlays require declination coordinates in this plot.',call.=FALSE) } else {
      for (i in 1:obj$n) {
        if (length(obj$epoch)==1) {
          abline(v=obj$decs[i], col=obj$col[i], lwd=obj$lwd[i], lty=obj$lty[i])
          if (show.obj.label) { text(obj$decs[i], .95*par('usr')[4], colnames(obj$decs)[i], col=obj$col[i], pos=4, offset=0.2, cex=0.7) }
        } else {
          xp <- c(obj$decs[1:2,i], rev(obj$decs[1:2,i]))
          yp <- rep(par('usr')[3:4], each=2)
          polygon(xp, yp, border=obj$col[i], col=.with_alpha(obj$col[i],.3),
                  lty=if (is.null(obj$lty)) 1 else obj$lty[i],
                  lwd=if (is.null(obj$lwd)) 1 else obj$lwd[i])
          if (show.obj.label) { text(mean(obj$decs[1:2,i]), .95*par('usr')[4], colnames(obj$decs)[i], col=obj$col[i], pos=4, offset=0.2, cex=0.7) }
        }
      }
    }
  }
}
