# Set the opacity while preserving the RGB channels and any colour names.
.with_alpha <- function(col, alpha) {
  channels <- grDevices::col2rgb(col, alpha=TRUE) / 255
  grDevices::rgb(channels[1,], channels[2,], channels[3,], alpha=alpha,
                 names=colnames(channels))
}

# At a fixed horizon altitude, declination is monotonic in cos(azimuth).
# Intersect in cosine space so a partly circumpolar interval ends at the meridian.
.compass_intervals <- function(decs, latitude, altitude) {
  rad <- pi/180
  scale <- cos(latitude*rad)*cos(altitude*rad)
  if (abs(scale) < 1e-14) return(NULL)
  limits <- (sin(range(decs)*rad) - sin(latitude*rad)*sin(altitude*rad))/scale
  if (limits[2] <= -1 || limits[1] >= 1) return(NULL)
  east <- rev(acos(pmax(-1, pmin(1, limits))))/rad
  rbind(rise=east, set=360-rev(east))
}

#' Plot orientation measurements on a compass grid
#' @description Inspect the directions of structural axes on a compass and compare them with candidate
#' rise or set bearings.
#' @details Measurements are drawn as spokes or rim ticks. North is at the top and
#' bearings increase clockwise; repeated bearings can overplot. A data frame supplies the
#' True.Azimuth column. For optional targets the function converts zero apparent altitude to geometric altitude
#' using the session atmosphere. It solves for rising and setting bearings on a level
#' horizon and shows the part of each declination interval that reaches it.
#' @seealso \code{\link{plotBars}}, \code{\link{orbit}}, \code{\link{sketchSky}}, \code{\link{planisphere}}
#' @return Draws measured bearings and optional target directions on the active graphics device. North is at the top and azimuth increases clockwise.
#' @param az Azimuths in degrees clockwise from north. Can be omitted if \emph{obj} is given.
#' @param style Analytical display style, either \code{"publication-colour"} (default) or
#' \code{"publication-mono"}. Monochrome uses grayscale and preserves line styles.
#' @param marks Draw measured orientations as equal-length \code{"spokes"} (default) or short \code{"rim"} ticks.
#' @param obj.position Draw target reference lines \code{"inside"} the circle or just \code{"outside"} its rim.
#' In both positions labels lie outside the circle. The first line names the target and the second gives rise or set.
#' Range endpoints share a label. Non-crossing targets have no reference bearing.
#' @param show.legend With TRUE, show the count and enabled target labels. FALSE centres the bare compass,
#' including its cardinal scale, measured bearings and target lines.
#' @param col Single color or color palette to use for plotting measurements.
#' @param lwd Line width for measured bearings. Defaults to 1.5.
#' @param lty Line type to plot measurements. Defaults to 1.
#' @param obj A \emph{skyscapeR.object} from \code{\link{skyObject}}, or
#' a list of these objects. Each group retains its own epochs and styles
#' for displaying the azimuths of celestial objects. Note that this assumes a single
#' location per group and a flat horizon of zero degrees. For a year range,
#' the target colour fills the interval of rising or setting bearings between
#' the stored declination limits. A zero-width interval is drawn as a line.
#' @param show.obj.labels Logical to control whether to display celestial objects names.
#' Defaults to TRUE.
#' @section Validation:
#' W01 checked compass bearings and the inner and outer edges of celestial year-range sectors against their spherical reference bearings. See the \href{https://stonexsky.github.io/skyscapeR/manual/#w01-plot-geometry}{W01 check in the technical manual}.
#' @export
#' @import utils stats
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
plotAzimuth <- function(az, obj, show.obj.labels=TRUE, col=NULL, lty=1, lwd=1.5,
                        marks=c("spokes","rim"), obj.position=c("inside","outside"),
                        show.legend=TRUE, style="publication-colour") {
  marks <- match.arg(marks)
  obj.position <- match.arg(obj.position)
  theme <- .plot_style(style)
  atmosphere <- skyscapeR.vars()[c("refraction","atm","temp")]
  groups <- if (missing(obj)) list() else .sky_groups(obj)
  if (!missing(az) && is.data.frame(az) && !"True.Azimuth" %in% names(az))
    stop("The measurement table needs a True.Azimuth column.",call.=FALSE)
  measured <- if (missing(az)) numeric() else if (is.data.frame(az)) az$True.Azimuth else az
  if (!is.numeric(measured) && length(measured)) stop("Supply numeric azimuths or a True.Azimuth column.")
  targets <- list()
  for (group.index in seq_along(groups)) {
    g <- groups[[group.index]]
    site <- if (inherits(g$loc,"skyscapeR.horizon")) g$loc$metadata$georef else g$loc
    if (!is.numeric(site) || length(site)<2L || any(!is.finite(head(site,3L))) || abs(site[1])>90)
      stop("Supply a site in each skyObject used on a compass.",call.=FALSE)
    elevation <- if (length(site) >= 3L) site[3] else 0
    pressure <- atmosphere$atm
    if (pressure == 0) pressure <- 1013.25*(1-.0065*elevation/288)^5.255
    altitude <- if (atmosphere$refraction) swephR::swe_refrac_extended(
      0, elevation, pressure, atmosphere$temp, .0065, 1L)$return else 0
    for (i in seq_len(g$n)) {
      bounds <- .compass_intervals(g$decs[,i], site[1], altitude)
      if (is.null(bounds)) next
      targets[[length(targets)+1L]] <- data.frame(
        lo=bounds[,1], hi=bounds[,2], az=rowMeans(bounds), event=rownames(bounds),
        name=.target_name(.target_id(g,i)), group=group.index,
        epoch=paste(vapply(g$epoch,function(t) as.character(BC.AD(t)),character(1)),collapse=" to "),
        col=as.character(g$col[i]), lty=as.character(g$lty[i]), lwd=g$lwd[i],
        stringsAsFactors=FALSE)
    }
  }
  refs <- do.call(rbind,targets)
  if (NROW(refs)) {
    duplicate.names <- names(which(vapply(split(refs$group,refs$name),function(ix) length(unique(ix))>1,logical(1))))
    use <- refs$name %in% duplicate.names
    refs$name[use] <- paste(refs$name[use],refs$epoch[use])
  }
  old <- par(c("mai","xpd","pty")); on.exit(par(old))
  par(mai=c(.08,.07,if (show.legend) .35 else .08,.07), pty="m", xpd=FALSE)
  plot.new()
  labelled <- show.legend && show.obj.labels && NROW(refs)>0L
  radius <- 1.24
  if (labelled) {
    # Measure text in physical units before establishing the equal-angle frame.
    fraction <- min(.23,max(strwidth(refs$name,units="inches",cex=.68))/par("pin")[1])
    radius <- 1.26/(1-2*fraction)
  }
  lim <- c(-radius,radius)
  plot.window(lim,c(-1.26,1.26),asp=1,xaxs="i",yaxs="i")
  inside <- obj.position=="inside"
  inner <- if (inside) .12 else 1.01
  outer <- if (inside) .99 else 1.16
  if (NROW(refs)) {
    cols <- .plot_colour(refs$col,theme)
    for (i in which(refs$hi > refs$lo)) {
      arc <- seq(refs$lo[i], refs$hi[i], length.out=ceiling(refs$hi[i]-refs$lo[i])+1L)
      polygon(rbind(.polar_xy(arc,outer), .polar_xy(rev(arc),inner)),
              col=.with_alpha(cols[i],.34), border=NA)
    }
  }
  lines(.polar_xy(0:360),col=theme$muted,lwd=.8)
  at <- seq(0,345,15)
  a <- .polar_xy(at,ifelse(at%%45==0,.93,.97)); b <- .polar_xy(at)
  segments(a[,1],a[,2],b[,1],b[,2],col=theme$muted,lwd=.7)
  text(.polar_xy(seq(0,315,45),1.13),c("N","NE","E","SE","S","SW","W","NW"),
       cex=rep(c(.84,.68),4),col=theme$ink)
  if (length(measured)) {
    a <- .polar_xy(measured,if (marks=="rim") .73 else 0); b <- .polar_xy(measured)
    segments(a[,1],a[,2],b[,1],b[,2],col=.plot_colour(if (is.null(col)) theme$data else col,theme),lty=lty,lwd=lwd)
  }
  points(0,0,pch=16,cex=.36,col=theme$ink)
  if (NROW(refs)) {
    line.types <- .line_types(refs$lty)
    bounds <- c(refs$lo,refs$hi)
    index <- rep(seq_len(nrow(refs)),2L)
    keep <- c(rep(TRUE,nrow(refs)),refs$hi != refs$lo)
    a <- .polar_xy(bounds[keep],inner)
    b <- .polar_xy(bounds[keep],outer)
    index <- index[keep]
    segments(a[,1],a[,2],b[,1],b[,2],col=cols[index],lty=line.types[index],lwd=refs$lwd[index])
    b <- .polar_xy(refs$az,outer)
    if (labelled) {
      # Point each label to the midpoint of its range.
      side <- ifelse(refs$event=="rise",1,-1)
      for (s in c(-1,1)) {
        ix <- which(side==s)
        if (!length(ix)) next
        label <- paste(refs$name[ix],refs$event[ix],sep="\n")
        labs <- unique(label)
        centres <- vapply(labs,function(k) mean(b[ix[label==k],2]),numeric(1))
        cx <- min(.68, .79*diff(par("usr")[1:2])/2/max(strwidth(labs),.01),
                  1.8/max(.01,2.5*strheight("M")*length(labs)))
        h <- strheight("M")*cx
        pos <- .label_rows(centres,-1.10,1.10,3.3*h)
        # The text anchor on either side is beyond the circle for every bearing.
        label.x <- s*1.23
        for (j in seq_along(labs)) {
          k <- ix[label==labs[j]]
          segments(b[k,1],b[k,2],label.x-s*.025,pos[j],col=cols[k],lwd=.65)
          text(label.x,pos[j],labs[j],adj=c(if(s>0) 0 else 1,.5),cex=cx,col=cols[k[1]])
        }
      }
    }
  }
  if (show.legend) {
    count <- sum(is.finite(measured))
    heading <- if (count) sprintf("%s measured bearing%s",count,if(count==1) "" else "s") else "Celestial bearings"
    text(0,1.36,heading,cex=.88,font=2,col=theme$ink,xpd=NA)
    if (length(groups)) {
      sites <- unique(vapply(groups,function(g) {
        if (inherits(g$loc,"skyscapeR.horizon")) .site_caption(g$loc$metadata$georef,g$loc$metadata$name)
        else .site_caption(g$loc)
      },character(1)))
      epochs <- unique(vapply(groups,function(g) paste(vapply(g$epoch,function(y) as.character(BC.AD(y)),
                                                character(1)),collapse=" to "),character(1)))
      lab <- c(if(length(sites)==1L) sites,if(length(epochs)==1L) epochs)
      if (length(lab)) text(0,1.24,paste("References:",paste(lab,collapse="  |  ")),
                            cex=.63,col=theme$muted,xpd=NA)
    }
  }
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
#' @param style Analytical display style, either \code{"publication-colour"} (default) or
#' \code{"publication-mono"}. Monochrome uses grayscale and preserves line styles.
#' @param geometry Use \code{"interval"} for sightline intervals with a centre point, or
#' \code{"bar"} for filled or outlined rectangles. Explicitly supplying the older \code{shade}
#' setting selects rectangles unless \code{geometry} is supplied.
#' @param show.legend Display target annotations when \code{show.obj.label} is also TRUE.
#' @param col Measurement colour. NULL uses the current style.
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
#' @param show.obj.label Logical. Label celestial overlays on the declination plot. Defaults to TRUE.
#' @section Validation:
#' W01 checked the endpoints and row labels of an azimuth interval spanning north. See the \href{https://stonexsky.github.io/skyscapeR/manual/#w01-plot-geometry}{W01 check in the technical manual}.
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
#' @section Coordinates and interpretation:
#' \code{unit} labels the measurement axis; it does not transform azimuths into
#' declinations. Celestial overlay values are declinations. Use them with
#' declination measurements, or convert measurements first using the appropriate
#' site and horizon. Bars show the supplied symmetric angular bounds, not fitted
#' confidence intervals. Use \code{\link{az.pdf}}, \code{\link{coordtrans}} and
#' \code{\link{spd}} when the full measurement distributions are needed.
plotBars <- function(val, unc, names, unit="Declination", obj, show.obj.label=TRUE,
                     xrange, yrange, sort=FALSE, shade=TRUE, mark=TRUE, col=NULL,
                     geometry=if (missing(shade)) "interval" else "bar",
                     show.legend=TRUE, style="publication-colour") {
  theme <- .plot_style(style)
  geometry <- match.arg(geometry,c("interval","bar"))
  groups <- if (missing(obj)) list() else .sky_groups(obj)
  if (!length(val) || !all(is.finite(val)) || !length(unc) %in% c(1,length(val)) ||
      !all(is.finite(unc)) || any(unc<0)) stop("Supply finite values and one nonnegative width per value, or a common width.")
  if (unit=="Declination" && any(abs(val)>90)) unit <- "Azimuth"
  unc <- rep_len(unc,length(val))
  labels <- if (missing(names)) as.character(seq_along(val)) else as.character(names)
  if (length(labels)!=length(val)) stop("Supply one name per value.")
  cols <- rep_len(.plot_colour(if (is.null(col)) theme$data else col,theme), length(val))
  idx <- if (sort) order(val,decreasing=TRUE) else seq_along(val)
  val <- val[idx]; unc <- unc[idx]; labels <- labels[idx]; cols <- cols[idx]
  if (missing(xrange)) xrange <- range(val-unc,val+unc)+c(-5,5)
  if (missing(yrange)) yrange <- c(.5,length(val)+.5)
  old <- par(c("mai","mgp","pty")); on.exit(par(old))
  par(mai=c(.62,.18,.18,.16),mgp=c(1.8,.55,0),pty="m")
  width <- max(strwidth(labels,units="inches"),.01)
  size <- min(.82,.35*par("fin")[1]/width)
  par(mai=c(.62,width*size+.14,.18,.16))
  plot.default(NA_real_,NA_real_,type="n",xlim=xrange,ylim=yrange,
               xlab=paste0(if (unit=="Azimuth") "True azimuth" else unit," (\u00b0)"),ylab="",
               axes=FALSE,xaxs="i",yaxs="i")
  xt <- .plot_ticks(xrange)
  abline(v=xt,col=theme$grid,lwd=.65)
  .plot_target_bands(groups,unit,show.legend && show.obj.label,theme)
  at <- seq_along(val)
  if (geometry=="interval") {
    segments(val-unc,at,val+unc,at,col=cols,lwd=2)
    if (mark) points(val,at,pch=16,col=cols,cex=.64)
  } else {
    rect(val-unc,at-.23,val+unc,at+.23,col=if(shade) .with_alpha(cols,.36) else NA,
         border=if(shade) NA else cols)
    if (mark) segments(val,at-.23,val,at+.23,col=cols,lwd=1.4)
  }
  .plot_axis(1,xt,theme=theme)
  rows <- at[at>=min(yrange) & at<=max(yrange)]
  axis(2,at=rows,labels=labels[rows],las=1,lwd=0,lwd.ticks=0,gap.axis=-1,
       cex.axis=size,col.axis=theme$ink)
  invisible(NULL)
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
#' @param style Frame and silhouette style. A sketch supplies its stored style by default. Without a
#' sketch the default is \code{"publication-colour"}. Other choices are \code{"natural"} and
#' \code{"publication-mono"}. Supplied panorama pixels and sky objects retain their own rendering.
#' @param show.legend Display the key for celestial tracks. The angular axes are always drawn.
#' @param col.ground Silhouette colour. NULL uses a dark muted foreground for natural skies, pale grey-green for publication colour or grey for monochrome.
#' @details A horizon with an attached DEM, PeakFinder or photographic PNG displays the image by default. Add show.skyline = TRUE to inspect its numerical profile, or use.image = FALSE to draw the silhouette alone. With xrange supplied, azimuth labels default to degrees. A partial photograph initially opens on its recorded sector. An image-backed horizon or a supplied sky uses a 2:1 plot box. Copy the PNG with the saved horizon and use relinkHor() after moving it. With automatic altitude limits, a 30-degree azimuth sector shows 15 degrees of altitude.
#' @param sky A colour string or a background from \code{\link{sketchSky}} or \code{\link{stellariumSky}}.
#' @param show.skyline Draw the numerical skyline. By default the line appears only without a panorama PNG.
#' @param fastLoad Use a sampled PNG for quick plots. New photographic imports store a preview. For older horizons the first preview loads the image. The numerical skyline and celestial tracks retain their precision.
#' @param use.image Use an attached external panorama PNG when available. Default TRUE.
#' Missing or changed files warn and fall back to the numerical silhouette.
#' @param col.sky Earlier name for a solid sky colour. Prefer sky.
#' @param ... Additional arguments to be passed to \emph{plot}.
#' @rdname plot.skyscapeR.horizon
#' @section Validation:
#' W01 checked the skyline ordinate and azimuth viewport on flat and raised profiles, including natural and publication sky overlays. See the \href{https://stonexsky.github.io/skyscapeR/manual/#w01-plot-geometry}{W01 check in the technical manual}.
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
#' @section Layers and coordinate consistency:
#' A reusable sky is calculated for its own stored location and instant. Pair it
#' with the horizon used to create it; drawing a saved sky does not recalculate
#' its positions for another observer. Celestial reference tracks from
#' \code{\link{skyObject}} use the plotted horizon's location and session
#' atmosphere. Keep that atmosphere consistent with the stored sketch when
#' comparing a track with a disc or point symbol.
#'
#' The foreground masks sky layers geometrically. It does not make the sky object
#' an extinction or visibility calculation through the terrain. Image calibration
#' and numerical skyline agreement remain the responsibility of the horizon
#' source. Azimuth-altitude plotting is a cylindrical display, with increasing
#' horizontal distortion towards the zenith. Retain consistent ranges when
#' comparing styles, epochs or constellation traditions.
plot.skyscapeR.horizon <- function(x, show.az=(!missing(xrange) || !missing(xlim)), xrange, yrange, obj, col.ground=NULL, col.sky, use.image=TRUE, ..., show.skyline=NULL, sky=NULL, xlim, ylim, fastLoad=FALSE, show.legend=TRUE, style=NULL) {
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
  if (is.null(style)) style <- if (inherits(sky,"skyscapeR.sketch")) sky$style else "publication-colour"
  if (identical(style,"publication")) style <- "publication-mono"
  style <- match.arg(style,c("natural","publication-colour","publication-mono"))
  theme <- .plot_style(if (style=="natural") "publication-colour" else style)
  if (is.null(col.ground)) col.ground <- switch(style,natural="#313D37","publication-colour"="#A5B2A8","publication-mono"="#B2B2B2")
  col.ground <- .plot_colour(col.ground,theme)
  force(show.az)
  if (missing(xrange)) {
    coverage <- x$metadata$source$coverage
    xrange <- if (identical(x$metadata$source$provider,"Panorama") && diff(coverage)<360-1e-8)
      coverage else c(0,360)
  }
  if (xrange[2]<xrange[1]) xrange[2] <- xrange[2]+360

  old.style <- par("mgp")
  on.exit(par(mgp=old.style), add=TRUE)
  par(mar=c(3.4,3.7,if(show.legend && length(groups)) 2.8 else .65,.85), mgp=c(2.0,.55,0))
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
      lines(orb$az+360*turn,orb$alt,col=.plot_colour(obj$col[i],theme),lty=obj$lty[i],lwd=obj$lwd[i])
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
  .plot_x_baseline(lwd=.5)

  .plot_axis(2,.plot_ticks(yrange),theme=theme)
  mtext("Altitude (\u00b0)",side=2,line=2.25,cex=.86)
  mtext(if(show.az) "True azimuth (\u00b0)" else "True azimuth",side=1,line=1.95,cex=.86)
  if (show.legend && length(groups)) {
    labs <- unlist(lapply(groups,function(g) .target_name(colnames(g$decs))),use.names=FALSE)
    cols <- unlist(lapply(groups,function(g) .plot_colour(g$col,theme)),use.names=FALSE)
    types <- unlist(lapply(groups,function(g) .line_types(g$lty)),use.names=FALSE)
    .density_key(labs,cols,types,list(),theme)
  }
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
#' probability mass. Every disjoint component is listed in the HDR readout. Standalone curves have a fixed label column.
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
#' @param style Analytical display style, either \code{"publication-colour"} (default) or
#' \code{"publication-mono"}. Monochrome uses grayscale and preserves line styles.
#' @param show.legend Show the HDR mass above stacked densities. The linked diagnostic always identifies
#' its component panels and reports the numerical interval or intervals.
#' @param hdr High density region to highlight. Defaults to 0.954
#' @param show.az Deprecated and ignored, with a warning when supplied. Provide both
#' \code{az} and \code{hor} for the combined panels, or omit both for distributions alone.
#' @param xlim Earlier name for xrange.
#' @param xrange Range of x-axis to plot.
#' @param col Fill colour for highest-density regions. NULL uses the current style.
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
#' @section Validation:
#' In W01 the shaded components matched the two intervals returned by hpdi() for a bimodal density. See the \href{https://stonexsky.github.io/skyscapeR/manual/#w01-plot-geometry}{W01 check in the technical manual}.
#' @export
plot.skyscapeR.pdf <- function(x, index, hdr=0.954, show.az, xrange, col=NULL, ..., az, hor, xlim, show.legend=TRUE, style="publication-colour") {
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

  theme <- .plot_style(style)
  fills <- .plot_colour(if (is.null(col)) .with_alpha(theme$data,.35) else col,theme)
  if (combined) {
    curve <- .unwrap_north_pdf(az$data[[index]])
    angle <- curve$x
    out <- x$data[[index]]
    xl <- if (missing(xrange)) range(out$x) else xrange
    yl <- range(angle)
    hd <- hpdi(out,hdr)
    profile <- .coord_horizon(hh)
    bearing <- .profile_bearings(angle,profile$az)
    alt <- approx(profile$az,profile$alt,xout=bearing)$y
    unc <- approx(profile$az,profile$alt.unc,xout=bearing)$y
    dec <- vapply(-2:2,function(k) az2dec(angle,hh,alt+k*unc,refraction=settings$refraction,
                                           atm=settings$atm,temp=settings$temp),numeric(length(angle)))
    old <- par(c("mai","mgp","xpd")); on.exit(par(old))
    par(mgp=c(1.65,.42,0),xpd=FALSE)
    panel <- .plot_panels()
    panel(.105,.40,.252,.82)
    plot(curve$y,angle,type="n",xlim=c(0,max(curve$y)*1.05),ylim=yl,axes=FALSE,
         xlab="",ylab="True azimuth (\u00b0)",xaxs="i",yaxs="i",...)
    polygon(c(curve$y,numeric(length(angle))),c(angle,rev(angle)),col=.with_alpha(theme$data,.22),border=NA)
    lines(curve$y,angle,col=theme$data,lwd=1.3)
    yt <- .plot_ticks(yl,4); .plot_axis(2,yt,yt%%360,theme=theme,cex.axis=.76)
    mtext("Source PDF",side=3,line=.3,cex=.74,adj=0)
    panel(.295,.40,.955,.82)
    plot(NA_real_,NA_real_,xlim=xl,ylim=yl,xlab="",ylab="",axes=FALSE,xaxs="i",yaxs="i",...)
    abline(v=.plot_ticks(xl),h=yt,col=theme$grid,lwd=.65)
    polygon(c(dec[,1],rev(dec[,5])),c(angle,rev(angle)),col=.with_alpha(theme$data,.15),border=NA)
    polygon(c(dec[,2],rev(dec[,4])),c(angle,rev(angle)),col=.with_alpha(theme$data,.32),border=NA)
    lines(dec[,3],angle,col=theme$ink,lwd=1.2)
    mtext("Azimuth to declination",side=3,line=.3,adj=0,cex=.76)
    panel(.295,.14,.955,.34)
    plot(NA_real_,NA_real_,xlim=xl,ylim=c(0,1.08*max(out$y)),xlab="Declination (\u00b0)",
         ylab="",axes=FALSE,xaxs="i",yaxs="i",...)
    .draw_pdf_curve(out,hd,fills[1],theme$data)
    .plot_axis(1,.plot_ticks(xl),theme=theme,cex.axis=.8)
    panel(.035,.075,.26,.345)
    plot.new();plot.window(c(0,1),c(0,1),xaxs="i",yaxs="i")
    lab <- paste(c(paste0(format(100*hdr),"% HDR"),sprintf("%.2f\u00b0 to %.2f\u00b0",hd[,1],hd[,2])),collapse="\n")
    size <- min(.8,.96/max(.01,strwidth(lab)),.95/max(.01,strheight(lab)))
    text(.02,.55,lab,adj=c(0,.5),cex=size,col=theme$ink)
    panel(.075,.92,.955,.998)
    plot.new();plot.window(c(0,1),c(0,1),xaxs="i",yaxs="i")
    lab <- x$metadata$name[index]
    if (length(lab) && !is.na(lab)) text(0,.50,lab,adj=c(0,.5),font=2,cex=min(.95,.96/max(.01,strwidth(lab,font=2))))
    return(invisible(NULL))
  }
  coord <- x$metadata$coord
  limits <- if (missing(xrange)) {if(coord=="az") c(0,360) else c(-90,90)} else xrange
  if (coord=="az" && limits[2]<limits[1]) limits[2] <- limits[2]+360
  labels <- x$metadata$name[index]
  if (!length(labels)) labels <- as.character(index)
  labels[is.na(labels)] <- ""
  par(mai=c(.55,.2,if(show.legend) .38 else .15,.15),mgp=c(1.8,.5,0))
  w <- max(strwidth(labels,units="inches"),.01)
  size <- min(.82,.33*par("fin")[1]/w,par("pin")[2]/(length(index)*par("csi")))
  mai <- par("mai");mai[2] <- w*size+.12;par(mai=mai)
  plot(NA_real_,NA_real_,type="n",xlim=limits,ylim=c(0,length(index)),axes=FALSE,
       xlab=.angle_title(coord),ylab="",xaxs="i",yaxs="i",...)
  abline(v=.plot_ticks(limits),col=theme$grid,lwd=.65)
  for (i in seq_along(index)) {
    d <- x$data[[index[i]]]; hd <- hpdi(d,hdr)
    # Use the same peak for the outline and every shaded component of this row.
    height <- max(d$y)
    if (!is.finite(height) || height<=0) next
    d$y <- .76*d$y/height
    base <- i-1+.09
    segments(min(limits),base,max(limits),base,col=theme$grid,lwd=.7)
    for (off in .curve_offsets(d$x,coord,limits))
      .draw_pdf_curve(d,hd,rep_len(fills,length(index))[i],theme$data,base,off)
  }
  axis(2,at=seq_along(index)-.53,labels=labels,las=1,cex.axis=size,
       lwd=0,lwd.ticks=0,gap.axis=-1,col.axis=theme$ink)
  .density_axis(coord,limits,theme)
  if (show.legend) mtext(sprintf("%g%% highest-density region",100*hdr),side=3,line=.3,adj=0,cex=.78,col=theme$muted)
  invisible(NULL)
}

.draw_pdf_curve <- function(d,hd,fill,line,base=0,offset=0) {
  gaps <- diff(d$x)
  step <- if (any(gaps>0)) median(gaps[gaps>0]) else Inf
  cuts <- gaps>5*step
  for (rows in split(seq_len(nrow(d)),cumsum(c(TRUE,cuts)))) {
    x <- d$x[rows]; y <- d$y[rows]
    for (k in seq_len(nrow(hd))) {
      keep <- which(x>=hd[k,1] & x<=hd[k,2])
      if (length(keep)) polygon(c(x[keep],rev(x[keep]))+offset,
          c(y[keep]+base,rep(base,length(keep))),col=fill,border=NA)
    }
    lines(x+offset,y+base,col=line,lwd=1.3)
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
#' @param col Observed density colour. NULL uses the current style.
#' @param style Analytical display style, either \code{"publication-colour"} (default) or
#' \code{"publication-mono"}. Monochrome uses grayscale and preserves line styles.
#' @param show.legend Show the observed-curve and celestial-target key when target overlays are supplied.
#' @param shading Logical to control whether to color the entire distribution. Defaults to TRUE
#' @param ... Additional arguments to be passed to \emph{plot}.
#' @param obj A skyObject or list of target groups for a density or significance plot. Declination plots show single epochs as lines and epoch ranges as bands. Azimuth plots use each group's site latitude to mark rising and setting boundary bearings at a geometric level horizon.
#' @rdname plot.skyscapeR.spd
#' @section Validation:
#' In W01 the drawn line matched the stored coordinates and heights of the summed density. See the \href{https://stonexsky.github.io/skyscapeR/manual/#w01-plot-geometry}{W01 check in the technical manual}.
#' @export
plot.skyscapeR.spd <- function(x, xrange, yrange, title=NULL, col=NULL, shading=TRUE, ..., obj=NULL, xlim, ylim,
                               show.legend=TRUE, style="publication-colour") {
  if (!missing(xlim)) {
    if (!missing(xrange)) stop("Supply only one of xrange and xlim.")
    xrange <- xlim
  }
  if (!missing(ylim)) {
    if (!missing(yrange)) stop("Supply only one of yrange and ylim.")
    yrange <- ylim
  }
  coord <- x$metadata$coord
  d <- x$data
  limits <- .density_limits(d,coord,if (missing(xrange)) NULL else xrange)
  if (missing(yrange)) yrange <- c(0,max(d$y,1e-12,na.rm=TRUE)*1.04)
  theme <- .plot_style(style)
  cc <- .plot_colour(if (is.null(col)) theme$data else col,theme)
  key <- show.legend && length(obj)>0
  par(mar=c(3.4,4,if (key) 2.7 else if (is.null(title)) .7 else 2.0,.8),mgp=c(2.2,.6,0))
  plot(NA_real_,NA_real_,type="n",main=title,xlab=.angle_title(coord),ylab="Density",
       xlim=limits,ylim=yrange,axes=FALSE,xaxs="i",yaxs="i",...)
  yt <- .plot_ticks(yrange); abline(h=yt,col=theme$grid,lwd=.65)
  targets <- .density_targets(obj,coord,limits,theme)
  for (off in .curve_offsets(d$x,coord,limits)) {
    at <- d$x+off
    if (shading) polygon(c(at,rev(at)),c(d$y,numeric(length(at))),
                         col=.with_alpha(cc,.30),border=NA)
    lines(at,d$y,col=cc,lwd=1.65)
  }
  .density_axis(coord,limits,theme); .plot_axis(2,yt,theme=theme)
  if (key) .density_key("Summed density",cc,1,targets,theme)
  invisible(NULL)
}

# Common coordinates, target styles and legends for orientation densities.
.density_limits <- function(d, coord, xrange=NULL) {
  if (is.null(xrange)) {
    ok <- is.finite(d$y) & d$y >= 1e-12
    if (!any(ok)) ok <- is.finite(d$x)
    xrange <- range(d$x[ok])
  }
  if (length(xrange)!=2 || any(!is.finite(xrange))) stop("Supply two finite angular limits.")
  if (coord=="az" && xrange[2]<xrange[1]) xrange[2] <- xrange[2]+360
  if (diff(xrange)==0) xrange <- xrange+c(-.5,.5)
  xrange
}

# Draw separate periodic copies; never connect two pieces across a gap.
.curve_offsets <- function(x, coord, limits) {
  if (coord!="az") return(0)
  lo <- ceiling((min(limits)-max(x))/360)
  hi <- floor((max(limits)-min(x))/360)
  if (lo>hi) return(numeric())
  360*seq(lo,hi)
}

.density_axis <- function(coord, limits, theme, pos=NULL) {
  at <- .plot_ticks(limits)
  lab <- if (coord=="az") at%%360 else at
  .plot_axis(1,at,labels=lab,pos=pos,theme=theme)
}

.density_targets <- function(obj, coord, limits, theme) {
  groups <- .sky_groups(obj)
  keys <- list()
  for (g in groups) {
    if (coord=="az") {
      site <- if (inherits(g$loc,"skyscapeR.horizon")) as.numeric(g$loc$metadata$georef) else g$loc
      if (!is.numeric(site) || length(site)<2L || !is.finite(site[1]))
        stop("Supply a site in each skyObject for an azimuth plot.",call.=FALSE)
    }
    for (i in seq_len(g$n)) {
      cc <- .plot_colour(g$col[i],theme)
      dec <- g$decs[,i]
      band <- coord=="dec" && length(g$epoch)>1L
      if (coord=="az") {
        cosA <- sin(dec*pi/180)/cos(site[1]*pi/180)
        angle <- acos(cosA[is.finite(cosA) & abs(cosA)<=1])*180/pi
        at <- unique(c(angle,360-angle))
        for (v in at) for (off in .curve_offsets(v,coord,limits))
          abline(v=v+off,col=cc,lty=g$lty[i],lwd=g$lwd[i])
      } else if (band) {
        rect(min(dec),par("usr")[3],max(dec),par("usr")[4],col=.with_alpha(cc,.19),
             border=cc,lty=g$lty[i],lwd=.65)
      } else abline(v=dec,col=cc,lty=g$lty[i],lwd=g$lwd[i])
      id <- .target_id(g,i)
      ep <- paste(vapply(g$epoch,function(y) as.character(BC.AD(y)),character(1)),collapse=" to ")
      keys[[length(keys)+1L]] <- list(name=paste(.target_name(id),ep,sep="  "),col=cc,lty=g$lty[i])
    }
  }
  keys
}

.density_key <- function(labels, colours, types, targets, theme, fills=NULL) {
  if (length(targets)) {
    labels <- c(labels,vapply(targets,`[[`,character(1),"name"))
    colours <- c(colours,vapply(targets,`[[`,character(1),"col"))
    types <- c(as.character(types),vapply(targets,function(k) as.character(k$lty),character(1)))
    # Convert mixed numeric and named line types into legal string styles.
    types <- .line_types(types)
  }
  if (is.null(fills)) fills <- rep(NA,length(labels)) else fills <- rep_len(c(fills,rep(NA,length(targets))),length(labels))
  if (!length(labels)) return(invisible())
  area <- par("usr")
  size <- min(.72,.90*par("pin")[1]/max(.01,max(strwidth(labels,units="inches"))+.35))
  legend(area[1],area[4]+yinch(.035),legend=labels,col=colours,lty=types,lwd=1.45,
         pch=ifelse(is.na(fills),NA_integer_,15L),pt.cex=1.2,
         ncol=max(1,min(length(labels),floor(par("pin")[1]/max(.65,max(strwidth(labels,units="inches",cex=size))+.34)))),
         xjust=0,yjust=0,xpd=NA,bty="n",cex=size,text.col=theme$ink,x.intersp=.5,y.intersp=.9)
}

.line_types <- function(x) {
  types <- as.character(x); ind <- types %in% as.character(0:6)
  types[ind] <- c("blank","solid","dashed","dotted","dotdash","longdash","twodash")[as.integer(types[ind])+1]
  types
}

#' Plot an orientation significance test
#' @description Locate the orientations at which the observed summed density departs from the simulated
#' comparison envelope. New tests use an envelope pooled from the observation and simulations.
#' @details The observed curve and envelope use the azimuth or declination grid stored in the test
#' result. The two-tailed display shades between its lower and upper pointwise bounds; the
#' one-tailed display shades up to the upper bound from zero. Optional annotations show the
#' stored global p-value and the local intervals, coloured by deficit or excess and
#' labelled by p-value symbols. The observed and envelope curves jointly set the visible vertical range. The lower panel identifies adjusted and legacy unadjusted local results. The plot reads the fitted result without repeating the
#' simulations. New local values use maximum-cluster adjustment under the complete null;
#' cluster endpoints remain descriptive, not confidence limits. Legacy fits retain their
#' original unadjusted values. Interpret local regions alongside the whole-curve test.
#' @return Draws the empirical density and pointwise comparison envelope. Optional annotations display the global p-value and local regions.
#' @param x A \emph{skyscapeR.sigTest} object created with \code{\link{randomTest}}
#' @param xlim Earlier name for xrange.
#' @param xrange Range of x-axis to plot.
#' @param title Title to add to the plot.
#' @param col Observed density colour. NULL uses the current style.
#' @param style Analytical display style, either \code{"publication-colour"} (default) or
#' \code{"publication-mono"}. Monochrome uses grayscale and preserves line styles.
#' @param show.legend Display the density/envelope key and any local significance symbol key.
#' \code{show.pval} independently controls the numeric global p-value.
#' @param show.pval Logical to control whether to print the global p-value. Default is TRUE.
#' @param show.local Draw a separate panel on the same angular scale, with Excess and Deficit rows. Default FALSE. Symbols are ***, **, *, + and an open circle, with upper cutoffs 0.001, 0.01, 0.05 and 0.10. The open circle is centred at the same height as the stars and marks values above 0.10.
#' @param pal Deficit and excess colours in that order. NULL uses the current style.
#' @param ... Additional arguments to be passed to \emph{plot}.
#' @inheritParams plot.skyscapeR.spd
#' @rdname plot.skyscapeR.sigTest
#' @section Validation:
#' The tested public fits plotted the stored local interval midpoints and significance symbols. A wrapped interval has start greater than end. A plot spanning the full turn displays its two visible segments with the same adjusted probability. See \href{https://stonexsky.github.io/skyscapeR/manual/#s03-scoring-checks}{S03 validation in the technical manual}.
#'
#' W01 checked clipped local intervals, the pointwise envelope and the global probability label against the stored fit. See the \href{https://stonexsky.github.io/skyscapeR/manual/#w01-plot-geometry}{W01 check in the technical manual}.
#' @export
#' @seealso \code{\link{randomTest}}, \code{\link{modelTest}}, \code{\link{print.skyscapeR.sigTest}}
#' @examples
#' az <- az.pdf(az = c(85, 95, 112), unc = 3, verbose = FALSE)
#' set.seed(42)
#' fit <- randomTest(az, nsims = 10, ncores = 1, .res = 2, verbose = FALSE)
#' plot(fit)
plot.skyscapeR.sigTest <- function(x, xrange, title=NULL, col=NULL, show.pval=TRUE, show.local=FALSE,
                                   pal=NULL, ..., obj=NULL, xlim, show.legend=TRUE, style="publication-colour") {
  if (!missing(xlim)) {
    if (!missing(xrange)) stop("Supply only one of xrange and xlim.")
    xrange <- xlim
  }
  d <- x$data$empirical
  null <- x$data$null.hyp
  coord <- x$metadata$coord
  limits <- .density_limits(d,coord,if (missing(xrange)) NULL else xrange)
  theme <- .plot_style(style)
  cc <- .plot_colour(if (is.null(col)) theme$data else col,theme)
  pal <- .plot_colour(if (is.null(pal)) c(theme$lower,theme$target) else pal,theme)
  if (length(pal)!=2L) stop("Supply deficit and excess colours in pal.")
  low <- if (x$metadata$tails==1L) numeric(length(d$x)) else null$CE.lower
  offsets <- .curve_offsets(d$x,coord,limits)
  peak <- 0
  for (off in offsets) {
    at <- d$x+off
    inside <- at>=min(limits) & at<=max(limits)
    for (y in list(d$y,null$CE.upper,low)) {
      edges <- if (length(at)>1L) approx(at,y,xout=limits,ties="ordered")$y else y
      peak <- max(peak,y[inside],edges,na.rm=TRUE)
    }
  }
  top <- max(1e-12,1.05*peak)
  ht <- par("fin")[2]
  bottom <- if (show.local) min(ht*.44,if (show.legend) 1.52 else 1.23) else min(ht*.23,.62)
  upper <- min(ht*.29,if (show.legend) .78 else if (show.pval || !is.null(title)) .49 else .14)
  par(mai=c(bottom,min(.72,par("fin")[1]*.20),upper,.13),mgp=c(1.9,.55,0))
  plot(NA_real_,NA_real_,type="n",xlim=limits,ylim=c(0,top),axes=FALSE,
       xlab=if (show.local) "" else .angle_title(coord),ylab="Density",xaxs="i",yaxs="i",...)
  yt <- .plot_ticks(c(0,top)); abline(h=yt,col=theme$grid,lwd=.65)
  for (off in offsets) {
    at <- d$x+off
    polygon(c(at,rev(at)),c(null$CE.upper,rev(low)),col=theme$envelope,border=NA)
  }
  targets <- .density_targets(obj,coord,limits,theme)
  for (off in offsets) lines(d$x+off,d$y,col=cc,lwd=1.65)
  .plot_axis(2,yt,theme=theme)
  if (show.legend) .density_key(c("Observed",if(length(x$metadata$conf)==1L) sprintf("%g%% pointwise envelope",100*x$metadata$conf) else "Pointwise envelope"),
        c(cc,theme$envelope),c(1,0),targets,theme,fills=c(NA,theme$envelope))
  if (!is.null(title)) mtext(title,side=3,line=if(show.legend) 2.1 else .65,adj=0,cex=.93)
  if (show.pval) {
    p <- x$metadata$global.pval
    p <- if(p==0) paste("<",format(signif(1/(x$metadata$nsims+1),3))) else paste("=",format.pval(p,digits=2))
    mtext(paste("Global p",p),side=3,adj=1,line=if(show.legend) 2.0 else .55,cex=.86,font=2)
  }
  if (!show.local) .density_axis(coord,limits,theme) else {
    # Lower panel uses the same x transform and its own space below the density plot.
    old <- par(xpd=NA); on.exit(par(old),add=TRUE)
    local <- x$metadata$local.pval
    y <- -yinch(bottom*c(.27,.48))
    adjusted <- identical(x$metadata$local.inference,"maximum-cluster-v1")
    text(limits[1],-yinch(bottom*.085),if (adjusted) "Adjusted local significance" else "Local significance (unadjusted)",
         adj=c(0,.5),cex=.73,font=2,col=theme$ink)
    segments(limits[1],y,limits[2],y,col=theme$grid,lwd=.65)
    text(limits[1]-xinch(.055),y,c("Excess","Deficit"),adj=c(1,.5),cex=.73,col=theme$muted)
    if (NROW(local)) for (i in seq_len(nrow(local))) {
      excess <- local[[1]][i]=="+"
      if (is.na(excess) || !is.finite(suppressWarnings(as.numeric(local[[2]][i]))) ||
          !is.finite(suppressWarnings(as.numeric(local[[3]][i]))) || is.na(local[[4]][i])) next
      row <- if (excess) y[1] else y[2]
      bounds <- as.numeric(c(local[[2]][i],local[[3]][i]))
      if (coord=="az" && bounds[2]<bounds[1]) bounds[2] <- bounds[2]+360
      for (off in .curve_offsets(bounds,coord,limits)) {
        l <- max(min(limits),bounds[1]+off); r <- min(max(limits),bounds[2]+off)
        if (l>r) next
        segments(l,row,r,row,col=pal[1+excess],lwd=3,lend="butt")
        .local_mark(mean(c(l,r)),row+yinch(min(.088,bottom*.083)),pval2stars(local[[4]][i]),theme$ink)
      }
    }
    .density_axis(coord,limits,theme,pos=-yinch(bottom*.60))
    text(mean(limits),-yinch(bottom*(if (show.legend) .83 else .91)),.angle_title(coord),cex=.82)
    if (show.legend) .local_key(limits,-yinch(bottom*.936),theme)

  }
  invisible(NULL)
}

.local_mark <- function(x,y,symbol,col,scale=1) {
  # Physical units keep every glyph centred at precisely the same height.
  if (symbol=="\u25cb") {
    theta <- seq(0,2*pi,length.out=41)
    lines(x+xinch(scale*1.45/72)*cos(theta),y+yinch(scale*1.45/72)*sin(theta),col=col,lwd=.65)
    return(invisible(NULL))
  }
  plus <- symbol=="+"
  n <- if (plus) 1L else nchar(symbol)
  angles <- if (plus) c(0,90) else c(90,30,150)
  for (i in seq_len(n)) {
    cx <- x+xinch((i-(n+1)/2)*scale*5.1/72)
    dx <- xinch(scale*2/72)*cos(angles*pi/180); dy <- yinch(scale*2/72)*sin(angles*pi/180)
    segments(cx-dx,y-dy,cx+dx,y+dy,col=col,lwd=.85,lend="round")
  }
}

.local_key <- function(limits,y,theme) {
  sym <- c("***","**","*","+","\u25cb")
  lab <- expression(phantom(p) <= .001, phantom(p) <= .01, phantom(p) <= .05,
                    phantom(p) <= .10, phantom(p) > .10)
  size <- min(.72,par("pin")[1]/5.45)
  icon <- xinch(.17)*size
  w <- strwidth(lab,cex=size)+icon+xinch(.10)*size
  start <- mean(limits)-(sum(w)+strwidth("Local p   ",cex=size))/2
  text(start,y,"Local p",adj=c(0,.5),cex=size,col=theme$muted)
  left <- start+strwidth("Local p   ",cex=size)
  for (i in seq_along(sym)) {
    .local_mark(left+icon/2,y,sym[i],theme$muted,size)
    text(left+icon,y,lab[i],adj=c(0,.5),cex=size,col=theme$muted)
    left <- left+w[i]
  }
}

#' Plot seasonal stellar visibility
#' @description Compare the seasonal visibility of one or more stars at the modelled site and epoch. Each result shows the stellar name, designation and epoch. To compare one star across epochs, compute a phase object for each year and pass the second object to plot as its second argument.
#' @details Read each colour as a daily model classification: visible rising and setting, rising
#' only, setting only, visible without a qualifying near-horizon passage, or invisible. The
#' function plots the stored phase sequence on a 365-day scale and can label seasons or
#' display a legend. A single result defaults to a circular view, advancing clockwise
#' from December solstice. Further results use concentric rings with one leader per displayed layer;
#' optional sorting uses the stored apparent magnitude or seasonal onset. Changing the
#' display does not recalculate the visibility thresholds. Use \code{plot(first, list(second, third))} to compare results.
#' @return Draws on the active graphics device and returns NULL invisibly.
#' @param x A star phase result.
#' @param y A second object of \emph{skyscapeR.starphases} format or
#' a list of \emph{skyscapeR.starphases} objects for plotting (see Examples).
#' @param show.grid Draw monthly dividers through each ring or strip. Defaults to TRUE.
#' Set FALSE to omit the monthly divisions and retain the seasonal axis.
#' @param show.legend Show a legend for seasons present among the plotted
#' stars (default TRUE). For circular plots FALSE also hides titles, site/epoch text and all layer or event callouts. The seasonal scale stays visible. Linear plots keep their data and event labels.
#' @param style Analytical display style, either \code{"publication-colour"} (default) or
#' \code{"publication-mono"}. Monochrome uses grayscale and preserves line styles.
#' @param pal Five colours in the order RS (rise and set), R (rise only), S (set only),
#' V (curtailed passage) and I (arising and lying hidden). The default uses five colours from the plot style. Monochrome adds textures to distinguish the classes. Supply all five even when only some
#' seasons occur.
#' @param ... Additional arguments to be passed to \emph{plot}.
#' @param show.labels Logical. Mark annual events for a single result in either layout. Defaults to TRUE.
#' @param layout \code{"auto"} and \code{"circular"} draw a single ring or concentric rings
#' in input (or requested sort) order from the outside inward. Leaders connect the layer identifiers
#' to their rings. Select \code{"linear"} for strips. Supply this argument by name.
#' @param sort.by Optional ordering by catalogue magnitude (\code{"mag"}), first R-only
#' sample (\code{"rise"}) or first S-only sample (\code{"set"}) on the stored
#' cycle. This is not sorting by the dated heliacal event. Missing classes are
#' placed last; omission retains input order.
#' @rdname plot.skyscapeR.starphases
#' @section Validation:
#' W01 checked season boundaries, event positions, quarterly landmarks and star, site and epoch labels in linear and circular layouts. See the \href{https://stonexsky.github.io/skyscapeR/manual/#w01-plot-geometry}{W01 check in the technical manual}.
#' @export
#' @seealso \code{\link{star.phases}}
#' @examples
#' \dontrun{
#' if (requireNamespace("swephRdata", quietly = TRUE)) {
#'   # Plot the seasonality of Aldebaran for 4001 BCE:
#'   s1 <- star.phases('Aldebaran',-4000, c(35,-8, 200))
#'   plot(s1) # circular by default
#'   plot(s1, layout="linear")
#'   # Plot seasonality of several stars
#'   s2 <- star.phases('Alnilam',-4000, c(35,-8, 200))
#'   s3 <- star.phases('Elnath',-4000, c(35,-8, 200))
#'   plot(s1, s2)
#'   plot(s1, list(s2,s3))
#'   plot(s1, show.legend=FALSE, show.grid=FALSE)
#' }
#' }
#' @section Seasonal axis and comparison:
#' The axis is a 365-sample synthetic solar-longitude cycle anchored at December
#' solstice. Equinoxes and solstices mark quarter-cycle positions. The monthly
#' dividers mark twelve equal intervals of 30 degrees in solar longitude,
#' beginning at December solstice. For calendar dates use the underlying result's
#' \code{metadata$events$date.julian}, as described in \code{\link{star.phases}}.
#'
#' Objects may represent different stars, epochs, sites or model assumptions.
#' Shared star, epoch and site appear once. Different epochs of the same star label individual rings or rows. Mixed comparisons identify both star and year. Hold the site,
#' extinction and other inputs fixed for a controlled epoch comparison, or
#' state their differences in the figure caption. Names on the supplied list can
#' override the displayed star names. The legend contains only the union of
#' season classes present in the plotted objects.
#'
#' Event annotations use transitions between adjacent samples, including the last
#' and first samples of the cycle. A season continuing through December solstice
#' has no event at that boundary. Entry into and exit from I mark acronycal
#' setting and heliacal rising; entry into and exit from V mark acronycal rising
#' and heliacal setting. Constant classifications have no transition labels.
#' A single ring labels each event by name on two lines just outside the rim.
#' Its short leader meets the ring at the event's model day. A single ring shows the star, epoch and site in its centre. Rings for one star at several epochs show the star, designation and shared site in the centre, with a year for each ring. The circular layout
#' labels the solstices and equinoxes outside its rim. Both layouts mark the event's model day,
#' the last sample before entry into I or V, or the first after exit.
#' The linear seasonal scale sits immediately below the lowest strip. Each quarter
#' label aligns with its solstice or equinox across every row. Rising-event labels
#' sit close above a single strip and setting-event labels close below it, leaving
#' room for the seasonal scale where needed. Both layouts use the same
#' sequence and colours; the circular layout makes the cycle boundary explicit.
#' Labels wrap and shrink to fit; crowded comparisons need a larger device.
plot.skyscapeR.starphases <- function(x, y, show.labels=TRUE, show.legend=TRUE, show.grid=TRUE, sort.by, pal, ..., layout = c("auto", "linear", "circular"), style="publication-colour") {
  layout <- match.arg(layout)
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
  theme <- .plot_style(style)
  if (missing(pal)) pal <- if (theme$mono) c("#575757","#ABABAB","#BDBDBD","#DADADA","#EEEEEE") else
    c("#3A7777","#77A9CF","#D99061","#B3A1C9","#E0E4E3")
  if (length(pal)!=5L) stop("Supply five season colours in pal.")
  pal <- .plot_colour(pal,theme)
  layout <- if (layout=="auto") "circular" else layout
  info <- .phase_labels(x)
  if (layout=="circular") .plot_phase_circle(x,info,show.labels,show.legend,show.grid,pal,theme,...)
  else .plot_phase_rows(x,info,show.labels,show.legend,show.grid,pal,theme,...)
  invisible(NULL)
}

.phase_labels <- function(objects) {
  value <- function(s) if (!length(s) || is.na(s[1])) "" else as.character(s[1])
  names0 <- vapply(objects,function(z) value(z$metadata$star$name),character(1))
  bayer <- vapply(objects,function(z) value(z$metadata$star$Bayer),character(1))
  names0[!nzchar(names0)] <- bayer[!nzchar(names0)]
  ids <- ifelse(nzchar(bayer),bayer,names0)
  original <- ifelse(nzchar(bayer) & bayer!=names0,paste(names0,bayer,sep="  "),names0)
  use <- rep(FALSE,length(objects))
  if (!is.null(names(objects))) {
    use <- !is.na(names(objects)) & nzchar(names(objects))
    names0[use] <- names(objects)[use]
  }
  star <- ifelse(nzchar(bayer) & bayer!=names0,paste(names0,bayer,sep="  "),names0)
  years <- vapply(objects,function(z) if (length(z$metadata$year)==1L)
      value(BC.AD(z$metadata$year)) else "",character(1))
  sites <- vapply(objects,function(z) .site_caption(z$metadata$loc,z$metadata$site.name),character(1))
  same.star <- length(unique(ids))==1L
  same.year <- length(unique(years))==1L
  same.site <- length(unique(sites))==1L
  row <- if (same.star && !same.year) years else if (same.year) star else paste(star,years,sep="  ")
  if (same.star && !same.year && any(use)) row[use] <- paste(names0[use],years[use],sep="  ")
  if (!same.site) row <- paste(row,sites,sep="  ")
  heading <- c(if (same.star) original[1], if(same.site) sites[1], if(same.year) years[1])
  list(rows=unname(row),epoch.comparison=same.star && !same.year,
       centre=c(original[1],if (same.year) years[1],if (same.site) sites[1]),heading=paste(heading[nzchar(heading)],collapse="  |  "))
}

.phase_runs <- function(z) {
  classes <- match(as.character(z$data$phase),c("RS","R","S","V","I"))
  if (!length(classes) || anyNA(classes)) stop("Unrecognised stellar phase classification.")
  run <- rle(classes)
  ends <- cumsum(run$lengths)
  data.frame(start=c(0,head(ends,-1)),end=ends,code=run$values)
}

.phase_month_edges <- (0:12) * 365 / 12

.phase_plot_events <- function(z) {
  events <- .phase_transitions(z$data$phase)
  events[order(events$day), , drop=FALSE]
}

.phase_fill <- function(x,y,k,pal,theme) {
  polygon(x,y,col=pal[k],border=NA)
  if (theme$mono && k %in% 2:4) {
    polygon(x,y,density=c(0,14,14,9,0)[k],angle=c(0,48,-48,0,0)[k],col="#555555",border=NA)
  }
}

.phase_key <- function(runs,pal,theme,offset=.12) {
  k <- sort(unique(unlist(lapply(runs,`[[`,"code"))))
  labels <- c("Rise and set","Rise only","Set only","Curtailed passage","Hidden")[k]
  sz <- min(.77,par("fin")[1]/max(1.1,sum(strwidth(labels,units="inches"))+length(k)*.35))
  bounds <- par("usr")
  legend(mean(bounds[1:2]),bounds[3]-yinch(offset),xjust=.5,yjust=1,xpd=NA,bty="n",
         legend=labels,fill=pal[k],horiz=TRUE,cex=sz,text.col=theme$ink,
         density=if(theme$mono) c(NA,14,14,9,NA)[k] else NULL,
         angle=c(0,48,-48,0,0)[k],border=theme$muted,x.intersp=.6)
}

.plot_phase_circle <- function(objects,info,show.labels,show.legend,show.grid,pal,theme,...) {
  n <- length(objects)
  old <- par(c("mai","pty","xpd")); on.exit(par(old))
  par(mai=if(show.legend) c(.51,.08,.38,.08) else rep(.08,4),
      pty="m",xpd=FALSE)
  bare <- !show.legend
  multi <- n>1L
  events <- if (!multi && show.labels) .phase_plot_events(objects[[1]])
  direct <- !bare && NROW(events)>0L
  bounds <- if (multi && !bare) c(-1.26,2.40) else c(-1.26,1.26)
  ybounds <- c(-1.26,1.26)
  if (direct) {
    names2 <- sub(" ","\n",events$event,fixed=TRUE)
    # Reserve space on either side in the actual device's font and dimensions.
    size <- min(.78, .23*par("pin")[1]/max(strwidth(names2,units="inches")))
    margin <- c(max(strwidth(names2,cex=size,units="inches")),
                max(strheight(names2,cex=size,units="inches"))) + .11
    limits <- pmax(1.24,1/pmax(.25,1-2*margin/par("pin")))
    bounds <- c(-1,1)*limits[1]
    ybounds <- c(-1,1)*limits[2]
  }
  args <- modifyList(list(x=NA_real_,y=NA_real_,xlim=bounds,ylim=ybounds,
                     asp=1,axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i"),list(...))
  if (bare) args[c("main","sub","xlab","ylab")] <- list("","","","")
  do.call(plot.default,args)
  if (multi && !info$epoch.comparison && show.legend && is.null(args$main) && nzchar(info$heading))
    mtext(info$heading,side=3,line=.32,cex=min(.9,par("pin")[1]/max(.01,strwidth(info$heading,units="inches",font=2))),font=2)
  runs <- lapply(objects,.phase_runs)
  edges <- seq(1,.49,length.out=n+1L); pad <- min(.006,.12*abs(diff(edges)[1]))
  for (j in seq_along(runs)) {
    r <- runs[[j]]
    for (k in seq_len(nrow(r))) {
      day <- seq(r$start[k],r$end[k],length.out=max(3L,r$end[k]-r$start[k]+1L))
      p <- rbind(.polar_xy(day*360/365,edges[j]-pad),.polar_xy(rev(day)*360/365,edges[j+1]+pad))
      .phase_fill(p[,1],p[,2],r$code[k],pal,theme)
    }
  }
  months <- if (show.grid) 0:11 else c(0,3,6,9)
  angles <- .phase_month_edges[months+1L] * 360/365
  a <- .polar_xy(angles,tail(edges,1)); b <- .polar_xy(angles,edges[1])
  segments(a[,1],a[,2],b[,1],b[,2],col="white",lwd=ifelse(months%%3==0,1.4,.9))
  if ((!multi || info$epoch.comparison) && show.legend) {
    caption <- info$centre[nzchar(info$centre)]
    if (length(caption)) {
      fonts <- 1L+(caption==info$centre[1])
      widths <- mapply(function(z,f) strwidth(z,font=f),caption,fonts)
      size.centre <- min(.91,.80/max(widths),.49/(length(caption)*strheight("M")*1.5))
      step <- 1.6*strheight("M",cex=size.centre)
      yy <- (mean(seq_along(caption))-seq_along(caption))*step
      text(0,yy,labels=caption,cex=size.centre,font=fonts,col=theme$ink)
    }
  }
  p <- .polar_xy(0:3*90,1.105)
  if (multi && show.legend) p[2,2] <- p[2,2]-.18
  for (i in 1:4) text(p[i,1],p[i,2],c("Dec. sol.","Mar. eq.","Jun. sol.","Sep. eq.")[i],
                       srt=if(multi && show.legend) 0 else c(0,270,0,90)[i],cex=.70,col=theme$muted)
  if (show.legend && multi) {
    angle <- if (n==2) c(48,120) else seq(30,140,length.out=n)
    xy <- .polar_xy(angle,(head(edges,-1)+tail(edges,-1))/2)
    y <- seq(.92,-.83,length.out=n)
    size <- min(.82,.86*abs(diff(y)[1])/max(.01,strheight("M")),
                (par("usr")[2]-1.47)/max(.01,strwidth(info$rows)))
    for (j in seq_len(n)) {
      lines(c(xy[j,1],1.21,1.40),c(xy[j,2],y[j],y[j]),lwd=.7,col=theme$muted)
      points(xy[j,,drop=FALSE],pch=16,cex=.43,col=theme$ink)
      text(1.45,y[j],info$rows[j],adj=c(0,.5),cex=size,col=theme$ink)
    }
  }
  if (NROW(events)) {
    rim <- .polar_xy(events$day*360/365)
    if (direct) .phase_ring_events(rim,names2,size,theme)
    else segments(rim[,1],rim[,2],1.035*rim[,1],1.035*rim[,2],col=theme$ink,lwd=.95)
  }
  if (show.legend) .phase_key(runs,pal,theme)
  invisible(NULL)
}

.phase_ring_events <- function(rim,labels,size,theme) {
  # Separate neighbours along each side without moving their anchors on the ring.
  gap <- xinch(.06)
  height <- max(strheight(labels,cex=size))
  for (side in c(-1,1)) {
    i <- which(if (side==1) rim[,1]>=0 else rim[,1]<0)
    if (!length(i)) next
    target <- rim[i,2]*(1+gap+height/2)
    bound <- 1+gap+height/2
    # Leave space for the rotated equinox label on each side of the ring.
    equinox <- if (side==1) "Mar. eq." else "Sep. eq."
    clearance <- strwidth(equinox,cex=.70)/2 + height/2 + gap/3
    y <- target
    for (sign in c(-1,1)) {
      group <- if (sign==1) target>=0 else target<0
      y[group] <- sign*.label_rows(sign*target[group],clearance,bound,height+yinch(.04))
    }
    # The entire two-line label should clear the curved edge, even near a solstice.
    near <- pmax(0,abs(y)-height/2)
    edge <- sqrt(pmax(0,1-near^2))
    # Clear the horizontal solstice names at the top and bottom.
    solstice <- ifelse(y>=0,"Dec. sol.","Jun. sol.")
    near.solstice <- abs(abs(y)-1.105)<(height+strheight("M",cex=.70))/2+gap/3
    reserved <- ifelse(near.solstice,strwidth(solstice,cex=.70)/2,0)
    x <- side*(pmax(edge,abs(rim[i,1]),reserved)+gap)
    segments(rim[i,1],rim[i,2],x-side*gap/3,y,col=theme$muted,lwd=.7)
    text(x,y,labels=labels[i],adj=c(if(side==1) 0 else 1,.5),cex=size,col=theme$ink)
  }
}

.plot_phase_rows <- function(objects,info,show.labels,show.legend,show.grid,pal,theme,...) {
  n <- length(objects)
  single <- n==1L
  h <- if (single) .21 else .32
  bottom.edge <- if (single) -h else 1-h
  runs <- lapply(objects,.phase_runs)
  old <- par(c("mai","mgp","xpd","pty")); on.exit(par(old))
  row.labels <- if(n==1) "" else info$rows
  w <- max(.01,strwidth(row.labels,units="inches"))
  size <- min(.82,.29*par("fin")[1]/w,par("fin")[2]/max(1,n)*2.9)
  bottom.margin <- if (single) .12 else .40
  par(mai=c(bottom.margin+if(show.legend) .34 else 0,.12+w*size,.43,.2),
      mgp=c(1.85,.5,0),xpd=FALSE,pty="m")
  ylim <- if (single) c(-.96,.96) else c(bottom.edge,n+.55)
  args <- modifyList(list(x=NA_real_,y=NA_real_,xlim=c(0,365),ylim=ylim,type="n",
                 xaxs="i",yaxs="i",axes=FALSE,xlab="",ylab=""),list(...))
  do.call(plot.default,args)
  if (is.null(args$main) && nzchar(info$heading))
    mtext(info$heading,side=3,adj=0,line=.25,cex=min(.88,par("pin")[1]/max(.01,strwidth(info$heading,units="inches"))))
  for (j in seq_along(runs)) {
    y <- if (single) 0 else n-j+1
    r <- runs[[j]]
    for (k in seq_len(nrow(r))) .phase_fill(c(r$start[k],r$end[k],r$end[k],r$start[k]),
                                      y+c(-h,-h,h,h),r$code[k],pal,theme)
    if (show.grid) segments(.phase_month_edges,y-h,.phase_month_edges,y+h,
                            col="white",lwd=ifelse((0:12)%%3==0,1.4,.9))
    if (nzchar(row.labels[j])) text(-xinch(.055),y,row.labels[j],adj=c(1,.5),cex=size,xpd=NA)
  }
  # Keep the shared scale a fixed physical distance from the coloured strip.
  at <- .phase_month_edges[c(1,4,7,10,13)]
  quarter <- c("Dec. solstice","Mar. equinox","June solstice","Sept. equinox","Dec. solstice")
  if (max(strwidth(quarter,cex=.72))*1.6 > 365/4)
    quarter <- c("Dec.\nsol.","Mar.\neq.","Jun.\nsol.","Sep.\neq.","Dec.\nsol.")
  quarter.width <- strwidth(quarter,cex=.72)
  left <- pmax(0,pmin(365-quarter.width,at-quarter.width/2))
  right <- left+quarter.width
  top <- bottom.edge-yinch(.055)
  bottom <- top-max(strheight(quarter,cex=.72))-yinch(.025)
  visible <- at>=min(par("usr")[1:2]) & at<=max(par("usr")[1:2])
  segments(at[visible],bottom.edge,at[visible],bottom.edge-yinch(.032),col=theme$muted,lwd=.65,xpd=NA)
  text(left[visible],top,labels=quarter[visible],adj=c(0,1),cex=.72,col=theme$muted,xpd=NA)
  if (single && show.labels) {
    ev <- .phase_plot_events(objects[[1]])
    for (side in c(-1,1)) {
      rows <- which(if (side==1) endsWith(ev$event,"rising") else endsWith(ev$event,"setting"))
      if (!length(rows)) next
      pos <- ev$day[rows]; labs <- ev$event[rows]
      size <- min(.78, .92*365/max(1,sum(strwidth(labs))+.03*365*length(rows)))
      width <- max(strwidth(labs,cex=size))
      xx <- .label_rows(pos,width/2,365-width/2,width+2)
      start <- side*h
      end <- rep(start+side*yinch(.055),length(pos))
      if (side==-1) {
        # A setting name needs a second line of space only under a quarter label.
        half <- strwidth(labs,cex=size)/2+xinch(.02)
        for (i in seq_along(pos)) {
          overlaps <- visible & left<=max(pos[i],xx[i]+half[i]) &
                               right>=min(pos[i],xx[i]-half[i])
          if (any(overlaps)) end[i] <- min(end[i],bottom-yinch(.025))
        }
      }
      # Interrupt a setting leader only where it passes through a quarter label.
      from <- rep(start,length(pos))
      from.x <- pos
      if (side==-1) {
        x.top <- pos+(xx-pos)*(top-start)/(end-start)
        x.bottom <- pos+(xx-pos)*(bottom-start)/(end-start)
        for (i in seq_along(pos)) {
          if (end[i]<bottom && any(visible & left<=max(x.top[i],x.bottom[i]) & right>=min(x.top[i],x.bottom[i]))) {
            segments(pos[i],start,x.top[i],top,col=theme$muted,lwd=.7)
            from.x[i] <- x.bottom[i]
            from[i] <- bottom
          }
        }
      }
      segments(from.x,from,xx,end,col=theme$muted,lwd=.7)
      text(xx,end+side*yinch(.03),labs,adj=c(.5,if (side==1) 0 else 1),cex=size)
    }
  }
  if (show.legend) .phase_key(runs,pal,theme,offset=if(single) .07 else .35)
  invisible(NULL)
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
#' @section Validation:
#' The tested public fits printed the stored global probability and local interval endpoints. A wrapped interval has start greater than end. A plot spanning the full turn displays its two visible segments with the same adjusted probability. See \href{https://stonexsky.github.io/skyscapeR/manual/#s03-scoring-checks}{S03 validation in the technical manual}.
#'
#' W01 checked the reported wrapped interval endpoints and adjusted probabilities against a known fit. See the \href{https://stonexsky.github.io/skyscapeR/manual/#w01-plot-geometry}{W01 check in the technical manual}.
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

.plot_target_bands <- function(groups,unit,show.obj.label,theme=.plot_style("publication-colour")) {
  for (g in groups) {
    if (unit=="Azimuth") {
      warning("Celestial overlays require declination coordinates in this plot.",call.=FALSE)
      next
    }
    for (i in seq_len(g$n)) {
      dec <- g$decs[,i]; cc <- .plot_colour(g$col[i],theme)
      if (length(g$epoch)==1) abline(v=dec,col=cc,lwd=g$lwd[i],lty=g$lty[i])
      else rect(min(dec),par("usr")[3],max(dec),par("usr")[4],
                col=.with_alpha(cc,.17),border=cc,lty=g$lty[i])
      if (show.obj.label) text(mean(dec),par("usr")[4]-.16,.target_name(.target_id(g,i)),
                               pos=4,offset=.12,cex=.65,col=cc)
    }
  }
}
