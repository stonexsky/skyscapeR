#' Import a sky-masked photographic panorama
#' @description Extract a numerical skyline from a transparent-sky panorama and attach the calibrated photograph.
#' @details Choose controls according to what you know:
#' \itemize{
#' \item If the image limits are already calibrated, leave \code{calibrate = FALSE} and omit \code{points}.
#' \item To rotate an already level panorama to true north, click one landmark, enter its true azimuth, leave altitude blank and choose Finish. This changes only the azimuth origin.
#' \item One altitude alone shifts the altitude origin. One landmark with both coordinates shifts both origins. Neither case estimates tilt.
#' \item Two altitude controls at different bearings can raise one side and lower the other. For tilt around a full turn, use at least three well-spaced bearings. These fits retain the supplied vertical scale. Two altitude controls at the same bearing but different pixel heights can instead calibrate vertical scale.
#' \item Multiple known azimuths at distinct horizontal positions can calibrate the horizontal scale of a partial panorama. A full turn keeps its 360-degree span.
#' }
#' The picker samples the overview. Choose Zoom and click two corners to display the selected region at native pixel detail, then choose Add control point. Full-turn plots cut the calibrated view at north. In scripted controls use \code{NA_real_} for each unknown angle. Levelling adjusts the near-horizon band. If an altitude shift takes image rows beyond the poles, only whole rows within the physical altitude range are kept. A one-point adjustment records the calibrated bounds and uses the original PNG. Keep it with the saved object. A tilt adjustment writes a levelled PNG in the working directory. Use a panorama with nominally linear angular axes and transparent sky above terrain in every column. A partial panorama supports altitude lookup within its photographed sector. Stellarium export requires a full turn.
#' @param file PNG with sky removed to transparency. Keep this file alongside saved horizon objects.
#' @param loc WGS84 latitude, longitude and optional elevation in metres.
#' @param az.lim True bearings of the left and right edges, increasing clockwise. A sector may cross north, e.g. c(330,30).
#' @param alt.lim Altitudes of the bottom and top edges in degrees. Set these explicitly for a cropped panorama.
#' @param name Site name.
#' @param calibrate Open the interactive graphics picker. Choose Zoom and click two corners to enlarge a region. Add a control by clicking a landmark and entering its true azimuth, altitude or both. Finish to apply the fit.
#' @param points Optional data frame with numeric x, y, az and alt. Pixel coordinates start at the upper-left image edge. Use NA for an unknown angle.
#' @param alpha.threshold Minimum opacity counted as terrain, from above zero to one.
#' @param alt.unc Assumed skyline altitude uncertainty in degrees.
#' @return A skyscapeR.horizon object with a numerical profile, external PNG reference, angular bounds and control points. Use saveRDS() to save the object. Copy the displayed panorama PNG with the RDS and use \code{\link{relinkHor}} to reconnect a moved horizon or list.
#' @examples
#' \dontrun{
#' # Known bearing at pixel (1240, 380) in a level, cropped panorama
#' anchor <- data.frame(x = 1240, y = 380, az = 96.4, alt = NA_real_)
#' h <- importPanorama("ridge.png", c(51.75, -9.37, 150),
#'                     az.lim = c(70, 120), alt.lim = c(-5, 12), points = anchor)
#' plot(h, show.skyline = TRUE)
#' }
#' @seealso \code{\link{hor2alt}}, \code{\link{exportPanorama}}
#' @export
importPanorama <- function(file, loc, az.lim=c(0,360), alt.lim=c(-90,90),
                           name=tools::file_path_sans_ext(basename(file)), calibrate=FALSE,
                           points=NULL, alpha.threshold=0.5, alt.unc=0.5) {
  if (!requireNamespace("png", quietly=TRUE)) stop("Install the png package first.")
  pair <- function(x) is.numeric(x) && length(x)==2L && all(is.finite(x))
  if (!pair(az.lim) || !pair(alt.lim) || diff(alt.lim)<=0 ||
      any(abs(alt.lim)>90)) stop("Supply finite angular bounds, with altitude increasing.")
  if (az.lim[2]<az.lim[1]) az.lim[2] <- az.lim[2]+360
  if (diff(az.lim)<=0 || diff(az.lim)>360) stop("Azimuth must span up to one turn clockwise.")
  if (!is.numeric(loc) || !length(loc) %in% 2:3 || any(!is.finite(loc)) ||
      abs(loc[1])>90 || abs(loc[2])>180) stop("Supply latitude, longitude and optional elevation.")
  if (length(loc)==2L) loc <- c(loc,NA_real_)
  if (length(alpha.threshold)!=1L || !is.finite(alpha.threshold) ||
      alpha.threshold<=0 || alpha.threshold>1 || length(alt.unc)!=1L ||
      !is.finite(alt.unc) || alt.unc<0) stop("Invalid alpha threshold or altitude uncertainty.")
  file <- normalizePath(file,mustWork=TRUE)
  img <- png::readPNG(file,native=TRUE)
  w <- ncol(img); h <- nrow(img)
  if (w<2 || h<2) stop("Panorama must have at least two rows and columns.")
  if (calibrate) points <- .panorama_picker(img,points)
  fit <- function(limits, pixel, value, size, bearing=FALSE) {
    ok <- is.finite(value)
    if (!any(ok)) return(limits)
    t <- pixel[ok]/size
    pred <- limits[1]+t*diff(limits)
    error <- value[ok]-pred
    if (bearing) error <- (error+180)%%360-180
    if (length(unique(t))==1L || (bearing && diff(limits)==360))
      return(limits+mean(error))
    coeff <- stats::lm.fit(cbind(1,t),pred+error)$coefficients
    unname(coeff[1]+c(0,1)*coeff[2])
  }
  source.rows <- NULL
  original <- file
  tilt <- NULL
  if (!is.null(points)) {
    if (!is.data.frame(points) || !all(c("x","y","az","alt") %in% names(points)) ||
        !all(vapply(points[c("x","y","az","alt")],is.numeric,logical(1))) ||
        any(!is.finite(points$x) | points$x<0 | points$x>w |
            !is.finite(points$y) | points$y<0 | points$y>h) ||
        any(is.na(points$az)&is.na(points$alt)) ||
        any(abs(points$alt)>90,na.rm=TRUE) || any(is.infinite(points$az)))
      stop("Invalid panorama control points.")
    az.lim <- fit(az.lim,points$x,points$az,w,TRUE)
    ok <- is.finite(points$alt)
    # Near the horizon a tilted vertical axis has a sinusoidal altitude error.
    # Centre the angular terms so a common altitude offset remains independent.
    if (length(unique(points$x[ok])) >= 2L) {
      phase <- function(x) (x/w-.5)*diff(az.lim)*pi/180
      basis <- function(x) cbind(sin(phase(x)),cos(phase(x)))
      B <- basis(points$x[ok])
      centre <- colMeans(B)
      A <- sweep(B,2,centre)
      dec <- svd(A)
      keep <- dec$d > max(1e-8,dec$d[1]*.01)
      raw <- alt.lim[2]-points$y[ok]/h*diff(alt.lim)
      error <- points$alt[ok]-raw
      if (any(keep)) {
        slope <- dec$v[,keep,drop=FALSE] %*%
          (crossprod(dec$u[,keep,drop=FALSE],error-mean(error))/dec$d[keep])
        shift <- as.vector(mean(error)+sweep(basis(seq_len(w)-.5),2,centre)%*%slope)
        tilt <- list(offset=mean(error),coefficients=as.vector(slope),centre=centre)
        step <- diff(alt.lim)/h
        # Preserve native row-major pixels. Pad with transparency and clip at the poles.
        top <- min(90,alt.lim[2]+ceiling(max(shift)/step)*step)
        bottom <- max(-90,alt.lim[1]+floor(min(shift)/step)*step)
        rows <- max(2L,round((top-bottom)/step))
        warped <- integer(rows*w)
        for (j in seq_len(w)) {
          target <- top-(seq_len(rows)-.5)*(top-bottom)/rows
          src <- floor((alt.lim[2]+shift[j]-target)/step)+1L
          valid <- src>=1 & src<=h
          warped[(which(valid)-1)*w+j] <- img[(src[valid]-1)*w+j]
        }
        dim(warped) <- c(rows,w)
        class(warped) <- "nativeRaster"
        img <- warped
        h <- rows
        alt.lim <- c(bottom,top)
      }
    }
    if (is.null(tilt)) alt.lim <- fit(alt.lim,h-points$y,points$alt,h)
  }
  if (!all(is.finite(az.lim)) || diff(az.lim)<=0 || diff(az.lim)>360+1e-8)
    stop("Azimuth controls must give a clockwise span of at most 360 degrees. Check the bearings and az.lim.")
  if (!all(is.finite(alt.lim)) || diff(alt.lim)<=0)
    stop("Altitude must increase from bottom to top. Check the altitude controls and alt.lim.")
  # Keep whole pixel rows on their fitted angular grid when a shift crosses a pole.
  if (any(abs(alt.lim)>90)) {
    step <- diff(alt.lim)/h
    above <- max(0,ceiling((alt.lim[2]-90)/step-1e-10))
    below <- max(0,ceiling((-90-alt.lim[1])/step-1e-10))
    if (h-above-below<2)
      stop("Fewer than two image rows fall between -90 and +90 degrees. Check alt.lim and the altitude controls.")
    if (above+below>0) {
      source.rows <- c(above+1,h-below)
      h <- h-above-below
    }
    alt.lim <- c(max(-90,alt.lim[1]+below*step),min(90,alt.lim[2]-above*step))
  }
  # Use the upper edge of the first terrain pixel, measured from the image top.
  first <- if (is.null(source.rows)) 1 else source.rows[1]
  edge <- rep(NA_integer_,w)
  pending <- seq_len(w)
  for (row in seq_len(h)) {
    pixel <- img[(first+row-2)*as.double(w)+pending]
    alpha <- bitwShiftR(pixel,24L)
    alpha[is.na(alpha)] <- 128L  # Packed 0x80000000 is R's NA integer.
    found <- alpha>=255*alpha.threshold
    if (any(found)) {
      edge[pending[found]] <- row-1L
      pending <- pending[!found]
      if (!length(pending)) break
    }
  }
  if (anyNA(edge) || any(edge==0L)) stop("Each column needs transparent sky above visible terrain. Crop or repair the mask first.")
  az.lim <- az.lim-360*floor(az.lim[1]/360)
  a <- az.lim[1]+(seq_len(w)-.5)*diff(az.lim)/w
  z <- alt.lim[2]-edge*diff(alt.lim)/h
  data <- data.frame(az=a,alt=z,alt.unc=alt.unc)
  if (abs(diff(az.lim)-360)<1e-8) data$az <- data$az%%360 else {
    data <- rbind(transform(data[1,],az=az.lim[1]),data,transform(data[w,],az=az.lim[2]))
  }
  data <- do.call(rbind,lapply(-1:1,function(i) {
    shifted <- data
    shifted$az <- data$az+360*i
    shifted
  }))
  if (!is.null(tilt)) {
    file <- tempfile("levelled-panorama-",tmpdir=getwd(),fileext=".png")
    png::writePNG(img,file)
  }
  image <- list(path=file,md5=unname(tools::md5sum(file)),azimuth=az.lim,altitude=alt.lim,
    width=w,height=h,source.rows=source.rows,projection="equirectangular",attribution="User panorama",license="User supplied")
  image$preview <- .panorama_window(img,c(1,w),c(first,first+h-1),1024)
  structure(list(metadata=list(name=name,georef=matrix(loc,1,dimnames=list(NULL,c("Lat","Lon","Elev"))),
    image=image,source=list(provider="Panorama",original=original,levelling=tilt,points=points,coverage=az.lim)),
    data=data[order(data$az),]),class="skyscapeR.horizon")
}

.panorama_window <- function(img, cols, rows, limit=Inf) {
  sample.axis <- function(bounds) {
    n <- diff(bounds)+1
    count <- min(n,limit)
    bounds[1]+floor((seq_len(count)-.5)*n/count)
  }
  x <- sample.axis(cols)
  y <- sample.axis(rows)
  pixels <- integer(as.double(length(x))*length(y))
  dest <- seq_along(x)
  for (j in seq_along(y))
    pixels[(j-1)*as.double(length(x))+dest] <- img[(y[j]-1)*as.double(ncol(img))+x]
  dim(pixels) <- c(length(y),length(x))
  class(pixels) <- "nativeRaster"
  pixels
}

.panorama_picker <- function(img,points) {
  if (!interactive()) stop("Interactive calibration needs an R graphics window. Supply points for scripted use.")
  w <- ncol(img); h <- nrow(img); x <- c(0,w); y <- c(0,h)
  old <- graphics::par(no.readonly=TRUE); on.exit(graphics::par(old))
  repeat {
    graphics::plot.new(); graphics::plot.window(x,y,asp=1,xaxs="i",yaxs="i")
    # Sample the visible window at display resolution, keeping clicks in source pixels.
    size <- graphics::par("pin")*grDevices::dev.size("px")/grDevices::dev.size("in")
    left <- floor(x[1]); right <- ceiling(x[2])
    bottom <- floor(y[1]); top <- ceiling(y[2])
    overview <- left==0 && right==w && bottom==0 && top==h
    limit <- if (overview) max(1,ceiling(max(size))) else Inf
    preview <- .panorama_window(img,c(left+1,right),c(h-top+1,h-bottom),limit)
    graphics::rasterImage(preview,left,bottom,right,top,interpolate=FALSE)
    if (NROW(points)) graphics::points(points$x,h-points$y,pch=3,col="red")
    action <- utils::menu(c("Add control point","Zoom to two corners","Reset zoom","Finish"))
    if (action %in% c(0,4)) break
    if (action==3) {x <- c(0,w); y <- c(0,h); next}
    p <- graphics::locator(if (action==2) 2 else 1)
    if (!length(p$x)) next
    if (action==2) {
      if (length(p$x)==2 && diff(range(p$x))>0 && diff(range(p$y))>0) {
        xnew <- pmax(0,pmin(w,range(p$x))); ynew <- pmax(0,pmin(h,range(p$y)))
        if (diff(xnew)>0 && diff(ynew)>0) {x <- xnew; y <- ynew}
      }
    } else {
      if (p$x<0 || p$x>w || p$y<0 || p$y>h) next
      az <- suppressWarnings(as.numeric(readline("True azimuth (blank if unknown): ")))
      alt <- suppressWarnings(as.numeric(readline("Altitude (blank if unknown): ")))
      if (is.finite(az) || is.finite(alt))
        points <- rbind(points,data.frame(x=p$x,y=h-p$y,az=az,alt=alt))
    }
  }
  points
}

.read_panorama_pixels <- function(im) {
  src <- png::readPNG(im$path,native=TRUE)
  rows <- im$source.rows
  if (is.null(rows)) return(src)
  w <- as.double(ncol(src))
  pixels <- src[seq.int((rows[1]-1)*w+1,rows[2]*w)]
  dim(pixels) <- c(diff(rows)+1,w)
  class(pixels) <- "nativeRaster"
  pixels
}

.draw_photo <- function(hor, im, col.ground, show.skyline, fastLoad=FALSE) {
  coverage <- hor$metadata$source$coverage
  bounds <- graphics::par("usr")
  pixels <- if (!is.null(im)) {
    if (fastLoad) .horizon_preview(im) else .read_panorama_pixels(im)
  }
  turns <- seq(floor((bounds[1]-coverage[2])/360),ceiling((bounds[2]-coverage[1])/360))
  for (k in turns) {
    x <- coverage+360*k
    if (x[2]<=bounds[1] || x[1]>=bounds[2]) next
    if (!is.null(pixels)) {
      # Cut each visible piece on source pixel boundaries, including at north.
      w <- ncol(pixels); h <- nrow(pixels)
      dx <- diff(x)/w; dy <- diff(im$altitude)/h
      left <- max(0,floor((bounds[1]-x[1])/dx))
      right <- min(w,ceiling((bounds[2]-x[1])/dx))
      top <- max(0,floor((im$altitude[2]-bounds[4])/dy))
      bottom <- min(h,ceiling((im$altitude[2]-bounds[3])/dy))
      if (right>left && bottom>top) {
        piece <- if (left==0 && right==w && top==0 && bottom==h) pixels else
          .panorama_window(pixels,c(left+1,right),c(top+1,bottom))
        graphics::rasterImage(piece,x[1]+left*dx,im$altitude[2]-bottom*dy,
                              x[1]+right*dx,im$altitude[2]-top*dy,interpolate=FALSE)
      }
    }
    if (is.null(pixels) || isTRUE(show.skyline)) {
      az <- seq(x[1],x[2],length.out=hor$metadata$image$width+1)
      alt <- hor2alt(hor,az)
      if (is.null(pixels)) graphics::polygon(c(az,rev(az)),c(alt,rep(-90,length(az))),col=col.ground,border=NA)
      if (is.null(pixels) && !identical(show.skyline,FALSE) || isTRUE(show.skyline))
        graphics::lines(az,alt,lwd=1.5)
    }
  }
  invisible(!is.null(pixels))
}
