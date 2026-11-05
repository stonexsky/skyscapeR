#' Export a horizon as an equirectangular PNG
#' @description Write a full-sphere 2:1 RGBA texture for planetaria or image editors.
#' @param hor A skyscapeR.horizon object.
#' @param file Destination PNG filename. Choose a new filename for each export.
#' @param width Even output width in pixels. Height is half the width. Default 4096.
#' @param use.image Use the attached terrain or photographic texture when available.
#' Set FALSE to export the numerical skyline as a silhouette.
#' @param col.ground Fill colour for a numerical silhouette.
#' @details The left edge is true north. Azimuth increases clockwise through east,
#' south and west. The top is the zenith, the middle is zero altitude and the
#' bottom is the nadir. Image exports preserve the calibrated angular bounds and
#' sample source pixels by nearest neighbour. Sky and areas outside photographic
#' coverage are transparent. For a numerical skyline, the ground extends to the
#' nadir over the measured azimuth sector. A missing or modified image triggers
#' a warning and the numerical silhouette is exported instead.
#' @return The absolute PNG path, invisibly.
#' @examples
#' \dontrun{
#' exportPanorama(hor, "site.png")
#' exportPanorama(hor, "site-silhouette.png", width = 2048, use.image = FALSE)
#' }
#' @seealso \code{\link{exportStellarium}}, \code{\link{importPanorama}}
#' @export
exportPanorama <- function(hor, file, width=4096, use.image=TRUE, col.ground="darkgreen") {
  if (!inherits(hor,"skyscapeR.horizon")) stop("Supply a skyscapeR.horizon object.")
  if (!is.character(file) || length(file)!=1L || is.na(file) || !grepl("\\.png$",file,ignore.case=TRUE))
    stop("Supply a destination filename ending in .png.")
  if (file.exists(file)) stop("Destination already exists: ",file)
  if (!is.numeric(width) || length(width)!=1L || !is.finite(width) || width<2 ||
      width%%2!=0) stop("Width must be a positive even number of pixels.")
  if (!requireNamespace("png",quietly=TRUE)) stop("Install png to export a panorama.")
  height <- width/2
  az <- (seq_len(width)-.5)*360/width
  alt <- 90-(seq_len(height)-.5)*180/height
  im <- if (use.image) .horizon_image(hor)
  out <- integer(as.double(width)*height)
  if (!is.null(im)) {
    src <- .read_panorama_pixels(im)
    bearing <- (az-im$azimuth[1])%%360
    cols <- which(bearing<diff(im$azimuth))
    ix <- 1+floor(bearing[cols]*ncol(src)/diff(im$azimuth))
    iy <- 1+floor((im$altitude[2]-alt)*nrow(src)/diff(im$altitude))
    for (i in which(iy>=1 & iy<=nrow(src)))
      out[(i-1)*width+cols] <- src[(iy[i]-1)*ncol(src)+ix]
  } else {
    if (identical(hor$metadata$source$provider,"Panorama")) {
      skyline <- hor2alt(hor,az)
    } else {
      data <- hor$data
      order.central <- order(!(data$az>=0 & data$az<360))
      x <- (data$az%%360)[order.central]
      keep <- !duplicated(x)
      x <- x[keep]; y <- data$alt[order.central][keep]
      skyline <- stats::approx(c(x-360,x,x+360),rep(y,3),xout=az)$y
    }
    rgba <- grDevices::col2rgb(col.ground,alpha=TRUE)
    if (ncol(rgba)!=1L) stop("Supply one ground colour.")
    pixel <- readBin(as.raw(rgba),integer(),n=1L,size=4,endian="little")
    for (i in seq_len(height)) out[(i-1)*width+which(skyline>=alt[i])] <- pixel
  }
  dim(out) <- c(height,width)
  class(out) <- "nativeRaster"
  tmp <- tempfile(tmpdir=dirname(file),fileext=".png")
  on.exit(unlink(tmp))
  png::writePNG(out,tmp)
  if (!file.copy(tmp,file,overwrite=FALSE)) stop("Could not save panorama: ",file)
  invisible(normalizePath(file))
}
