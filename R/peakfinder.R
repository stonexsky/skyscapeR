#' Import a PeakFinder horizon and optional landscape drawing
#' @description Bring a terrain-derived horizon into orientation and visibility
#' calculations, optionally retaining foreground ridge lines for visual interpretation.
#' @details The numerical horizon describes angular obstruction at the observer's
#' location. A separate vector drawing supplies foreground detail but does not
#' explicitly record its vertical angular calibration. Both exports are therefore
#' needed to create a calibrated panorama PNG. Calculations always use the numerical
#' horizon, never the PNG. Terrain data do not include every tree or building.
#'
#' In PeakFinder choose Export, Horizon for Stellarium and save the ZIP. For an
#' image also choose Export, SVG for printing, Silhouettes (Vector SVG), from 0
#' to 360 degrees, preferably picture format 16:2. Use exactly the same viewpoint
#' and observer height for both. Contours (Bitmap SVG) is not supported. The PNG
#' is a white terrain silhouette with dark foreground ridge lines and transparent
#' sky, without labels or shading. Only PeakFinder's untransformed M/L vector
#' path format is supported; edited or future formats may require adaptation.
#'
#' Azimuth ruler positions calibrate horizontal coordinates. The upper envelope
#' of the vector paths is fitted to y = intercept + slope * tan(altitude) using
#' the ZIP horizon. Values at 0.01 degrees are excluded as potentially clipped.
#' Calibration needs at least 20 matches spanning 180 degrees and 0.5 degrees
#' of altitude relief. RMS and maximum angular residuals must not exceed 0.05
#' and 0.2 degrees. These are compatibility checks, not estimates of terrain
#' accuracy. For low-relief, clipped profiles the importer can instead measure the alpha edge of the eight calibrated, north-first Stellarium textures. It checks coverage, concurrence with unclipped heights, and agreement with the SVG before using the texture-derived numerical horizon. Cropped or ambiguous tiles cannot supply this calibration. Calibration failure,
#' unsupported SVGs or mismatched sites produce a warning and a horizon-only object.
#'
#' Vector vertices are transformed to linear azimuth and altitude coordinates
#' and rendered to a PNG. Its lower extent is -90 degrees (solid terrain below
#' the drawing), preventing transparent holes when exported to Stellarium. The
#' default width is 4096 pixels, about 0.088 degrees per pixel horizontally;
#' the numerical horizon retains its own sampling. No image array is stored in
#' the returned object. PNGs are read on demand as compact native rasters.
#'
#' The default image directory is the working directory. A supplied directory
#' is useful for project files. The object records an absolute PNG path and MD5
#' checksum. Saving it with saveRDS does not bundle the PNG. When moving a project,
#' copy the PNG too and reconnect it with \code{\link{relinkHor}}. Missing or
#' modified images trigger a warning and ordinary horizon plotting/export.
#' Repeating an import with the same files and rendering settings reuses the PNG filename.
#'
#' PeakFinder exports can report an elevation of 0 m plus eye height even at a
#' mountain site. A warning flags this; use elevation to supply a verified absolute
#' observer elevation. This override does not recalculate the downloaded horizon.
#' The default altitude uncertainty of 0.5 degrees is an explicit working assumption,
#' not a measured PeakFinder error estimate. Supply an appropriate value for research.
#' @param file Path to a PeakFinder Stellarium ZIP export.
#' @param svg Optional path to the matching full-circle Vector SVG export.
#' @param image.dir Directory for generated PNGs; defaults to the working directory.
#' @param width PNG width in pixels, an integer from 512 to 8192. Default 4096.
#' @param alt.unc Assumed horizon altitude uncertainty in degrees, a nonnegative scalar.
#' @param elevation Optional verified absolute observer elevation in metres.
#' @return A skyscapeR.horizon containing a periodically extended numerical horizon,
#' site and source metadata, and, on successful calibration, metadata$image with
#' the external PNG path, checksum, angular bounds and calibration diagnostics.
#' @references \href{https://www.peakfinder.com/}{PeakFinder}, exported drawings
#' attributed to PeakFinder under CC BY 4.0. Geometry conventions follow the
#' \href{https://stellarium.org/guide/}{Stellarium User Guide}.
#' @seealso \code{\link{createHor}}, \code{\link{downloadHWT}},
#' \code{\link{plot.skyscapeR.horizon}}, \code{\link{exportStellarium}}
#' @export
#' @examples
#' \dontrun{
#' hor <- importPeakFinder("landscape.zip") # warns that no PNG was created
#' hor <- importPeakFinder("landscape.zip", svg="panorama.svg",
#'                        image.dir="project-data/horizons", elevation=850)
#' plot(hor)
#' exportStellarium(hor, name="PeakFinder example")
#' }
importPeakFinder <- function(file, svg=NULL, elevation=NULL, alt.unc=0.5, width=4096, image.dir=getwd()) {
  scalar <- function(x) is.numeric(x) && length(x)==1L && is.finite(x)
  if (!scalar(alt.unc) || alt.unc<0) stop("alt.unc must be a nonnegative number.")
  if (!scalar(width) || width!=as.integer(width) || width<512 || width>8192)
    stop("width must be an integer between 512 and 8192.")
  if (!is.null(elevation) && !scalar(elevation)) stop("elevation must be finite metres.")
  if (!is.character(file) || length(file)!=1L || is.na(file) || !file.exists(file))
    stop("Provide an existing PeakFinder Stellarium ZIP file.")
  members <- utils::unzip(file,list=TRUE)
  # Never extract arbitrary archive paths or follow archive-provided filenames.
  if (any(grepl("(^/|^[A-Za-z]:|\\\\|(^|/)\\.\\.(/|$))",members$Name)))
    stop("Unsafe paths in ZIP archive.")
  ini <- members$Name[basename(members$Name)=="landscape.ini"]
  txt <- members$Name[basename(members$Name)=="horizon.txt"]
  if (length(ini)!=1L || length(txt)!=1L || dirname(ini)!=dirname(txt))
    stop("ZIP must contain one landscape.ini and matching horizon.txt.")
  if (any(members$Length[members$Name %in% c(ini,txt)]>1e7)) stop("Oversized horizon metadata.")
  read.member <- function(n) {
    con <- unz(file,n,open="rb"); on.exit(close(con))
    size <- members$Length[match(n,members$Name)]
    text <- rawToChar(readBin(con,"raw",n=size))
    Encoding(text) <- "UTF-8"
    strsplit(text,"\n",fixed=TRUE)[[1L]]
  }
  cfg <- read.member(ini); raw <- read.member(txt)
  get <- function(key, default=NULL) {
    hit <- grep(paste0("^\\s*",key,"\\s*="),cfg,value=TRUE)
    if (!length(hit)) return(default)
    if (length(hit)!=1L) stop("Ambiguous INI field: ",key)
    trimws(sub("^[^=]*=","",hit))
  }
  if (!any(grepl("peakfinder",c(cfg,raw),ignore.case=TRUE))) stop("Not a PeakFinder export.")
  if (!identical(get("polygonal_horizon_list","horizon.txt"),"horizon.txt") ||
      !identical(get("polygonal_horizon_list_mode","azDeg_altDeg"),"azDeg_altDeg") ||
      !identical(as.numeric(get("polygonal_angle_rotatez","0")),0))
    stop("Unsupported horizon filename, coordinate mode or rotation.")
  geo <- suppressWarnings(as.numeric(c(get("latitude"),get("longitude"),get("altitude"))))
  if (length(geo)!=3L || any(!is.finite(geo)) || abs(geo[1])>90 || abs(geo[2])>180)
    stop("Invalid PeakFinder location metadata.")
  original.geo <- geo
  if (!is.null(elevation)) geo[3] <- elevation else if (any(grepl(" 0m \\+",raw)))
    warning("PeakFinder reports zero ground elevation plus eye height; supply a verified elevation.",call.=FALSE)
  samples <- utils::read.table(text=paste(raw,collapse="\n"),comment.char="#")
  if (ncol(samples)!=2L || nrow(samples)<4L || !all(vapply(samples,is.numeric,logical(1))) ||
      any(!is.finite(as.matrix(samples))) || any(abs(samples[[2]])>=90)) stop("Invalid horizon samples.")
  names(samples) <- c("az","alt"); samples$az <- samples$az %% 360
  samples <- samples[order(samples$az),]
  if (anyDuplicated(samples$az)) {
    for (a in unique(samples$az[duplicated(samples$az)]))
      if (diff(range(samples$alt[samples$az==a]))>0.001) stop("Conflicting duplicate horizon bearings.")
    samples <- samples[!duplicated(samples$az),]
  }
  if (max(diff(c(samples$az,samples$az[1]+360)))>90) stop("Horizon does not cover a full circle.")
  clipped <- sum(abs(samples$alt-0.01)<1e-6)>=3L
  if (clipped && is.null(svg))
    warning("Horizon contains repeated 0.01-degree values; negative horizons may have been clipped by PeakFinder.",call.=FALSE)
  # Periodic copies preserve sparse source vertices without allocating a dense grid.
  periodic <- function(samples) do.call(rbind,lapply(c(-360,0,360),function(turn)
    data.frame(az=samples$az+turn,alt=samples$alt,alt.unc=alt.unc)))
  hor <- structure(list(metadata=list(name=get("name","PeakFinder"),
    georef=matrix(geo,nrow=1L,dimnames=list(NULL,c("Lat","Lon","Elev"))),
    source=list(provider="PeakFinder",archive=normalizePath(file),zip.md5=unname(tools::md5sum(file)),
      exported.georef=original.geo,altitude.uncertainty="user-supplied assumption",
      attribution="PeakFinder, www.peakfinder.com",license=get("license","See original export"))),
    data=periodic(samples)),class="skyscapeR.horizon")
  if (is.null(svg)) {
    warning("No vector SVG supplied: imported numerical horizon only; no panorama PNG was created.",call.=FALSE)
    return(hor)
  }
  tryCatch({
    method <- "horizon.txt"
    reference <- samples
    enough <- samples$alt[abs(samples$alt-.01)>1e-6]
    if (clipped && (length(enough)<20L || diff(range(enough))<.5)) {
      reference <- .pf_tile_profile(file,members,dirname(ini),get)
      method <- "calibrated texture alpha"
      # Independently check any unclipped elevations in the polygonal export.
      ok <- abs(samples$alt-.01)>1e-6
      if (any(ok)) {
        curve <- periodic(reference)
        delta <- stats::approx(curve$az,curve$alt,xout=samples$az[ok])$y-samples$alt[ok]
        if (max(abs(delta))>.1) stop("Texture edges disagree with the unclipped numerical horizon.")
      }
    }
    image <- .pf_image(svg,reference,original.geo,image.dir,width,hor$metadata$source$zip.md5,
                       texture=method!="horizon.txt")
    image$calibration$reference <- method
    hor$metadata$image <- image
    if (method!="horizon.txt") {
      hor$data <- periodic(reference)
      hor$metadata$source$profile <- method
      .sky_message("The exported horizon is clipped. Recovered the numerical horizon and calibrated the SVG using the Stellarium texture edges.")
    } else if (clipped) warning("Some exported horizon values may be clipped at 0.01 degrees.",call.=FALSE)
  },error=function(e) warning("Imported horizon only; no calibrated PNG: ",conditionMessage(e),
    if (clipped) " The numerical horizon may be clipped at 0.01 degrees.",call.=FALSE))
  hor
}

# Decode one tile at a time. Only accept the documented full-tile north-first layout.
.pf_tile_profile <- function(file,members,folder,get) {
  keys <- c(type="old_style",calibrated="true",nbside="8",nbsidetex="8",nb_decor_repeat="1")
  if (!all(vapply(names(keys),function(k) identical(get(k),unname(keys[k])),logical(1))) ||
      as.numeric(get("decor_angle_rotatez"))!= -90 || get("tan_mode","false")!="false")
    stop("Unsupported Stellarium texture calibration.")
  bottom <- as.numeric(get("decor_angle_shift")); span <- as.numeric(get("decor_alt_angle"))
  if (length(bottom)!=1L || length(span)!=1L || !all(is.finite(c(bottom,span))) ||
      span<=0 || bottom< -90 || bottom+span>90) stop("Invalid texture altitude limits.")
  tiles <- lapply(0:7,function(i) {
    key <- paste0("tex",i); name <- paste0("peakfinder",i,".png")
    if (!identical(get(key),name) || !identical(get(paste0("side",i)),paste0(key,":0:0:1:1")))
      stop("Unsupported tile order or texture crop.")
    member <- if (folder==".") name else paste(folder,name,sep="/")
    size <- members$Length[members$Name==member]
    if (length(size)!=1L || size<33 || size>1e8) stop("Missing or oversized texture tile.")
    con <- unz(file,member,open="rb"); on.exit(close(con))
    bytes <- readBin(con,"raw",n=size)
    dims <- readBin(bytes[17:24],integer(),2L,size=4L,endian="big")
    if (any(dims<2L | dims>4096L) || bytes[26]!=as.raw(6)) stop("Unsupported RGBA tile dimensions or format.")
    im <- png::readPNG(bytes,native=TRUE)
    cols <- unique(round(seq(1,ncol(im),length.out=min(512L,ncol(im)))))
    edge <- vapply(cols,function(j) {
      # nativeRaster stores complete scanlines consecutively.
      alpha <- bitwShiftR(im[seq.int(j,length(im),by=ncol(im))],24L)
      ground <- which(alpha>=128L)
      if (!length(ground) || ground[1]==1L || any(alpha[ground[1]:length(alpha)]<128L))
        stop("Texture edge is cropped or ambiguous.")
      ground[1]-1L
    },integer(1))
    data.frame(az=45*(i+(cols-.5)/ncol(im)),alt=bottom+span*(1-edge/nrow(im)))
  })
  do.call(rbind,tiles)
}

.pf_image <- function(svg,samples,geo,image.dir,width,zip.md5,texture=FALSE) {
  if (!requireNamespace("xml2",quietly=TRUE)) stop("Install xml2 to import vector SVG drawings.")
  if (!is.character(svg) || length(svg)!=1L || is.na(svg) || !file.exists(svg) ||
      file.info(svg)$size>2e7) stop("Provide a vector SVG file smaller than 20 MB.")
  doc <- xml2::read_xml(svg,options="NONET")
  if (length(xml2::xml_find_all(doc,"//*[local-name()='image']")))
    stop("Use Silhouettes (Vector SVG), not Contours (Bitmap SVG).")
  desc <- xml2::xml_text(xml2::xml_find_first(doc,"//*[local-name()='desc']"))
  coord <- regmatches(desc,regexec("([0-9.]+)[^0-9]+([NS]),? *([0-9.]+)[^0-9]+([EW])",desc))[[1L]]
  if (length(coord)!=5L) stop("SVG lacks readable viewpoint coordinates.")
  pos <- as.numeric(coord[c(2,4)]) * ifelse(coord[c(3,5)] %in% c("S","W"),-1,1)
  if (any(abs(pos-geo[1:2])>0.00002)) stop("SVG and ZIP viewpoints differ.")
  elev <- regmatches(desc,regexec("([-+0-9.]+)m *\\+([-+0-9.]+)m",desc))[[1L]]
  if (length(elev)!=3L || abs(sum(as.numeric(elev[2:3]))-geo[3])>0.1)
    stop("SVG and ZIP observer elevations differ or cannot be verified.")
  ruler <- xml2::xml_find_all(doc,"//*[@id='ruler']/*[local-name()='text']")
  bearing <- suppressWarnings(as.numeric(gsub("[^0-9.-]","",xml2::xml_text(ruler))))
  xpos <- suppressWarnings(as.numeric(xml2::xml_attr(ruler,"x")))
  if (length(bearing)<3L || any(!is.finite(c(bearing,xpos))) ||
      min(bearing)!=0 || abs(max(bearing)+stats::median(diff(bearing))-360)>1e-6) stop("SVG must have a full 0-to-360-degree ruler.")
  xf <- stats::lm(xpos ~ bearing); xb <- stats::coef(xf)
  box <- as.numeric(strsplit(xml2::xml_attr(doc,"viewBox")," +")[[1L]])
  if (length(box)!=4L || any(!is.finite(box)) || abs(box[3]-2*xb[1]-360*xb[2])>2)
    stop("SVG drawing width does not describe a full circle.")
  if (xb[2]<=0 || max(abs(stats::residuals(xf)))>1) stop("Unsupported SVG azimuth ruler.")
  nodes <- xml2::xml_find_all(doc,"//*[starts-with(@id,'section')]/*[local-name()='path']")
  if (!length(nodes)) stop("No PeakFinder vector ridge paths found.")
  if (length(xml2::xml_find_all(nodes,"ancestor-or-self::*[@transform]")))
    stop("Transformed SVG paths are not supported.")
  paths <- lapply(xml2::xml_attr(nodes,"d"),function(d) {
    if (grepl("[A-DF-KN-Za-df-z]",d) || !grepl("^M",d)) stop("Unsupported SVG path commands.")
    nums <- regmatches(d,gregexpr("[-+]?(?:[0-9]*\\.)?[0-9]+(?:[eE][-+]?[0-9]+)?",d,perl=TRUE))[[1L]]
    if (length(nums)<4L || length(nums)%%2L) stop("Invalid vector path.")
    xy <- matrix(as.numeric(nums),ncol=2L,byrow=TRUE)
    if (any(!is.finite(xy)) || any(diff(xy[,1])<0)) stop("Nonmonotone vector ridge path.")
    xy
  })
  x <- xb[1]+xb[2]*samples$az; y <- rep(Inf,length(x))
  for (p in paths) {
    if (length(unique(p[,1]))<2L) next
    ok <- x>=p[1,1] & x<=p[nrow(p),1]
    y[ok] <- pmin(y[ok],stats::approx(p[,1],p[,2],xout=x[ok],ties=min)$y)
  }
  good <- is.finite(y) & (texture | abs(samples$alt-0.01)>1e-6)
  if (sum(good)<20L || diff(range(samples$alt[good]))<(if (texture) .1 else .5) ||
      360-max(diff(c(sort(samples$az[good]),min(samples$az[good])+360)))<180)
    stop("Insufficient unclipped horizon relief or coverage for calibration.")
  if (texture && (mean(good)<.98 || stats::cor(y[good],samples$alt[good])> -.98))
    stop("The SVG ridges do not match the calibrated texture edges.")
  talt <- tan(samples$alt[good]*pi/180)
  fit <- stats::lm(y[good] ~ talt); b <- unname(stats::coef(fit))
  if (any(!is.finite(b)) || b[2]>=0) stop("Invalid vertical calibration.")
  error <- atan((y[good]-b[1])/b[2])*180/pi-samples$alt[good]
  rms <- sqrt(mean(error^2)); worst <- max(abs(error))
  if (rms>0.05 || worst>0.2) stop("Horizon and SVG do not agree within calibration tolerances.")
  # Verify interleaved holdouts to catch fits driven by a few vertices.
  k <- seq_along(talt)%%2L==0L
  f2 <- stats::lm(y[good][k] ~ talt[k]); bb <- unname(stats::coef(f2))
  hold <- atan((y[good][!k]-bb[1])/bb[2])*180/pi-samples$alt[good][!k]
  if (any(!is.finite(hold)) || max(abs(hold))>0.2) stop("Calibration holdout check failed.")
  lines <- lapply(paths,function(p) {
    # Subdivide before the nonlinear vertical conversion to preserve curved segments.
    xx <- sort(unique(c(p[,1],seq(p[1,1],p[nrow(p),1],length.out=
      max(2L,ceiling((p[nrow(p),1]-p[1,1])/xb[2]/0.1)+1L)))))
    if (length(unique(p[,1]))<2L) return(cbind((p[,1]-xb[1])/xb[2],atan((p[,2]-b[1])/b[2])*180/pi))
    yy <- stats::approx(p[,1],p[,2],xout=xx,ties=min)$y
    cbind((xx-xb[1])/xb[2],atan((yy-b[1])/b[2])*180/pi)
  })
  top <- min(89,max(samples$alt)+1)
  height <- as.integer(ceiling(width*(top+90)/360))
  if (!is.character(image.dir) || length(image.dir)!=1L || is.na(image.dir)) stop("Invalid image.dir.")
  dir.create(image.dir,recursive=TRUE,showWarnings=FALSE)
  dest <- file.path(normalizePath(image.dir,mustWork=TRUE),paste0("peakfinder-",zip.md5,"-",
    unname(tools::md5sum(svg)),"-",width,"-v3.png"))
  tmp <- tempfile(tmpdir=image.dir,fileext=".png"); on.exit(unlink(tmp),add=TRUE)
  previous <- grDevices::dev.cur()
  if (identical(Sys.info()[["sysname"]], "Darwin") && capabilities("aqua")) {
    grDevices::png(tmp, width=width, height=height, bg="transparent", type="quartz")
  } else {
    grDevices::png(tmp, width=width, height=height, bg="transparent")
  }
  if (identical(previous, grDevices::dev.cur()))
    stop("Could not open a PNG device. Check getOption('bitmapType') in R.")
  dev <- grDevices::dev.cur()
  on.exit(if (dev %in% grDevices::dev.list()) grDevices::dev.off(dev),add=TRUE)
  graphics::par(mar=rep(0,4),xaxs="i",yaxs="i")
  graphics::plot.new(); graphics::plot.window(c(0,360),c(-90,top))
  az <- c(samples$az-360,samples$az,samples$az+360)
  alt <- rep(samples$alt,3)
  keep <- c(0,samples$az[samples$az>0 & samples$az<360],360)
  alt <- stats::approx(az,alt,xout=keep)$y; az <- keep
  graphics::polygon(c(az,rev(az)),c(alt,rep(-90,length(az))),col="white",border=NA)
  for (p in lines) graphics::lines(p[,1],p[,2],col="#333333",lwd=0.5)
  graphics::lines(az,alt,col="#191919",lwd=0.7)
  grDevices::dev.off(dev)
  if (file.exists(dest) && identical(unname(tools::md5sum(tmp)),unname(tools::md5sum(dest)))) {
    unlink(tmp)
  } else if (!file.rename(tmp,dest)) {
    stop("Could not save panorama PNG (existing files are not removed on failure).")
  }
  list(path=dest,md5=unname(tools::md5sum(dest)),azimuth=c(0,360),altitude=c(-90,top),
       width=width,height=height,projection="equirectangular",style="vector line drawing",
       attribution="PeakFinder, www.peakfinder.com; rendered by skyscapeR",license="CC BY 4.0",
       calibration=list(rms.deg=rms,max.deg=worst,n=sum(good),x=unname(xb),y=b))
}

.horizon_image <- function(hor) {
  im <- hor$metadata$image
  if (is.null(im)) return(NULL)
  valid <- is.list(im) && is.character(im$path) && length(im$path)==1L &&
    !is.na(im$path) && file.exists(im$path) &&
    identical(unname(tools::md5sum(im$path)),im$md5) &&
    identical(im$projection,"equirectangular") && (identical(hor$metadata$source$provider,"Panorama") || identical(im$azimuth,c(0,360))) &&
    is.numeric(im$altitude) && length(im$altitude)==2L && all(is.finite(im$altitude)) &&
    im$altitude[1]>=-90 && im$altitude[2]>-90 && im$altitude[2]<=90
  if (!isTRUE(valid)) {
    warning("Panorama PNG is missing, changed or invalid; using numerical horizon only.",call.=FALSE)
    return(NULL)
  }
  im
}

.horizon_preview <- local({
  last <- NULL
  function(im) {
    if (!is.null(im$preview)) return(im$preview)
    key <- list(im$md5,im$source.rows)
    if (!identical(last$key,key)) {
      pixels <- .read_panorama_pixels(im)
      last <<- list(key=key,image=.panorama_window(pixels,c(1,ncol(pixels)),c(1,nrow(pixels)),1024))
    }
    last$image
  }
})

.draw_horizon <- function(hor,col.ground="darkgreen",use.image=TRUE,show.skyline=NULL,fastLoad=FALSE) {
  im <- if (use.image) .horizon_image(hor) else NULL
  if (identical(hor$metadata$source$provider,"Panorama")) return(.draw_photo(hor,im,col.ground,show.skyline,fastLoad))
  raster <- NULL
  if (!is.null(im)) {
    raster <- tryCatch({
      if (!requireNamespace("png",quietly=TRUE)) stop("Install png to plot panorama images.")
      path <- im$path
      usr <- graphics::par("usr")
      windows <- list()
      for (turn in seq(floor(usr[1]/360),floor(usr[2]/360))) {
        if (fastLoad || im$width<=8192) {
          if (!length(windows)) pixels <- if (fastLoad) .horizon_preview(im) else png::readPNG(path,native=TRUE)
          windows[[length(windows)+1L]] <- list(pixels,360*turn,im$altitude[1],360*(turn+1),im$altitude[2])
        } else {
          grid <- suppressWarnings(terra::rast(path))
          terra::ext(grid) <- c(360*turn,360*(turn+1),im$altitude)
          box <- c(max(usr[1],360*turn),min(usr[2],360*(turn+1)),
                   max(usr[3],im$altitude[1]),min(usr[4],im$altitude[2]))
          if (box[2]<=box[1] || box[4]<=box[3]) next
          # A virtual window avoids materialising the full-resolution panorama.
          terra::window(grid) <- terra::ext(box)
          bounds <- as.vector(terra::ext(grid))
          if (terra::ncell(grid)>8e6)
            grid <- terra::spatSample(grid,8e6,method="regular",as.raster=TRUE)
          preview <- tempfile(fileext=".png")
          on.exit(unlink(c(preview,paste0(preview,".aux.xml"))),add=TRUE)
          terra::writeRaster(grid,preview,datatype="INT1U",NAflag=NA)
          windows[[length(windows)+1L]] <- list(png::readPNG(preview,native=TRUE),
            bounds[1],bounds[3],bounds[2],bounds[4])
        }
      }
      windows
    },error=function(e) { warning(conditionMessage(e)," Using numerical horizon only.",call.=FALSE); NULL })
  }
  bounds <- graphics::par("usr")[1:2]
  az <- sort(unique(c(bounds,hor$data$az[hor$data$az>bounds[1] & hor$data$az<bounds[2]])))
  alt <- stats::approx(hor$data$az,hor$data$alt,xout=az,rule=2)$y
  if (is.null(raster)) graphics::polygon(c(az,rev(az)),c(alt,rep(min(-90,graphics::par("usr")[3]),length(az))),
                    col=col.ground,border=NA)
  if (!is.null(raster)) {
    for (piece in raster)
      graphics::rasterImage(piece[[1]],piece[[2]],piece[[3]],piece[[4]],piece[[5]],
                            interpolate=abs(diff(bounds))>=360)
  }
  if (isTRUE(show.skyline) || (is.null(show.skyline) && is.null(raster))) graphics::lines(az,alt,lwd=1.5)
  if (!is.null(raster) && !identical(im$attribution,"User-supplied DEM")) graphics::mtext(paste0(im$attribution," (",im$license,")"),side=3,adj=0,cex=0.5)
  invisible(!is.null(raster))
}
