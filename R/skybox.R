#' Capture and step through the Stellarium sky above a horizon
#' @description Capture an equirectangular sky panorama at the observing site and display it behind a horizon.
#' @param hor Site horizon with latitude, longitude and elevation in metres.
#' @param jd Numeric UT Julian day. Supply either jd or time.
#' @param time A date and clock time such as "2026/12/21 16:00", or a POSIXct or POSIXlt instant.
#' @param timezone For character time, a zone name or fixed UTC offset as in time2jd(). Defaults to the package setting.
#' @param calendar For character time, the Gregorian or Julian calendar as in time2jd(). Defaults to the package setting.
#' @param display Named logical switches: atmosphere, landscape, fog, cardinals, compass.marks (numbered horizon bearings), azimuthal.grid, equatorial.grid, constellation.lines, constellation.labels (names) and constellation.art. Use galactic.equator to toggle the Galactic equator independently of the equatorial grid. Use ecliptic and equator for the lines of date, poles for the celestial poles and zenith for the Stellarium zenith/nadir markers. Older Stellarium versions may require the Compass Marks plugin to be loaded to show the bearings. Landscape and fog default to FALSE so the site horizon supplies the foreground. Other omitted switches use the current Stellarium setting.
#' @param width Even screenshot width in pixels, up to 8192. The default \code{NULL} uses the size configured in Stellarium. Try 4096 if a larger capture fails. For a narrow sector, \code{\link{sketchSky}} gives scalable vector detail.
#' @param screenshot.dir Existing local directory where Stellarium writes screenshots.
#' @param port Local Remote Control server port.
#' @param timeout Maximum wait for the capture in seconds.
#' @param sky A background returned by sketchSky(), stellariumSky() or skyStep().
#' @param minutes Signed elapsed minutes from the captured instant. Fractions are allowed.
#' @param ... Arguments passed to plot(hor, sky=...), such as xrange.
#' @details Enable Remote Control under Configuration > Plugins in Stellarium Desktop, select Load at startup and restart if needed. Enable the server in the plugin settings. Use the local server without a password and install the R packages httr2 and png. Close dialogs and use a normal rectangular viewport. Point screenshot.dir at the directory configured in Stellarium.
#'  
#' The script freezes time and takes one full-sphere Cylinder Fill screenshot. The centre is south, the seam is north, and the horizon is halfway up. Choose width for more pixels per degree. R keeps the native pixels for each plotted sector. For a 30-degree sector a width of 8192 supplies about 683 source pixels. The physical plot size and graphics device determine the displayed resolution. If the script times out, check the screenshot path and that your Stellarium version offers Cylinder Fill.
#'  
#' Capture restores the requested display switches, location, time and time rate after taking the picture. The GUI is shown again. Each capture leaves one PNG in the screenshot directory. Replot an existing sky to inspect another bearing without recapturing.
#'  
#' To move through time, pass either kind of sky to \code{skyStep()}. For a natural \code{\link{sketchSky}} background, each step retains the Sun and Moon display settings and recalculates illumination from the enabled sources. Positive minutes advance the clock and negative minutes rewind it. Assign the result back to the sky variable to continue from the new instant. Repeat \code{xrange} and any other plot limits to keep the framing. A sketch is recalculated locally with its style, display switches, constellations and magnitude limit. Older sketches without saved catalogue or calendar settings use the defaults. A Stellarium step takes a fresh screenshot using the saved capture settings.
#' @return stellariumSky() returns a skyscapeR.sky object with the native raster, PNG path, Julian day, observer location and capture settings. skyStep() invisibly returns the updated sky of the same kind.
#' @examples
#' \dontrun{
#' hor <- malta_horizons[[1]]
#' sector <- c(220, 240)
#' sky <- stellariumSky(hor, time = "2026/12/21 16:00", timezone = "UTC",
#'                     display = c(atmosphere = TRUE, cardinals = FALSE))
#' plot(hor, sky = sky, xrange = sector)
#' sky <- skyStep(hor, sky, minutes = 2, xrange = sector)
#' sky <- skyStep(hor, sky, minutes = -0.5, xrange = sector)
#' }
#' @seealso \code{\link{time2jd}}, \code{\link{plot.skyscapeR.horizon}}
#' @export
stellariumSky <- function(hor, jd=NULL, screenshot.dir=getOption("skyscapeR.screenshot.dir",
    "~/Pictures/Stellarium"), port=8090L, timeout=60, time=NULL, timezone, calendar,
    display=c(landscape=FALSE, fog=FALSE), width=NULL) {
  if (!inherits(hor,"skyscapeR.horizon")) stop("Supply a site horizon.")
  if (!requireNamespace("httr2",quietly=TRUE) || !requireNamespace("png",quietly=TRUE))
    stop("Install httr2 and png to capture a Stellarium sky.")
  loc <- as.numeric(hor$metadata$georef)
  if (length(loc)!=3 || any(!is.finite(loc))) stop("The horizon must record latitude, longitude and elevation.")
  jd <- .sky_instant(jd,time,timezone,calendar)
  if (!is.numeric(timeout) || length(timeout)!=1L || !is.finite(timeout) || timeout<=0)
    stop("Supply a positive timeout in seconds.")
  if (length(port)!=1L || !is.numeric(port) || !is.finite(port) || port!=floor(port) ||
      port<1 || port>65535) stop("Supply a local port between 1 and 65535.")
  if (!is.null(width) && (!is.numeric(width) || length(width)!=1L || !is.finite(width) ||
      width<2 || width>.Machine$integer.max || width%%2!=0))
    stop("width must be a positive even number of pixels.")
  if (!is.null(width) && width>8192)
    stop("Requested width exceeds the 8192-pixel capture limit. Use width=8192 or sketchSky() for a narrow sector.")
  dir <- normalizePath(screenshot.dir,mustWork=TRUE)
  req <- function(path) httr2::req_timeout(httr2::request(
    sprintf("http://127.0.0.1:%d/api/scripts/%s",as.integer(port),path)),timeout)
  prefix <- basename(tempfile("skyscapeR-"))
  view <- NULL
  if (!is.null(width)) {
    url <- sprintf("http://127.0.0.1:%d/api/stelproperty/list", as.integer(port))
    properties <- names(httr2::resp_body_json(httr2::req_perform(
      httr2::req_timeout(httr2::request(url),timeout))))
    for (module in c("StelMainView","MainView")) {
      ids <- paste0(module, ".", c("flagUseCustomScreenshotSize","customScreenshotWidth","customScreenshotHeight"))
      if (all(ids %in% properties)) view <- module
    }
    if (is.null(view)) stop("This Stellarium server does not expose custom screenshot dimensions. Use width=NULL and set the size in Stellarium.")
  }
  script <- .stellarium_script(loc, jd, prefix, display, width=width, view=view)
  if (isTRUE(httr2::resp_body_json(httr2::req_perform(req("status")))$scriptIsRunning))
    stop("Stellarium is running another script.")
  files <- file.path(dir,paste0(prefix,".png"))
  reply <- httr2::resp_body_string(httr2::req_perform(httr2::req_body_form(req("direct"),code=script)))
  if (!startsWith(tolower(trimws(reply)),"ok")) stop("Stellarium declined capture: ",reply)
  start <- Sys.time()
  idle.since <- NULL
  repeat {
    running <- tryCatch(httr2::resp_body_json(httr2::req_perform(req("status")))$scriptIsRunning,
      error=function(e) stop("Lost contact with Stellarium during capture. Check whether Desktop is still running, then retry with a smaller width. ",conditionMessage(e),call.=FALSE))
    if (!is.logical(running) || length(running)!=1L || is.na(running))
      stop("Stellarium returned an unrecognised script status.")
    if (!running) {
      if (file.exists(files)) {
        image <- tryCatch(png::readPNG(files,native=TRUE), error=function(e) NULL)
        if (!is.null(image)) break
      }
      if (is.null(idle.since)) idle.since <- Sys.time()
    } else idle.since <- NULL
    if (as.numeric(difftime(Sys.time(),start,units="secs")) >= timeout ||
        (!is.null(idle.since) && as.numeric(difftime(Sys.time(),idle.since,units="secs")) >= 3)) {
      if (running) stop("Stellarium is still running its script after ",timeout," seconds.")
      if (file.exists(files)) stop("Stellarium finished but R could not read the PNG: ",files)
      stop("Stellarium is idle and R found no capture at ",files,
           ". Check its configured screenshot folder against screenshot.dir. If they match, check the script console for an error.")
    }
    Sys.sleep(.25)
  }
  if (!is.null(width) && !all(dim(image)==c(width/2,width)))
    stop(sprintf("Requested %d x %d pixels, but Stellarium saved %d x %d. Try a smaller width. PNG: %s",
      width,width/2,ncol(image),nrow(image),files))
  structure(list(image=image, files=files, jd=jd, loc=loc,
    settings=list(screenshot.dir=dir, port=port, timeout=timeout, display=display, width=width)),
    class="skyscapeR.sky")
}

#' @rdname stellariumSky
#' @export
skyStep <- function(hor, sky, minutes=1, ...) {
  if (!inherits(sky,"skyscapeR.sky")) stop("Start with sketchSky() or stellariumSky().")
  if (length(minutes)!=1L || !is.numeric(minutes) || !is.finite(minutes))
    stop("Supply one signed step in minutes.")
  args <- c(list(hor=hor, jd=sky$jd + minutes/1440), sky$settings)
  if (inherits(sky,"skyscapeR.sketch")) {
    args$style <- if (identical(sky$style,"publication")) "publication-mono" else sky$style
    args$display <- sky$display
    if (!is.null(args[["const"]])) {
      args$display <- as.list(args$display)
      args$display$constellations <- args[["const"]]
      args[["const"]] <- NULL
    }
    next.sky <- do.call(sketchSky,args)
  } else {
    next.sky <- do.call(stellariumSky,args)
  }
  graphics::plot(hor, sky=next.sky, ...)
  invisible(next.sky)
}

.sky_instant <- function(jd,time,timezone,calendar) {
  if (!is.null(time)) {
    if (!is.null(jd)) stop("Supply either time or jd.")
    jd <- if (inherits(time, "POSIXt")) as.numeric(as.POSIXct(time))/86400 + 2440587.5 else
      time2jd(time, timezone, calendar)
  }
  if (!is.numeric(jd) || length(jd)!=1L || !is.finite(jd)) stop("Supply one instant using time or jd.")
  jd
}

.stellarium_script <- function(loc, jd, prefix, display, path=system.file(
    "scripts", "skyscape-panorama.ssc", package="skyscapeR", mustWork=TRUE), width=NULL, view=NULL) {
  methods <- c(landscape="LandscapeMgr.FlagLandscape", fog="LandscapeMgr.FlagFog",
    atmosphere="LandscapeMgr.FlagAtmosphere", cardinals="LandscapeMgr.FlagCardinalsPoints",
    ecliptic="GridLinesMgr.FlagEclipticLine", equator="GridLinesMgr.FlagEquatorLine",
    galactic.equator="GridLinesMgr.FlagGalacticEquatorLine",
    poles="GridLinesMgr.FlagCelestialPoles", zenith="GridLinesMgr.FlagZenithNadir",
    azimuthal.grid="GridLinesMgr.FlagAzimuthalGrid", equatorial.grid="GridLinesMgr.FlagEquatorGrid",
    constellation.lines="ConstellationMgr.FlagLines", constellation.labels="ConstellationMgr.FlagLabels",
    constellation.art="ConstellationMgr.FlagArt", compass.marks="captureCompass.Marks")
  if (is.list(display) && length(display) && all(vapply(display,function(x)
      is.logical(x) && length(x)==1L && !is.na(x),logical(1)))) display <- unlist(display)
  if (length(display) && (!is.logical(display) || anyNA(display) || is.null(names(display)) ||
      anyDuplicated(names(display)) || !all(names(display) %in% names(methods))))
    stop("display must contain uniquely named TRUE or FALSE switches. See ?stellariumSky.")
  flags <- c(landscape=FALSE, fog=FALSE)
  if (length(display)) flags[names(display)] <- display
  getter <- sub(".", ".get", methods[names(flags)], fixed=TRUE)
  setter <- sub(".", ".set", methods[names(flags)], fixed=TRUE)
  if ("cardinals" %in% names(flags)) {
    i <- match("cardinals", names(flags))
    # Stellarium renamed the cardinal-point slots. Select a matching pair before capture.
    getter[i] <- "getCardinals"
    setter[i] <- "setCardinals"
  }
  fmt <- function(x) format(x,digits=16,scientific=FALSE,decimal.mark=".",trim=TRUE)
  script <- paste(readLines(path),collapse="\n")
  if (!is.null(width)) {
    ids <- paste0(view, ".", c("customScreenshotWidth","customScreenshotHeight","flagUseCustomScreenshotSize"))
    sizes <- c(fmt(width),fmt(width/2),"true")
    for (i in seq_along(ids)) {
      script <- sub("@SAVE@", sprintf('var size%d=core.getStelProperty("%s");\n@SAVE@',i,ids[i]),script,fixed=TRUE)
      script <- sub("@APPLY@", sprintf('core.setStelProperty("%s",%s);\n@APPLY@',ids[i],sizes[i]),script,fixed=TRUE)
      script <- sub("@RESTORE@", sprintf('core.setStelProperty("%s",size%d);\n@RESTORE@',ids[i],i),script,fixed=TRUE)
    }
  }
  for (key in c("LAT","LON","ELEV","JD","PREFIX","SAVE","APPLY","RESTORE")) {
    val <- switch(key,LAT=fmt(loc[1]),LON=fmt(loc[2]),ELEV=fmt(round(loc[3])),JD=fmt(jd),PREFIX=prefix,
      SAVE=paste(sprintf("var saved%d=%s();",seq_along(flags),getter),collapse="\n"),
      APPLY=paste(sprintf("%s(%s);",setter,tolower(flags)),collapse="\n"),
      RESTORE=paste(sprintf("%s(saved%d);",setter,seq_along(flags)),collapse="\n"))
    script <- gsub(paste0("@",key,"@"),val,script,fixed=TRUE)
  }
  script
}

.draw_stellarium_sky <- function(sky) {
  if (inherits(sky, "skyscapeR.sky")) {
    bounds <- graphics::par("usr")
    span <- abs(diff(bounds[1:2]))
    available <- ncol(sky$image)*span/360
    panel <- graphics::par("pin")[1]*grDevices::dev.size("px")[1]/grDevices::dev.size("in")[1]
    if (span < 180 && is.finite(panel) && available < panel)
      warning(sprintf("This sector contains %.0f Stellarium pixels across a %.0f-pixel plot. %s",
        available,panel,if (ncol(sky$image)<8192) "Increase capture width (up to 8192) or use sketchSky()." else "Use sketchSky() for a sharp narrow view."),call.=FALSE)
    for (turn in seq(floor(min(bounds[1:2])/360), ceiling(max(bounds[1:2])/360)-1)) {
      img <- sky$image
      w <- ncol(img)
      h <- nrow(img)
      x <- pmax(0,pmin(w,(sort(bounds[1:2])/360-turn)*w))
      y <- pmax(0,pmin(h,(90-rev(sort(bounds[3:4])))/180*h))
      left <- floor(x[1]); right <- ceiling(x[2])
      top <- floor(y[1]); bottom <- ceiling(y[2])
      if (left>=right || top>=bottom) next
      if (left>0 || right<w || top>0 || bottom<h) {
        cols <- seq(left+1,right)
        rows <- seq(top,bottom-1)
        pixels <- integer(length(cols)*length(rows))
        offsets <- seq_along(cols)
        for (i in seq_along(rows))
          pixels[offsets+(i-1)*length(cols)] <- img[cols+rows[i]*w]
        img <- pixels
        dim(img) <- c(length(rows),length(cols))
        class(img) <- "nativeRaster"
      }
      graphics::rasterImage(img,360*(turn+left/w),90-bottom*180/h,
        360*(turn+right/w),90-top*180/h,interpolate=FALSE)
    }
    return(invisible(NULL))
  }
  if (!inherits(sky,"skyscapeR.skybox")) stop("Supply a capture from stellariumSky().")
  bounds <- graphics::par("usr")
  size <- pmax(2L, pmin(4096L,round(graphics::par("pin")*grDevices::dev.size("px")/grDevices::dev.size("in"))))
  az <- (bounds[1]+(seq_len(size[1])-.5)*diff(bounds[1:2])/size[1])*pi/180
  alt <- (bounds[4]-(seq_len(size[2])-.5)*diff(bounds[3:4])/size[2])*pi/180
  dims <- dim(sky$tiles[[1]])
  out <- integer(prod(size))
  for (i in seq_along(alt)) {
    east <- cos(alt[i])*sin(az); north <- cos(alt[i])*cos(az); up <- sin(alt[i])
    face <- ifelse(abs(north)>=abs(east),ifelse(north>=0,1L,3L),ifelse(east>=0,2L,4L))
    face[abs(up)>pmax(abs(north),abs(east))] <- if (up>=0) 5L else 6L
    for (f in unique(face)) {
      j <- which(face==f); e <- east[j]; n <- north[j]
      depth <- switch(f,n,e,-n,-e,up,-up)
      u <- switch(f,e,-n,-e,n,-n,-n)/depth
      v <- switch(f,up,up,up,up,-e,e)/depth
      x <- pmin(dims[2],pmax(1,floor((dims[2]+u*dims[1])/2)+1))
      y <- pmin(dims[1],pmax(1,floor((1-v)*dims[1]/2)+1))
      out[(i-1)*size[1]+j] <- sky$tiles[[f]][(y-1)*dims[2]+x]
    }
  }
  dim(out) <- rev(size); class(out) <- "nativeRaster"
  graphics::rasterImage(out,bounds[1],bounds[3],bounds[2],bounds[4],interpolate=FALSE)
}
