# Private per-process state for DEM cluster workers.
.dem_worker_state <- new.env(parent=emptyenv())

# Each worker keeps a bounded set of native-grid tiles across neighbouring rays.
# A one-cell border supplies all four bilinear neighbours at tile seams.
.dem_tile_cache <- function(dem) {
  cache <- new.env(parent=emptyenv())
  cache$tiles <- new.env(parent=emptyenv())
  cache$order <- character()
  cache$extent <- as.vector(terra::ext(dem))
  cache$res <- terra::res(dem)
  cache$nc <- ncol(dem)
  cache$nr <- nrow(dem)
  cache
}

.dem_tile_heights <- function(dem,xy,cache) {
  x <- (xy[,1]-cache$extent[1])/cache$res[1]
  y <- (cache$extent[4]-xy[,2])/cache$res[2]
  col <- floor(x-.5)
  row <- floor(y-.5)
  inside <- which(col>=0 & row>=0 & col+1<cache$nc & row+1<cache$nr)
  out <- rep(NA_real_,length(x))
  fallback <- which((x>=0 & x<=cache$nc & y>=0 & y<=cache$nr) &
                    (col<0 | row<0 | col+1>=cache$nc | row+1>=cache$nr))
  tile.cols <- ceiling(cache$nc/256)
  tile <- (row[inside]%/%256)*tile.cols+col[inside]%/%256
  for (idx in split(inside,tile)) {
    left <- col[idx[1]]%/%256*256
    top <- row[idx[1]]%/%256*256
    key <- paste(left,top)
    width <- min(257,cache$nc-left)
    values <- cache$tiles[[key]]
    if (is.null(values)) {
      values <- terra::readValues(dem,row=top+1,nrows=min(257,cache$nr-top),
                                  col=left+1,ncols=width,mat=FALSE)
      if (length(cache$order)>=128L) {
        rm(list=cache$order[1L],envir=cache$tiles)
        cache$order <- cache$order[-1L]
      }
      cache$tiles[[key]] <- values
    }
    cache$order <- c(cache$order[cache$order!=key],key)
    i <- (row[idx]-top)*width+col[idx]-left+1
    a <- values[i]; b <- values[i+1]; c <- values[i+width]; d <- values[i+width+1]
    dx <- x[idx]-.5-col[idx]; dy <- y[idx]-.5-row[idx]
    out[idx] <- (1-dy)*((1-dx)*a+dx*b)+dy*((1-dx)*c+dx*d)
    fallback <- c(fallback,idx[!is.finite(out[idx])])
  }
  if (length(fallback)) out[fallback] <- .dem_native_heights(dem,xy[fallback,,drop=FALSE])
  out
}

.dem_native_heights <- function(dem,native) {
  cell <- terra::extract(dem,native,method="simple")[[1L]]
  height <- terra::extract(dem,native,method="bilinear")[[1L]]
  height[!is.finite(cell)] <- NA_real_
  height
}

# Transform the sample locations to the native DEM grid. The local azimuthal
# equidistant projection preserves geodesic distances and bearings from its origin.
.dem_heights <- function(dem, xy, local.crs, site=NULL, allow.missing=FALSE) {
  native <- terra::project(xy, local.crs, terra::crs(dem))
  if (any(!is.finite(native))) stop("Cannot project the requested DEM sample locations.")
  cache <- attr(dem,"skyscape.tile.cache",exact=TRUE)
  if (is.null(site) && !is.null(cache)) return(.dem_tile_heights(dem,native,cache))
  cell <- terra::extract(dem, native, method="simple")[[1L]]
  height <- terra::extract(dem, native, method="bilinear")[[1L]]
  bad <- which(!is.finite(cell) | !is.finite(height))
  if (allow.missing && length(bad)) return(NA_real_)
  if (is.null(site) && length(bad)) {
    height[bad] <- NA_real_
    return(height)
  }
  if (length(bad)) {
    i <- bad[1L]
    outside <- is.na(terra::cellFromXY(dem,native[i,,drop=FALSE]))
    reason <- if (outside) "outside the raster extent" else if (!is.finite(cell[i])) "in a cell with missing elevation" else "in a cell without a finite bilinear elevation"
    where <- if (is.null(site)) sprintf("sample %.0f m from the observer",sqrt(sum(xy[i,]^2))) else
      sprintf("observer (latitude %.6f, longitude %.6f)",site[1],site[2])
    stop(sprintf("DEM %s is %s (DEM x=%.3f, y=%.3f). %s",where,reason,native[i,1],native[i,2],
      if (is.null(site)) "Check DEM coverage along the ray." else "Check loc order (latitude, longitude), DEM CRS, footprint and NoData at the observer."),call.=FALSE)
  }
  height
}

.dem_style <- function(terrain.palette, land_cols, sea_cols, sea_level, n_sea, az, alt, ambient, haze, palette_n) {
  scalar <- function(x) is.numeric(x) && length(x)==1L && is.finite(x)
  colours <- function(x) {
    if (!is.character(x) || !length(x) || anyNA(x)) return(FALSE)
    tryCatch(all(grDevices::col2rgb(x,alpha=TRUE)[4,]==255),error=function(e) FALSE)
  }
  if (is.null(terrain.palette)) stop("Supply a named terrain.palette.")
  preset <- demPalettes(terrain.palette)
  if (is.null(land_cols)) land_cols <- preset$land_cols
  if (is.null(sea_cols)) sea_cols <- preset$sea_cols
  if (!colours(land_cols) || !colours(sea_cols)) stop("Supply opaque land_cols and sea_cols.")
  if (!is.null(sea_level) && !scalar(sea_level)) stop("sea_level must be finite metres or NULL.")
  if (!scalar(palette_n) || palette_n!=floor(palette_n) || palette_n<4 || palette_n>65536)
    stop("palette_n must be a whole number from 4 to 65536.")
  if (!scalar(n_sea) || n_sea!=floor(n_sea) || n_sea<2 || n_sea>65534 || (!is.null(sea_level) && n_sea>palette_n-2))
    stop("n_sea must reserve at least two colours for sea and two for land.")
  if (!scalar(az) || !scalar(alt) || alt<0 || alt>90 ||
      !scalar(ambient) || ambient<0 || ambient>1 || !scalar(haze) || haze<0 || haze>1)
    stop("Supply finite lighting angles, altitude from 0 to 90, and ambient and haze from 0 to 1.")
  list(terrain.palette=terrain.palette,land_cols=land_cols,sea_cols=sea_cols,sea_level=sea_level,n_sea=as.integer(n_sea),palette_n=as.integer(palette_n),
       light.azimuth=az%%360,light.altitude=alt,ambient=ambient,haze=haze)
}

.dem_ramp <- function(colours,n) {
  ramp <- if (length(colours)==1L) rep(colours,n) else
    grDevices::colorRampPalette(colours,space="Lab")(n)
  ramp[c(1L,n)] <- colours[c(1L,length(colours))]
  grDevices::col2rgb(ramp)
}

.dem_colours <- function(z, limits, style, land, sea) {
  index <- function(z,low,high,n) {
    if (high<=low) return(rep.int(1L,length(z)))
    1L+floor((n-1L)*pmax(0,pmin(1,(z-low)/(high-low))))
  }
  level <- style$sea_level
  low <- if (is.null(level)) limits[1] else level
  if (identical(style$terrain.palette,"temperate")) low <- max(low,limits[1])
  base <- land[,index(z,low,limits[2],ncol(land)),drop=FALSE]
  if (!is.null(level)) {
    submerged <- z<=level
    water.index <- if (limits[1]>=level) rep.int(ncol(sea),sum(submerged)) else
      index(z[submerged],limits[1],level,ncol(sea))
    base[,submerged] <- sea[,water.index,drop=FALSE]
  }
  base
}

# Dot product of the surface normal with the selected illustrative light.
.dem_light <- function(east,north,along,d,z,observer,sin.az,cos.az,sun,style) {
  normal <- sqrt(1+east^2+north^2)
  solar <- sun[3]-sun[1]*east-sun[2]*north
  if (!identical(style$lighting,"relief")) return(solar/normal)
  dz <- observer-z
  range <- sqrt(d^2+dz^2)
  horizontal <- d/range
  vertical <- dz/range
  viewer <- along*horizontal+vertical
  alignment <- sun[3]*vertical-horizontal*(sun[1]*sin.az+sun[2]*cos.az)
  norm <- sqrt(.75^2+.25^2+2*.75*.25*alignment)
  pmax(.15,abs(.75*viewer+.25*solar)/(normal*norm))
}

# Paint the first visible terrain along each line of sight. A nearer sample
# hides all later samples below its elevation angle.
.dem_column <- function(angle, height, distance, light, rows, coverage, terrain.range, style, land, sea, feature.colour=NULL) {
  # All fully covered pixels below the nearest sample have identical RGBA.
  total <- length(rows)
  if (total>1L && coverage[total]==1 && rows[total]<angle[1L]) {
    constant <- which(rows<angle[1L] & coverage==1)
    last <- constant[1L]
    rows <- rows[seq_len(last)]
    coverage <- coverage[seq_len(last)]
  }
  index <- findInterval(rows, cummax(angle)) + 1L
  # At an interpolated horizon, extend the last ridge colour to the profile.
  index[index > length(angle)] <- which.max(angle)
  previous <- pmax(1L,index-1L)
  rise <- angle[index]-angle[previous]
  weight <- rep(1,length(rows))
  climbs <- rise>0
  weight[climbs] <- pmax(0,pmin(1,(rows[climbs]-angle[previous[climbs]])/rise[climbs]))
  interpolate <- function(value) value[previous]+weight*(value[index]-value[previous])
  shade <- style$ambient+(1-style$ambient)*pmax(0,interpolate(light))
  base <- .dem_colours(interpolate(height),terrain.range,style,land,sea)
  if (!is.null(feature.colour)) {
    colours <- feature.colour[index]
    marked <- !is.na(colours)
    if (any(marked)) base[,marked] <- grDevices::col2rgb(colours[marked])
  }
  haze <- style$haze*(1-exp(-interpolate(distance)/20000))
  channel <- function(i,sky) round(base[i,]*shade*(1-haze)+sky*haze)
  red <- channel(1L,184)
  green <- channel(2L,199)
  blue <- channel(3L,210)
  rgba <- red + 256*green + 65536*blue + 16777216*round(255*coverage)
  colour <- as.integer(ifelse(rgba >= 2^31,rgba-2^32,rgba))
  if (total>length(colour)) colour <- c(colour,rep.int(tail(colour,1),total-length(colour)))
  colour
}

# Probe actual points on each image ray. A distant suffix can be discarded
# only if even the highest elevation in the DEM stays below a nearer point.
# Keep two neighbours for interpolation and surface normals. A gap at that
# boundary falls back to the full ray to preserve its next valid neighbour.
.dem_screen_block <- function(dem,local.crs,az,distance,observer,radius,patat,ceiling.angle) {
  nr <- length(distance)
  xy <- function(d) cbind(as.vector(outer(d,sin(az))),as.vector(outer(d,cos(az))))
  probe <- unique(c(1L,round(exp(seq(0,log(nr),length.out=64L)))))
  h <- matrix(.dem_heights(dem,xy(distance[probe]),local.crs),nrow=length(probe))
  d <- distance[probe]
  angles <- if (!is.null(patat)) .patat_angle(h,observer,d,radius,patat) else if (is.finite(radius)) {
    atan2(h-observer-(radius+h)*(2*sin(d/radius/2)^2),(radius+h)*sin(d/radius))*180/pi
  } else atan2(h-observer,d)*180/pi
  angles[!is.finite(angles)] <- -Inf
  # Every ray, including the two shading neighbours, must clear the bound.
  lower <- apply(angles,1L,min)
  candidate <- which(cummax(lower)[findInterval(seq_len(nr),probe)] > ceiling.angle+1e-9)
  end <- nr
  if (length(candidate)) end <- min(nr,candidate[1L]+2L)
  z <- matrix(.dem_heights(dem,xy(distance[seq_len(end)]),local.crs),nrow=end)
  if (end < nr && any(!is.finite(z[(end-1L):end,]))) {
    rest <- matrix(.dem_heights(dem,xy(distance[(end+1L):nr]),local.crs),nrow=nr-end)
    return(rbind(z,rest))
  }
  z
}

.dem_chunk <- function(dem, local.crs, bearings, distance, observer, radius, profile=NULL, terrain.range=c(Inf,-Inf), style=NULL, features=list(), ids.all=seq_along(bearings),patat=NULL,ceiling.angle=NULL) {
  if (!terra::inMemory(dem)) {
    terra::readStart(dem)
    on.exit(terra::readStop(dem),add=TRUE)
  }
  if (length(features)) return(.hf_trace(dem,local.crs,bearings,distance,observer,radius,profile,terrain.range,style,features,ids.all,patat))
  render <- !is.null(profile)
  count <- length(bearings)
  nr <- length(distance)
  # Use the numerical skyline only to decide whether probing is worthwhile.
  # All decisions to discard terrain use exact points on the image rays.
  if (!is.null(ceiling.angle) && max(profile[ids.all]) <= ceiling.angle[ceiling(.85*nr)])
    ceiling.angle <- NULL
  radians <- pi/180
  incomplete <- logical(length(ids.all))
  highest <- numeric(length(ids.all))
  # One small group of rays plus two neighbours for surface normals.
  block <- max(1L, min(count, floor(50000/nr)))
  if (render) {
    land.count <- if (is.null(style$sea_level)) style$palette_n else style$palette_n-style$n_sea
    land <- .dem_ramp(style$land_cols,land.count)
    sea <- if (!is.null(style$sea_level)) .dem_ramp(style$sea_cols,style$n_sea)
    altitude <- style$light.altitude*radians
    bearing <- style$light.azimuth*radians
    sun <- c(cos(altitude)*sin(bearing),cos(altitude)*cos(bearing),sin(altitude))
    # nativeRaster uses row-major pixel storage, unlike an ordinary R matrix.
    pixels <- matrix(0L, length(ids.all), count/2L)
    rows <- 90 - (seq_len(ncol(pixels))-.5)*180/ncol(pixels)
  }
  if (is.null(patat) && is.finite(radius)) {
    central <- distance/radius
    bend <- 2*sin(central/2)^2
    arc <- sin(central)
  }
  for (first in seq.int(ids.all[1L],tail(ids.all,1L),by=block)) {
    ids <- first:min(tail(ids.all,1L),first+block-1L)
    columns <- if (render) c(first-1L,ids,tail(ids,1)+1L) else ids
    az <- bearings[(columns-1L) %% count + 1L]*radians
    # Adjacent blocks share two columns needed for lateral shading.
    reuse <- render && first!=ids.all[1L]
    if (!is.null(ceiling.angle)) {
      z <- .dem_screen_block(dem,local.crs,az,distance,observer,radius,patat,ceiling.angle)
    } else {
      new.az <- if (reuse) az[-(1:2)] else az
      xy <- cbind(as.vector(outer(distance,sin(new.az))),as.vector(outer(distance,cos(new.az))))
      fresh <- matrix(.dem_heights(dem,xy,local.crs),nrow=nr)
      z <- if (reuse) cbind(overlap,fresh) else fresh
      if (render) overlap <- z[,c(ncol(z)-1L,ncol(z)),drop=FALSE]
    }
    for (j in seq_along(ids)) {
      k <- j+as.integer(render)
      keep <- which(is.finite(z[,k]))
      if (!length(keep)) stop("No terrain samples inside the DEM on a bearing. Reduce distance.step or enlarge the DEM.",call.=FALSE)
      elevation <- z[keep,k]
      d <- distance[keep]
      incomplete[ids[j]-ids.all[1]+1L] <- length(keep)<nr
      if (!is.null(patat)) {
        angle <- .patat_angle(elevation,observer,d,radius,patat)
      } else {
      if (is.finite(radius)) {
        up <- elevation-observer - (radius+elevation)*bend[keep]
        forward <- (radius+elevation)*arc[keep]
      } else {
        up <- elevation-observer
        forward <- d
      }
      angle <- atan2(up,forward)/radians
      }
      if (!render) {
        highest[ids[j]-ids.all[1]+1L] <- max(angle)
        exposed <- c(TRUE,diff(cummax(angle))>0)
        terrain.range <- range(terrain.range,elevation[exposed],finite=TRUE)
      }
      if (render) {
        slope <- diff(elevation)/diff(d)
        along <- if (!length(slope)) 0 else c(slope[1L],(head(slope,-1L)+tail(slope,-1L))/2,tail(slope,1L))
        across <- (z[keep,k+1L]-z[keep,k-1L])/(2*d*sin(2*pi/count))
        across[!is.finite(across)] <- 0
        east <- along*sin(az[k])+across*cos(az[k])
        north <- along*cos(az[k])-across*sin(az[k])
        light <- .dem_light(east,north,along,d,elevation,observer,sin(az[k]),cos(az[k]),sun,style)
        limit <- profile[ids[j]]
        coverage <- pmax(0,pmin(1,.5+(limit-rows)*count/360))
        ground <- which(coverage>0)
        pixels[ids[j]-ids.all[1]+1L,ground] <- .dem_column(angle,elevation,d,light,rows[ground],coverage[ground],terrain.range,style,land,sea)
      }
    }
  }
  if (!render) return(list(altitude=highest,terrain.range=terrain.range,incomplete=incomplete))
  pixels
}

.dem_trace <- function(dem, local.crs, bearings, distance, observer, radius,
                       profile=NULL,terrain.range=c(Inf,-Inf),style=NULL,features=list(),
                       cluster=NULL,progress=FALSE,tile.dir=NULL,patat=NULL,ceiling.angle=NULL) {
  progress <- progress && isTRUE(getOption("skyscapeR.verbose",TRUE))
  n <- length(bearings)
  args <- list(local.crs=local.crs,bearings=bearings,distance=distance,observer=observer,
               radius=radius,profile=profile,terrain.range=terrain.range,style=style,features=features,patat=patat,ceiling.angle=ceiling.angle)
  workers <- max(1L,length(cluster))
  max.columns <- 128L
  if (!is.null(tile.dir)) {
    # Limit simultaneous tile sources in the final scanline pass. Bound each
    # pixel buffer to 64 MiB even at the maximum panorama width.
    max.columns <- max(128L,min(ceiling(n/64),floor(64*1024^2/(4*n/2))))
  }
  chunksize <- max(1L,min(max.columns,ceiling(n/(workers*4))))
  ids <- seq_len(n)
  # Trace the first percent locally so the bar can advance before any
  # parallel batch returns. Keep those real results in the final profile.
  starter <- if (progress && is.null(profile)) min(n,ceiling(n/100)) else 0L
  ids <- ids[ids>starter]
  if (progress && is.null(profile)) chunksize <- min(chunksize,max(1L,ceiling(n/(50*workers))))
  groups <- ceiling((ids-starter)/chunksize)
  if (!is.null(tile.dir)) {
    # Return small initial tiles promptly, then amortise disk writes on larger tiles.
    warm <- min(8L,chunksize)
    boundary <- min(n,warm*workers)
    groups <- ifelse(ids<=boundary,ceiling(ids/warm),workers+ceiling((ids-boundary)/chunksize))
  }
  chunks <- split(ids,groups)
  if (progress) {
    bar <- utils::txtProgressBar(max=n,style=3)
    on.exit(close(bar),add=TRUE)
  }
  render <- !is.null(profile)
  if (!render) highest <- numeric(n) else if (is.null(tile.dir)) pixels <- matrix(0L,n,n/2L)
  incomplete <- logical(n)
  if (starter) {
    for (j in seq_len(starter)) {
      part <- do.call(.dem_chunk,c(list(dem=dem,ids.all=j),args))
      highest[j] <- part$altitude
      incomplete[j] <- part$incomplete
      terrain.range <- range(terrain.range,part$terrain.range,finite=TRUE)
      utils::setTxtProgressBar(bar,j)
    }
  }
  sources <- character(length(chunks))
  if (!is.null(cluster)) parallel::clusterCall(cluster,function(args, tile.dir) {
    .dem_worker_state$args <- args
    .dem_worker_state$tile.dir <- tile.dir
    NULL
  },args,tile.dir)
  for (first in seq_along(chunks)[(seq_along(chunks)-1L) %% workers==0L]) {
    batch <- chunks[first:min(first+workers-1L,length(chunks))]
    parts <- if (is.null(cluster)) {
      list(do.call(.dem_chunk,c(list(dem=dem,ids.all=batch[[1L]]),args)))
    } else parallel::parLapply(cluster,batch,.dem_worker)
    for (i in seq_along(parts)) {
      if (render && !is.null(tile.dir)) {
        ids <- batch[[i]]
        tile <- sprintf("tile-%d.png",ids[1])
        if (is.null(cluster)) .dem_write_tile(parts[[i]],file.path(tile.dir,tile))
        sources[first+i-1L] <- sprintf(paste0('<SimpleSource><SourceFilename relativeToVRT="1">%s</SourceFilename>',
          '<SourceBand>BAND</SourceBand><SrcRect xOff="0" yOff="0" xSize="%d" ySize="%d"/><DstRect xOff="%d" yOff="0" xSize="%d" ySize="%d"/></SimpleSource>'),
          tile,length(ids),n/2,ids[1]-1L,length(ids),n/2)
      } else if (render) pixels[batch[[i]],] <- parts[[i]] else {
        incomplete[batch[[i]]] <- parts[[i]]$incomplete
        highest[batch[[i]]] <- parts[[i]]$altitude
        terrain.range <- range(terrain.range,parts[[i]]$terrain.range,finite=TRUE)
      }
    }
    if (progress) {
      done <- max(tail(batch,1)[[1L]])
      utils::setTxtProgressBar(bar,done)
    }
  }
  if (!render) return(list(altitude=highest,terrain.range=terrain.range,incomplete=incomplete))
  if (!is.null(tile.dir)) {
    vrt <- file.path(tile.dir,"image.vrt")
    bands <- vapply(1:4,function(b) paste0('<VRTRasterBand dataType="Byte" band="',b,
      '"><ColorInterp>',c("Red","Green","Blue","Alpha")[b],'</ColorInterp>',
      paste(sub("BAND",as.character(b),sources,fixed=TRUE),collapse=""),'</VRTRasterBand>'),character(1))
    writeLines(c(sprintf('<VRTDataset rasterXSize="%d" rasterYSize="%d"><GeoTransform>0,1,0,0,0,-1</GeoTransform>',n,n/2),bands,"</VRTDataset>"),vrt)
    return(vrt)
  }
  dim(pixels) <- rev(dim(pixels))
  class(pixels) <- "nativeRaster"
  attr(pixels,"channels") <- 4L
  pixels
}

.dem_write_tile <- function(pixels,path) {
  dim(pixels) <- rev(dim(pixels))
  class(pixels) <- "nativeRaster"
  attr(pixels,"channels") <- 4L
  png::writePNG(pixels,path)
  invisible(NULL)
}

.dem_worker <- function(ids) {
  pixels <- do.call(.dem_chunk,c(list(dem=.dem_worker_state$terrain,ids.all=ids),.dem_worker_state$args))
  directory <- .dem_worker_state$tile.dir
  if (is.null(directory)) return(pixels)
  .dem_write_tile(pixels,file.path(directory,sprintf("tile-%d.png",ids[1])))
  NULL
}

.dem_image <- function(dem, local.crs, distance, observer, radius, hor, image.dir,
                       width, attribution, license, style, features=list(), cluster=NULL, progress=FALSE,patat=NULL,verbose=FALSE) {
  if (!requireNamespace("png",quietly=TRUE)) stop("Install png to render a terrain panorama.")
  if (!dir.exists(image.dir)) dir.create(image.dir,recursive=TRUE,showWarnings=FALSE)
  directory <- normalizePath(image.dir,mustWork=TRUE)
  az <- (seq_len(width)-.5)*360/width
  limit <- stats::approx(hor$data$az,hor$data$alt,xout=az)$y
  tile.dir <- NULL
  if (width>8192) {
    tile.dir <- tempfile("panorama-tiles-")
    dir.create(tile.dir)
    on.exit(unlink(tile.dir,recursive=TRUE),add=TRUE)
  }
  if (!terra::inMemory(dem) && as.double(width)*length(distance)>1e6) {
    attr(dem,"skyscape.tile.cache") <- .dem_tile_cache(dem)
    if (!is.null(cluster)) parallel::clusterCall(cluster,function() {
      terrain <- .dem_worker_state$terrain
      attr(terrain,"skyscape.tile.cache") <- .dem_tile_cache(terrain)
      .dem_worker_state$terrain <- terrain
      NULL
    })
  }
  ceiling.angle <- NULL
  # Amortise the raster maximum scan and the per-ray probes on large jobs.
  # Feature heights can exceed the DEM ceiling and require the feature path.
  if (!length(features) && length(distance)>=256L &&
      as.double(width)*length(distance)>max(1e6,terra::ncell(dem)/8)) {
    if (verbose) .sky_message("Reading the DEM maximum for occlusion checks.")
    maximum <- terra::global(dem,"max",na.rm=TRUE)[1,1]
    if (is.finite(maximum) && (!is.finite(radius) || radius+observer>0)) {
      upper <- if (!is.null(patat)) .patat_angle(maximum,observer,distance,radius,patat) else if (is.finite(radius)) {
        atan2(maximum-observer-(radius+maximum)*2*sin(distance/radius/2)^2,
              (radius+maximum)*sin(distance/radius))*180/pi
      } else atan2(maximum-observer,distance)*180/pi
      ceiling.angle <- rev(cummax(rev(upper)))
      if (verbose) .sky_message("Rendering terrain columns with occlusion checks.")
    }
  }
  pixels <- .dem_trace(dem,local.crs,az,distance,observer,radius,limit,hor$metadata$source$terrain.range,style,features,cluster,progress,tile.dir,patat,ceiling.angle)
  tmp <- tempfile("dem-",tmpdir=directory,fileext=".png")
  on.exit(unlink(c(tmp,paste0(tmp,".aux.xml"))),add=TRUE)
  if (is.null(tile.dir)) png::writePNG(pixels,tmp) else {
    if (isTRUE(verbose)) .sky_message("Assembling the terrain tiles and writing the final PNG.")
    terra::writeRaster(terra::rast(pixels),tmp,filetype="PNG",datatype="INT1U",NAflag=NA,
                      memmax=.25,gdal="ZLEVEL=1")
  }
  checksum <- unname(tools::md5sum(tmp))
  stem <- gsub("[^[:alnum:]_-]+","-",trimws(hor$metadata$name))
  stem <- substr(gsub("^-+|-+$","",stem),1,80)
  if (!length(stem) || !nzchar(stem)) stem <- "horizon"
  dest <- file.path(directory,paste0(stem,".png"))
  if (file.exists(dest) && !identical(unname(tools::md5sum(dest)),checksum))
    dest <- file.path(directory,paste0(stem,"-",checksum,".png"))
  if (!file.exists(dest) || !identical(unname(tools::md5sum(dest)),checksum)) {
    if (!file.rename(tmp,dest) &&
        !(file.exists(dest) && identical(unname(tools::md5sum(dest)),checksum)))
      stop("Could not save terrain panorama in ",directory)
  }
  list(path=dest,md5=checksum,azimuth=c(0,360),altitude=c(-90,90),width=width,
       height=width/2,distance.step=distance[1L],projection="equirectangular",style="shaded DEM terrain",
       attribution=attribution,license=license,render.settings=style)
}

.dem_horizon <- function(dem, loc, name, alt.unc, max.distance, eye.height,
                         az.step, distance.step, curvature, refraction.k, panorama,
                         image.dir, width, attribution, license, style, features=list(), loc.crs=NULL,ncores=NULL,verbose=TRUE,progress=FALSE,atm,temp,tempgrad,panorama.distance.step="adaptive", panorama.acuity=1) {
  if (!requireNamespace("terra",quietly=TRUE)) stop("Install terra to construct a horizon from a DEM.")
  scalar <- function(x) is.numeric(x) && length(x)==1L && is.finite(x)
  flag <- function(x) is.logical(x) && length(x)==1L && !is.na(x)
  label <- function(x) is.character(x) && length(x)==1L && !is.na(x) && nzchar(x) && !grepl("[\r\n]",x)
  if (!flag(verbose) || !flag(progress)) stop("verbose and progress must be TRUE or FALSE.")
  if (!is.null(ncores) && (!scalar(ncores) || ncores<1 || ncores!=floor(ncores))) stop("ncores must be a positive integer.")
  if (inherits(loc,"SpatVector")) {
    if (!is.null(loc.crs) || terra::geomtype(loc)!="points" || nrow(loc)!=1L ||
        nrow(terra::crds(loc))!=1L || !nzchar(terra::crs(loc)))
      stop("Supply one observer point with a CRS, without loc.crs.")
    loc.crs <- terra::crs(loc)
    loc <- as.numeric(terra::crds(loc))
  }
  if (!is.null(loc.crs)) {
    if (!is.numeric(loc) || !length(loc) %in% 2:3 || any(!is.finite(loc)))
      stop("Supply loc as x, y and optional ground elevation for loc.crs.")
    ll <- terra::project(matrix(loc[1:2],1L),loc.crs,"OGC:CRS84")
    loc <- c(ll[1,2],ll[1,1],loc[-(1:2)])
  }
  if (!is.numeric(loc) || !length(loc) %in% 2:3 || any(!is.finite(loc)) ||
      abs(loc[1L])>=90 || abs(loc[2L])>180) stop("Supply loc as latitude, longitude and optional ground elevation.")
  if (!scalar(max.distance) || max.distance<=0) stop("Supply a positive max.distance in metres.")
  if (!scalar(eye.height) || eye.height<0) stop("eye.height must be nonnegative metres.")
  if (!scalar(alt.unc) || alt.unc<0) stop("Supply a nonnegative scalar DEM horizon uncertainty in degrees.")
  if (!scalar(az.step) || az.step<.01 || az.step>45) stop("az.step must be between 0.01 and 45 degrees.")
  if (!flag(curvature) || !flag(panorama)) stop("curvature and panorama must be TRUE or FALSE.")
  if (!is.null(refraction.k) && (!scalar(refraction.k) || refraction.k>=1 || (!curvature && refraction.k!=0)))
    stop("refraction.k must be less than one and requires curvature when nonzero.")
  if (!label(attribution) || !label(license)) stop("Supply single-line DEM attribution and license text.")
  if (panorama && (!scalar(width) || width<32 || width>81920 || width%%2!=0 || !label(image.dir)))
    stop("Supply an image directory and an even width from 32 to 81920 pixels.")
  if (verbose) .sky_message(sprintf("Preparing DEM horizon: search radius %g m (%g km).",max.distance,max.distance/1000))
  source <- "supplied SpatRaster"
  if (is.character(dem) && length(dem)==1L && file.exists(dem)) {
    source <- basename(dem)
    dem <- terra::rast(dem)
  }
  if (!inherits(dem,"SpatRaster") || terra::nlyr(dem)!=1L || !nzchar(terra::crs(dem)) ||
      isTRUE(terra::is.factor(dem))) stop("Supply a single continuous DEM band with a known CRS and elevations in metres.")
  units <- tolower(terra::units(dem))
  if (length(units) && any(units %in% c("ft","feet","foot","us survey foot")))
    stop("Convert DEM elevations to metres before calculating the horizon.")
  if (source=="supplied SpatRaster") {
    files <- terra::sources(dem)
    if (length(files)==1L && nzchar(files)) source <- basename(files)
  }
  local.crs <- sprintf("+proj=aeqd +lat_0=%.12f +lon_0=%.12f +datum=WGS84 +units=m",loc[1],loc[2])
  centre <- matrix(c(0,0),nrow=1L)
  sampled <- .dem_heights(dem,centre,local.crs,site=loc,allow.missing=length(loc)==3L)
  ground <- if (length(loc)==3L) loc[3L] else sampled
  adaptive <- panorama && identical(panorama.distance.step,"adaptive")
  if (is.null(distance.step) || adaptive) {
    xy <- as.numeric(terra::project(centre,local.crs,terra::crs(dem)))
    resolution <- terra::res(dem)
    adjacent <- rbind(xy+c(resolution[1],0),xy+c(0,resolution[2]))
    metres <- terra::project(adjacent,terra::crs(dem),local.crs)
    cell.size <- min(sqrt(rowSums(metres^2)))
    if (is.null(distance.step)) distance.step <- cell.size/2
  }
  if (!scalar(distance.step) || distance.step<=0) stop("distance.step must be positive metres.")
  n <- max(2,ceiling(max.distance/distance.step))
  if (n>250000) stop("Too many samples per ray. Increase distance.step or reduce max.distance.")
  distance <- seq_len(n)*(max.distance/n)
  image.distance <- distance
  if (adaptive) {
    if (!scalar(panorama.acuity) || panorama.acuity<=0 || panorama.acuity>=5400)
      stop("panorama.acuity must be a positive angle below 5400 arcminutes.")
    angle <- max(panorama.acuity/60,360/width)
    image.distance <- .dem_adaptive_distances(max.distance,cell.size,angle)
  } else if (panorama && !is.null(panorama.distance.step)) {
    if (!scalar(panorama.distance.step) || panorama.distance.step<=0)
      stop("panorama.distance.step must be positive metres or 'adaptive'.")
    count.image <- max(2,ceiling(max.distance/panorama.distance.step))
    if (count.image>250000) stop("Too many panorama samples per ray. Increase panorama.distance.step.")
    image.distance <- seq_len(count.image)*(max.distance/count.image)
  }
  observer <- ground+eye.height
  patat <- NULL
  if (curvature && is.null(refraction.k)) {
    if (!scalar(atm) || atm<=0 || !scalar(temp) || temp<= -273.15 || !scalar(tempgrad))
      stop("Patat refraction requires positive pressure, temperature above absolute zero and a finite temperature gradient.")
    lat <- loc[1]*pi/180
    radius <- 6378137*6356752.314245/sqrt((6378137*sin(lat))^2+(6356752.314245*cos(lat))^2)
    patat <- .252*atm*(34.2+tempgrad)/(temp+273.15)^2
  } else radius <- if (curvature) 6371008.8/(1-refraction.k) else Inf
  if (max.distance>=pi*min(6371008.8,radius)/2)
    stop("The requested distance exceeds the regional horizon model.")
  count <- ceiling(360/az.step)
  bearings <- (seq_len(count)-1)*360/count
  prepared <- .hf_prepare(features,local.crs,max.distance)
  if (length(prepared$items)) {
    bearings <- sort(unique(c(bearings,prepared$bearings)))
    if (panorama) {
      image.az <- (seq_len(width)-.5)*360/width
      missed <- vapply(prepared$items,function(f) {
        theta <- atan2(f$edges[,1],f$edges[,2])*180/pi %% 360
        centre <- atan2(mean(f$edges[,1]),mean(f$edges[,2]))*180/pi %% 360
        arc <- (theta-centre+180) %% 360-180
        pixel <- (image.az-centre+180) %% 360-180
        !any(pixel>=min(arc) & pixel<=max(arc))
      },logical(1))
      if (any(missed)) warning("Some feature footprints are narrower than the PNG sampling: ",
        paste(unique(vapply(prepared$items[missed],`[[`,character(1),"label")),collapse=", "),
        ". Increase width or move the viewpoint closer.",call.=FALSE)
    }
  }
  numerical.features <- Filter(function(f) f$height!=0,prepared$items)
  edge.az <- unique(c(bearings,if (panorama) (seq_len(width)-.5)*360/width,seq(0,359.9,by=.1)))*pi/180
  edge <- terra::project(max.distance*cbind(sin(edge.az),cos(edge.az)),local.crs,terra::crs(dem))
  partial <- anyNA(terra::cellFromXY(dem,edge))
  if (partial) warning("The DEM does not cover the entire requested viewshed. Each bearing uses the available terrain. Check coverage before interpreting the horizon.",call.=FALSE)
  workers <- .worker_count(ncores,length(bearings),verbose)
  cluster <- NULL
  if (workers>1L) {
    # Pack a file reference so each worker can read native raster cells locally.
    if (terra::inMemory(dem)) {
      path <- tempfile(fileext=".tif")
      on.exit(unlink(path),add=TRUE)
      dem <- terra::writeRaster(dem,path,datatype="FLT8S")
    }
    packed <- terra::wrap(dem,proxy=TRUE)
    cluster <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cluster),add=TRUE,after=FALSE)
    parallel::clusterCall(cluster,.libPaths,.libPaths())
    parallel::clusterCall(cluster,function(p) {
      .dem_worker_state$terrain <- terra::unwrap(p)
      NULL
    },packed)
  }
  if (verbose) .sky_message(sprintf("Tracing %d horizon bearings at %.2f m distance steps using %d worker(s).",length(bearings),max.distance/n,workers))
  trace <- .dem_trace(dem,local.crs,bearings,distance,observer,radius,features=numerical.features,cluster=cluster,progress=progress,patat=patat)
  if (any(trace$incomplete) && !partial) warning("DEM coverage has gaps. Missing terrain samples were skipped; check incomplete.bearings before interpreting the horizon.",call.=FALSE)
  hor <- createHor(bearings,trace$altitude,alt.unc=alt.unc,loc=c(loc[1:2],observer),name=name)
  hor$metadata$source <- list(provider="DEM",dataset=source,crs=terra::crs(dem),
    terrain.range=trace$terrain.range,ground.elevation=ground,sampled.ground=sampled,eye.height=eye.height,
    partial.coverage=partial || any(trace$incomplete),incomplete.bearings=bearings[trace$incomplete],max.distance=max.distance,distance.step=max.distance/n,az.step=360/count,
    curvature=curvature,refraction.k=refraction.k,patat=if (!is.null(patat)) list(radius=radius,atm=atm,temp=temp,tempgrad=tempgrad),attribution=attribution,license=license)
  if (length(prepared$summary)) hor$metadata$features <- list(
    elevation="additive",overlap="maximum positive plus minimum negative",layers=prepared$summary)
  if (adaptive && verbose) .sky_message(sprintf("Adaptive panorama: %.2f m near the observer, up to %.2f m on this ray (%d samples).",
    image.distance[1L],max(diff(c(0,image.distance))),length(image.distance)))
  if (panorama && !adaptive && verbose) .sky_message(sprintf("Rendering %d x %d terrain panorama at %.2f m distance steps.",width,width/2,image.distance[1L]))
  if (panorama) hor$metadata$image <- tryCatch(
    .dem_image(dem,local.crs,image.distance,observer,radius,hor,image.dir,width,attribution,license,style,prepared$items,cluster,progress,patat,verbose),
    error=function(e) { warning("DEM horizon created, but panorama failed: ",conditionMessage(e),call.=FALSE); NULL })
  if (adaptive && !is.null(hor$metadata$image)) {
    hor$metadata$image$sampling <- list(method="adaptive",cell.size=cell.size,
      acuity.arcmin=panorama.acuity,effective.arcmin=angle*60,
      near.step=image.distance[1L],max.step=max(diff(c(0,image.distance))),samples=length(image.distance))
  }
  if (verbose) .sky_message("Horizon complete.")
  hor
}

.patat_angle <- function(elevation,observer,d,radius,coefficient) {
  (atan((elevation-observer)/d-d/(2*radius))+coefficient*d/radius)*180/pi
}

# Use the transverse angular footprint as a proxy for radial terrain detail.
.dem_adaptive_distances <- function(radius,cell,angle) {
  if (!is.finite(cell) || cell<=0) stop("Cannot determine the DEM cell size at the observer.")
  near <- min(cell,radius/2)
  limit <- min(250000,ceiling(radius/near)+1)
  out <- numeric(limit)
  footprint <- tan(angle*pi/180)
  d <- 0
  for (i in seq_len(limit)) {
    step <- max(near,d*footprint)
    d <- min(radius,d+step)
    if (radius-d<=sqrt(.Machine$double.eps)*near) d <- radius
    out[i] <- d
    if (d>=radius) return(out[seq_len(i)])
  }
  stop("Too many adaptive panorama samples. Supply a numeric panorama.distance.step or reduce max.distance.")
}
