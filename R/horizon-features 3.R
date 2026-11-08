#' Reconstruct or highlight archaeological structures in a DEM horizon
#' @description Explore how stones, poles and earthworks change the visible
#' landscape and the horizon used in an orientation analysis. Describe one vector
#' layer of structures, then pass one or several layers to \code{createHor}.
#' @details Reconstructions are explicit hypotheses about additional relief above
#' the supplied terrain. Heights are added to the existing DEM, including any
#' surviving remains already represented there. Ditch depths are subtracted.
#' No original ground surface or historical height is inferred. Zero height or
#' depth highlights a footprint without changing elevation.
#'
#' Every feature and multipart component is used. Points become circular
#' footprints with the supplied diameter. Lines become strips with half the width
#' on each side and round ends. Polygons use their whole area by default for
#' stones, poles, mounds and platforms, and their boundary for banks, ditches and
#' walls. Polygon holes remain holes in area mode; boundary mode includes their
#' edges. Dimensions are metres in a local azimuthal equidistant projection.
#'
#' Flat profiles extrude the footprint with vertical sides and a top following
#' the underlying DEM. Bell profiles use sin(pi * t / 2)^2, where t is distance
#' inward from the footprint edge divided by the maximum inward distance of its
#' component. For buffered points and lines the normalising distance is their
#' radius or half width. Polygon mound components use an estimated largest
#' inscribed radius, refined from an interior grid. Irregular mounds follow their
#' footprint, may have multiple summits, and holes affect the surface. This is a
#' compact rounded profile with zero slope at its edge and crest, not a Gaussian.
#'
#' Overlaps use the greatest positive and most negative displacement at each
#' location, added together. Positive heights do not accumulate. Later features
#' (then later rows within a layer) determine colour where footprints overlap.
#' NULL colour restores the terrain palette there. Colours receive terrain
#' shading and haze and are only drawn where the modified surface is visible.
#'
#' The renderer intersects footprint boundaries with individual rays, adding
#' samples at vertical faces and across curved profiles without refining the
#' entire DEM. Circles and round buffers use 128-sided approximations. Numerical
#' rays also include feature edge and central bearings. PNG columns remain a
#' finite angular grid; a structure narrower than a pixel can be missed. Increase
#' width and refine az.step and distance.step to assess convergence. The returned
#' horizon still has the package's 0.01-degree grid. These are terrain-attached
#' solids with no overhangs, openings, leaning stones, cast shadows or textures.
#' The observer elevation is unchanged by reconstruction, so supply loc's ground
#' elevation explicitly if the observer should stand on a reconstructed feature.
#'
#' Save feature objects separately for reproducibility. They contain packed vector
#' geometry and settings, but no DEM or PNG. The horizon stores compact settings,
#' feature counts and a geometry checksum, with its PNG linked externally as usual.
#' @param geometry A vector filename (including a shapefile or GeoPackage), a
#' terra SpatVector, or an object accepted by terra::vect, with a known CRS.
#' A layer must contain one geometry family: points, lines or polygons.
#' @param type Structure type: stone, pole, bank, ditch, mound, platform or wall.
#' Stones and poles accept points or polygons; mounds and platforms require
#' polygons; banks, ditches and walls accept lines or polygons.
#' @param height Additional height in metres, nonnegative. Default 0. Ignored
#' for ditches, which use depth. Positive height and depth cannot be combined.
#' @param depth Nonnegative ditch depth in metres. Default 0; only for ditches.
#' @param width Total width in metres for lines or polygon boundaries.
#' @param diameter Diameter in metres for point footprints.
#' @param profile Either bell or flat. Defaults to bell for banks, ditches and
#' mounds, and flat for stones, poles, platforms and walls.
#' @param use For polygons, area or boundary. NULL chooses the type's default.
#' @param colour One opaque R colour, or NULL for the terrain palette.
#' @param name Descriptive label for the layer. Defaults to the structure type.
#' @param height_field,depth_field,width_field,diameter_field,colour_field Optional
#' attribute column names supplying per-row settings. Do not also supply the
#' corresponding common value. Numeric fields must be numeric, finite and meet
#' the same bounds as common values. Colour fields must contain valid colours.
#' @return A serialisable skyscapeR.horFeature object. Pass it or a list of these
#' objects as features to createHor with a DEM.
#' @seealso \code{\link{createHor}}, \code{\link{demPalettes}},
#' \code{\link{exportStellarium}}, \code{\link{plot.skyscapeR.horizon}}
#' @export
#' @examples
#' if (requireNamespace("terra", quietly=TRUE)) {
#'   points <- terra::vect(data.frame(x=c(-1.75,-1.751), y=c(51.19,51.191),
#'     h=c(2,3)), geom=c("x","y"), crs="EPSG:4326")
#'   stones <- horFeature(points, "stone", diameter=0.8,
#'     height_field="h", colour="grey50")
#' }
#' \dontrun{
#' banks <- horFeature("banks.shp", "bank", width=6, height=1.2,
#'   colour="#8A9B54")
#' ditches <- horFeature("ditches.shp", "ditch", width=4, depth=0.8,
#'   colour="#A18467")
#' hor <- createHor(dem="site-dtm.tif", loc=c(51.19,-1.75),
#'   max.distance=5000, features=list(stones, banks, ditches),
#'   image.dir="horizons")
#' plot(hor)
#' }
horFeature <- function(geometry, type=c("stone","pole","bank","ditch","mound","platform","wall"),
                       height=0, height_field=NULL, depth=0, depth_field=NULL, width=NULL,
                       width_field=NULL, diameter=NULL, diameter_field=NULL, profile=NULL,
                       use=NULL, colour=NULL, colour_field=NULL, name=NULL) {
  if (!requireNamespace("terra",quietly=TRUE)) stop("Install terra to describe horizon features.")
  type <- match.arg(type)
  if (!inherits(geometry,"SpatVector")) geometry <- terra::vect(geometry)
  if (!nrow(geometry) || !nzchar(terra::crs(geometry)) ||
      any(!terra::is.valid(geometry))) stop("Supply nonempty valid feature geometry with a known CRS.")
  coordinates <- terra::geom(geometry)
  if (any(!is.finite(coordinates[,c("x","y"),drop=FALSE])) ||
      length(unique(coordinates[,"geom"]))!=nrow(geometry))
    stop("Every feature row must contain finite, nonempty geometry.")
  family <- terra::geomtype(geometry)
  allowed <- switch(type,stone=,pole=c("points","polygons"),mound=,platform="polygons",c("lines","polygons"))
  if (!family %in% allowed) stop(type," requires ",paste(allowed,collapse=" or ")," geometry.")
  if (is.null(use)) use <- if (type %in% c("bank","ditch","wall")) "boundary" else "area"
  use <- match.arg(use,c("area","boundary"))
  if (is.null(profile)) profile <- if(type %in% c("bank","ditch","mound")) "bell" else "flat"
  profile <- match.arg(profile,c("bell","flat"))
  if (is.null(name)) name <- type
  if (!is.character(name) || length(name)!=1L || is.na(name) || !nzchar(name)) stop("Supply a feature name.")
  read <- function(value,field,label,explicit,positive=FALSE,optional=FALSE) {
    if (!is.null(field)) {
      if (explicit) stop("Supply either ",label," or ",label,"_field.")
      if (!is.character(field) || length(field)!=1L || is.na(field) || !field %in% names(geometry))
        stop("Unknown ",label,"_field.")
      value <- as.data.frame(geometry)[[field]]
    } else {
      if (optional && is.null(value)) return(rep(NA_real_,nrow(geometry)))
      if (length(value)!=1L) stop(label," must be a common value or use ",label,"_field.")
      value <- rep(value,nrow(geometry))
    }
    if (label=="colour") {
      ok <- vapply(value,function(x) is.character(x) && !is.na(x) &&
        tryCatch(grDevices::col2rgb(x,alpha=TRUE)[4L]==255,error=function(e) FALSE),logical(1))
    } else {
      ok <- if(is.numeric(value)) is.finite(value) & if(positive) value>0 else value>=0 else rep(FALSE,length(value))
    }
    if (any(!ok)) stop("Invalid ",label," in feature row(s) ",paste(which(!ok),collapse=", "),".")
    value
  }
  h <- read(height,height_field,"height",!missing(height))
  d <- read(depth,depth_field,"depth",!missing(depth))
  if ((type=="ditch" && any(h!=0)) || (type!="ditch" && any(d!=0))) stop("Use depth for ditches and height for other structures.")
  w <- read(width,width_field,"width",!missing(width),TRUE,TRUE)
  diam <- read(diameter,diameter_field,"diameter",!missing(diameter),TRUE,TRUE)
  if ((family=="lines" || (family=="polygons" && use=="boundary")) && anyNA(w)) stop("Supply width for lines and polygon boundaries.")
  if (family=="points" && anyNA(diam)) stop("Supply diameter for points.")
  col <- if(is.null(colour) && is.null(colour_field)) rep(NA_character_,nrow(geometry)) else
    read(colour,colour_field,"colour",!missing(colour))
  settings <- data.frame(row=seq_len(nrow(geometry)),height=if(type=="ditch") -d else h,
                         width=w,diameter=diam,colour=col,stringsAsFactors=FALSE)
  structure(list(geometry=terra::wrap(geometry),type=type,name=name,use=use,profile=profile,
                 settings=settings),class="skyscapeR.horFeature")
}

# Each component is a small polygon, represented by oriented segments. Even/odd
# containment preserves holes. All calculations below use local metric coordinates.
.hf_edges <- function(g) {
  rings <- split(seq_len(nrow(g)),interaction(g[,"part"],g[,"hole"],drop=TRUE))
  do.call(rbind,lapply(rings,function(i) {
    xy <- g[i,c("x","y"),drop=FALSE]
    if (!all(xy[1,]==xy[nrow(xy),])) xy <- rbind(xy,xy[1,])
    cbind(head(xy,-1),tail(xy,-1))
  }))
}

.hf_inside_distance <- function(xy,e) {
  x <- xy[,1]; y <- xy[,2]; inside <- rep(FALSE,length(x)); distance <- rep(Inf,length(x))
  for (i in seq_len(nrow(e))) {
    a <- e[i,1]; b <- e[i,2]; dx <- e[i,3]-a; dy <- e[i,4]-b
    if (dy!=0) inside <- xor(inside,((b>y)!=(b+dy>y)) & (x < a+(y-b)*dx/dy))
    den <- dx*dx+dy*dy
    t <- if(den>0) pmax(0,pmin(1,((x-a)*dx+(y-b)*dy)/den)) else 0
    distance <- pmin(distance,sqrt((x-a-t*dx)^2+(y-b-t*dy)^2))
  }
  list(inside=inside | distance<1e-9,distance=distance)
}

.hf_radius <- function(e,box) {
  # Grid seeds plus bounded local optimisation, independently for each component.
  xy <- as.matrix(expand.grid(seq(box[1],box[2],length.out=25),seq(box[3],box[4],length.out=25)))
  v <- .hf_inside_distance(xy,e); score <- ifelse(v$inside,v$distance,0)
  ids <- head(order(score,decreasing=TRUE),8)
  objective <- function(p) { v <- .hf_inside_distance(matrix(p,nrow=1),e); if(v$inside) -v$distance else v$distance }
  results <- vapply(ids,function(i) {
    o <- stats::optim(xy[i,],objective,method="Nelder-Mead",control=list(maxit=300,reltol=1e-9))
    max(score[i],-o$value)
  },numeric(1))
  r <- max(results)
  if (!is.finite(r) || r<=0) stop("Cannot determine a mound's interior width. Check its geometry.")
  r
}

.hf_prepare <- function(features,crs,max.distance) {
  if (is.null(features) || identical(features,list())) return(list(items=list(),summary=list(),bearings=numeric()))
  if (inherits(features,"skyscapeR.horFeature")) features <- list(features)
  if (!is.list(features) || !length(features) || !all(vapply(features,inherits,logical(1),"skyscapeR.horFeature")))
    stop("features must be a horFeature object or a list of horFeature objects.")
  items <- summary <- list(); bearings <- numeric()
  for (layer in seq_along(features)) {
    f <- features[[layer]]; v <- terra::project(terra::unwrap(f$geometry),crs)
    tmp <- tempfile(); on.exit(unlink(tmp),add=TRUE)
    saveRDS(f,tmp,version=2); checksum <- unname(tools::md5sum(tmp)); unlink(tmp)
    summary[[layer]] <- list(name=f$name,type=f$type,use=f$use,profile=f$profile,
                            count=nrow(v),settings=f$settings,checksum=checksum)
    for (i in seq_len(nrow(v))) {
      s <- f$settings[i,]; shape <- v[i,]; family <- terra::geomtype(shape)
      original <- terra::geom(shape)
      original.theta <- atan2(original[,"x"],original[,"y"])*180/pi %% 360
      radius <- NULL
      if (family=="points") {
        radius <- s$diameter/2
        # Older terra versions buffer only the first coordinate of a multipoint.
        shape <- terra::vect(original[,c("x","y"),drop=FALSE],crs=crs)
        shape <- terra::buffer(shape,width=radius,quadsegs=32)
      } else if (family=="lines" || f$use=="boundary") {
        radius <- s$width/2
        if (family=="polygons") shape <- terra::as.lines(shape)
        shape <- terra::buffer(shape,width=radius,quadsegs=32)
      }
      g <- terra::geom(shape)
      for (component in split(seq_len(nrow(g)),interaction(g[,"geom"],g[,"part"],drop=TRUE))) {
        edges <- .hf_edges(g[component,,drop=FALSE])
        box <- c(range(edges[,c(1,3)]),range(edges[,c(2,4)]))
        nearest <- sqrt(max(box[1],0,-box[2])^2+max(box[3],0,-box[4])^2)
        if (nearest>max.distance) next
        item <- list(edges=edges,box=box,height=s$height,colour=s$colour,profile=f$profile,
                     radius=if(f$profile=="bell") if(is.null(radius)) .hf_radius(edges,box) else radius else 1,
                     label=paste(f$name,"row",i),layer=layer,row=i)
        items[[length(items)+1L]] <- item
        if (s$height!=0) bearings <- c(bearings,original.theta)
        # Include edges just outside and inside the angular footprint, plus centre.
        theta <- atan2(edges[,1],edges[,2])*180/pi %% 360
        ordered <- sort(unique(theta)); gaps <- diff(c(ordered,ordered[1]+360)); start <- which.max(gaps)
        arc <- (theta-ordered[start]) %% 360; lo <- min(arc[arc>0],360); hi <- 360
        if (length(ordered)==1L) { lo <- 0; hi <- 0 }
        centre <- (ordered[start]+(lo+hi)/2) %% 360
        if (s$height!=0) bearings <- c(bearings,(ordered[start]+c(lo-1e-6,lo+1e-6,hi-1e-6,hi+1e-6)) %% 360,centre)
      }
    }
  }
  list(items=items,summary=summary,bearings=sort(unique(bearings)))
}

.hf_values <- function(xy,items) {
  positive <- negative <- numeric(nrow(xy)); colour <- rep(NA_character_,nrow(xy))
  for (f in items) {
    b <- f$box
    ids <- which(xy[,1]>=b[1]-1e-8 & xy[,1]<=b[2]+1e-8 & xy[,2]>=b[3]-1e-8 & xy[,2]<=b[4]+1e-8)
    if (!length(ids)) next
    v <- .hf_inside_distance(xy[ids,,drop=FALSE],f$edges); ids <- ids[v$inside]
    if (!length(ids)) next
    t <- pmin(1,v$distance[v$inside]/f$radius)
    h <- f$height * if(f$profile=="bell") sin(pi*t/2)^2 else 1
    positive[ids] <- pmax(positive[ids],h); negative[ids] <- pmin(negative[ids],h)
    colour[ids] <- f$colour
  }
  list(height=positive+negative,colour=colour)
}

.hf_ray_samples <- function(bearing,distance,items) {
  direction <- c(sin(bearing*pi/180),cos(bearing*pi/180)); extra <- numeric(); active <- list()
  cross <- function(x,y) x[,1]*y[2]-x[,2]*y[1]
  for (f in items) {
    b <- f$box; near <- 0; far <- max(distance)
    for (axis in 1:2) {
      limits <- b[if(axis==1) 1:2 else 3:4]
      if (abs(direction[axis])<1e-14) {
        if (limits[1]>0 || limits[2]<0) { far <- -1; break }
      } else {
        bounds <- sort(limits/direction[axis]); near <- max(near,bounds[1]); far <- min(far,bounds[2])
      }
    }
    if (far<near) next
    e <- f$edges; a <- e[,1:2,drop=FALSE]; delta <- e[,3:4,drop=FALSE]-a
    den <- cross(delta,direction)
    good <- abs(den)>1e-12
    u <- rep(NA_real_,nrow(e)); u[good] <- -cross(a,direction)[good]/den[good]
    intersections <- a+delta*u
    d <- as.vector(intersections %*% direction)
    hits <- sort(unique(d[good & is.finite(u) & u>=-1e-9 & u<=1+1e-9 & d>=0 & d<=max(distance)]))
    breaks <- sort(unique(c(0,hits,max(distance))))
    mid <- (head(breaks,-1)+tail(breaks,-1))/2
    inside <- .hf_inside_distance(outer(mid,direction),e)$inside
    if (!any(inside)) next
    active[[length(active)+1L]] <- f
    for (j in which(inside)) {
      lo <- breaks[j]; hi <- breaks[j+1L]
      eps <- min(1e-5,(hi-lo)/1000)
      extra <- c(extra,lo-eps,lo,lo+eps,hi-eps,hi,hi+eps,
                 if(f$profile=="bell") c(seq(lo,hi,length.out=129),
                   lo+seq(0,min(hi-lo,2*f$radius),length.out=129),
                   hi-seq(0,min(hi-lo,2*f$radius),length.out=129)) else (lo+hi)/2)
    }
  }
  samples <- sort(unique(c(distance,extra[extra>0 & extra<=max(distance)])))
  list(distance=samples,items=active,direction=direction)
}

# Feature-aware tracing keeps only one ray and two offset samples in memory.
.hf_trace <- function(dem,local.crs,bearings,distance,observer,radius,profile,
                      terrain.range,style,features,ids.all=seq_along(bearings),patat=NULL) {
  render <- !is.null(profile); highest <- numeric(length(ids.all)); radians <- pi/180
  if (render) {
    count <- length(bearings)
    land.count <- if(is.null(style$sea_level)) style$palette_n else style$palette_n-style$n_sea
    land <- .dem_ramp(style$land_cols,land.count)
    sea <- if(!is.null(style$sea_level)) .dem_ramp(style$sea_cols,style$n_sea)
    a <- style$light.altitude*radians; b <- style$light.azimuth*radians
    sun <- c(cos(a)*sin(b),cos(a)*cos(b),sin(a))
    pixels <- matrix(0L,length(ids.all),count/2L)
    rows <- 90-(seq_len(ncol(pixels))-.5)*360/count
  }
  incomplete <- logical(length(ids.all))
  for (j in ids.all) {
    ray <- .hf_ray_samples(bearings[j],distance,features)
    d <- ray$distance; xy <- outer(d,ray$direction)
    values <- .hf_values(xy,ray$items)
    z <- .dem_heights(dem,xy,local.crs)+values$height
    incomplete[j-ids.all[1]+1L] <- any(!is.finite(z))
    keep <- which(is.finite(z))
    if (!length(keep)) stop("No terrain samples inside the DEM on a bearing. Reduce distance.step or enlarge the DEM.",call.=FALSE)
    z <- z[keep]; d <- d[keep]; xy <- xy[keep,,drop=FALSE]
    values$colour <- values$colour[keep]
    if (!is.null(patat)) {
      angle <- .patat_angle(z,observer,d,radius,patat)
    } else {
    if (is.finite(radius)) {
      up <- z-observer-2*(radius+z)*sin(d/radius/2)^2
      forward <- (radius+z)*sin(d/radius)
    } else { up <- z-observer; forward <- d }
    angle <- atan2(up,forward)/radians
    }
    highest[j-ids.all[1]+1L] <- max(angle)
    if (!render) {
      exposed <- c(TRUE,diff(cummax(angle))>0)
      terrain.range <- range(terrain.range,z[exposed],finite=TRUE)
    } else {
      slope <- diff(z)/diff(d)
      along <- if (!length(slope)) 0 else c(slope[1L],(head(slope,-1L)+tail(slope,-1L))/2,tail(slope,1L))
      # A centimetre offset estimates the local surface normal without a global grid.
      side <- c(ray$direction[2],-ray$direction[1])*.01
      left <- sweep(xy,2,side,"-"); right <- sweep(xy,2,side,"+")
      pair <- rbind(left,right)
      side.z <- .dem_heights(dem,pair,local.crs)+.hf_values(pair,features)$height
      across <- (tail(side.z,length(d))-head(side.z,length(d)))/.02
      across[!is.finite(across)] <- 0
      east <- along*ray$direction[1]+across*ray$direction[2]
      north <- along*ray$direction[2]-across*ray$direction[1]
      light <- .dem_light(east,north,along,d,z,observer,ray$direction[1],ray$direction[2],sun,style)
      # Trace each image ray directly, so interpolation cannot erase a narrow solid.
      coverage <- pmax(0,pmin(1,.5+(highest[j-ids.all[1]+1L]-rows)*count/360))
      ground <- which(coverage>0)
      pixels[j-ids.all[1]+1L,ground] <- .dem_column(angle,z,d,light,rows[ground],coverage[ground],
                                     terrain.range,style,land,sea,values$colour)
    }
  }
  if (!render) return(list(altitude=highest,terrain.range=terrain.range,incomplete=incomplete))
  pixels
}
