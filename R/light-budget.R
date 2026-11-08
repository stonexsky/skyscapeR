# One unit of integrated flux is a V=0 source above the atmosphere.
# Raster intensities use that flux per square degree; area is on the sphere.
.source_flux <- function(magnitude, altitude, k) {
  10^(-.4*magnitude) * .band_transmission(altitude,k)
}

.mw_radiance <- local({
  levels <- NULL
  function() {
    if (is.null(levels)) {
      # Provisional relative intensities; anchor their area-weighted mean over
      # the published four-degree Carina field (van der Kruit 2024, Table 3).
      ratio <- c(0,1,0,1.35,2.1)
      field <- expand.grid(l=seq(284.95,289.05,.1),b=seq(-4.05,.05,.1))
      inside <- .cluster_distance(287,-2,field$l,field$b) <= 2
      l <- field$l[inside]*pi/180
      b <- field$b[inside]*pi/180
      id <- .galactic_classes(cos(b)*cos(l),cos(b)*sin(l),sin(b),diag(3))
      levels <<- ratio * (3600^2 * 10^(-.4*20.7)) /
        stats::weighted.mean(ratio[id+1L],cos(b))
    }
    levels
  }
})

# Fixed-angular quadrature. The same sky samples are integrated into the
# requested field, including at the north seam or across the zenith.
.light_samples <- function(az,alt,radius,flux,kind,fraction=1,angle=0,step=.03,cluster.power=3) {
  n <- max(64,min(65536,ceiling(4*(radius/step)^2)))
  q <- (seq_len(n)-.5)/n
  phi <- seq_len(n)*(pi*(3-sqrt(5)))
  # Cluster profile proportional to (1-r^2/R^2)^cluster.power. Disc samples are uniform in projected area.
  r <- radius*switch(kind,cluster=sqrt(1-(1-q)^(1/(cluster.power+1))),moon=sqrt(q))
  x <- r*cos(phi)
  y <- r*sin(phi)
  if (kind=="moon") {
    u <- fraction*x+(1-fraction)*sqrt(pmax(0,radius^2-y^2))
    x <- u*cos(angle)-y*sin(angle)
    y <- u*sin(angle)+y*cos(angle)
  }
  rad <- pi/180
  d <- sqrt(x*x+y*y)*rad
  bearing <- atan2(x,y)
  a <- alt*rad
  h <- asin(pmax(-1,pmin(1,sin(a)*cos(d)+cos(a)*sin(d)*cos(bearing))))
  lon <- az+atan2(sin(bearing)*sin(d)*cos(a),cos(d)-sin(a)*sin(h))/rad
  cbind(az=lon%%360,alt=h/rad,flux=rep(flux/n,n))
}

# Linear raster before display compression. Each cell contains average surface
# brightness; multiplying by its solid angle recovers integrated source flux.
.shared_radiance <- function(sky,box,nx,ny,colour=FALSE,point.block=1L) {
  dx <- diff(box[1:2])/nx
  dy <- diff(box[3:4])/ny
  top <- box[4]-(seq_len(ny)-1)*dy
  rad <- pi/180
  row.area <- dx*(sin(pmin(90,pmax(-90,top))*rad)-
                    sin(pmin(90,pmax(-90,top-dy))*rad))/rad
  az <- rep(box[1]+(seq_len(nx)-.5)*dx,each=ny)
  alt <- rep(top-dy/2,nx)
  # K&S nanoLamberts to V=0 units per square degree.
  per.nL <- 3600^2/(34.08*exp(20.7233))
  sky.light <- .sky_luminance(sky$light,az,alt)*per.nL
  cells <- length(alt)
  mw <- numeric(cells)
  if (!is.null(sky$milky.way.frame)) {
    atm <- c(sky$atmosphere,list(jd=sky$jd,geopos=sky$loc[c(2,1,3)]))
    geom <- rep(.geometric_altitude(top-dy/2,atm),nx)*rad
    ids <- .galactic_classes(cos(geom)*sin(az*rad),cos(geom)*cos(az*rad),
                             sin(geom),sky$milky.way.frame)
    mw <- .mw_radiance()[ids+1L]*.band_transmission(alt,sky$light$k)
  }
  rays <- list()
  resolved <- numeric(cells)
  chroma <- if (colour) matrix(0,cells,2)
  tint <- function(col) {
    rgb <- grDevices::col2rgb(col)[,1]/255
    rgb <- ifelse(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)^2.4)
    rgb[c(1,3)]/sum(rgb*c(.2126,.7152,.0722))-1
  }
  budget <- sky$light.budget
  moon <- sky$bodies$moon
  add <- function(az,alt,radius,flux,kind,fraction=1,angle=0,col="white",cluster.power=3) {
    if (!is.finite(az+alt+flux+radius) || flux<=0 || alt+radius<box[3] ||
        alt-radius>box[4] || alt+radius<0 || fraction<=0) return()
    if (!is.null(moon) && kind %in% c("star","planet") &&
        .cluster_distance(moon[1],moon[2],az,alt)<sky$diameters["moon"]/2) return()
    weights <- if(colour) tint(col) else c(0,0)
    if (kind=="star") {
      # Bilinear deposition in sampling pixels. A star has no resolved limb.
      # Reflect an adjacent centre across a pole if the stencil straddles it.
      delta.az <- ((az-mean(box[1:2])+180)%%360)-180
      px <- dx*point.block
      py <- dy*point.block
      xx <- (diff(box[1:2])/2+delta.az)/px-.5
      yy <- (box[4]-alt)/py-.5
      grid <- expand.grid(x=floor(xx)+0:1,y=floor(yy)+0:1)
      w <- (1-abs(grid$x-xx))*(1-abs(grid$y-yy))
      h <- box[4]-(grid$y+.5)*py
      a <- box[1]+(grid$x+.5)*px
      if (point.block>1L) {
        # Subdivide each display sample at constant surface brightness. A high-DPI
        # export then preserves the apparent size and flux of unresolved stars.
        offsets <- seq_len(point.block)-(point.block+1)/2
        sub <- expand.grid(x=offsets,y=offsets,i=seq_along(w))
        h <- h[sub$i]+sub$y*dy
        a <- a[sub$i]+sub$x*dx
        area <- abs(cos(h*rad))
        w <- w[sub$i]*area/stats::ave(area,sub$i,FUN=sum)
      }
      pole <- abs(h)>90
      a[pole] <- a[pole]+180
      h[pole] <- sign(h[pole])*(180-abs(h[pole]))
      # Redistribute a sample below the geometric horizon into the first
      # visible row when its stellar centre is above the horizon.
      h <- abs(h)
      rays[[length(rays)+1L]] <<- cbind(a%%360,h,flux*w,weights[1],weights[2])
      return()
    }
    if (kind=="planet") kind <- "moon"
    span <- if(abs(alt)+radius>=90) 180 else asin(sin(radius*rad)/cos(alt*rad))/rad
    offset <- (az-box[1]+span)%%360
    if (offset>diff(box[1:2])+2*span) return()
    step <- max(.001,min(dy,dx*max(.01,abs(cos(rad*alt))))/2)
    # For resolved images evaluate a continuous profile. This avoids sampling
    # speckles when a zoom fills the device with a single disc or cluster.
    if (radius > 5*max(dy,dx*abs(cos(alt*rad))) && (kind!="moon" || fraction>.04)) {
      j <- which(abs(rep(top-dy/2,nx)-alt)<=radius+dy &
                 abs((rep(box[1]+(seq_len(nx)-.5)*dx,each=ny)-az+180)%%360-180)<=span+dx)
      if (!length(j)) return()
      xx <- box[1]+(floor((j-1)/ny)+.5)*dx
      yy <- top[(j-1L)%%ny+1L]-dy/2
      value <- numeric(length(j))
      for (ox in c(-.25,.25)) for (oy in c(-.25,.25)) {
        hh <- (yy+oy*dy)*rad
        da <- (xx+ox*dx-az)*rad
        east <- cos(hh)*sin(da)/rad
        up <- (sin(hh)*cos(alt*rad)-cos(hh)*sin(alt*rad)*cos(da))/rad
        rr <- (east^2+up^2)/radius^2
        profile <- switch(kind,
          cluster=(cluster.power+1)*pmax(0,1-rr)^cluster.power,
          moon={
            x <- east*cos(angle)+up*sin(angle)
            y <- -east*sin(angle)+up*cos(angle)
            (x >= (1-2*fraction)*sqrt(pmax(0,radius^2-y^2)))/fraction
          })
        value <- value+ifelse(rr<=1 & abs(hh)<=pi/2 & hh>=0,profile,0)/4
      }
      norm <- pi*radius^2
      centred <- (az-box[1])%%360
      whole <- (diff(box[1:2])>=360 || (centred>=span && centred<=diff(box[1:2])-span)) &&
        alt-radius>=max(0,box[3]) && alt+radius<=box[4]
      sampled <- sum(value*row.area[(j-1L)%%ny+1L])
      if (whole && sampled>0) norm <- sampled
      added <- value*(flux/norm)
      resolved[j] <<- resolved[j]+added
      if (colour) {
        chroma[j,1] <<- chroma[j,1]+added*weights[1]
        chroma[j,2] <<- chroma[j,2]+added*weights[2]
      }
      return()
    }
    rays[[length(rays)+1L]] <<- cbind(.light_samples(az,alt,radius,flux,kind,fraction,angle,step,cluster.power),weights[1],weights[2])
  }
  for(i in seq_len(NROW(sky$stars))) {
    col <- if (length(sky$star.colours)) sky$star.colours[i] else "white"
    add(sky$stars[i,1],sky$stars[i,2],0,budget$stars[i],"star",col=col)
  }
  for(i in seq_len(NROW(sky$clusters))) {
    cl <- sky$clusters[i,]
    add(cl$az,cl$alt,cl$radius,cl$residual.flux,"cluster",
        cluster.power=if (cl$identifier=="M 45") 1 else 3)
  }
  for (id in names(budget$bodies)) {
    p <- sky$bodies[[id]]
    if (id=="moon") {
      radius <- sky$diameters[id]/2
      covered <- .cluster_distance(p[1],p[2],az,alt)<=radius
      mw[covered] <- 0
      add(p[1],p[2],radius,budget$bodies[id],"moon",sky$lunar$fraction,sky$lunar$angle)
    } else {
      disc <- budget$planet.discs[[id]]
      if (is.null(disc)) {
        # Shared sketches saved before planet dimensions were retained.
        ph <- .sketch_pheno(sky$jd,id,sky$loc)
        solar <- body.position(sky$jd,obj="sun",loc=sky$loc,verbose=FALSE)$horizontal
        disc <- list(fraction=ph[2],angle=.lunar_light_angle(p,unlist(solar)))
        radius <- ph[4]/2
      } else radius <- sky$diameters[id]/2
      palette <- c(mercury="#DAD7D3",venus="#FFF9E8",mars="#FFAF88",
                   jupiter="#F5DDC1",saturn="#EDD9A0")
      add(p[1],p[2],radius,budget$bodies[id],"planet",disc$fraction,disc$angle,palette[[id]])
    }
  }
  direct <- resolved
  if (length(rays)) {
    r <- do.call(rbind,rays)
    # Repeat only where a requested plot spans more than one complete turn.
    for (turn in 0:floor(diff(box[1:2])/360)) {
      x <- floor(((r[,1]-box[1])%%360+360*turn)/dx)
      y <- floor((box[4]-r[,2])/dy)
      keep <- x<nx & y>=0 & y<ny & r[,2]>=0
      if (!any(keep)) next
      index <- 1L+y[keep]+ny*x[keep]
      values <- if (colour) cbind(r[keep,3],r[keep,3]*r[keep,4],r[keep,3]*r[keep,5]) else r[keep,3]
      summed <- rowsum(values,index,reorder=FALSE)
      j <- as.integer(rownames(summed))
      area <- row.area[(j-1L)%%ny+1L]
      direct[j] <- direct[j]+summed[,1]/area
      if (colour) chroma[j,] <- chroma[j,]+summed[,2:3,drop=FALSE]/area
    }
  }
  if (colour) {
    # Existing sky palette for atmospheric light. Diffuse Galactic and cluster
    # light retain a neutral tint in the absence of spectral calibration.
    sky.col <- .sky_colour(sky$light,az,alt,opaque=TRUE)
    pal <- unique(sky.col)
    tint.sky <- t(vapply(pal,tint,numeric(2)))
    chroma <- chroma+sky.light*tint.sky[match(sky.col,pal),,drop=FALSE]
  }
  list(chroma=chroma,radiance=matrix(sky.light+mw+direct,ny,nx),direct=matrix(direct,ny,nx),
       milky.way=matrix(mw,ny,nx),area=row.area)
}

.draw_shared_light <- function(sky,box) {
  dpi <- grDevices::dev.size("px")/grDevices::dev.size("in")
  block <- max(1L,min(3L,round(min(dpi)/120)))
  base <- pmin(floor(c(2400,1500)/block),pmax(96,ceiling(graphics::par("pin")*120)))
  size <- base*block
  field <- .shared_radiance(sky,box,size[1],size[2],colour=TRUE,point.block=block)
  pixels <- .shared_pixels(field,sky)
  graphics::rasterImage(pixels,box[1],box[3],box[2],box[4],interpolate=FALSE)
}

.shared_pixels <- function(field,sky) {
  # Preserve the previous neutral sRGB tone curve in linear display luminance.
  adaptation <- max(1,.sky_luminance(sky$light,0,90)/60)^.75
  lum <- as.vector(field$radiance)
  grey <- pmin(1,asinh(lum/(.04*adaptation))/asinh(2500))
  y <- ifelse(grey<=.04045,grey/12.92,((grey+.055)/1.055)^2.4)
  rb <- field$chroma*(y/pmax(lum,.Machine$double.eps))
  delta <- cbind(rb[,1],-(.2126*rb[,1]+.0722*rb[,2])/.7152,rb[,2])
  # Desaturate towards neutral where individual channels reach the gamut.
  strength <- pmin(1,(1-y)/pmax(delta[,1],delta[,2],delta[,3],1e-15),
                   y/pmax(-delta[,1],-delta[,2],-delta[,3],1e-15))
  rgb <- pmax(y+delta*strength,0)
  rgb <- ifelse(rgb<=.0031308,12.92*rgb,1.055*rgb^(1/2.4)-.055)
  matrix(grDevices::rgb(pmin(1,rgb[,1]),pmin(1,rgb[,2]),pmin(1,rgb[,3])),nrow(field$radiance))
}

#' Inspect the linear light field of a natural sky
#'
#' Samples a natural sketch before display compression or masking by local
#' terrain. One integrated flux unit is a magnitude-zero V-band source.
#' @param sky A natural sketch from \code{\link{sketchSky}}.
#' @param xrange Increasing true azimuth limits in degrees. A sector may cross
#' north, for example c(340,380).
#' @param yrange Increasing altitude limits within -90 to 90 degrees.
#' @param size Positive integer column and row counts, up to 1600 by 1000.
#' @return A list containing azimuth and altitude pixel centres, total radiance,
#' direct light (stars, residual cluster light, planets and lunar disc), and
#' Milky Way radiance. Matrices are in magnitude-zero flux units per square
#' degree, ordered from high to low altitude and increasing azimuth. Multiply
#' by area (square degrees per pixel by row) and sum to obtain flux. The
#' difference between total and the two source matrices is the atmospheric
#' background. The Sun contributes to the background; its disc is a plot overlay.
#' @details Direct light includes stars, residual cluster light, planetary discs
#' and the illuminated lunar disc. The Milky Way uses the accepted tracing and
#' a provisional Carina surface-brightness calibration. Stars remain unresolved
#' as the view is enlarged. See \code{\link{sketchSky}} for sources, assumptions
#' and interpretation.
#' @export
skyRadiance <- function(sky,xrange=c(0,360),yrange=c(0,90),size=c(720,180)) {
  if (!inherits(sky,"skyscapeR.sketch") || is.null(sky$light.budget))
    stop("Supply a current natural sketch from sketchSky().")
  valid <- function(x) is.numeric(x) && length(x)==2L && all(is.finite(x))
  if (!valid(xrange) || diff(xrange)<=0 || diff(xrange)>360 ||
      !valid(yrange) || diff(yrange)<=0 || any(abs(yrange)>90))
    stop("Supply increasing azimuth and altitude ranges covering at most one sky.")
  if (!valid(size) || any(size!=floor(size)) || any(size<1) || any(size>c(1600,1000)))
    stop("Supply integer size between c(1,1) and c(1600,1000).")
  result <- .shared_radiance(sky,c(xrange,yrange),size[1],size[2])
  result$chroma <- NULL
  result$azimuth <- xrange[1]+(seq_len(size[1])-.5)*diff(xrange)/size[1]
  result$altitude <- yrange[2]-(seq_len(size[2])-.5)*diff(yrange)/size[2]
  result
}
