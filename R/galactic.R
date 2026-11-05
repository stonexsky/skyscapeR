# Use two non-collinear directions as a basis for the Swiss Ephemeris
# frame rotation. Both frames are evaluated at the SAME instant, so the
# catalogue stars' proper motions do not impart motion to the diffuse band.
.galactic_paths <- function(jd, outline, equator) {
  normalise <- function(x) x / sqrt(sum(x^2))
  basis <- function(j2000) {
    direction <- function(id) {
      # Swiss, barycentric true Cartesian equatorial coordinates.
      flags <- 2L + 16L + 2048L + 4096L + 16384L + if (j2000) 32L else 0L
      pos <- swephR::swe_fixstar2_ut(id, jd, flags)
      if (pos$return < 0L || bitwAnd(pos$return, 7L) != 2L)
        .ephemeris_failure(jd, id, pos$serr)
      normalise(pos$xx[1:3])
    }
    x <- direction(",alf Lyr")
    y <- direction(",alf Vir")
    y <- normalise(y - sum(x*y)*x)
    z <- x[c(2,3,1)]*y[c(3,1,2)] - x[c(3,1,2)]*y[c(2,3,1)]
    cbind(x,y,z)
  }
  # Rows are the Galactic axes in J2000 equatorial coordinates.
  axes <- rbind(c(-.0548755604,-.8734370902,-.4838350155),
                c( .4941094279,-.4448296300, .7469822445),
                c(-.8676661490,-.1980763734, .4559837762))
  transform <- axes %*% basis(TRUE) %*% t(basis(FALSE))
  path <- function(l,b) {
    l <- l*pi/180
    b <- b*pi/180
    v <- cbind(cos(b)*cos(l),cos(b)*sin(l),sin(b)+0*l) %*% transform
    cbind(ra=(atan2(v[,2],v[,1])*180/pi)%%360,
          dec=atan2(v[,3],sqrt(v[,1]^2+v[,2]^2))*180/pi)
  }
  lines <- list()
  if (outline) {
    file <- system.file("extdata","milky-way.csv",package="skyscapeR",mustWork=TRUE)
    edges <- utils::read.csv(file)
    for (id in c("south","north")) {
      edge <- edges[edges$edge == id, ]
      lines[[id]] <- path(edge$longitude,edge$latitude)
    }
  }
  if (equator) lines$galactic.equator <- path(seq(0,360,by=.5),0)
  lines
}

.band_transmission <- function(altitude,k) {
  10^(-0.4*.stellar_extinction(pmax(0,altitude),k))
}

# Slice the two closed sky curves at their vertices. Within each longitude
# interval the crossings are straight segments and their order is constant.
# A known interior direction establishes whether a polar cap is in the band.
.draw_galactic_band <- function(south, north, polar.limit=NULL, visibility=1, k=0) {
  if (is.null(south) || is.null(north)) return(invisible(NULL))
  delta <- function(x) (x+180)%%360-180
  segments <- function(p) {
    a <- head(p[,1],-1)
    b <- a+delta(diff(p[,1]))
    y <- head(p[,2],-1)
    z <- tail(p[,2],-1)
    i <- a != b
    data.frame(a=a[i],b=b[i],y=y[i],slope=(z[i]-y[i])/(b[i]-a[i]))
  }
  edges <- rbind(segments(south),segments(north))
  crossings <- function(x) {
    d <- delta(x-edges$a)
    end <- edges$b-edges$a
    i <- which((d >= 0 & d < end) | (d < 0 & d >= end))
    heights <- edges$y[i]+d[i]*edges$slope[i]
    o <- order(heights)
    cbind(height=heights[o], slope=edges$slope[i[o]])
  }
  # The short great-circle arc between corresponding edges crosses the
  # interior. The first vertex of each contour is at Galactic longitude zero.
  unit <- function(p) {
    p <- p*pi/180
    cbind(cos(p[,2])*cos(p[,1]),cos(p[,2])*sin(p[,1]),sin(p[,2]))
  }
  mid <- unit(south[1,,drop=FALSE])+unit(north[1,,drop=FALSE])
  i <- 1L
  probe.x <- atan2(mid[i,2],mid[i,1])*180/pi
  probe.y <- atan2(mid[i,3],sqrt(sum(mid[i,1:2]^2)))*180/pi
  bottom.inside <- sum(crossings(probe.x)[,1] < probe.y) %% 2 == 0
  cuts <- sort(unique(c(0:360,south[,1]%%360,north[,1]%%360)))
  strips <- lapply(seq_len(length(cuts)-1L),function(j) {
    x <- cuts[j:(j+1L)]
    half <- diff(x)/2
    if (half < 1e-10) return(NULL)
    hits <- crossings(mean(x))
    if (bottom.inside) hits <- rbind(c(-90,0),hits)
    if (nrow(hits) %% 2L) hits <- rbind(hits,c(90,0))
    if (!nrow(hits)) return(NULL)
    left <- hits[,1]-half*hits[,2]
    right <- hits[,1]+half*hits[,2]
    lo <- seq(1L,nrow(hits),2L)
    hi <- lo+1L
    cbind(x=rep(c(x,rev(x),NA_real_),length(lo)),
          y=as.vector(rbind(left[lo],right[lo],right[hi],left[hi],NA_real_)))
  })
  band <- do.call(rbind,strips)
  band <- band[-nrow(band),,drop=FALSE]  # separators belong between subpaths
  fill <- grDevices::adjustcolor("grey55",alpha.f=.35*visibility)
  if (!is.null(polar.limit)) {
    r <- pmin(polar.limit,90-band[,2])
    a <- band[,1]*pi/180
    graphics::polypath(r*sin(a),r*cos(a),col=fill,border=NA,rule="winding")
  } else {
    box <- graphics::par("usr")
    paint <- function(colour) {
      for (turn in seq(floor(min(box[1:2])/360),floor(max(box[1:2])/360)))
        graphics::polypath(band[,1]+360*turn,band[,2],col=colour,border=NA,rule="winding")
    }
    if (k == 0) return(paint(fill))
    # Clip the existing contour so gaps and the azimuth seam survive the shading.
    breaks <- sort(unique(c(box[3:4], -90:90)))
    breaks <- breaks[breaks >= box[3] & breaks <= box[4]]
    on.exit(graphics::clip(box[1],box[2],box[3],box[4]),add=TRUE)
    for (i in seq_len(length(breaks)-1L)) {
      heights <- breaks[i+0:1]
      graphics::clip(box[1],box[2],heights[1],heights[2])
      transmission <- .band_transmission(mean(heights),k)
      paint(grDevices::adjustcolor("grey55",alpha.f=.35*visibility*transmission))
    }
  }
}
