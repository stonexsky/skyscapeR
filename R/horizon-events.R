.panorama_sector <- function(hor) {
  if (!identical(hor$metadata$source$provider, "Panorama")) return(NULL)
  coverage <- hor$metadata$source$coverage
  if (length(coverage)!=2L || any(!is.finite(coverage)) || diff(coverage)<=0)
    stop("Panorama has invalid azimuth coverage.")
  if (diff(coverage)<360-1e-8) coverage else NULL
}

.require_full_skyline <- function(hor) {
  if (!is.null(.panorama_sector(hor)))
    stop("A full 360-degree skyline is required; the photographed sector has incomplete coverage.", call.=FALSE)
}

# Terrain vertices and interior extrema supplement the time grid. This helper
# accepts a smooth astronomical position function for independent analytic tests.
.horizon_profile <- function(az,altitude) {
  if(length(az)<2L || length(az)!=length(altitude) ||
     any(!is.finite(c(az,altitude))) || any(diff(az)<=0) || min(az)>0 || max(az)<360)
    stop("Horizon must contain finite, increasing azimuths covering 0 to 360 degrees.")
  horizon <- stats::approxfun(az,altitude,rule=1)
  if (abs(horizon(0)-horizon(360))>1e-7)
    stop('Horizon altitudes at 0 and 360 degrees must agree.')
  # Remove redundant points on straight terrain segments, retaining corners.
  slopes <- diff(altitude)/diff(az)
  keep <- c(TRUE,abs(diff(slopes))>1e-9,TRUE)
  corners <- az[keep];corners <- corners[corners>=0 & corners<=360]
  attr(horizon,"range") <- range(altitude[az>=0 & az<=360],horizon(c(0,360)))
  attr(horizon,"corners") <- corners
  horizon
}

.horizon_crossings <- function(position,jd0,az,altitude,.step=1/1440,
                               .duration=2,.all=FALSE,horizon=.horizon_profile(az,altitude)) {
  corners <- attr(horizon,"corners")
  extended <- sort(unique(c(corners-360,corners,corners+360)))
  wrap <- function(x) (x+180)%%360-180
  residual <- function(t) {
    p <- position(t)
    p[2]-horizon(p[1]%%360)
  }
  tol <- .01/86400
  grid <- sort(unique(c(seq(jd0,jd0+.duration,by=.step/4),jd0+.duration)))
  pos <- vapply(grid,position,numeric(2))
  if(any(!is.finite(pos))) stop('Nonfinite position during horizon crossing search.')
  # Include azimuth turning points detected on the subgrid. Unwrapping locally
  # avoids treating north (360/0) as a discontinuity in a continuous horizon.
  da <- wrap(diff(pos[1,]))
  turns <- which(head(da,-1)*tail(da,-1)<0)+1L
  extra <- numeric()
  for(i in turns) {
    ref <- pos[1,i]
    f <- function(t) ref+wrap(position(t)[1]-ref)
    origin <- grid[i-1]
    ex <- stats::optimize(function(seconds) f(origin+seconds/86400),
                         c(0,(grid[i+1]-origin)*86400),
                         maximum=da[i-1]>0,tol=.001)
    extremum <- origin+if(da[i-1]>0) ex$maximum/86400 else ex$minimum/86400
    extra <- c(extra,extremum)
  }
  grid <- sort(unique(c(grid,extra)))
  result <- c(NA_real_,NA_real_)
  crossings <- numeric()
  # Recursively split rapid azimuth motion before locating terrain vertices.
  segment <- function(a,b,depth=0L) {
    pa <- position(a);pb <- position(b)
    delta <- wrap(pb[1]-pa[1])
    middle <- wrap(position((a+b)/2)[1]-pa[1])
    reverses <- middle<min(0,delta)-1e-8 || middle>max(0,delta)+1e-8
    if((abs(delta)>10 || reverses) && b-a>tol && depth<24L) {
      m <- (a+b)/2
      return(sort(unique(c(segment(a,m,depth+1L),segment(m,b,depth+1L)))))
    }
    knots <- numeric()
    if(abs(delta)>1e-12) {
      # Binary lookup avoids scanning the entire DEM profile for each bracket.
      ends <- sort(c(pa[1],pa[1]+delta))
      indices <- findInterval(ends,extended)
      selected <- numeric()
      if(indices[2]>indices[1]) {
        targets <- extended[seq.int(indices[1]+1L,indices[2])]
        selected <- targets[targets<ends[2]]-pa[1]
      }
      for(target in selected) {
        f <- function(t) wrap(position(t)[1]-pa[1])-target
        root <- stats::uniroot(f,c(a,b),tol=tol)$root
        if(abs(f(root))<1e-3) knots <- c(knots,root)
      }
    }
    sort(unique(c(a,knots,b)))
  }
  for(i in seq_len(length(grid)-1L)) {
    knots <- segment(grid[i],grid[i+1L])
    for(j in seq_len(length(knots)-1L)) {
      a <- knots[j];b <- knots[j+1L]
      y <- c(residual(a),residual(b))
      if(any(!is.finite(y))) stop('Nonfinite altitude during horizon crossing search.')
      times <- c(a,b)
      # Between corners, a smooth path can enter and leave without an endpoint
      # sign change. Search its interior extremum, including grazing passages.
      if(y[1]*y[2]>0) {
        direction <- if(y[1]<0) 1 else -1
        # On each smooth short segment, test endpoint slopes before the
        # unimodal extremum search. Monotone segments need no optimization.
        probe <- min((b-a)/4,.001/86400)
        left <- direction*(residual(a+probe)-y[1])
        right <- direction*(y[2]-residual(b-probe))
        if(left>0 && right<0) {
          f <- function(t) direction*residual(t)
          ex <- stats::optimize(function(seconds) f(a+seconds/86400),
                                c(0,(b-a)*86400),maximum=TRUE,tol=.001)
          if(ex$objective>0) times <- c(a,a+ex$maximum/86400,b)
        }
      }
      yy <- vapply(times,residual,numeric(1))
      for(k in seq_len(length(times)-1L)) {
        direction <- if(yy[k]<=0 && yy[k+1]>0) 1L else
          if(yy[k]>=0 && yy[k+1]<0) 2L else 0L
        if(direction && (.all || is.na(result[direction]))) {
          root <- stats::uniroot(residual,times[k+0:1],tol=tol)$root
          # A touch at an exactly sampled zero is not a rise or set.
          if(yy[k]==0 || yy[k+1]==0) {
            sides <- c(residual(root-tol),residual(root+tol))
            if(prod(sides)>=0) next
          }
          if(.all) crossings <- c(crossings,root) else result[direction] <- root
        }
      }
    }
    if(!.all && all(is.finite(result))) break
  }
  if(.all) sort(unique(crossings)) else result
}

# Find the first crossings in a bounded two-day search after local midnight.
.horizon_events <- function(obj,jd0,hor,refraction,atm,temp,.step=1/1440) {
  loc <- as.numeric(hor$metadata$georef)
  az <- hor$data$az;altitude <- hor$data$alt
  swephR::swe_set_topo(loc[2],loc[1],loc[3])
  body <- checkbody(obj)
  position <- function(jd) {
    eq <- .file_position(jd,body,2048L+32768L)
    if (eq$return < 0L) stop(eq$serr)
    horizontal <- swephR::swe_azalt(jd,1L,c(loc[2],loc[1],loc[3]),
                                   atm,temp,eq$xx[1:3])$xaz
    c((horizontal[1]+180)%%360,horizontal[if(refraction) 3L else 2L])
  }
  .horizon_crossings(position,jd0,az,altitude,.step)
}
