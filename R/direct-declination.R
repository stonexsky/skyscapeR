# Integrate the spherical transformation when the four-sigma altitude
# endpoints approach a celestial pole or have coincident declinations.
.direct_declination <- function(curve, profile, horizon, bounds, prec, res,
                                output.res, refraction, atm, temp) {
  rad <- pi / 180
  lat <- horizon$metadata$georef[1] * rad
  # Resolve the input PDF even when the requested bearing grid is coarse.
  az <- seq(min(curve$x), max(curve$x), length.out =
              ceiling(diff(range(curve$x)) / min(prec, 0.025)) + 1)
  weight <- approx(curve$x, curve$y, az)$y * (az[2] - az[1])
  weight[c(1,length(az))] <- weight[c(1,length(az))] / 2
  bearing <- .profile_bearings(az, profile$az)
  height <- approx(profile$az, profile$alt, bearing)$y
  sigma <- approx(profile$az, profile$alt.unc, bearing)$y
  nodes <- qnorm(((1:256) - 0.5) / 256)
  if (refraction) {
    limits <- range(height + min(nodes)*sigma, height + max(nodes)*sigma)
    heights <- seq(floor(limits[1]*100)/100, ceiling(limits[2]*100)/100+.01, by=.01)
    correction <- vapply(heights, function(h)
      swephR::swe_refrac_extended(h, horizon$metadata$georef[3], atm, temp, 0, 0)$return - h, numeric(1))
  }
  step <- min(res, output.res, .025)
  grid <- seq(bounds[1], bounds[2], length.out=ceiling(diff(bounds)/step)+1)
  step <- grid[2]-grid[1]
  mass <- numeric(length(grid))
  for (z in nodes) {
    h <- height + z*sigma
    if (refraction) h <- h - approx(heights, correction, h)$y
    h <- h * rad
    dec <- asin(pmax(-1,pmin(1,sin(lat)*sin(h)+cos(lat)*cos(h)*cos(az*rad)))) / rad
    keep <- dec >= bounds[1] & dec <= bounds[2]
    if (!any(keep)) next
    u <- (dec[keep]-bounds[1])/step
    i <- pmin(floor(u)+1L,length(grid)-1L)
    frac <- pmin(1,pmax(0,u-i+1))
    # Deposit probability between neighbouring nodes, preserving the integral.
    contribution <- rowsum(c(weight[keep]*(1-frac),weight[keep]*frac),
                           c(i,i+1L),reorder=FALSE)
    ids <- as.integer(rownames(contribution))
    mass[ids] <- mass[ids] + contribution[,1]
  }
  # Integrate the mass into output bins. Bin edges halfway between output
  # nodes keep narrow features when the caller requests a coarse output grid.
  x <- seq(bounds[1],bounds[2],by=output.res)
  if (length(x)<2L) stop("Requested output grid has fewer than two points.")
  edges <- c(x[1],head(x,-1)+output.res/2,tail(x,1))
  cumulative <- approx(c(grid[1]-step/2,grid+step/2),c(0,cumsum(mass)),edges,rule=2)$y
  cumulative[1] <- 0
  if (tail(x,1) == tail(grid,1)) cumulative[length(cumulative)] <- sum(mass)
  data.frame(x=x,y=diff(cumulative)/(length(nodes)*diff(edges)))
}
