# Naked-eye point-source threshold from the two Hecht regimes (Schaefer 1986).
.hecht_limit <- function(brightness) {
  dark <- brightness < 10^3.17
  coefficient <- ifelse(dark, 10^-9.8, 10^-8.35)
  scale <- ifelse(dark, 10^-1.9, 10^-5.9)
  -16.57 - 2.5 * log10(coefficient * (1 + sqrt(scale * brightness))^2)
}

# Combine a terrain-clear fixed-star track with the Sun on actual civil days.
.schaefer_calendar <- function(year, star, loc, track, k, limit, atm, temp, refraction) {
  first <- swephR::swe_julday(year, 1L, 1L, 0, 1L)
  end <- swephR::swe_julday(year + 1L, 1L, 1L, 0, 1L)
  days <- seq(first - 1, end)
  result <- data.frame(local.jd=days, morning=FALSE, evening=FALSE)
  open <- which(track$open)
  if (!length(open)) return(result)
  # Interior samples resolve even a short clearance between terrain crossings.
  widths <- diff(track$time)[open]
  times <- sort(c(track$time[open] + widths/4, track$time[open] + widths*3/4))
  times <- times[times < track$time[1] + 360/360.985647]
  rad <- pi/180
  eq <- as.numeric(star$coord[2:3])
  pos <- vapply(times, function(jd) ff(jd,loc,atm,temp,eq), numeric(3))
  alt <- pos[if (refraction) 3L else 2L,]
  az <- (pos[1,] + 180) * rad
  sid <- vapply(times, swephR::swe_sidtime, numeric(1)) * 15
  west <- sin((sid + loc[2] - eq[1]) * rad) > 0
  # At a geographic pole the track has no rising or setting branch.
  if (abs(loc[1]) == 90 || abs(eq[2]) == 90) return(result)
  z <- (90 - alt) * rad
  extinguished <- star$app.mag + .stellar_extinction(pmax(alt,0), k)
  dark <- limit + k >= .hecht_limit(10^3.17)
  C <- if (dark) 10^-9.8 else 10^-8.35
  K <- if (dark) 10^-1.9 else 10^-5.9
  background <- (sqrt(10^(-.4*(16.57+limit+k))/C)-1)^2/K
  midnight <- days - loc[2]/360
  nodes <- seq(midnight[1], tail(midnight,1)+1, by=.5)
  sun <- matrix(.file_position(nodes,0L,2048L)$xx,nrow=length(nodes))
  ra <- c(sun[1,1], sun[1,1]+cumsum((diff(sun[,1])+180)%%360-180))
  dec <- sun[,2]
  lst <- vapply(midnight,swephR::swe_sidtime,numeric(1))*15
  lat <- loc[1]*rad
  for (i in seq_along(days)) {
    t <- ((sid-lst[i])%%360)/360.985647
    twice <- which(t < 1-360/360.985647)
    ix <- c(seq_along(t),twice)
    t <- c(t,t[twice]+360/360.985647)
    ix0 <- 2L*i-1L + as.integer(t>=.5)
    mix <- (2*t)%%1
    hs <- (sid[ix]+loc[2]-ra[ix0]-mix*(ra[ix0+1L]-ra[ix0]))*rad
    ds <- (dec[ix0]+mix*(dec[ix0+1L]-dec[ix0]))*rad
    h <- asin(pmin(1,pmax(-1,sin(lat)*sin(ds)+cos(lat)*cos(ds)*cos(hs))))
    solar.az <- atan2(-cos(ds)*sin(hs), sin(ds)*cos(lat)-cos(ds)*cos(hs)*sin(lat))
    delta <- abs((solar.az-az[ix]+pi)%%(2*pi)-pi)
    log.b <- 4.75-delta*z[ix]/3 + h*(12+8.21*z[ix])+2.86*z[ix]
    b <- background + k/.2*pmax(0,10^log.b-118)
    seen <- h<0 & .hecht_limit(b)>=extinguished[ix]
    result$morning[i] <- any(seen & !west[ix] & t<.5)
    result$evening[i] <- any(seen & west[ix] & t>=.5)
  }
  result
}
