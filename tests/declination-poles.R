library(skyscapeR)
skyscapeR.vars(refraction=FALSE)
# Positive density at 2.02 degrees has distinct full-precision endpoints.
az <- c(0,3,7,90,180,270,357,360)
hor <- createHor(az=az,alt=c(43,46,15,2,0,3,15,43),loc=c(45,0,0),
                 alt.unc=.5,verbose=FALSE)
p <- az.pdf(3,1,pdf="epanechnikov",verbose=FALSE)
quantiles <- function(h,p,step) {
  d <- coordtrans(p,h,.prec=step,verbose=FALSE)$data[[1]]
  mass <- c(0,cumsum(diff(d$x)*(head(d$y,-1)+tail(d$y,-1))/2))
  stopifnot(all(is.finite(d$y)), all(d$y>=0), abs(tail(mass,1)-1)<1e-9)
  approx(mass,d$x,c(.025,.5,.975),ties="ordered")$y
}
stopifnot(max(abs(quantiles(hor,p,1)-c(83.933,87.729,88.403)))<.08)
exact <- createHor(az=c(0,360),alt=c(50,50),loc=c(50,0,0),alt.unc=.5,verbose=FALSE)
q <- quantiles(exact,az.pdf(0,2,verbose=FALSE),.1)
stopifnot(max(abs(q[c(1,3)]-c(87.071,89.819)))<.05)
