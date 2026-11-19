library(skyscapeR)
g <- skyscapeR:::.milky_way_grid()
stopifnot(is.raw(g), identical(dim(g),c(3600L,1201L)),identical(sort(unique(as.integer(g))),0:4))
# Exercise the actual coordinate lookup over every stored cell, in bounded chunks.
for (j in seq_len(ncol(g))) {
 l <- seq(-180,179.9,.1)*pi/180;b <- (-60+(j-1)*.1)*pi/180
 got <- skyscapeR:::.galactic_classes(cos(b)*cos(l),cos(b)*sin(l),rep(sin(b),length(l)),diag(3))
 stopifnot(identical(got,as.integer(g[,j])))
}
l <- c(-180,180,540)*pi/180;b <- .1
stopifnot(length(unique(skyscapeR:::.galactic_classes(cos(b)*cos(l),cos(b)*sin(l),rep(sin(b),3),diag(3))))==1L)
# The visually selected Emu neighbour remains in the dark class on this grid.
l <- -78.25*pi/180;b <- -.75*pi/180
stopifnot(skyscapeR:::.galactic_classes(cos(b)*cos(l),cos(b)*sin(l),sin(b),diag(3))==2L)
# Verify calibration remains finite and retains its area-weighted anchor.
f <- expand.grid(l=seq(284.95,289.05,.1),b=seq(-4.05,.05,.1))
keep <- skyscapeR:::.cluster_distance(287,-2,f$l,f$b)<=2
l <- f$l[keep]*pi/180;b <- f$b[keep]*pi/180
id <- skyscapeR:::.galactic_classes(cos(b)*cos(l),cos(b)*sin(l),sin(b),diag(3))
v <- skyscapeR:::.mw_radiance();stopifnot(all(is.finite(v)),all(v>=0))
mag <- -2.5*log10(weighted.mean(v[id+1L],cos(b))/3600^2)
stopifnot(abs(mag-20.7)<1e-10)
# Transmission decreases with increasing extinction, and increases with altitude.
h <- c(0,2,5,10,20,45,90)
a <- skyscapeR:::.band_transmission(h,.2);b <- skyscapeR:::.band_transmission(h,.3)
stopifnot(all(diff(a)>0),all(a>b),all(a>0 & a<1))
cat('Passed: all 4,323,600 cells, longitude seam, selected patch, Carina anchor and extinction ordering.\n')
