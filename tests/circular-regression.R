library(skyscapeR)
p <- az.pdf(0,2,.res=.005,verbose=FALSE)$data[[1]]
stopifnot(abs(p$y[p$x==0]-stats::dnorm(0,0,2))<1e-12,
          sum(abs(p$x-359.995)<1e-9)==1)
d <- dec.pdf(c(-10,12),1,verbose=FALSE)
stopifnot(identical(subset(d,c(2,2,1))$metadata$dec,c(12,12,-10)))
set.seed(90213)
f <- randomTest(az.pdf(c(359,2),3,.res=.5,verbose=FALSE),
                nsims=15,conf=.8,.res=.5,ncores=1,verbose=FALSE)
r <- f$metadata$local.pval
stopifnot(nrow(r)==1L,r$start==349.5,r$end==10,
          r$p.value==3/16,f$metadata$global.pval==7/16)
