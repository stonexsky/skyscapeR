library(skyscapeR)
# Route 3 must reproduce individual transformations, including family/shape.
h1 <- createHor(az=c(0,90,180,270,360), alt=c(1,4,1,2,1),
                alt.unc=.5, loc=c(-11,-77,0))
h2 <- createHor(az=c(0,90,180,270,360), alt=c(2,1,3,2,2),
                alt.unc=.4, loc=c(35,-8,100))
horizons <- list(h1,h2)
centres <- c(75,115)
fixed <- lapply(centres,function(x) data.frame(x=x,w=1))
for (family in c('normal','uniform','epanechnikov','triangular',
                 'tapered_uniform','uniform_normal')) {
  extra <- switch(family, tapered_uniform=list(taper=.5),
                  uniform_normal=list(survey.sd=.5), list())
  a <- do.call(az.pdf,c(list(az=centres,unc=c(1,2),pdf=family,
                           .res=.5,verbose=FALSE),extra))
  fit <- modelTest(a,fixed,hor=horizons,nsims=3,ncores=1,.res=.5,
                   .transform.res=.1,.prec=.5,.output.res=.1,verbose=FALSE)
  expected <- spd(coordtrans(a,horizons,.res=.1,.prec=.5,.output.res=.1,
                             verbose=FALSE),.res=.5,xrange=c(-90,90))
  stopifnot(identical(fit$data$empirical,expected$data),
            identical(fit$data$empirical$y,fit$data$null.hyp$CE.upper),
            fit$metadata$coord=='dec',fit$metadata$input.coord=='az',
            fit$metadata$route=='az-to-dec',fit$metadata$global.pval==1)
}
# Uniform null: compare modelTest with explicit draws and coordtrans calls.
a <- az.pdf(centres,c(1,2),.res=.5,verbose=FALSE)
m <- data.frame(x=seq(0,359.5,by=.5),w=1)
set.seed(651)
draws <- replicate(5, c(m$x[sample.int(nrow(m),1,prob=m$w)],
                        m$x[sample.int(nrow(m),1,prob=m$w)]),simplify=FALSE)
manual <- t(vapply(draws,function(x) {
  sim <- az.pdf(x,c(1,2),.res=.5,verbose=FALSE)
  spd(coordtrans(sim,horizons,.prec=.5,.output.res=.1,verbose=FALSE),
      xrange=c(-90,90),.res=.5,normalise=TRUE)$data$y
},numeric(361)))
observed <- spd(coordtrans(a,horizons,.prec=.5,.output.res=.1,verbose=FALSE),
                xrange=c(-90,90),.res=.5,normalise=TRUE)$data
reference <- skyscapeR:::.simulation_inference(observed,5,list(1:5),character(),manual,.95,2)
set.seed(651)
fit <- modelTest(a,m,hor=horizons,nsims=5,ncores=1,.res=.5,
                 .prec=.5,.output.res=.1,normalise=TRUE,verbose=FALSE)
stopifnot(identical(fit$data$null.hyp,reference$envelope),
          identical(fit$metadata$global.pval,reference$global),
          identical(fit$metadata$local.pval,reference$local))
# A shared profile is equivalent to a repeated ordered profile list.
set.seed(9)
single <- modelTest(a,m,hor=h1,nsims=3,ncores=1,.res=1,verbose=FALSE)
set.seed(9)
repeated <- modelTest(a,m,hor=list(h1,h1),nsims=3,ncores=1,.res=1,verbose=FALSE)
stopifnot(identical(single$data,repeated$data))
# Invalid route combinations and unmatched horizon lists fail explicitly.
d <- dec.pdf(c(-20,10),1,verbose=FALSE)
stopifnot(inherits(try(modelTest(d,m,hor=h1,verbose=FALSE),silent=TRUE),'try-error'),
          inherits(try(modelTest(a,m,hor=list(h1),verbose=FALSE),silent=TRUE),'try-error'),
          inherits(try(modelTest(a,data.frame(x=-10,w=1),hor=h1,verbose=FALSE),silent=TRUE),'try-error'))

# Preserve asymmetric triangular proportions when replacing azimuth widths.
a <- az.pdf(c(75,115),pdf='triangular',lower=c(73,112),upper=c(76,119),
            .res=.5,verbose=FALSE)
fixed <- list(data.frame(x=80,w=1),data.frame(x=120,w=1))
fit <- modelTest(a,fixed,unc=c(3,4),hor=horizons,nsims=3,ncores=1,
                 .res=.5,verbose=FALSE)
sim <- skyscapeR:::.measurement_pdf(c(80,120),'az',c(3,4),'triangular',
                                    a$metadata$shape,c('a','b'),.5,1e-5,FALSE)
y <- spd(coordtrans(sim,horizons,verbose=FALSE),xrange=c(-90,90),.res=.5)$data$y
stopifnot(isTRUE(all.equal((fit$data$empirical$y+3*y)/4,fit$data$null.hyp$CE.mean)))
