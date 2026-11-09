library(skyscapeR)

# Known weights and masses for the field-guide example, including clipped
# declinations and a modal (rather than mean-centred) asymmetric reconstruction.
density <- function(...) az.pdf(..., .res=0.01, .cutoff=0, verbose=FALSE)$data[[1]]
at <- function(curve, x) approx(curve$x, curve$y, x)$y
near <- function(x, y, tol=1e-10) stopifnot(all(abs(x-y) < tol))
ep <- density(90, 2, "epanechnikov")
soft <- density(90, 2, "tapered_uniform", taper=0.5)
conv <- density(90, 2, "uniform_normal", survey.sd=0.2)
tri <- density(90, pdf="triangular", lower=89, upper=93)
near(at(ep, c(88,89,90,91,92)), c(0,9/32,3/8,9/32,0))
near(at(soft, c(87.5,87.75,88,90,92,92.25,92.5)), c(0,1/9,2/9,2/9,2/9,1/9,0))
near(at(tri, c(89,89.5,90,91,93)), c(0,.25,.5,1/3,0))
near(at(conv, c(87.8,88,88.2,90)), c(pnorm(-1),.5,pnorm(1),1)/4)
for (curve in list(ep,soft,conv,tri)) near(sum(curve$y)*.01,1,1e-5)
near(at(dec.pdf(90,2,"epanechnikov",.cutoff=0,verbose=FALSE)$data[[1]], 90), 3/8)
stopifnot(identical(density(0,2,"uniform")$y, density(0,2,"uniform_normal",survey.sd=0)$y),
          identical(density(0,2,"uniform")$y, density(0,2,"tapered_uniform",taper=0)$y))
north <- density(0, pdf="triangular", lower=359, upper=3)
near(at(north,c(359.5,0,1,2)), c(.25,.5,1/3,1/6))
one <- density(0, pdf="triangular", lower=0, upper=200)
near(at(one,c(0,190,210,359)), c(.01,.0005,0,0))
# A very broad wrapped survey error remains normalised at north and the antipode.
wide <- density(0, 20, "uniform_normal", survey.sd=400)
near(sum(wide$y)*.01,1)
for (x in c(0,180,359)) {
  ref <- sum(vapply(-12:12,function(k)
    (pnorm((x+360*k+20)/400)-pnorm((x+360*k-20)/400))/40,numeric(1)))
  near(at(wide,x), ref)
}
fails <- function(expr) inherits(try(expr,silent=TRUE), "try-error")
stopifnot(fails(az.pdf(90,2,"tapered_uniform")),
          fails(az.pdf(90,2,"normal",taper=1)), fails(az.pdf(90,2,"uniform_normal",survey.sd=-1)),
          fails(az.pdf(90,pdf="triangular",lower=93,upper=95)),
          fails(az.pdf(90,2,"triangular",lower=89,upper=93)),
          fails(dec.pdf(20,pdf="triangular",lower=19,upper=18)))
# When each candidate is its observed bearing, both PDFs must match exactly,
# even after reordering observations with distinct shape parameters.
for (coord in c("az","dec")) {
  constructor <- if (coord=="az") az.pdf else dec.pdf
  assumed <- constructor(c(10,30),2,"uniform_normal",verbose=FALSE)
  stated <- constructor(c(10,30),2,"uniform_normal",survey.sd=.5,verbose=FALSE)
  stopifnot(identical(assumed, stated), identical(assumed$metadata$shape$survey.sd, c(.5,.5)),
            fails(constructor(10,2,"normal",survey.sd=.5)))
  for (opts in list(list(unc=c(2,3),pdf="epanechnikov"),
                    list(unc=c(2,3),pdf="uniform_normal"),
                    list(unc=c(2,3),pdf="tapered_uniform",taper=c(.5,1)),
                    list(unc=c(2,3),pdf="uniform_normal",survey.sd=c(.2,.4)),
                    list(pdf="triangular",lower=c(8,27),upper=c(15,31)))) {
    obs <- do.call(constructor,c(list(c(10,30),.res=.1,verbose=FALSE),opts))
    obs <- subset(obs,c(2,1,2))
    model <- lapply(c(30,10,30),function(a) data.frame(x=a,w=1))
    fit <- modelTest(obs,model,nsims=3,ncores=1,.res=.1,verbose=FALSE)
    stopifnot(identical(fit$data$empirical$y,fit$data$null.hyp$CE.upper),
              fit$metadata$global.pval==1)
    if (coord=="az") {
      set.seed(719)
      null <- randomTest(obs,nsims=3,ncores=1,.res=1,verbose=FALSE)
      stopifnot(identical(null$metadata$shape,obs$metadata$shape),
                is.finite(null$metadata$global.pval))
    }
  }
}
