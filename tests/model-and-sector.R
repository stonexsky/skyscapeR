library(skyscapeR)
for (coord in c("az","dec")) {
  constructor <- if (coord=="az") az.pdf else dec.pdf
  input <- constructor(c(10,30),c(3,5),.res=.5,.cutoff=.002,verbose=FALSE)
  limits <- if (coord=="az") c(0,360) else c(-90,90)
  # A deterministic site-specific null reproduces the measurements exactly.
  model <- list(data.frame(x=10,w=1),data.frame(x=30,w=1))
  fit <- modelTest(input,model,nsims=3,ncores=1,.res=1,verbose=FALSE)
  stopifnot(identical(fit$data$empirical$y,fit$data$null.hyp$CE.upper),
            fit$metadata$global.pval==1)
}
if (requireNamespace("png",quietly=TRUE)) {
  f <- tempfile(fileext=".png")
  img <- array(0,c(20,60,4))
  img[11:20,,4] <- 1
  png::writePNG(img,f)
  sector <- importPanorama(f,c(50,0,0),az.lim=c(330,30),alt.lim=c(-10,10))
  good <- az.pdf(359,3,pdf="uniform",.res=.5,verbose=FALSE)
  bad <- az.pdf(29,3,pdf="uniform",.res=.5,verbose=FALSE)
  stopifnot(is.na(hor2alt(sector,180)),
            inherits(coordtrans(good,sector,verbose=FALSE),"skyscapeR.pdf"),
            inherits(try(coordtrans(bad,sector,verbose=FALSE),silent=TRUE),"try-error"),
            inherits(try(riseset(loc=sector,date="2000/01/02"),silent=TRUE),"try-error"))
  unlink(f)
}
