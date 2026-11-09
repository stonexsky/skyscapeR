library(skyscapeR)
# A non-default source and sum grid must be reproduced for a null draw
# identical to the observed bearings, even when the shared skyline is used.
az <- c(30,120)
p <- az.pdf(az,c(22,32),.res=30,.cutoff=1e-4,verbose=FALSE)
h <- createHor(az=c(0,90,180,270,360),alt=c(2,8,2,5,2),
               loc=c(46,0,0),alt.unc=1,verbose=FALSE)
settings <- list(family="normal",unc=c(22,32),source.res=30,source.cutoff=1e-4,
                 coord="dec",quick=TRUE,xrange=c(-90,90),resolution=2,normalise=FALSE,
                 transform=list(hor=h,refraction=FALSE,atm=1013.25,temp=15,.res=2))
skyscapeR.vars(refraction=FALSE)
empirical <- coordtrans(spd(p,.res=2),h,xrange=c(-90,90),.res=2,verbose=FALSE)
stopifnot(identical(empirical$data$y,skyscapeR:::.random_curve(az,settings)))
