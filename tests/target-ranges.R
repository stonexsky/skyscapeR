library(skyscapeR)
# Named local limits use the observer and every epoch even for short ranges.
site<-c(50,0,0)
x<-skyObject("June Solstice",2000:2001,loc=site)
stopifnot(isTRUE(all.equal(unname(x$decs[,1]),range(jS(2000:2001,site,verbose=FALSE)))))
stopifnot(isTRUE(all.equal(zenith(20,year=2000),20)),is.null(antizenith(75,year=2000)))
