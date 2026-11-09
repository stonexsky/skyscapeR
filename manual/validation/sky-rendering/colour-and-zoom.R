library(skyscapeR)
hor <- createHor(0:359,rep(0,360),c(-30,149,0),name="Southern sky")
skyscapeR.vars(refraction=FALSE)
s <- sketchSky(hor,jd=2461268.916666667,display=list(planets=TRUE))
p <- s$bodies$venus
m <- s$bodies$moon
print(skyscapeR:::.sketch_pheno(s$jd,"venus",s$loc)[c(2,4)]*c(1,3600))
print(s$diameters["moon"]*60)
winter <- sketchSky(hor,time="2026/01/01 12:30",timezone="UTC",display=list(moon=FALSE,planets=TRUE))

old <- packageVersion("skyscapeR")=="2.0.9001"
S <- s
S$light.budget$stars[] <- 0
S$light.budget$bodies <- S$light.budget$bodies["venus"]
S$clusters <- NULL
for(w in c(5,.3,.04)){
 x <- p[1]+c(-.5,.5)*w
 y <- p[2]+c(-.25,.25)*w
 f <- skyRadiance(S,x,y,c(600,300))
 field <- f$radiance
 K <- .04*max(1,skyscapeR:::.sky_luminance(S$light,0,90)/60)^.75
 tone <- function(L) pmin(asinh(L/K)/asinh(2500),1)
 excess <- tone(field)-tone(field-f$direct)
 hit <- which(excess>max(excess)*.5,arr.ind=TRUE)
 radial <- skyscapeR:::.cluster_distance(p[1],p[2],f$azimuth[hit[,2]],f$altitude[hit[,1]])
 # Full transverse illuminated width perpendicular to local Sun direction.
 angle <- skyscapeR:::.lunar_light_angle(p,S$bodies$sun)
 across <- -(f$azimuth[hit[,2]]-p[1])*cos(p[2]*pi/180)*sin(angle)+(f$altitude[hit[,1]]-p[2])*cos(angle)
 width <- diff(range(across))*3600
 cat(sprintf("Venus field %.3f deg, half-contrast width %.2f arcsec, flux %.6f\n",w,width,sum(f$direct*f$area)))
 if(!old && w==.04) stopifnot(abs(width-s$diameters["venus"]*3600)<.7)
}
a <- winter$stars["alf CMa",1:2]
winter$light.budget$stars[setdiff(names(winter$light.budget$stars),"alf CMa")] <- 0
winter$light.budget$bodies <- numeric()
winter$clusters <- NULL
for(w in c(5,.3,.04)){
 f <- skyRadiance(winter,a[1]+c(-.5,.5)*w,a[2]+c(-.25,.25)*w,c(600,300))
 j<-which(f$direct>0,arr.ind=TRUE)
 cat(sprintf("Sirius field %.3f deg, footprint %d x %d pixels, flux %.6f\n",w,diff(range(j[,2]))+1,diff(range(j[,1]))+1,sum(f$area*f$direct)))
 if(!old) stopifnot(nrow(j)<=4)
}
if(!old){
 field <- skyscapeR:::.shared_radiance(s,c(p[1]+c(-3,3),p[2]+c(-1.5,1.5)),240,120,colour=TRUE)
 mono <- skyscapeR:::.shared_radiance(s,c(p[1]+c(-3,3),p[2]+c(-1.5,1.5)),240,120)
 stopifnot(identical(mono$radiance,field$radiance),identical(mono$direct,field$direct))
 testColour <- function(s,box){
 f <- skyscapeR:::.shared_radiance(s,box,240,120,colour=TRUE)
 rgb <- grDevices::col2rgb(as.vector(skyscapeR:::.shared_pixels(f,s)))/255
 linear <- ifelse(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)^2.4)
 actual <- drop(c(.2126,.7152,.0722)%*%linear)
 K <- .04*max(1,skyscapeR:::.sky_luminance(s$light,0,90)/60)^.75
 grey <- pmin(1,asinh(as.vector(f$radiance)/K)/asinh(2500))
 recovered <- ifelse(actual<=.0031308,actual*12.92,1.055*actual^(1/2.4)-.055)
 err <- max(abs(recovered-grey))*255
 cat("display luminance discrepancy (8-bit levels):",err,"\n")
 stopifnot(err<.52)
 }
 testColour(s,c(185,285,0,50))
 testColour(sketchSky(hor,jd=s$jd-.08,display=list(planets=TRUE)),c(160,300,0,45))
}
