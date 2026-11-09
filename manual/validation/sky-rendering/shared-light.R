library(skyscapeR)
skyscapeR.vars(refraction=FALSE)
create <- function(id,month,lat.offset=20){
  jd <- swephR::swe_julday(2026,month,1,0,1L)
  eq <- swephR::swe_fixstar2_ut(paste0(",",id),jd,2L+2048L)$xx[1:2]
  lon <- (eq[1]-15*swephR::swe_sidtime(jd)+180)%%360-180
  h <- createHor(0:359,rep(0,360),c(eq[2]+lat.offset,lon,0),name=id)
  s <- sketchSky(h,jd=jd,style="natural",max.mag=5.5,
                 display=list(moon=FALSE,planets=TRUE))
  list(h=h,s=s)
}
integral <- function(field) sum(field$direct*field$area)
f <- skyscapeR:::.shared_radiance
a <- create("alf CMa",1)
b <- create("M 45",12)
cl <- subset(b$s$clusters,identifier=="M 45")
d <- skyscapeR:::.cluster_distance(cl$az,cl$alt,b$s$stars[,1],b$s$stars[,2])
ii <- which(d <= cl$radius)
sf <- a$s$light.budget$stars["alf CMa"]
cf <- sum(b$s$light.budget$stars[ii])+cl$residual.flux
print(c(Sirius=sf,cluster=cf,Sirius.over.cluster=unname(sf/cf)))
single <- a$s
single$light.budget$stars[] <- 0
single$light.budget$stars["alf CMa"] <- sf
single$light.budget$bodies <- numeric()
single$clusters <- NULL
p <- single$stars["alf CMa",]
for(width in c(90,45,5)){
  field <- f(single,c(p[1]+c(-.5,.5)*width,p[2]+c(-.25,.25)*width),240,120)
  print(c(width=width,recovery=integral(field)/sf))
  stopifnot(abs(integral(field)/sf-1)<1e-9)
}
target <- b$s
target$light.budget$stars[-ii] <- 0
target$light.budget$bodies <- numeric()
target$clusters <- cl
field <- f(target,c(cl$az+c(-7,7),cl$alt+c(-2,2)),240,120)
stopifnot(abs(integral(field)/cf-1)<1e-9,abs(sf/cf-16.7)<.7)
jd <- swephR::swe_julday(2026,3,1,0,1L)
D <- pi/180
v <- c(cos(-2*D)*cos(287*D),cos(-2*D)*sin(287*D),sin(-2*D))%*%skyscapeR:::.galactic_frame(jd)
h <- createHor(0:359,rep(0,360),c(asin(v[3])/D,(atan2(v[2],v[1])/D-15*swephR::swe_sidtime(jd)+180)%%360-180,0),name="Carina")
s <- sketchSky(h,jd=jd,display=list(moon=FALSE))
field <- f(s,c(0,360,88,90),720,100)
car <- sum(field$milky.way*field$area)/(720*sum(field$area))
print(c(Carina=car,Sirius.over.Carina=unname(sf/car),cluster.over.Carina=cf/car))
stopifnot(abs(sf/car-56)<2,abs(cf/car-3.4)<.15)
# Actual Moon phases and visible planets; seek two different illuminated phases.
found <- numeric()
checked <- character()
for(day in 0:29){
  jd <- swephR::swe_julday(2026,2,1,0,1L)+day
  eq <- swephR::swe_calc_ut(jd,1L,2048L+2L)$xx[1:2]
  lon <- (eq[1]-15*swephR::swe_sidtime(jd)+180)%%360-180
  h <- createHor(0:359,rep(0,360),c(eq[2]+30,lon,0),name="Lunar test")
  obj <- sketchSky(h,jd=jd,display=list(stars=FALSE,planets=TRUE))
  if(obj$bodies$moon[2]<10) next
  for(id in names(obj$light.budget$bodies)){
    p<-obj$bodies[[id]]
    if(p[2]<5) next
    checked <- c(checked,id)
    one<-obj;one$light.budget$bodies<-obj$light.budget$bodies[id]
    field<-f(one,c(p[1]+c(-4,4),p[2]+c(-2,2)),96,48)
    stopifnot(abs(integral(field)/one$light.budget$bodies[[id]]-1)<1e-8)
  }
  found<-c(found,obj$lunar$fraction)
}
print(range(found))
print(table(checked))
stopifnot(all(c("moon","mercury","venus","mars","jupiter","saturn") %in% checked))
stopifnot(diff(range(found))>.55)
# The public diagnostic reproduces the linear integrator; test the north seam
# and a sufficiently tight zoom to exercise resolved rendering.
seam <- single
seam$stars["alf CMa",1] <- 359.95
for(span in c(4,2)){
  r <- skyRadiance(seam,c(360-span/2,360+span/2),c(69,71),c(300,200))
  stopifnot(abs(integral(r)/sf-1)<1e-9)
}
off <- sketchSky(h,jd=jd,display=list(moon=FALSE,stars=FALSE))
stopifnot(is.null(off$light$moon),is.null(off$bodies$moon),length(off$light.budget$bodies)==0)
zero <- skyRadiance(off,c(355,365),c(68,70),c(30,10))
stopifnot(all(zero$direct==0),all(zero$milky.way==0))
cat("All shared-light checks passed.\n")
