library(skyscapeR)
skyscapeR.vars(refraction=FALSE)
hor<-createHor(0:359,rep(0,360),c(51.18,-1.83,100),name="Wiltshire flat-horizon benchmark")
events<-do.call(rbind,lapply(c("Sirius","Aldebaran"),function(n){
 e<-star.phases(n,2026,hor,refraction=FALSE,k=.2)$metadata$events
 print(e);cbind(star=n,e)
}))
write.csv(events,"events.csv",row.names=FALSE)
