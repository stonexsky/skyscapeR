library(skyscapeR)
f<-skyscapeR:::.sky_luminance
P<-c(0,10,30,50,70,90)
# Koomen Table II, Maryland; none of these measurements enter the fit.
obs<-rbind(c(13.5,13,7.7,3.7,2.5,2),c(2,2,4,3.2,2,2),
           c(.6,.65,.3,.1,.075,.06),c(.09,.08,.07,.07,.05,.06))
for(i in 1:4){
 l<-list(sun=c(0,if(i<=2)-3 else -6),k=.2)
 az<-if(i%%2)0 else 180
 # compare the zenith only once per solar altitude
 for(model in c("old","new")){
   l$calibrated.twilight<-model=="new"
   p<-f(l,az,P)/(10.76391*pi*1e5)
   ii<-if(i%%2)1:6 else 1:5
   df<-data.frame(h=l$sun[2],az=az,alt=P[ii],observed=obs[i,ii],pred=p[ii],model=model)
   if(i==1 && model=="old") data<-df else data<-rbind(data,df)
 }
}
write.csv(data,"holdout.csv",row.names=FALSE)
print(tapply(abs(2.5*log10(data$pred/data$observed)),data$model,mean))
print(subset(data,h==-3 & alt%in%c(10,90)))
