library(skyscapeR)
skyscapeR.vars(refraction=FALSE)
loc<-c(51.18,-1.83,100); geo<-loc[c(2,1,3)]
hor<-createHor(0:359,rep(0,360),loc,name="Wiltshire benchmark")
ev<-read.csv("events.csv")
B0<-60*1e-5/pi
increment<-function(B) ifelse(B<.0708,(6.505e-4*B^.25-8.461e-4*sqrt(B))^2,
                                    (1.772e-4*B^.25+7.167e-5*sqrt(B))^2)
threshold<-function(nL) 6-2.5*log10(increment(nL*1e-5/pi)/increment(B0))
f<-skyscapeR:::.sky_luminance
toDate<-function(j) as.Date(j-2440587.5,origin="1970-01-01")
per.nl<-3600^2/(34.08*exp(20.7233))
records<-list()
for(i in seq_len(nrow(ev))){
 e<-ev[i,];morning<-e$event=="heliacal rising"
 ref<-time2jd(paste(e$date.julian,"00:00"),timezone="UTC",calendar="J")
 key<-if(e$star=="Sirius") ",alf CMa" else ",alf Tau"
 horflag<-2L+2048L+32768L
 first<-sketchSky(hor,jd=ref,display=list(moon=FALSE),max.mag=2)
 swephR::swe_set_topo(geo[1],geo[2],geo[3])
 M<-swephR::swe_fixstar2_mag(key)$mag
 for(d in -32:32){
   day<-ref+d
   # Local date convention is UTC here. Rising windows in morning, setting in evening.
   hours<-if(morning) seq(0,10,1/60) else seq(15,24-1/60,1/60)
   jd<-day+hours/24
   eq<-swephR::swe_fixstar2_ut(key,day+.5,horflag)$xx[1:3]
   altaz<-function(j,eq)swephR::swe_azalt(j,1L,geo,1013.25,15,eq)$xaz
   starpos<-t(vapply(jd,altaz,numeric(3),eq=eq))
   ix<-which(starpos[,2]>0 & starpos[,2]<=10)
   if(!length(ix))next
   dt<-jd[ix];starpos<-starpos[ix,,drop=FALSE]
   sun<-t(vapply(dt,function(j)altaz(j,swephR::swe_calc_ut(j,0L,horflag)$xx[1:3]),numeric(3)))
   rising<-(starpos[,1]+180)%%360<180
   keep<-which(rising==morning & sun[,2]<0)
   if(!length(keep)){
     records[[length(records)+1L]]<-data.frame(i,d,jd=day,margin=-Inf,gate=-Inf,best=NA)
     next
   }
   v<-eq*pi/180
   id<-skyscapeR:::.galactic_classes(cos(v[2])*cos(v[1]),cos(v[2])*sin(v[1]),sin(v[2]),t(skyscapeR:::.galactic_frame(day)))
   mw<-skyscapeR:::.mw_radiance()[id+1]*skyscapeR:::.band_transmission(starpos[keep,2],.2)/per.nl
   bg<-vapply(keep,function(k)f(list(sun=c((sun[k,1]+180)%%360,sun[k,2]),k=.2,calibrated.twilight=TRUE),(starpos[k,1]+180)%%360,starpos[k,2]),numeric(1))+mw
   margin<-threshold(bg)-M-skyscapeR:::.stellar_extinction(starpos[keep,3],.2)
   gated<-margin[sun[keep,2] < -6]
   records[[length(records)+1L]]<-data.frame(i,d,jd=day,margin=max(margin),gate=if(length(gated))max(gated) else -Inf,best=dt[keep[which.max(margin)]])
 }
}
z<-do.call(rbind,records);write.csv(z,"windows.csv",row.names=FALSE)
for(i in seq_len(nrow(ev))){
 q<-z[z$i==i,];morning<-ev$event[i]=="heliacal rising"
 print(ev[i,c("star","event")])
 for(col in c("margin","gate"))for(limit in c(5.5,6,6.5)){
   yes<-which(q[[col]]+limit-6>=0)
   if(!length(yes))next
   j<-if(morning)min(yes) else max(yes)
   cat(col,limit,as.character(toDate(q$jd[j])),q$d[j],"\n")
   stopifnot(abs(q$d[j])<32)
   if(col=="margin" && limit==6){
    # Re-evaluate the transition instant through the normal renderer.
    view<-sketchSky(hor,jd=q$best[j],display=list(moon=FALSE),max.mag=2)
    name<-if(ev$star[i]=="Sirius")"alf CMa" else "alf Tau"
    pos<-view$stars[name,]; rr<-skyRadiance(view,pos[1]+c(-.01,.01),pos[2]+c(-.01,.01),c(1,1))
    B<-(rr$radiance-rr$direct)/per.nl
    check<-threshold(B)+2.5*log10(view$light.budget$stars[name])
    cat("direct renderer margin",check,"scan",q$margin[j],"\n")
    stopifnot(abs(check-q$margin[j])<.03)
   }
 }
}
