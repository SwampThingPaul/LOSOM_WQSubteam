## 
## LOSOM WQ Subteam
## Nutrient Loads
##
##
## Code was compiled by Paul Julian
## contact info: paul.julian@floridadep.gov

## BAD ## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

#Libraries
library(AnalystHelper);#devtools::install_github("SwampThingPaul/AnalystHelper")
library(plyr)
library(reshape)

# Additional Analysis libraries
library(zyp)
library(mblm)
library(rkt)
library(car)
library(pedometrics)
library(MASS)

#Paths
wd="D:/Work/LOSOM_WQ"

paths=paste0(wd,c("/Plots/","/Export/","/Data/"))
#Folder.Maker(paths);#One and done. Creates folders in working directory.
plot.path=paths[1]
export.path=paths[2]
data.path=paths[3]

## Functions
qq.function=function(y){
  n=length(y)
  r=order(order(y))
  if(n>10){p=(r-1/2)/n}else{p=(r-3/8)/(n+1/4)}
  qqnorm.val=qnorm(p)
  return(qqnorm.val)
}

decimal.WY=function(date,WY.type="FL"){
  # calculates decimal water year (based on Florida WY) 
  # similar to lubridate::decimal_dates()
  
  require(AnalystHelper)
  require(lubridate)
  
  Y <- year(date)
  WY <- WY(date,WY.type=WY.type)
  
  if(WY.type=="FL"){
    start <- make_datetime(WY-1, 5L, 1L, tz = tz(date))
    end <- make_datetime(WY, 5L, 1L, tz = tz(date))
  }
  else if(WY.type=="Fed"){
    start <- make_datetime(WY-1, 10L, 1L, tz = tz(date))
    end <- make_datetime(WY, 10L, 1L, tz = tz(date))
  }else{NA}
  
  sofar <- as.numeric(difftime(date, start, units = "secs"))
  total <- as.numeric(difftime(end, start, units = "secs"))
  
  dec.dateWY <- WY + sofar/total
  return(dec.dateWY)
}

##
dates=as.Date(c("1998-05-01","2019-04-30"))
# CRE Hydro ---------------------------------------------------------------
# Discharge
q.dbkeys=data.frame(SITE=c("S79","S78","S77"),DBKEY=c("00865","DJ236","15635"))
q.cre.dat=data.frame()
for(i in 1:nrow(q.dbkeys)){
  tmp=DBHYDRO_daily(dates[1],dates[2],q.dbkeys$DBKEY[i])
  tmp$DBKEY=as.character(q.dbkeys$DBKEY[i])
  q.cre.dat=rbind(q.cre.dat,tmp)
  print(i)
}
q.cre.dat=merge(q.cre.dat,q.dbkeys,"DBKEY")
q.cre.dat$Date.EST=date.fun(q.cre.dat$Date)

q.cre.dat.xtab=cast(q.cre.dat,Date.EST~SITE,value="Data.Value",mean)
q.cre.dat.xtab$month=format(q.cre.dat.xtab$Date,"%m")
q.cre.dat.xtab$CY=format(q.cre.dat.xtab$Date,"%Y")
q.cre.dat.xtab$monCY=with(q.cre.dat.xtab,date.fun(paste(CY,month,"01",sep="-")))
q.cre.dat.xtab$WY=WY(q.cre.dat.xtab$Date.EST)

q.cre.dat.xtab$S77=with(q.cre.dat.xtab,ifelse(S77<0,0,S77))
q.cre.dat.xtab$S78=with(q.cre.dat.xtab,ifelse(S78<0,0,S78))
q.cre.dat.xtab$C43.basin.in=with(q.cre.dat.xtab,ifelse(S79<S77,0,S79-S77))
q.cre.dat.xtab$basin.ratio=with(q.cre.dat.xtab,ifelse(S79==0,NA,C43.basin.in/S79))

q.cre.dat.xtab.mon=ddply(q.cre.dat.xtab,c("monCY","WY"),summarise,S77=sum(S77,na.rm=T),S78=sum(S78,na.rm=T),S79=sum(S79,na.rm=T))
q.cre.dat.xtab.mon$C43Basin=with(q.cre.dat.xtab.mon,ifelse(S79<S77,0,S79-S77))
q.cre.dat.xtab.mon$basin.q.ratio=with(q.cre.dat.xtab.mon,C43Basin/S79)

plot(C43Basin~monCY,q.cre.dat.xtab.mon,type="l")
plot(basin.q.ratio~monCY,q.cre.dat.xtab.mon,type="l")
plot(S78~monCY,q.cre.dat.xtab.mon,type="l")
with(q.cre.dat.xtab.mon,lines(monCY,S79,col="red"))

# Stage data
stg.dbkeys=data.frame(SITE=c("S79_H","S79_H","S235_T","S235_T"),DBKEY=c("00864","AN786","15566","38259"))
stg.dbkeys=subset(stg.dbkeys,DBKEY!="00864");#not sure of the datum (NGVD vs NAVD)
stg.cre.dat=data.frame()
for(i in 1:nrow(stg.dbkeys)){
  tmp=DBHYDRO_daily(dates[1],dates[2],stg.dbkeys$DBKEY[i])
  tmp$DBKEY=as.character(stg.dbkeys$DBKEY[i])
  stg.cre.dat=rbind(stg.cre.dat,tmp)
  print(i)
}
stg.cre.dat=merge(stg.cre.dat,stg.dbkeys,"DBKEY")
stg.dat.da.xtab=cast(stg.cre.dat,Date~SITE,value="Data.Value",mean)
stg.dat.da.xtab$WY=WY(stg.dat.da.xtab$Date)
stg.dat.da.xtab$month=format(stg.dat.da.xtab$Date,"%m")
stg.dat.da.xtab$CY=format(stg.dat.da.xtab$Date,"%Y")
stg.dat.da.xtab$monCY=with(stg.dat.da.xtab,date.fun(paste(CY,month,"01",sep="-")))
stg.dat.da.xtab$grad=with(stg.dat.da.xtab,(S235_T-S79_H))
#stg.dat.da.xtab$grad=with(stg.dat.da.xtab,(S235_T-S79_H)/213356.3)
head(stg.dat.da.xtab)

plot(grad~Date,stg.dat.da.xtab,type="l")

stg.dat.da.xtab.mon=ddply(stg.dat.da.xtab,c("monCY","WY"),summarise,
                          S79_H=mean(S79_H,na.rm=T),
                          S235_T=mean(S235_T,na.rm=T),
                          grad=mean(grad,na.rm=T))

plot(grad~monCY,stg.dat.da.xtab.mon,type="l");mtext(side=3,"test")


cre.hydro.mon=merge(q.cre.dat.xtab.mon,stg.dat.da.xtab.mon,c("monCY","WY"))
plot(grad~monCY,cre.hydro.mon,type="l")



#### 
####
####
###

# CRE WQ ------------------------------------------------------------------
params=data.frame(Test.Number=c(18,21,80,20,25,23,61,179,100),param=c("NOx","TKN","TN","NH4","TP","SRP","Chla","Chla","TOC"))
wq.dat=DBHYDRO_WQ(dates[1],dates[2],c("S79","S77"),params$Test.Number)
wq.dat=merge(wq.dat,params,"Test.Number")
unique(wq.dat$Collection.Method)
wq.dat=subset(wq.dat,Collection.Method=="G")
wq.dat.xtab=cast(wq.dat,Station.ID+Date.EST~param,value="HalfMDL",mean)
wq.dat.xtab$DIN=with(wq.dat.xtab,NH4+NOx)
wq.dat.xtab$TN=with(wq.dat.xtab, TN_Combine(NOx,TKN,TN))

# Reversal Evaluation
wq.dat.xtab$TPReversal=with(wq.dat.xtab,ifelse(is.na(SRP)==T|is.na(TP)==T,0,ifelse(SRP>(TP*1.3),1,0)));# Reversals identified as 1 reversals consistent with TP rule evaluation
wq.dat.xtab$TNReversal=with(wq.dat.xtab,ifelse(is.na(DIN)==T|is.na(TN)==T,0,ifelse(DIN>(TN*1.3),1,0)));

sum(wq.dat.xtab$TNReversal,na.rm=T)
subset(wq.dat.xtab,TNReversal==T)
sum(wq.dat.xtab$TPReversal,na.rm=T)

par(family="serif",oma=c(1,1,1,1),mar=c(4,4,1,1))
layout(matrix(1:2,1,2,byrow=F))
plot(TN~DIN,wq.dat.xtab,ylab="TN (mg L\u207B\u00b9)",xlab="DIN (mg L\u207B\u00b9)",pch=21,bg=ifelse(TNReversal==1,"dodgerblue1",NA),col=adjustcolor("grey",0.8));abline(0,1,col="dodgerblue1")
plot(TP~SRP,wq.dat.xtab,ylab="TP (mg L\u207B\u00b9)",xlab="SRP (mg L\u207B\u00b9)",pch=21,bg=ifelse(TPReversal==1,"red",NA),col=adjustcolor("grey",0.8));abline(0,1,col="red")

dev.off()

####
wq.dat.xtab$WY=WY(wq.dat.xtab$Date.EST)
wq.dat.xtab$month=format(wq.dat.xtab$Date.EST,"%m")
wq.dat.xtab$CY=format(wq.dat.xtab$Date.EST,"%Y")
wq.dat.xtab$monCY=with(wq.dat.xtab,date.fun(paste(CY,month,"01",sep="-")))

wq.dat.xtab.mon=ddply(subset(wq.dat.xtab,Station.ID=="S79"),c("WY","monCY"),summarise,mean.TP=mean(TP,na.rm=T),mean.TN=mean(TN,na.rm=T),mean.TOC=mean(TOC,na.rm=T))
mean(wq.dat.xtab.mon$mean.TOC,na.rm=T)

plot(mean.TOC~format(wq.dat.xtab.mon$monCY,"%m"),wq.dat.xtab.mon)
wq.dat.xtab.mon.TOC=ddply(subset(wq.dat.xtab,Station.ID=="S79"),c("month"),summarise,POR.mean.TOC=mean(TOC,na.rm=T))
with(wq.dat.xtab.mon.TOC,points(month,POR.mean.TOC,pch=21,bg="Red"))

plot(mean.TN~format(wq.dat.xtab.mon$monCY,"%m"),wq.dat.xtab.mon)
plot(mean.TP~format(wq.dat.xtab.mon$monCY,"%m"),wq.dat.xtab.mon)

cre.hydro.wq.mon=merge(cre.hydro.mon,wq.dat.xtab.mon,c("WY","monCY"))
cre.hydro.wq.mon$hydro.season=FL.Hydroseason(cre.hydro.wq.mon$monCY)
cre.hydro.wq.mon$month=as.numeric(format(cre.hydro.wq.mon$monCY,"%m"))
cre.hydro.wq.mon=merge(cre.hydro.wq.mon,data.frame(month=c(5:12,1:4),month.plot=1:12),"month",all.x=T)
cre.hydro.wq.mon=cre.hydro.wq.mon[order(cre.hydro.wq.mon$monCY),]


#tiff(filename=paste0(plot.path,"S79_MonthlyTPTN.tiff"),width=4.75,height=6,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
#png(filename=paste0(plot.path,"png/S79_MonthlyTPTN.png"),width=4.75,height=6,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.75,0.75),oma=c(2,1,0.5,0.5));
layout(matrix(1:2,2,1))

ylim.val=c(0.05,0.35);ymaj=c(0.05,log.scale.fun(ylim.val,"major"),0.3);ymin=log.scale.fun(ylim.val,"minor")
xlim.val=c(1,12);xmaj=1:12;xmaj.lab=c(5:12,1:4)
plot(mean.TP~month.plot,cre.hydro.wq.mon,log="y",ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(cre.hydro.wq.mon,points(month.plot,mean.TP,pch=21,bg=adjustcolor("dodgerblue1",0.25),col="dodgerblue1",lwd=0.1,cex=1.25))
k=predict(loess(mean.TP~month.plot,cre.hydro.wq.mon),data.frame(month.plot=seq(1,12,length.out = 24)),se=T)
lines(seq(1,12,length.out = 24),k$fit,lwd=2,col="red")
lines(seq(1,12,length.out = 24),k$fit - qt(0.975,k$df)*k$se, lty=2,lwd=1,col="red")
lines(seq(1,12,length.out = 24),k$fit + qt(0.975,k$df)*k$se, lty=2,lwd=1,col="red")
axis_fun(1,line=-0.5,xmaj,xmaj,NA)
axis_fun(2,ymaj,ymin,format(ymaj*1000));box(lwd=1)
mtext(side=3,"S-79")

mtext(side=2,line=2.5,"Total Phosphorus (\u03BCg L\u207B\u00B9)")
legend("topright",c("Monthly Mean","LOESS \u00B1 95% CI"),
       pch=c(21,NA),
       lty=c(NA,1),lwd=c(0.01,1),
       col=c("dodgerblue1","red"),
       pt.bg=c(adjustcolor("dodgerblue1",0.25),NA),
       pt.cex=1,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)

ylim.val=c(0.5,3);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)#log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
plot(mean.TN~month.plot,cre.hydro.wq.mon,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(cre.hydro.wq.mon,points(month.plot,mean.TN,pch=21,bg=adjustcolor("dodgerblue1",0.25),col="dodgerblue1",lwd=0.1,cex=1.25))
k=predict(loess(mean.TN~month.plot,cre.hydro.wq.mon),data.frame(month.plot=seq(1,12,length.out = 24)),se=T)
lines(seq(1,12,length.out = 24),k$fit,lwd=2,col="red")
lines(seq(1,12,length.out = 24),k$fit - qt(0.975,k$df)*k$se, lty=2,lwd=1,col="red")
lines(seq(1,12,length.out = 24),k$fit + qt(0.975,k$df)*k$se, lty=2,lwd=1,col="red")
axis_fun(1,line=-0.5,xmaj,xmaj,xmaj.lab)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=1,line=1.75,"Month")
mtext(side=2,line=2.5,"Total Nitrogen (mg L\u207B\u00B9)")
dev.off()

## trend analysis
# Monthly
wq.dat.xtab.mon$decWY=decimal.WY(wq.dat.xtab.mon$monCY)
wq.dat.xtab.mon$month=as.numeric(format(wq.dat.xtab.mon$monCY,"%m"))
wq.dat.xtab.mon$CY=as.numeric(format(wq.dat.xtab.mon$monCY,"%Y"))

WYs=WY(dates[1]):WY(dates[2])
months=rep(1:12,length(WYs))
wq.dat.xtab.mon2=merge(data.frame(WY=sort(rep(WYs,12)),month=months,fill=1),wq.dat.xtab.mon,c("WY","month"),all.x=T)

#TP
#Autocorrelation
acf(wq.dat.xtab.mon$mean.TP)
pacf(wq.dat.xtab.mon$mean.TP)
#test=data.frame(pacf=ar.yw(wq.dat.xtab.mon$mean.TP,floor(10 * (log10(length(wq.dat.xtab.mon$mean.TP)))))$partialacf)
#points(1:22,test$pacf)
# 95% confidence (i.e. dashed blue lines) qnorm((1+0.95)/2)/sqrt(n.used)

acf.dat.rslt=data.frame()
for(h in 0:24){
  #demean
  #x=sweep(as.matrix(wq.dat.xtab.mon$mean.TP),2,colMeans(as.matrix(wq.dat.xtab.mon$mean.TP),na.rm=T))
  lagged=lag(as.zoo(wq.dat.xtab.mon$mean.TP),-h,na.pad=T)
  tmp.dat=as.zoo(wq.dat.xtab.mon$mean.TP)
  stat=with(data.frame(lag=lagged,dat=tmp.dat),cor.test(lag,dat,method="pearson"))
  acf.dat.rslt=rbind(acf.dat.rslt,data.frame(lag=h,estimate=as.numeric(stat$estimate),pval=stat$p.value))
}
acf.dat.rslt
acf(wq.dat.xtab.mon$mean.TP)
points(acf.dat.rslt$lag,acf.dat.rslt$estimate)

tseries::adf.test(wq.dat.xtab.mon$mean.TP)
#the null hypothesis is reject and accept alternative(i.e. stationary)
#https://nwfsc-timeseries.github.io/atsa-labs/sec-boxjenkins-aug-dickey-fuller.html
#http://rstudio-pubs-static.s3.amazonaws.com/273431_12485645a6b743cf8e42731415b8003c.html
#https://otexts.com/fpp2/

ylim.val=c(-0.50,1.1);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,24);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
#tiff(filename=paste0(plot.path,"S79_TPACF.tiff"),width=5,height=3.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
#png(filename=paste0(plot.path,"png/S79_TPACF.png"),width=5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.75,0.75),oma=c(2,1,0.5,0.5));
plot(estimate~lag,acf.dat.rslt,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
ci.val=qnorm((1+0.95)/2)/sqrt(length(wq.dat.xtab.mon$mean.TP))
#abline(h=c(ci.val,-ci.val),lty=2,lwd=1.5,col="blue")
polygon(c(-1,25,25,-1),c(ci.val,ci.val,-ci.val,-ci.val),col=adjustcolor("grey",0.5),border=0)
with(acf.dat.rslt,segments(lag,0,lag,estimate,lwd=1.5,lty=2))
with(acf.dat.rslt,points(lag,estimate,pch=21,bg=ifelse(pval<0.05,"indianred1","dodgerblue1"),lwd=0.01))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,expression(paste("ACF ",italic("r")["Pearson"])))
mtext(side=1,line=1.5,"Lag")
legend("topright",legend=c("\u03C1 < 0.05","\u03C1 > 0.05","95% CI"),
       pch=c(21,21,22),pt.bg=c("indianred1","dodgerblue1",adjustcolor("grey",0.5)),col=c("black","black",NA),
       lty=NA,lwd=c(0.1,0.1,0),pt.cex=1.5,cex=0.7,ncol=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

ylim.val=c(0.05,0.35);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")#by.y=0.05;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=ylim.val;xmaj=log.scale.fun(xlim.val,"major");xmin=log.scale.fun(xlim.val,"minor")#by.x=0.05;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
#tiff(filename=paste0(plot.path,"S79_lagplots.tiff"),width=7,height=2.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
#png(filename=paste0(plot.path,"png/S79_lagplots.png"),width=7,height=2.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.75,0.75),oma=c(3,1.5,0.5,0.5));
layout(matrix(1:4,1,4))
for(h in 0:3){
  tmp.dat=data.frame(dat=as.zoo(wq.dat.xtab.mon$mean.TP),lag=lag(as.zoo(wq.dat.xtab.mon$mean.TP),-h,na.pad=T))
  plot(dat~lag,tmp.dat,ylim=ylim.val,xlim=xlim.val,log="xy",type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  abline(0,1,lty=2)
  with(tmp.dat,points(lag,dat,pch=21,bg=adjustcolor("grey",0.25),lwd=0.01,cex=1.25))
  #LOESS model
  x.val=seq(min(tmp.dat$lag,na.rm=T),max(tmp.dat$lag,na.rm=T),length.out=100)
  k=predict(loess(dat~lag,tmp.dat),data.frame(lag=x.val),se=T)
  lines(x.val,k$fit,lwd=2,col="red")
  lines(x.val,k$fit - qt(0.975,k$df)*k$se, lty=2,lwd=1.25,col="red")
  lines(x.val,k$fit + qt(0.975,k$df)*k$se, lty=2,lwd=1.25,col="red")
  if(h==0){legend("topleft",legend="LOESS \u00B1 95% CI",lty=1,lwd=1,col="red",cex=0.75,ncol=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)}
  if(h==0){axis_fun(2,ymaj,ymin,ymaj*1000,cex=1)}else{axis_fun(2,ymaj,ymin,NA)}
  if(h==0){mtext(side=2,line=2,"TP (\u03BCg L\u207B\u00B9)",cex=1)}
  axis_fun(1,line=-0.5,xmaj,xmin,xmaj*1000,cex=1);box(lwd=1)
  mtext(side=3,line=0.2,paste0("Lag ",h),cex=0.75)
}
mtext(side=1,line=1,outer=T,"Lagged TP (\u03BCg L\u207B\u00B9)",cex=1)
dev.off()


#Trend
with(wq.dat.xtab.mon,cor.test(mean.TP,decWY,method="kendall"))
#with(subset(wq.dat.xtab.mon,is.na(mean.TP)==F),zyp.sen(mean.TP~decWY))
with(subset(wq.dat.xtab.mon,is.na(mean.TP)==F),mblm(mean.TP~decWY,repeated = T))
plot(mean.TP~decWY,wq.dat.xtab.mon,las=1,type="b")
abline(with(subset(wq.dat.xtab.mon,is.na(mean.TP)==F),mblm(mean.TP~decWY,repeated = T)),col="red")

#seasonal kendall
wq.dat.xtab.mon$hydro.season=FL.Hydroseason(wq.dat.xtab.mon$monCY)
with(subset(wq.dat.xtab.mon,is.na(mean.TP)==F),rkt(decimal_date(monCY),mean.TP,month))

#TN
#Autocorrelation
acf(subset(wq.dat.xtab.mon,is.na(mean.TN)==F)$mean.TN)
pacf(subset(wq.dat.xtab.mon,is.na(mean.TN)==F)$mean.TN)

acf.TN.dat.rslt=data.frame()
for(h in 0:24){
  #demean
  #x=sweep(as.matrix(wq.dat.xtab.mon$mean.TP),2,colMeans(as.matrix(wq.dat.xtab.mon$mean.TP),na.rm=T))
  lagged=lag(as.zoo(na.omit(wq.dat.xtab.mon$mean.TN)),-h,na.pad=T)
  tmp.dat=as.zoo(na.omit(wq.dat.xtab.mon$mean.TN))
  stat=with(data.frame(lag=lagged,dat=tmp.dat),cor.test(lag,dat,method="pearson"))
  acf.TN.dat.rslt=rbind(acf.TN.dat.rslt,data.frame(lag=h,estimate=as.numeric(stat$estimate),pval=stat$p.value))
}
acf.TN.dat.rslt
acf(na.omit(wq.dat.xtab.mon$mean.TN))
points(acf.TN.dat.rslt$lag,acf.TN.dat.rslt$estimate)

tseries::adf.test(na.omit(wq.dat.xtab.mon$mean.TN))

ylim.val=c(-0.05,1.1);by.y=0.2;ymaj=seq(max(ylim.val[1],0),ylim.val[2],by.y);ymin=seq(max(ylim.val[1],0),ylim.val[2],by.y/2)
xlim.val=c(0,24);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
#tiff(filename=paste0(plot.path,"S79_TNACF.tiff"),width=5,height=3.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
#png(filename=paste0(plot.path,"png/S79_TPACF.png"),width=5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.75,0.75),oma=c(2,1,0.5,0.5));
plot(estimate~lag,acf.TN.dat.rslt,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
ci.val=qnorm((1+0.95)/2)/sqrt(length(na.omit(wq.dat.xtab.mon$mean.TN)))
#abline(h=c(ci.val,-ci.val),lty=2,lwd=1.5,col="blue")
polygon(c(-1,25,25,-1),c(ci.val,ci.val,-ci.val,-ci.val),col=adjustcolor("grey",0.5),border=0)
with(acf.TN.dat.rslt,segments(lag,0,lag,estimate,lwd=1.5,lty=2))
with(acf.TN.dat.rslt,points(lag,estimate,pch=21,bg=ifelse(pval<0.05,"indianred1","dodgerblue1"),lwd=0.01))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,expression(paste("ACF ",italic("r")["Pearson"])))
mtext(side=1,line=1.5,"Lag")
legend("topright",legend=c("\u03C1 < 0.05","\u03C1 > 0.05","95% CI"),
       pch=c(21,21,22),pt.bg=c("indianred1","dodgerblue1",adjustcolor("grey",0.5)),col=c("black","black",NA),
       lty=NA,lwd=c(0.1,0.1,0),pt.cex=1.5,cex=0.7,ncol=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

ylim.val=c(0.75,3);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=ylim.val;by.x=0.5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
#tiff(filename=paste0(plot.path,"S79_TN_lagplots.tiff"),width=7,height=2.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
#png(filename=paste0(plot.path,"png/S79_TN_lagplots.png"),width=7,height=2.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2.5,0.75,0.75),oma=c(3,1.5,0.5,0.5));
layout(matrix(1:4,1,4))
for(h in 0:3){
  tmp.dat=data.frame(dat=as.zoo(na.omit(wq.dat.xtab.mon$mean.TN)),lag=lag(as.zoo(na.omit(wq.dat.xtab.mon$mean.TN)),-h,na.pad=T))
  plot(dat~lag,tmp.dat,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  abline(0,1,lty=2)
  with(tmp.dat,points(lag,dat,pch=21,bg=adjustcolor("grey",0.25),lwd=0.01,cex=1.25))
  #LOESS model
  x.val=seq(min(tmp.dat$lag,na.rm=T),max(tmp.dat$lag,na.rm=T),length.out=100)
  k=predict(loess(dat~lag,tmp.dat),data.frame(lag=x.val),se=T)
  lines(x.val,k$fit,lwd=2,col="red")
  lines(x.val,k$fit - qt(0.975,k$df)*k$se, lty=2,lwd=1.25,col="red")
  lines(x.val,k$fit + qt(0.975,k$df)*k$se, lty=2,lwd=1.25,col="red")
  if(h==0){legend("topleft",legend="LOESS \u00B1 95% CI",lty=1,lwd=1,col="red",cex=0.75,ncol=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)}
  if(h==0){axis_fun(2,ymaj,ymin,ymaj,cex=1)}else{axis_fun(2,ymaj,ymin,NA)}
  if(h==0){mtext(side=2,line=2.5,"TN (mg L\u207B\u00B9)",cex=1)}
  axis_fun(1,line=-0.5,xmaj,xmin,xmaj,cex=1);box(lwd=1)
  mtext(side=3,line=0.2,paste0("Lag ",h),cex=0.75)
}
mtext(side=1,line=1,outer=T,"Lagged TN (mg L\u207B\u00B9)",cex=1)
dev.off()

#Trend
with(wq.dat.xtab.mon,cor.test(mean.TN,decWY,method="kendall"))
with(subset(wq.dat.xtab.mon,is.na(mean.TN)==F),mblm(mean.TN~decWY,repeated = T))
plot(mean.TN~decWY,wq.dat.xtab.mon,las=1,type="b")
abline(with(subset(wq.dat.xtab.mon,is.na(mean.TN)==F),mblm(mean.TN~decWY,repeated = T)),col="red")

#seasonal kendall
with(subset(wq.dat.xtab.mon,is.na(mean.TN)==F),rkt(decimal_date(monCY),mean.TN,month))

# Annual
wq.dat.xtab.WY=ddply(subset(wq.dat.xtab,Station.ID=="S79"),c("WY"),summarise,mean.TP=mean(TP,na.rm=T),N.TP=N(TP),mean.TN=mean(TN,na.rm=T),N.TN=N(TN))
with(subset(wq.dat.xtab.WY,N.TP>=4),cor.test(mean.TP,WY,method="kendall"))
with(subset(wq.dat.xtab.WY,is.na(mean.TP)==F&N.TP>=4),mblm(mean.TP~WY,repeated = T))
plot(mean.TP~WY,wq.dat.xtab.WY,las=1,type="b")
abline(with(subset(wq.dat.xtab.WY,is.na(mean.TP)==F&N.TP>=4),mblm(mean.TP~WY,repeated = T)),col="Red")

with(subset(wq.dat.xtab.WY,N.TN>=4),cor.test(mean.TN,WY,method="kendall"))
with(subset(wq.dat.xtab.WY,is.na(mean.TN)==F&N.TN>=4),mblm(mean.TN~WY,repeated = T))
plot(mean.TN~WY,wq.dat.xtab.WY,las=1,type="b")
abline(with(subset(wq.dat.xtab.WY,is.na(mean.TN)==F&N.TP>=4),mblm(mean.TN~WY,repeated = T)),col="Red")


#tiff(filename=paste0(plot.path,"S79_MonthlyTPTN_TS.tiff"),width=5,height=6,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,3,0.75,0.75),oma=c(2,1,0.5,0.5));
layout(matrix(1:2,2,1))
ylim.val=c(0.05,0.4);ymaj=c(0.05,log.scale.fun(ylim.val,"major"),0.4);ymin=log.scale.fun(ylim.val,"minor")
xlim.val=c(1999,2020);xmaj=seq(xlim.val[1],xlim.val[2],5);xmin=seq(xlim.val[1],xlim.val[2],1)
plot(mean.TP~decWY,wq.dat.xtab.mon,log="y",ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(wq.dat.xtab.mon,pt_line(decWY,mean.TP,2,adjustcolor("dodgerblue1",0.5),1,21,adjustcolor("dodgerblue1",0.5),cex=1))
mblm.mod=with(subset(wq.dat.xtab.mon,is.na(mean.TP)==F),mblm(mean.TP~decWY,repeated = T))
mod.pred=predict(mblm.mod,data.frame(decWY=seq(1999,2020,1)),interval="confidence")
shaded.range(seq(1999,2020,1),mod.pred[,2],mod.pred[,3],"black",lty=0)
lines(seq(1999,2020,1),mod.pred[,1],lty=2)
#axis_fun(1,line=-0.5,xmaj,xmin,paste(5,xmaj-1,sep=" - "))
axis_fun(1,line=-0.5,xmaj,xmin,NA)
axis_fun(2,ymaj,ymin,format(ymaj*1000));box(lwd=1)
mtext(side=3,"S-79")
mtext(side=2,line=2.5,"Total Phosphorus (\u03BCg L\u207B\u00B9)")
legend("topright",c("Monthly Mean","Thiel-Sen \u00B1 95% CI"),
       pch=c(NA,22),
       lty=c(2,2),lwd=c(0.01,1),
       col=c("dodgerblue1","grey"),
       pt.bg=c(adjustcolor("dodgerblue1",0.25),"grey"),
       pt.cex=1,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,text.col="white")
legend("topright",c("Monthly Mean","Thiel-Sen \u00B1 95% CI"),
       pch=c(21,NA),
       lty=c(NA,2),lwd=c(0.01,1),
       col=c("black","black"),
       pt.bg=c(adjustcolor("dodgerblue1",0.5),"grey"),
       pt.cex=1,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)

ylim.val=c(0.9,3);by.y=0.5;ymaj=seq(max(1,ylim.val[1]),ylim.val[2],by.y);ymin=seq(max(1,ylim.val[1]),ylim.val[2],by.y/2)#log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
plot(mean.TN~decWY,wq.dat.xtab.mon,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(wq.dat.xtab.mon,pt_line(decWY,mean.TN,2,adjustcolor("dodgerblue1",0.5),1,21,adjustcolor("dodgerblue1",0.5),cex=1))
mblm.mod=with(subset(wq.dat.xtab.mon,is.na(mean.TN)==F),mblm(mean.TN~decWY,repeated = T))
mod.pred=predict(mblm.mod,data.frame(decWY=seq(1999,2020,1)),interval="confidence")
shaded.range(seq(1999,2020,1),mod.pred[,2],mod.pred[,3],"black",lty=0)
lines(seq(1999,2020,1),mod.pred[,1],lty=2)
axis_fun(1,line=-0.5,xmaj,xmin,paste(5,xmaj-1,sep="/"))
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=1,line=1.75,"Date (Month - Year)")
mtext(side=2,line=2.5,"Total Nitrogen (mg L\u207B\u00B9)")
dev.off()

#tiff(filename=paste0(plot.path,"S79_MonthlyTP_hist.tiff"),width=3,height=6.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(0.25,3,0.75,2),oma=c(4,1,0.5,0.5));
layout(matrix(1:12,12,1))

xlim.val=c(0.06,0.35);by.x=0.05;xmaj=seq(min(0,xlim.val[1]),xlim.val[2],by.x);xmin=seq(min(0,xlim.val[1]),xlim.val[2],by.x/2)
ylim.val=c(0,6);by.y=3;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
for(i in 1:12){
  tmp.dat=subset(wq.dat.xtab.mon,month==i)$mean.TP
  w=0.01;#width
  startcat <- floor(min(tmp.dat,na.rm=TRUE)/w)*w
  breaks.val <- seq(startcat,max(tmp.dat,na.rm=TRUE)+w,w)
  hist(tmp.dat,main=NA,xlim=xlim.val,ylim=ylim.val,breaks=breaks.val,xlab=NA,ylab=NA,col="skyblue3",yaxs="i",axes=F)
  if(i==12){axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj))}else(axis_fun(1,xmaj,xmin,NA))
  axis_fun(2,ymaj,ymin,ymaj)
  box(lwd=1)
  mtext(month.abb[i],las=3,side=4)
  if(i==12){mtext(side=1,line=2.5,"Total Phosphorus (\u03BCg L\u207B\u00B9)")}
}
mtext(side=2,line=-0.5,outer=T,"Frequency")
dev.off()

#tiff(filename=paste0(plot.path,"S79_MonthlyTPWY.tiff"),width=3,height=6.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(0.25,3.5,0.75,1.75),oma=c(4,1,0.5,0.5));
layout(matrix(1:12,12,1))

xlim.val=c(1999,2019);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(0,0.35);by.y=0.15;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
for(i in 1:12){
  tmp.dat=subset(wq.dat.xtab.mon2,month==i)
  plot(mean.TP~WY,tmp.dat,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  #with(tmp.dat,points(WY,mean.TP,pch=21,bg="skyblue2",lwd=0.1,cex=1.25))
  with(tmp.dat,pt_line(WY,mean.TP,1,"skyblue2",1.5,21,"skyblue2",pt.lwd=0.1,cex=1.1))
  if(i==12){axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj))}else(axis_fun(1,xmaj,xmin,NA))
  axis_fun(2,ymaj,ymin,format(ymaj*1000))
  box(lwd=1)
  mtext(month.abb[i],las=3,side=4,cex=0.9,line=0.25)
  if(i==12){mtext(side=1,line=2.5,"Water Year")}
}
mtext(side=2,line=-0.5,outer=T,"Total Phosphorus (\u03BCg L\u207B\u00B9)")
dev.off()

#tiff(filename=paste0(plot.path,"S79_MonthlyTNWY.tiff"),width=3,height=6.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(0.25,3.5,0.75,1.75),oma=c(4,1,0.5,0.5));
layout(matrix(1:12,12,1))

xlim.val=c(1999,2019);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(0,3);by.y=1.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
for(i in 1:12){
  tmp.dat=subset(wq.dat.xtab.mon2,month==i)
  plot(mean.TN~WY,tmp.dat,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  #with(tmp.dat,points(WY,mean.TP,pch=21,bg="skyblue2",lwd=0.1,cex=1.25))
  with(tmp.dat,pt_line(WY,mean.TN,1,"skyblue2",1.5,21,"skyblue2",pt.lwd=0.1,cex=1.1))
  if(i==12){axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj))}else(axis_fun(1,xmaj,xmin,NA))
  axis_fun(2,ymaj,ymin,format(ymaj))
  box(lwd=1)
  mtext(month.abb[i],las=3,side=4,cex=0.9,line=0.25)
  if(i==12){mtext(side=1,line=2.5,"Water Year")}
}
mtext(side=2,line=-0.5,outer=T,"Total Nitrogen (mg L\u207B\u00B9)")
dev.off()



###
## TP Model
#vars=c("mean.TP","grad","C43Basin","S78","S79","S77","basin.q.ratio","hydro.season")
vars=c("mean.TP","grad","C43Basin","S77","basin.q.ratio","hydro.season")
S79.TP.mod=lm(log(mean.TP)~.,na.omit(cre.hydro.wq.mon[,vars]))
layout(matrix(1:4,2,2));plot(S79.TP.mod)
shapiro.test(residuals(S79.TP.mod));hist(residuals(S79.TP.mod))
vif(S79.TP.mod)
summary(S79.TP.mod)

#AIC model
S79.TP.mod.sw=stepAIC(S79.TP.mod,direction="both",trace=F)
S79.TP.mod.sw$anova

#tiff(filename=paste0(plot.path,"S79_TPModel_diag.tiff"),width=6,height=5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
#png(filename=paste0(plot.path,"png/S79_TPModel_diag2.png"),width=6,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,3,1,0.5),oma=c(2,1,0.25,0.75));
layout(matrix(1:4,2,2))

cols="grey"#rainbow(length(S79.TP.mod$fitted.values))
ylim.val=c(-0.6,0.6);by.y=0.4;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(-2.6,-1.2);by.x=0.4;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(S79.TP.mod$fitted.values,S79.TP.mod$residuals,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
abline(h=0)
with(S79.TP.mod,points(fitted.values,residuals,pch=21,bg=adjustcolor(cols,0.5),col="grey",lwd=0.1,cex=1.25))
fit.vals=S79.TP.mod$fitted.values
res.vals=S79.TP.mod$residuals
with(lowess(fit.vals,res.vals),lines(x,y,col="red",lwd=2))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,format(round(ymaj,2)));box(lwd=1)
mtext(side=2,line=2.5,"Residuals")
mtext(side=1,line=1.75,"Fitted Values")

#cols=viridisLite::viridis(length(S79.TP.mod$fitted.values))
ylim.val=c(-3,3);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(-3,3);by.x=1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
rstd=rstandard(S79.TP.mod)
qq.x=qq.function(S79.TP.mod$residuals)
plot(rstd~qq.x,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
points(qq.x,rstd,pch=21,bg=adjustcolor(cols,0.5),col="grey",lwd=0.1,cex=1.25)
abline(0,1,lty=3)
axis_fun(1,line=-0.5,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2.5,"Standardized Residuals")
mtext(side=1,line=1.75,"Theoretical Quantiles")

#cols=wesanderson::wes_palette("Zissou1",length(S79.TP.mod$fitted.values),"continuous")
ylim.val=c(0,2);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(-2.6,-1.2);by.x=0.4;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(sqrt(abs(rstd))~fit.vals,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
points(fit.vals,sqrt(abs(rstd)),pch=21,bg=adjustcolor(cols,0.5),col="grey",lwd=0.1,cex=1.25)
with(lowess(fit.vals,sqrt(abs(rstd))),lines(x,y,col="red",lwd=2))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2,expression(sqrt("Standardized Residuals")))
mtext(side=1,line=1.75,"Fitted Values")

#plot(S79.TP.mod,which=4)
#plot(cooks.distance(S79.TP.mod))
#plot(S79.TP.mod,which=5)
#cols=viridisLite::magma(length(S79.TP.mod$fitted.values))
ylim.val=c(-3,3);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,0.35);by.x=0.1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
lev=hatvalues(S79.TP.mod)
plot(rstd~lev,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
abline(h=0)
points(lev,rstd,pch=21,bg=adjustcolor(cols,0.5),col="grey",lwd=0.1,cex=1.25)
with(lowess(lev,rstd),lines(x,y,col="red",lwd=2))
inf=lm.influence(S79.TP.mod)
hh=seq(min(range(inf)[1],range(inf)[2]/100),xlim.val[2]+(xlim.val[2]*0.5),length.out=101)
hh=hh[hh>0]
crit.cook=0.5
cl.h=sqrt(crit.cook*length(coef(S79.TP.mod))*(1-hh)/hh)
lines(hh,cl.h,lty=2,col=2)
lines(hh,-cl.h,lty=2,col=2)
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj));axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2,"Standardized Residuals")
mtext(side=1,line=1.75,"Leverage")
dev.off()


shapiro.test(residuals(S79.TP.mod))
hist(residuals(S79.TP.mod))
#hist(residuals(S79.TP.mod),xlab="Model Residuals",main=NA,col="grey",las=1,ylim=c(0,40),yaxs="i",xaxt="n");box(lwd=1)
#mtext(side=2,line=2.5,"Frequency")
#mtext(side=1,line=1.75,"Model Residuals")
#axis_fun(1,line=-0.5,seq(-0.6,0.6,0.4),seq(-0.6,0.6,0.2),seq(-0.6,0.6,0.4))
#dev.off()

S79.TP.pred.mod=with(na.omit(cre.hydro.wq.mon[,vars]),exp(predict(S79.TP.mod,data.frame(grad=grad,C43Basin=C43Basin,S77=S77,basin.q.ratio=basin.q.ratio,hydro.season=hydro.season),interval="confidence")))

delta=na.omit(cre.hydro.wq.mon[,vars])$mean.TP-S79.TP.pred.mod[,1]
mod.d=density(na.omit(as.numeric(S79.TP.pred.mod[,1])))
samp.d=density(cre.hydro.wq.mon$mean.TP)
mod.compare=data.frame(dat=c(S79.TP.pred.mod[,1],cre.hydro.wq.mon$mean.TP),group=c(rep("mod",length(S79.TP.pred.mod[,1])),rep("obs",length(cre.hydro.wq.mon$mean.TP))))
kruskal.test(dat~group,mod.compare)

ylim.val=c(0.04,0.25);ymaj=c(0.05,log.scale.fun(ylim.val,"major"),0.5);ymin=log.scale.fun(ylim.val,"minor")
#tiff(filename=paste0(plot.path,"S79_TPModel_Obs.tiff"),width=3.5,height=3,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,3,1,0.5),oma=c(2,1,0.75,0.75));

boxplot(dat~group,mod.compare,ylim=ylim.val,axes=F,ylab=NA,xlab=NA,outline=F,col="grey",log="y")
axis_fun(1,1:2,1:2,c("Modelled","Observed"))
axis_fun(2,ymaj,ymin,format(ymaj*1000));box(lwd=1)
mtext(side=2,line=2.5,"Total Phosphorus (\u03BCg L\u207B\u00B9)")
dev.off()

#tiff(filename=paste0(plot.path,"S79_TPModel_compare.tiff"),width=6,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(2,3,1,0.5),oma=c(2,1,0.75,0.75));
layout(matrix(c(1,2,3,3),2,2,byrow=T))

ylim.val=c(0.04,0.5);ymaj=c(0.05,log.scale.fun(ylim.val,"major"),0.5);ymin=log.scale.fun(ylim.val,"minor")
xlim.val=ylim.val;xmaj=c(0.05,log.scale.fun(xlim.val,"major"),0.5);xmin=log.scale.fun(xlim.val,"minor")
plot(mean.TP~S79.TP.pred.mod[,1],na.omit(cre.hydro.wq.mon[,vars]),ylim=ylim.val,xlim=xlim.val,type="n",log="xy",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(na.omit(cre.hydro.wq.mon[,vars]),points(S79.TP.pred.mod[,1],mean.TP,pch=21,bg=adjustcolor("forestgreen",0.25),lwd=0.1,col="grey"))
abline(0,1,lty=2)
axis_fun(1,line=-0.5,xmaj,xmin,xmaj*1000);axis_fun(2,ymaj,ymin,ymaj*1000);box(lwd=1)
mtext(side=1,line=1.5,"Modelled TP (\u03BCg L\u207B\u00B9)",cex=0.9)
mtext(side=2,line=2.5,"Observed TP (\u03BCg L\u207B\u00B9)",cex=0.9)

ylim.val=c(0,25);by.y=5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0.05,0.5);xmaj=c(0.05,log.scale.fun(xlim.val,"major"),0.5);xmin=log.scale.fun(xlim.val,"minor")#by.x=0.05;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(mod.d,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA,main=NA,log="x",yaxs="i")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
polygon(mod.d,col=adjustcolor("indianred1",0.5),lwd=1,border="indianred1")
polygon(samp.d,col=adjustcolor("dodgerblue1",0.5),lwd=1,border="dodgerblue1")
axis_fun(1,line=-0.5,xmaj,xmin,xmaj*1000);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=1,line=1.5,"TP (\u03BCg L\u207B\u00B9)",cex=0.9)
mtext(side=2,line=2,"Density",cex=0.9)
legend("topright",c("Observed","Modelled"),
       pch=c(22,22),
       lty=c(NA,NA),lwd=c(1,1),
       col=c("dodgerblue1","indianred1"),
       pt.bg=adjustcolor(c("dodgerblue1","indianred1",0.25)),
       pt.cex=1.5,ncol=1,cex=0.9,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)

ylim.val=c(-0.15,0.15);by.y=0.1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=date.fun(dates);xmaj=seq(xlim.val[1],xlim.val[2],"4 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
plot(mean.TP~monCY,na.omit(cre.hydro.wq.mon[,c("monCY",vars)]),ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
segments(na.omit(cre.hydro.wq.mon[,c("monCY",vars)])$monCY,0,na.omit(cre.hydro.wq.mon[,c("monCY",vars)])$monCY,delta,lty=1,col="red")
points(na.omit(cre.hydro.wq.mon[,c("monCY",vars)])$monCY,delta,pch=21,bg="red",lwd=0.1)
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"));axis_fun(2,ymaj,ymin,format(round(ymaj,2)*1000));box(lwd=1)
mtext(side=2,line=2.5,"\u0394 TP (\u03BCg L\u207B\u00B9)",cex=0.9)
mtext(side=1,line=2,"Date (Month-Year)",cex=0.9)
mtext(side=3,line=-0.8,cex=0.5,"Observed - Predicted")
dev.off()

ylim.val=c(0.04,0.4);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
xlim.val=date.fun(dates);xmaj=seq(xlim.val[1],xlim.val[2],"4 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
#tiff(filename=paste0(plot.path,"S79_TPModel.tiff"),width=6,height=3.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
#png(filename=paste0(plot.path,"S79_TPModel.png"),width=6,height=3.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,3,0.25,0.5),oma=c(2,1.75,0.75,0.75));

plot(mean.TP~monCY,na.omit(cre.hydro.wq.mon[,c("monCY",vars)]),ylim=ylim.val,xlim=xlim.val,log="y",type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
#with(cre.hydro.wq.mon,pt_line(monCY,S79.TP.pred.mod,2,adjustcolor("indianred1",0.5),1,21,adjustcolor("indianred1",0.5),cex=1.25))
#shaded.range(cre.hydro.wq.mon$monCY,S79.TP.pred.mod[,2],S79.TP.pred.mod[,3],"indianred1",lty=1)
with(na.omit(cre.hydro.wq.mon[,c("monCY",vars)]),lines(monCY,S79.TP.pred.mod[,1],lty=1,col=adjustcolor("indianred1",0.5),lwd=2))
with(na.omit(cre.hydro.wq.mon[,c("monCY",vars)]),lines(monCY,S79.TP.pred.mod[,2],lty=2,col=adjustcolor("indianred1",0.5),lwd=1))
with(na.omit(cre.hydro.wq.mon[,c("monCY",vars)]),lines(monCY,S79.TP.pred.mod[,3],lty=2,col=adjustcolor("indianred1",0.5),lwd=1))
with(na.omit(cre.hydro.wq.mon[,c("monCY",vars)]),pt_line(monCY,mean.TP,2,adjustcolor("dodgerblue1",0.5),1,21,adjustcolor("dodgerblue1",0.5),cex=1))

axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"))
axis_fun(2,ymaj,ymin,ymaj*1000);box(lwd=1)
mtext(side=2,line=2.5,"Total Phosphorus (\u03BCg L\u207B\u00B9)")
mtext(side=1,line=2,"Date (Month-Year)")

legend("bottomleft",c("Observed","Modelled \u00B1 95% CI"),
       pch=c(NA,NA),
       lty=c(2,1),lwd=c(1,1),
       col=adjustcolor(c("dodgerblue1","indianred1",0.5)),
       pt.bg=adjustcolor(c("dodgerblue1","indianred1",0.5)),
       pt.cex=1.5,ncol=2,cex=0.9,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,text.col = "white")
legend("bottomleft",c("Observed","Modelled \u00B1 95% CI"),
       pch=c(21,NA),
       lty=c(NA,NA),lwd=c(0.1,1),
       col=adjustcolor(c("dodgerblue1","indianred1",0.5)),
       pt.bg=c(adjustcolor("dodgerblue1",0.5),NA),
       pt.cex=1.5,ncol=2,cex=0.9,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()


#tiff(filename=paste0(plot.path,"S79_TPParameters.tiff"),width=7,height=6,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(2,3.5,1.5,1.75),oma=c(2,1,0.5,0.5));
layout(matrix(1:6,3,2))
ylim.val=c(0.05,0.35);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
xlim.val=c(5.5,8.6);by.x=1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)

plot(mean.TP~grad,cre.hydro.wq.mon,log="y",ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(cre.hydro.wq.mon,points(grad,mean.TP,pch=21,bg="indianred1",lwd=0.1))
mod=lm(log(mean.TP)~grad,cre.hydro.wq.mon)
x.val=seq(min(cre.hydro.wq.mon$grad,na.rm=T),max(cre.hydro.wq.mon$grad,na.rm=T),length.out=100)
mod.pred=predict(mod,data.frame(grad=x.val),interval="confidence")
shaded.range(x.val,exp(mod.pred[,2]),exp(mod.pred[,3]),"grey",lty=1)
lines(x.val,exp(mod.pred[,1]),lty=2)
axis_fun(1,line=-0.5,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,ymaj*1000);box(lwd=1)
#mtext(side=1,line=1.75,"Surface Water Gradient (S235TW - S79HW; Ft NGVD29)",cex=0.8)
mtext(side=1,line=1.75,"Surface Water Gradient (Ft NGVD29)",cex=0.8)

xlim.val=c(0,25e4);by.x=5e4;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(mean.TP~C43Basin,cre.hydro.wq.mon,log="y",ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(cre.hydro.wq.mon,points(C43Basin,mean.TP,pch=21,bg="indianred1",lwd=0.1))
mod=lm(log(mean.TP)~C43Basin,cre.hydro.wq.mon)
x.val=seq(min(cre.hydro.wq.mon$C43Basin,na.rm=T),max(cre.hydro.wq.mon$C43Basin,na.rm=T),length.out=100)
mod.pred=predict(mod,data.frame(C43Basin=x.val),interval="confidence")
shaded.range(x.val,exp(mod.pred[,2]),exp(mod.pred[,3]),"grey",lty=1)
lines(x.val,exp(mod.pred[,1]),lty=2)
axis_fun(1,line=-0.5,xmaj,xmin,xmaj/1e4)
axis_fun(2,ymaj,ymin,ymaj*1000);box(lwd=1)
mtext(side=1,line=1.75,expression(paste("Q"["C43Basin"]," x10"^3,"; CFS")),cex=0.8)

xlim.val=c(0,20e4);by.x=5e4;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(mean.TP~S77,cre.hydro.wq.mon,log="y",ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(cre.hydro.wq.mon,points(S77,mean.TP,pch=21,bg="indianred1",lwd=0.1))
mod=lm(log(mean.TP)~S77,cre.hydro.wq.mon)
x.val=seq(min(cre.hydro.wq.mon$S77,na.rm=T),max(cre.hydro.wq.mon$S77,na.rm=T),length.out=100)
mod.pred=predict(mod,data.frame(S77=x.val),interval="confidence")
shaded.range(x.val,exp(mod.pred[,2]),exp(mod.pred[,3]),"grey",lty=1)
lines(x.val,exp(mod.pred[,1]),lty=2)
axis_fun(1,line=-0.5,xmaj,xmin,xmaj/1e4)
axis_fun(2,ymaj,ymin,ymaj*1000);box(lwd=1)
mtext(side=1,line=2,expression(paste("Q"["S77"]," x10"^3,"; CFS")),cex=0.8)

xlim.val=c(0,1);by.x=0.2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(mean.TP~basin.q.ratio,cre.hydro.wq.mon,log="y",ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(cre.hydro.wq.mon,points(basin.q.ratio,mean.TP,pch=21,bg="indianred1",lwd=0.1))
mod=lm(log(mean.TP)~basin.q.ratio,cre.hydro.wq.mon)
x.val=seq(min(cre.hydro.wq.mon$basin.q.ratio,na.rm=T),max(cre.hydro.wq.mon$basin.q.ratio,na.rm=T),length.out=100)
mod.pred=predict(mod,data.frame(basin.q.ratio=x.val),interval="confidence")
shaded.range(x.val,exp(mod.pred[,2]),exp(mod.pred[,3]),"grey",lty=1)
lines(x.val,exp(mod.pred[,1]),lty=2)
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj))
axis_fun(2,ymaj,ymin,ymaj*1000);box(lwd=1)
mtext(side=1,line=2,expression(paste("Q"["C43Basin"],": Q"["S79"])),cex=0.8)

tmp=ddply(cre.hydro.wq.mon,"hydro.season",summarise,mean.val=mean(mean.TP,na.rm=T),SE.val=SE(mean.TP))
ylim.val=c(0.05,0.35);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
#boxplot(mean.TP~hydro.season,cre.hydro.wq.mon)
x=barplot(tmp$mean.val,ylim=ylim.val,log="y",yaxt="n")
with(tmp,errorbars(x,mean.val,SE.val,col="black"))
axis_fun(2,ymaj,ymin,ymaj*1000)
axis_fun(1,x,x,c("Wet","Dry"));box(lwd=1)
mtext(side=1,line=2,"Hydrologic Season",cex=0.8)
mtext(side=2,line=-1,outer=T,"Total Phosphorus (\u03BCg L\u207B\u00B9)")

dev.off()

##
## TN Model
#vars=c("mean.TP","grad","C43Basin","S78","S79","S77","basin.q.ratio","hydro.season")
vars=c("mean.TN","grad","C43Basin","S77","basin.q.ratio","hydro.season")
S79.TN.mod=lm((1/mean.TN)~.,na.omit(cre.hydro.wq.mon[,vars]))
layout(matrix(1:4,2,2));plot(S79.TN.mod)
shapiro.test(residuals(S79.TN.mod));hist(residuals(S79.TN.mod))
vif(S79.TN.mod)
summary(S79.TN.mod)

#AIC model
S79.TN.mod.sw=stepAIC(S79.TN.mod,direction="both",trace=F)
S79.TN.mod.sw$anova
summary(S79.TN.mod.sw)
layout(matrix(1:4,2,2));plot(S79.TN.mod.sw)
shapiro.test(residuals(S79.TN.mod.sw));
vif(S79.TN.mod.sw)

#tiff(filename=paste0(plot.path,"S79_TNModel_diag.tiff"),width=6,height=5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(2,3,1,0.5),oma=c(2,1,0.25,0.75));
layout(matrix(1:4,2,2))

ylim.val=c(-0.3,0.3);by.y=0.4;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0.5,0.9);by.x=0.1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(S79.TN.mod.sw$fitted.values,S79.TN.mod.sw$residuals,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
abline(h=0)
with(S79.TN.mod.sw,points(fitted.values,residuals,pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1,cex=1.25))
fit.vals=S79.TN.mod.sw$fitted.values
res.vals=S79.TN.mod.sw$residuals
with(lowess(fit.vals,res.vals),lines(x,y,col="red",lwd=2))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,format(round(ymaj,2)));box(lwd=1)
mtext(side=2,line=2.5,"Residuals")
mtext(side=1,line=1.75,"Fitted Values")

ylim.val=c(-3,3);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(-3,3);by.x=1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
rstd=rstandard(S79.TN.mod.sw)
qq.x=qq.function(S79.TN.mod.sw$residuals)
plot(rstd~qq.x,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
points(qq.x,rstd,pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1,cex=1.25)
abline(0,1,lty=3)
axis_fun(1,line=-0.5,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2.5,"Standardized Residuals")
mtext(side=1,line=1.75,"Theoretical Quantiles")

ylim.val=c(0,2);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0.5,0.9);by.x=0.1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(sqrt(abs(rstd))~fit.vals,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
points(fit.vals,sqrt(abs(rstd)),pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1,cex=1.25)
with(lowess(fit.vals,sqrt(abs(rstd))),lines(x,y,col="red",lwd=2))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2,expression(sqrt("Standardized Residuals")))
mtext(side=1,line=1.75,"Fitted Values")

#plot(S79.TN.mod.sw,which=4)
#plot(cooks.distance(S79.TN.mod.sw))
#plot(S79.TN.mod.sw,which=5)

ylim.val=c(-3,3);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,0.35);by.x=0.1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
lev=hatvalues(S79.TN.mod.sw)
plot(rstd~lev,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
abline(h=0)
points(lev,rstd,pch=21,bg=adjustcolor("grey",0.5),col="grey",lwd=0.1,cex=1.25)
with(lowess(lev,rstd),lines(x,y,col="red",lwd=2))
inf=lm.influence(S79.TN.mod.sw)
hh=seq(min(range(inf)[1],range(inf)[2]/100),xlim.val[2]+(xlim.val[2]*0.5),length.out=101)
hh=hh[hh>0]
crit.cook=0.5
cl.h=sqrt(crit.cook*length(coef(S79.TN.mod.sw))*(1-hh)/hh)
lines(hh,cl.h,lty=2,col=2)
lines(hh,-cl.h,lty=2,col=2)
crit.cook=1
cl.h=sqrt(crit.cook*length(coef(S79.TN.mod.sw))*(1-hh)/hh)
lines(hh,cl.h,lty=2,col=2)
lines(hh,-cl.h,lty=2,col=2)
axis_fun(1,line=-0.5,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2,"Standardized Residuals")
mtext(side=1,line=1.75,"Leverage")
dev.off()

S79.TN.pred.mod=with(na.omit(cre.hydro.wq.mon[,vars]),1/(predict(S79.TN.mod.sw,data.frame(grad=grad,C43Basin=C43Basin),interval="confidence")))

delta=na.omit(cre.hydro.wq.mon[,vars])$mean.TN-S79.TN.pred.mod[,1]
mod.d=density(na.omit(as.numeric(S79.TN.pred.mod[,1])))
samp.d=density(na.omit(cre.hydro.wq.mon[,vars])$mean.TN)
mod.compare=data.frame(dat=c(S79.TN.pred.mod[,1],na.omit(cre.hydro.wq.mon[,vars])$mean.TN),group=c(rep("mod",length(S79.TN.pred.mod[,1])),rep("obs",length(na.omit(cre.hydro.wq.mon[,vars])$mean.TN))))
kruskal.test(dat~group,mod.compare)
boxplot(dat~group,mod.compare)

ylim.val=c(0.75,2);by.y=0.25;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
#tiff(filename=paste0(plot.path,"S79_TNModel_Obs.tiff"),width=3.5,height=3,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,3,1,0.5),oma=c(2,1,0.75,0.75));

boxplot(dat~group,mod.compare,ylim=ylim.val,axes=F,ylab=NA,xlab=NA,outline=F,col="grey")
axis_fun(1,1:2,1:2,c("Modelled","Observed"))
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"Total Nitrogen (mg L\u207B\u00B9)")
dev.off()


#tiff(filename=paste0(plot.path,"S79_TNModel_compare.tiff"),width=6,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(2,3,1,0.5),oma=c(2,1,0.75,0.75));
layout(matrix(c(1,2,3,3),2,2,byrow=T))

ylim.val=c(0.75,2.6);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=ylim.val;by.x=by.y;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(mean.TN~S79.TN.pred.mod[,1],na.omit(cre.hydro.wq.mon[,vars]),ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(na.omit(cre.hydro.wq.mon[,vars]),points(S79.TN.pred.mod[,1],mean.TN,pch=21,bg=adjustcolor("forestgreen",0.25),lwd=0.1,col="grey"))
abline(0,1,lty=2)
axis_fun(1,line=-0.5,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=1,line=1.5,"Modelled TN (mg L\u207B\u00B9)",cex=0.9)
mtext(side=2,line=2.5,"Observed TN (mg L\u207B\u00B9)",cex=0.9)

ylim.val=c(0,10);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0.75,2.6);by.x=0.5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(mod.d,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA,main=NA,yaxs="i")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
polygon(mod.d,col=adjustcolor("indianred1",0.5),lwd=1,border="indianred1")
polygon(samp.d,col=adjustcolor("dodgerblue1",0.5),lwd=1,border="dodgerblue1")
axis_fun(1,line=-0.5,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=1,line=1.5,"TN (mg L\u207B\u00B9)",cex=0.9)
mtext(side=2,line=2,"Density",cex=0.9)
legend("topright",c("Observed","Modelled"),
       pch=c(22,22),
       lty=c(NA,NA),lwd=c(1,1),
       col=c("dodgerblue1","indianred1"),
       pt.bg=adjustcolor(c("dodgerblue1","indianred1",0.25)),
       pt.cex=1.5,ncol=1,cex=0.9,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)

ylim.val=c(-0.5,1.25);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=date.fun(dates);xmaj=seq(xlim.val[1],xlim.val[2],"4 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
plot(mean.TN~monCY,na.omit(cre.hydro.wq.mon[,c("monCY",vars)]),ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
segments(na.omit(cre.hydro.wq.mon[,c("monCY",vars)])$monCY,0,na.omit(cre.hydro.wq.mon[,c("monCY",vars)])$monCY,delta,lty=1,col="red")
points(na.omit(cre.hydro.wq.mon[,c("monCY",vars)])$monCY,delta,pch=21,bg="red",lwd=0.1)
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"));axis_fun(2,ymaj,ymin,format(round(ymaj,2)));box(lwd=1)
mtext(side=2,line=2.5,"\u0394 TN (mg L\u207B\u00B9)",cex=0.9)
mtext(side=1,line=2,"Date (Month-Year)",cex=0.9)
mtext(side=3,line=-0.8,cex=0.5,"Observed - Predicted")
dev.off()

ylim.val=c(0.75,2.75);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=date.fun(dates);xmaj=seq(xlim.val[1],xlim.val[2],"4 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
#tiff(filename=paste0(plot.path,"S79_TNModel.tiff"),width=6,height=3.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
#png(filename=paste0(plot.path,"S79_TNModel.png"),width=6,height=3.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,2.5,0.25,0.5),oma=c(2,1.75,0.75,0.75));

plot(mean.TN~monCY,na.omit(cre.hydro.wq.mon[,c("monCY",vars)]),ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
#with(cre.hydro.wq.mon,pt_line(monCY,S79.TP.pred.mod,2,adjustcolor("indianred1",0.5),1,21,adjustcolor("indianred1",0.5),cex=1.25))
#shaded.range(cre.hydro.wq.mon$monCY,S79.TP.pred.mod[,2],S79.TP.pred.mod[,3],"indianred1",lty=1)
with(na.omit(cre.hydro.wq.mon[,c("monCY",vars)]),lines(monCY,S79.TN.pred.mod[,1],lty=1,col=adjustcolor("indianred1",0.5),lwd=2))
with(na.omit(cre.hydro.wq.mon[,c("monCY",vars)]),lines(monCY,S79.TN.pred.mod[,2],lty=2,col=adjustcolor("indianred1",0.5),lwd=1))
with(na.omit(cre.hydro.wq.mon[,c("monCY",vars)]),lines(monCY,S79.TN.pred.mod[,3],lty=2,col=adjustcolor("indianred1",0.5),lwd=1))
with(na.omit(cre.hydro.wq.mon[,c("monCY",vars)]),pt_line(monCY,mean.TN,2,adjustcolor("dodgerblue1",0.5),1,21,adjustcolor("dodgerblue1",0.5),cex=1))

axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"))
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2.75,"Total Nitrogen (mg L\u207B\u00B9)")
mtext(side=1,line=2,"Date (Month-Year)")

legend("bottomleft",c("Observed","Modelled \u00B1 95% CI"),
       pch=c(NA,NA),
       lty=c(2,1),lwd=c(1,1),
       col=adjustcolor(c("dodgerblue1","indianred1",0.5)),
       pt.bg=adjustcolor(c("dodgerblue1","indianred1",0.5)),
       pt.cex=1.5,ncol=2,cex=0.9,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,text.col = "white")
legend("bottomleft",c("Observed","Modelled \u00B1 95% CI"),
       pch=c(21,NA),
       lty=c(NA,NA),lwd=c(0.1,1),
       col=adjustcolor(c("dodgerblue1","indianred1",0.5)),
       pt.bg=c(adjustcolor("dodgerblue1",0.5),NA),
       pt.cex=1.5,ncol=2,cex=0.9,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

#tiff(filename=paste0(plot.path,"S79_TNParameters.tiff"),width=7,height=6,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(2,3.5,1.5,1.75),oma=c(2,1,0.5,0.5));
layout(matrix(1:6,3,2))
ylim.val=c(0.75,2.75);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(5.5,8.6);by.x=1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)

plot(mean.TN~grad,cre.hydro.wq.mon,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(cre.hydro.wq.mon,points(grad,mean.TN,pch=21,bg="indianred1",lwd=0.1))
mod=lm(log(mean.TN)~grad,cre.hydro.wq.mon)
x.val=seq(min(cre.hydro.wq.mon$grad,na.rm=T),max(cre.hydro.wq.mon$grad,na.rm=T),length.out=100)
mod.pred=predict(mod,data.frame(grad=x.val),interval="confidence")
shaded.range(x.val,exp(mod.pred[,2]),exp(mod.pred[,3]),"grey",lty=1)
lines(x.val,exp(mod.pred[,1]),lty=2)
axis_fun(1,line=-0.5,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
#mtext(side=1,line=1.75,"Surface Water Gradient (S235TW - S79HW; Ft NGVD29)",cex=0.8)
mtext(side=1,line=1.75,"Surface Water Gradient (Ft NGVD29)",cex=0.8)

xlim.val=c(0,25e4);by.x=5e4;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(mean.TN~C43Basin,cre.hydro.wq.mon,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(cre.hydro.wq.mon,points(C43Basin,mean.TN,pch=21,bg="indianred1",lwd=0.1))
mod=lm(log(mean.TN)~C43Basin,cre.hydro.wq.mon)
x.val=seq(min(cre.hydro.wq.mon$C43Basin,na.rm=T),max(cre.hydro.wq.mon$C43Basin,na.rm=T),length.out=100)
mod.pred=predict(mod,data.frame(C43Basin=x.val),interval="confidence")
shaded.range(x.val,exp(mod.pred[,2]),exp(mod.pred[,3]),"grey",lty=1)
lines(x.val,exp(mod.pred[,1]),lty=2)
axis_fun(1,line=-0.5,xmaj,xmin,xmaj/1e4)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=1,line=1.75,expression(paste("Q"["C43Basin"]," x10"^3,"; CFS")),cex=0.8)

xlim.val=c(0,20e4);by.x=5e4;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(mean.TN~S77,cre.hydro.wq.mon,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(cre.hydro.wq.mon,points(S77,mean.TN,pch=21,bg="indianred1",lwd=0.1))
mod=lm(log(mean.TN)~S77,cre.hydro.wq.mon)
x.val=seq(min(cre.hydro.wq.mon$S77,na.rm=T),max(cre.hydro.wq.mon$S77,na.rm=T),length.out=100)
mod.pred=predict(mod,data.frame(S77=x.val),interval="confidence")
shaded.range(x.val,exp(mod.pred[,2]),exp(mod.pred[,3]),"grey",lty=1)
lines(x.val,exp(mod.pred[,1]),lty=2)
axis_fun(1,line=-0.5,xmaj,xmin,xmaj/1e4)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=1,line=2,expression(paste("Q"["S77"]," x10"^3,"; CFS")),cex=0.8)

xlim.val=c(0,1);by.x=0.2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(mean.TN~basin.q.ratio,cre.hydro.wq.mon,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(cre.hydro.wq.mon,points(basin.q.ratio,mean.TN,pch=21,bg="indianred1",lwd=0.1))
mod=lm(log(mean.TN)~basin.q.ratio,cre.hydro.wq.mon)
x.val=seq(min(cre.hydro.wq.mon$basin.q.ratio,na.rm=T),max(cre.hydro.wq.mon$basin.q.ratio,na.rm=T),length.out=100)
mod.pred=predict(mod,data.frame(basin.q.ratio=x.val),interval="confidence")
shaded.range(x.val,exp(mod.pred[,2]),exp(mod.pred[,3]),"grey",lty=1)
lines(x.val,exp(mod.pred[,1]),lty=2)
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj))
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=1,line=2,expression(paste("Q"["C43Basin"],": Q"["S79"])),cex=0.8)

tmp=ddply(cre.hydro.wq.mon,"hydro.season",summarise,mean.val=mean(mean.TN,na.rm=T),SE.val=SE(mean.TN))
#ylim.val=c(0.05,0.35);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
#boxplot(mean.TN~hydro.season,cre.hydro.wq.mon)
x=barplot(tmp$mean.val,ylim=ylim.val,yaxt="n",xpd=F)
with(tmp,errorbars(x,mean.val,SE.val,col="black"))
axis_fun(2,ymaj,ymin,ymaj)
axis_fun(1,x,x,c("Wet","Dry"));box(lwd=1)
mtext(side=1,line=2,"Hydrologic Season",cex=0.8)
mtext(side=2,line=-0.5,outer=T,"Total Nitrogen (mg L\u207B\u00B9)")

dev.off()