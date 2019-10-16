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
cre.hydro.wq.mon$hydro.seaon=FL.Hydroseason(cre.hydro.wq.mon$monCY)
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
#TP
#Autocorrelation
acf(wq.dat.xtab.mon$mean.TP)
pacf(wq.dat.xtab.mon$mean.TP)
#lag 1
cor.test(wq.dat.xtab.mon$mean.TP[-nrow(wq.dat.xtab.mon)],wq.dat.xtab.mon$mean.TP[-1],method="spearman")
#lag 2
cor.test(wq.dat.xtab.mon$mean.TP[-((nrow(wq.dat.xtab.mon)-1):nrow(wq.dat.xtab.mon))],wq.dat.xtab.mon$mean.TP[-(1:2)],method="spearman")
acf(wq.dat.xtab.mon$mean.TP,plot=F)
acf(wq.dat.xtab.mon$mean.TP)

#Trend
with(wq.dat.xtab.mon,cor.test(mean.TP,decWY,method="kendall"))
#with(subset(wq.dat.xtab.mon,is.na(mean.TP)==F),zyp.sen(mean.TP~decWY))
with(subset(wq.dat.xtab.mon,is.na(mean.TP)==F),mblm(mean.TP~decWY,repeated = T))
plot(mean.TP~decWY,wq.dat.xtab.mon,las=1,type="b")
abline(with(subset(wq.dat.xtab.mon,is.na(mean.TP)==F),mblm(mean.TP~decWY,repeated = T)),col="red")

#seasonal kendall
with(wq.dat.xtab.mon,rkt(decimal_date(monCY),mean.TP,month))
#par(family="serif",mar=c(1,1,1,1),oma=c(0.5,0.5,0.5,0.5));
#layout(matrix(1:12,3,4))
#for(i in 1:12){
#  hist(subset(wq.dat.xtab.mon,month==i)$mean.TP,xlim=c(0.06,0.3))
#  mtext(side=3,month.name[i])
#}

#TN
#Autocorrelation
acf(subset(wq.dat.xtab.mon,is.na(mean.TN)==F)$mean.TN)
pacf(subset(wq.dat.xtab.mon,is.na(mean.TN)==F)$mean.TN)
#lag 1
cor.test(wq.dat.xtab.mon$mean.TN[-nrow(wq.dat.xtab.mon)],wq.dat.xtab.mon$mean.TN[-1],method="spearman")
#lag 2
cor.test(wq.dat.xtab.mon$mean.TN[-((nrow(wq.dat.xtab.mon)-1):nrow(wq.dat.xtab.mon))],wq.dat.xtab.mon$mean.TN[-(1:2)],method="spearman")
acf(subset(wq.dat.xtab.mon,is.na(mean.TN)==F)$mean.TN,plot=F)
acf(subset(wq.dat.xtab.mon,is.na(mean.TN)==F)$mean.TN)

#Trend
with(wq.dat.xtab.mon,cor.test(mean.TN,decWY,method="kendall"))
with(subset(wq.dat.xtab.mon,is.na(mean.TN)==F),mblm(mean.TN~decWY,repeated = T))
plot(mean.TN~decWY,wq.dat.xtab.mon,las=1,type="b")
abline(with(subset(wq.dat.xtab.mon,is.na(mean.TN)==F),mblm(mean.TN~decWY,repeated = T)),col="red")

#seasonal kendall
with(wq.dat.xtab.mon,rkt(decimal_date(monCY),mean.TN,month))

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


###
## TP Model