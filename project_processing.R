#Load the required libraries
library(reshape)
library(ggplot2)
library(grid)
library(zoo)
library(tseries)
library(vars)
library(survival)

#Load all of the data, manipulate to make manageable
#Modis
modisData <- read.csv("C:/Arc_Projects/env_bio/randomptsmodis2.csv",head=T)
mdata = melt(modisData,id=c("ID","X","Y"))
mdata$variable = substr(mdata$variable,start=2,stop=9)
colnames(mdata)[4]="juldate"
mdata$juldate = as.numeric(mdata$juldate)

#Corn - was only able to find monthly data
cornprice2 <- read.csv("C:\\Users\\patrick.lawrence.MSU\\Documents\\models\\corn_03_13.csv",head=T)
cornprice2sub = subset(cornprice2,yrfrac>2006)

#Wheat
wheatprices <- read.csv("C:\\Users\\patrick.lawrence.MSU\\Documents\\models\\Great Falls Historical Pricing\\combinedpricedata3.csv",head=T)
wheatpricesub = subset(wheatprices,year>2006)
wheatpricesub$dns_13_lo = wheatpricesub$dns_13_lo*1.01604691 #convert imperial to metric tons
wheatpricesub$totalfrac = wheatpricesub$year+wheatpricesub$yearfrac%%1

#Soils
scandata <- read.csv("C:\\Users\\patrick.lawrence.MSU\\Documents\\school\\Environmental_biophysics\\project_data\\SCAN\\all_scan_data_edit_final.txt",head=T)

#Replace -99.9 values with NA for all rows and columns
scandata <- as.data.frame(lapply(scandata, function(x){replace(x, x ==-99.9,NA)}))

#Calculate ET

#Rn comes in a J/s, so the conversion is (input*3600)/1000000
#Rn = MJ/m^2/hour
#G = MJ/m^2/hour
#u2 = mean hourly wind speed m/2 at 2 meters - comes in at mph, so conv is .44704*input
#es comes in at kPa, as does ea, which are saturated vapor pressure and actual vapor pressure
#ST1,ST2 are soil temp and come in as degrees C - which is fine
#ht diff between sensors is 6" = .1524 meters
#assume thermal conductivity of a loam to be .5 W/m/C - change later to reflect site differences?

ETstandardized = function(Rn,Temp,elev,u2,es,ea,ST1,ST2){  #veg height of 0.12 m, measurements at 2 m
  G = -0.5*(ST2 - ST1)/(.1524)*3600/1000000 #for now - W/m^2  to MJ/m^2/hr
  if (is.na(Rn)){
    Rn=0
  }
  Rns = (1-.23)*Rn  #ASCE standard
#   if (Rn>0){
#     rs = 50 #s/m - hourly surface resistance, daytime - nighttime is 200
#   } else {
#     rs = 200 #s/m  - hourly nighttime
#   }
  sigma = 2.042E-10 #MJ/m^2/h - stefan boltzmann
  if (Rn > 10){ #Accounting for radiometers, set clear-sky radiation to zero in the night
    fcd = 1.35*(Rn/Rso) - .35
  } else {
    fcd = 1 #simplification
  }
  Rnl = sigma*fcd*(.34-.14*sqrt(ea))*(Temp+273)^4
  Rn = Rns - Rnl
  Cn = 37 #mm/hr
  Cd = .24 #mm/hr
  #P = Pressurecalc(elev)*.1 #in kPa
  Press = 101.3*((293-.0065*elev)/293)^5.26  #ASCE specifications
  lambda = 2.45 #MJ/kg
  Cp = .001005 #MJ/kg/C
  #gamma = Cp*P/0.622*lambda    #psychrometric constant - kPa/degrees C
  gamma = .000665*Press  #ASCE specifications
  #s = .04145*exp(.06088*Temp)  #slope of vapor pressure curve in kPa/degrees C
  s = (2503*exp((1727*Temp)/(Temp+237.3)))/((T+237.3)^2)
  ET = (.408*s*(Rn - G) + gamma*(Cn/(Temp+273))*u2*(es - ea))/(s+gamma*(1+Cd*u2))
  return (ET)                                                              
}

Pressurecalc = function(elev){
  Ps = 1013.25 #hPa at z = 0
  H = 7600 #m
  pz = Ps*exp(-elev/H)
  return(pz)
}

site.elevations = list("2119"=1297.229,"2121"=859.536,"2117"=1129.589,"2118"=982.98,"2019"=826.008,
                       "581"=875.0808,"2120"=693.1152,"2020"=588.264,"2021"=499.872)
ETtimeseries = vector()
for (i in 1:length(scandata[,1])){    #
  Rn = scandata$SRADV.H.1..watt.[i]*3600/1000000
  Temp = scandata$TAVG.H.1..degC.[i]
  site = scandata$Site.Id[i]
  elev = get(as.character(site),site.elevations)
  u2 = mean(c(scandata$WSPDV.H.1.124..mph.[i]*.44704,scandata$WSPDV.H.1.116..mph.[i]*.44704,
            scandata$WSPDV.H.1.123..mph.[i]*.44704,scandata$WSPDV.H.1.137..mph.[i]*.44704,
            scandata$WSPDV.H.1.128..mph.[i]*.44704,scandata$WSPDV.H.1.120..mph.[i]*.44704,
            scandata$WSPDV.H.1.118..mph.[i]*.44704,scandata$WSPDV.H.1.127..mph.[i]*.44704,
            scandata$WSPDV.H.1.131..mph.[i]*.44704,scandata$WSPDV.H.1.30..mph.[i]*.44704),na.rm=T)
  es = scandata$SVPV.H.1..kPa.[i]
  if (is.na(es)){
    es = .6108*exp((17.27*Temp)/(Temp+237.3))
  }
  ea = scandata$PVPV.H.1..kPa.[i]
  if (is.na(ea)){
    ea = mean(c(scandata$RHUM.H.1..pct.[i],scandata$RHUM.I.1..pct.[i]),na.rm=T)*es/100
  }
  ST1 = scandata$STO.I.1..2..degC.[i]
  ST2 = scandata$STO.I.1..8..degC.[i]
  ETout = ETstandardized(Rn,Temp,elev,u2,es,ea,ST1,ST2)
  ETtimeseries = c(ETtimeseries,ETout)
}
scandata$ET = ETtimeseries



#Convert SCAN hourly data into daily data, collapse over sites
#Create julian dates
scandata$juldate = NA
for (i in 1:length(scandata[,1])){
  endval = strptime(scandata$Date[i],"%m/%d/%Y")$yday+1
  if (as.numeric(endval)<10){
    juldate = paste(substr(strptime(scandata$Date[i],"%m/%d/%Y"),start=1,stop=4),endval,sep="00")
  } else if (as.numeric(endval) >= 10 & as.numeric(endval) < 100){
    juldate = paste(substr(strptime(scandata$Date[i],"%m/%d/%Y"),start=1,stop=4),endval,sep="0")
  } else {
    juldate = paste(substr(strptime(scandata$Date[i],"%m/%d/%Y"),start=1,stop=4),endval,sep="")
  }
  scandata$juldate[i] = juldate
}
#too slow - did it in excel then bringing it in
scjuld = read.csv("C:\\Users\\patrick.lawrence.MSU\\Documents\\school\\Environmental_biophysics\\project_data\\SCAN\\scandata_juldate.csv",head=T)
scandata$juldate=as.numeric(scjuld$juldate)

SMS2=SMS4=SMS8=Tavg=Tmax=Tmin=Prcp=ET=vector()
for (date in unique(scandata$juldate)){
  datesubset = subset(scandata,juldate==date)
  SMS2 = c(SMS2,mean(datesubset$SMS.I.1..2..pct....loam.,na.rm=T))
  SMS4 = c(SMS4,mean(datesubset$SMS.I.1..4..pct....loam.,na.rm=T))
  SMS8 = c(SMS8,mean(datesubset$SMS.I.1..8..pct....loam.,na.rm=T))
  Tavg = c(Tavg,mean(datesubset$TAVG.H.1..degC.,na.rm=T))
  Tmax = c(Tmax,mean(datesubset$TMAX.H.1..degC.,na.rm=T))
  Tmin = c(Tmin,mean(datesubset$TMIN.H.1..degC.,na.rm=T))
  prcptally = vector()
  ETtally = vector()
  for (site in unique(datesubset$Site.Id)){
    sub = subset(datesubset,Site.Id==site)
    prcpsum = sum(sub$PREC.I.1..in.,na.rm=T)
    ETsum = sum(sub$ET,na.rm=T)
    prcptally = c(prcpsum,prcptally)
    ETtally = c(ETsum,ETtally)
  }
  Prcp = c(Prcp,mean(prcptally)) 
  ET = c(ET,mean(ETtally))
}
scancondense = as.data.frame(cbind(unique(scandata$juldate),SMS2,SMS4,SMS8,Tavg,Tmax,Tmin,Prcp,ET))
colnames(scancondense)[1]="juldate"

#Make more friendly column names - that ggplot can assimilate
colnames(scandata)[12:14]=c("SMS2","SMS4","SMS8")
colnames(scandata)[10]="Tavg"
colnames(scandata)[4]="Prcp"
colnames(scandata)[8:9]=c("Tmax","Tmin")

#Generate initial plots of the data
plotlist = list()
for (year in 2007:2012){ #:2012
  yearscansubset = subset(scandata,substr(juldate,start=1,stop=4)==year)
  yearscancondense = subset(scancondense,substr(juldate,start=1,stop=4)==year)
  yearmodis = subset(mdata,substr(juldate,start=1,stop=4)==year)
  yearmodiscondense = subset(meanmodis,substr(juldate,start=1,stop=4)==year)
  yearcorn = subset(cornprice2sub,substr(Juldate,start=1,stop=4)==year)
  yearwheat = wheatpricesub[wheatpricesub$year==year,]
  assign(paste("SMS2",year,sep=""),ggplot(yearscansubset,aes(x=juldate,y=SMS2))+geom_line()+
    geom_line(aes(juldate,SMS2,colour="red",size=3),data=yearscancondense)+xlab("Julian Day")+ylab("Soil moisture percent (@2 in.)")+ theme(legend.position="none"))
  assign(paste("SMS4",year,sep=""),ggplot(yearscansubset,aes(x=juldate,y=SMS4))+geom_line()+
    geom_line(aes(juldate,SMS4,colour="red",size=3),data=yearscancondense)+xlab("Julian Day")+ylab("Soil moisture percent (@4 in.)")+ theme(legend.position="none"))
  assign(paste("SMS8",year,sep=""),ggplot(yearscansubset,aes(x=juldate,y=SMS8))+geom_line()+
    geom_line(aes(juldate,SMS8,colour="red",size=3),data=yearscancondense)+xlab("Julian Day")+ylab("Soil moisture percent (@8 in.)")+ theme(legend.position="none"))
  assign(paste("Tavg",year,sep=""),ggplot(yearscansubset,aes(x=juldate,y=Tavg))+geom_line()+
    geom_line(aes(juldate,Tavg,colour="red",size=3),data=yearscancondense)+xlab("Julian Day")+ylab("Avg. Daily Temp (deg C)")+ theme(legend.position="none"))
  assign(paste("Tmax",year,sep=""),ggplot(yearscansubset,aes(x=juldate,y=Tmax))+geom_line()+
    geom_line(aes(juldate,Tmax,colour="red",size=3),data=yearscancondense)+xlab("Julian Day")+ylab("Max Temp (deg C)")+ theme(legend.position="none"))
  assign(paste("Tmin",year,sep=""),ggplot(yearscansubset,aes(x=juldate,y=Tmin))+geom_line()+
    geom_line(aes(juldate,Tmin,colour="red",size=3),data=yearscancondense)+xlab("Julian Day")+ylab("Min Temp (deg C)")+ theme(legend.position="none"))
  assign(paste("Prcp",year,sep=""),ggplot(yearscansubset,aes(x=juldate,y=Prcp))+geom_line()+
    geom_line(aes(juldate,Prcp,colour="red",size=3),data=yearscancondense)+xlab("Julian Day")+ylab("Avg Daily Prcp (hundredths of an inch)")+ theme(legend.position="none"))
  assign(paste("ET",year,sep=""),ggplot(yearscansubset,aes(x=juldate,y=ET))+geom_line()+
    geom_line(aes(juldate,ET,colour="red",size=3),data=yearscancondense)+xlab("Julian Day")+ylab("Avg Daily ET (mm/hr)")+ theme(legend.position="none"))
  assign(paste("modis",year,sep=""),ggplot(yearmodis,aes(x=juldate,y=value))+geom_line(aes(group=ID))+
    geom_line(aes(juldate,EVI,colour="red",size=3),data=yearmodiscondense)+xlab("Julian Day")+ylab("Avg Daily EVI")+ theme(legend.position="none"))
  assign(paste("corn",year,sep=""),ggplot(yearcorn,aes(x=Juldate,y=Price))+geom_line(aes(size=2))+xlab("Julian Day")+ylab("Monthly Corn Price ($/metric Ton)")+ theme(legend.position="none"))
  assign(paste("wheat",year,sep=""),ggplot(yearwheat,aes(x=yjday,y=dns_13_lo))+geom_line(aes(size=2))+xlab("Julian Day")+ylab("Daily Wheat Price ($/metric Ton)")+ theme(legend.position="none"))

  jpeg(filename = paste(year,"inputdata.jpg",sep=""),width=1000,height=1500,units="px",bg="white")
  multiplot(get(paste("SMS2",year,sep="")),get(paste("SMS4",year,sep="")),get(paste("SMS8",year,sep="")),get(paste("Tavg",year,sep="")),
            get(paste("Tmax",year,sep="")),get(paste("modis",year,sep="")),get(paste("Prcp",year,sep="")),get(paste("ET",year,sep="")),
            get(paste("corn",year,sep="")),get(paste("wheat",year,sep="")),cols=2)
  dev.off()
  
}

cornplot = ggplot(cornprice2,aes(yrfrac,Price))+geom_line()



p1 = ggplot(mdata,aes(variable,value,colour=ID))+geom_line(aes(group=ID)) 
p2 = ggplot(cornprices,aes(Yjday,Value))+geom_line(aes(group=Year))
p3 = ggplot(wheatpricesub,aes(totalfrac,hrw_12_lo))+geom_line()


p5 = ggplot(scandata,aes(yearfrac,SMSI1_2))+geom_line(aes(group=Site.Id))

multiplot(p1,p4,p3,p5,cols=1)+
grid.arrange(p1,p2)


#Average the modis data on a daily basis
meanmodis=data.frame(matrix(ncol=0,nrow=108))
i = 1
julout = vector()
mout = vector()
for (date in unique(mdata$juldate)){
  msubset = subset(mdata,juldate==date)
  julout = c(julout,date)
  mout = c(mout,mean(msubset$value,na.rm=T))
  #meanmodis$juldate[i] = date
  #meanmodis$EVI[i] = mean(msubset$value,na.rm=T)
  i = i+1
}
meanmodis$juldate=as.numeric(julout)
meanmodis$EVI=mout


#Average the modis data on a daily basis
meanmodis = data.frame(nrow=(ncol(scandata)-2),ncol=2)
for (i in 3:ncol(scandata)){
  meanmodis[i,1] = scandata[1,i]
  meanmodis[i,2] = mean(scandata[,i])
}




#Match up all of the disparate time series, interpolate if necessary
#first divide by years
for (year in 2007:2012){
  scansub = subset(scancondense,substr(juldate,start=1,stop=4)==year)
  scansub = subset(scansub,as.numeric(substr(juldate,start=5,stop=7))>59 & as.numeric(substr(juldate,start=5,stop=7))<306)
  modissub = subset(meanmodis,substr(juldate,start=1,stop=4)==year)
  modissub = subset(modissub,as.numeric(substr(juldate,start=5,stop=7))>59 & as.numeric(substr(juldate,start=5,stop=7))<306)
  cornsub = subset(cornprice2sub,substr(Juldate,start=1,stop=4)==year)
  cornsub = subset(cornsub,as.numeric(substr(Juldate,start=5,stop=7))>59 & as.numeric(substr(Juldate,start=5,stop=7))<306)
  wheatsub = wheatpricesub[wheatpricesub$year==year,]
  wheatsub = wheatsub[as.numeric(substr(wheatsub$yjday,start=5,stop=7))>59 & as.numeric(substr(wheatsub$yjday,start=5,stop=7))<306,]
  zs1 = zoo(scansub$SMS2,scansub$juldate)
  zs2 = zoo(scansub$SMS4,scansub$juldate)
  zs3 = zoo(scansub$SMS8,scansub$juldate)
  ztmean = zoo(scansub$Tavg,scansub$juldate)
  ztmax = zoo(scansub$Tmax,scansub$juldate)
  zp = zoo(scansub$Prcp,scansub$juldate)
  zet = zoo(scansub$ET,scansub$juldate)
  zmodis = zoo(modissub$EVI,modissub$juldate)
  zcorn = zoo(cornsub$Price,as.numeric(cornsub$Juldate))
  zwheat = zoo(wheatsub$dns_13_lo,as.numeric(unique(wheatsub$yjday)))    
  zout = merge(zs1,zs2,zs3,ztmean,ztmax,zp,zet,zmodis,zcorn,zwheat)
  zname = paste("z",year,sep="")
  assign(zname,zout)
}

library(zoo)
zs1 = zoo(scancondense$SMS2,scancondense$juldate)
zs2 = zoo(scancondense$SMS4,scancondense$juldate)
zs3 = zoo(scancondense$SMS8,scancondense$juldate)
ztmean = zoo(scancondense$Tavg,scancondense$juldate)
ztmax = zoo(scancondense$Tmax,scancondense$juldate)
zp = zoo(scancondense$Prcp,scancondense$juldate)
zet = zoo(scancondense$ET,scancondense$juldate)
zmodis = zoo(meanmodis$EVI,meanmodis$juldate)
zcorn = zoo(cornprice2sub$Price,as.numeric(cornprice2sub$Juldate))
zwheat = zoo(wheatpricesub$dns_13_lo,as.numeric(unique(wheatpricesub$yjday)))

#Merge
z = merge(zs1,zs2,zs3,ztmean,ztmax,zp,zet,zmodis,zcorn,zwheat)


#Interpolate the corn data from monthly to daily, modis to daily
#too lazy to not repeat the years!
z2007$zcorn = na.approx(z2007$zcorn)
z2007$zwheat = na.approx(z2007$zwheat)
z2007$zmodis = na.approx(z2007$zmodis,na.rm=F)
z2007$zet = na.approx(z2007$zet,na.rm=F)
z2007$zp = na.approx(z2007$zp,na.rm=F)
z2007$ztmax = na.approx(z2007$ztmax,na.rm=F)
z2007$ztmean = na.approx(z2007$ztmean,na.rm=F)
z2007$zs3 = na.approx(z2007$zs3,na.rm=F)
z2007$zs2 = na.approx(z2007$zs2,na.rm=F)
z2007$zs1 = na.approx(z2007$zs1,na.rm=F)
z2008$zcorn = na.approx(z2008$zcorn,na.rm=F)
z2008$zwheat = na.approx(z2008$zwheat,na.rm=F)
z2008$zmodis = na.approx(z2008$zmodis,na.rm=F)
z2008$zet = na.approx(z2008$zet,na.rm=F)
z2008$zp = na.approx(z2008$zp,na.rm=F)
z2008$ztmax = na.approx(z2008$ztmax,na.rm=F)
z2008$ztmean = na.approx(z2008$ztmean,na.rm=F)
z2008$zs3 = na.approx(z2008$zs3,na.rm=F)
z2008$zs2 = na.approx(z2008$zs2,na.rm=F)
z2008$zs1 = na.approx(z2008$zs1,na.rm=F)
z2009$zcorn = na.approx(z2009$zcorn,na.rm=F)
z2009$zwheat = na.approx(z2009$zwheat,na.rm=F)
z2009$zmodis = na.approx(z2009$zmodis,na.rm=F)
z2009$zet = na.approx(z2009$zet,na.rm=F)
z2009$zp = na.approx(z2009$zp,na.rm=F)
z2009$ztmax = na.approx(z2009$ztmax,na.rm=F)
z2009$ztmean = na.approx(z2009$ztmean,na.rm=F)
z2009$zs3 = na.approx(z2009$zs3,na.rm=F)
z2009$zs2 = na.approx(z2009$zs2,na.rm=F)
z2009$zs1 = na.approx(z2009$zs1,na.rm=F)
z2010$zcorn = na.approx(z2010$zcorn,na.rm=F)
z2010$zwheat = na.approx(z2010$zwheat,na.rm=F)
z2010$zmodis = na.approx(z2010$zmodis,na.rm=F)
z2010$zet = na.approx(z2010$zet,na.rm=F)
z2010$zp = na.approx(z2010$zp,na.rm=F)
z2010$ztmax = na.approx(z2010$ztmax,na.rm=F)
z2010$ztmean = na.approx(z2010$ztmean,na.rm=F)
z2010$zs3 = na.approx(z2010$zs3,na.rm=F)
z2010$zs2 = na.approx(z2010$zs2,na.rm=F)
z2010$zs1 = na.approx(z2010$zs1,na.rm=F)
z2011$zcorn = na.approx(z2011$zcorn,na.rm=F)
z2011$zwheat = na.approx(z2011$zwheat,na.rm=F)
z2011$zmodis = na.approx(z2011$zmodis,na.rm=F)
z2011$zet = na.approx(z2011$zet,na.rm=F)
z2011$zp = na.approx(z2011$zp,na.rm=F)
z2011$ztmax = na.approx(z2011$ztmax,na.rm=F)
z2011$ztmean = na.approx(z2011$ztmean,na.rm=F)
z2011$zs3 = na.approx(z2011$zs3,na.rm=F)
z2011$zs2 = na.approx(z2011$zs2,na.rm=F)
z2011$zs1 = na.approx(z2011$zs1,na.rm=F)
z2012$zcorn = na.approx(z2012$zcorn,na.rm=F)
z2012$zwheat = na.approx(z2012$zwheat,na.rm=F)
z2012$zmodis = na.approx(z2012$zmodis,na.rm=F)
z2012$zet = na.approx(z2012$zet,na.rm=F)
z2012$zp = na.approx(z2012$zp,na.rm=F)
z2012$ztmax = na.approx(z2012$ztmax,na.rm=F)
z2012$ztmean = na.approx(z2012$ztmean,na.rm=F)
z2012$zs3 = na.approx(z2012$zs3,na.rm=F)
z2012$zs2 = na.approx(z2012$zs2,na.rm=F)
z2012$zs1 = na.approx(z2012$zs1,na.rm=F)

#Get rid of leading and trailing NAs
z2007 = na.omit(z2007)
z2008 = na.omit(z2008)
z2009 = na.omit(z2009)
z2010 = na.omit(z2010)
z2011 = na.omit(z2011)
z2012 = na.omit(z2012)

#Test for stationarity of the time series
library(tseries)

adftable = data.frame(yearvariable=I(character()),statistic=numeric(),lagparam=numeric(),pval=numeric(),alternative=I(character()))
i = 1
varlist = c("zs1","zs2","zs3","ztmean","ztmax","zp","zet","zmodis","zwheat","zcorn")
for (year in 2007:2012){
  for (variable in varlist){
    outtest = adf.test(as.data.frame(get(paste("z",year,sep="")))[,variable])
    adftable[i,"yearvariable"]=I(paste(year,variable,sep=""))
    adftable[i,"statistic"]=outtest$statistic
    adftable[i,"lagparam"]=outtest$parameter
    adftable[i,"pval"]=outtest$p.value
    adftable[i,"alternative"]=outtest$alternative
    i=i+1
  }
}

#Need to first difference practically everything
dly <- diff(log(y))

#convert to tseries objects
namelist = vector()
for (year in 2007:2012){
  for (variable in varlist){
    assign(paste(variable,year,sep=""),ts(as.vector(get(paste("z",year,sep=""))[,variable])))
    namelist = c(namelist,get(paste(variable,year,sep="")))
  }
}

#bind them together
y2007=ts.union(zs12007,zs22007,zs32007,ztmean2007,ztmax2007,zp2007,zet2007,zmodis2007,zcorn2007,zwheat2007)

#Run first differences for all of the variables - they're all non-stationary!
y2007diff=ts.union(diff(zs12007),diff(zs22007),diff(zs32007),diff(ztmean2007),
               diff(ztmax2007),diff(zp2007),diff(zet2007),diff(zmodis2007),diff(zcorn2007),diff(zwheat2007))

y2007simple=ts.union(zp2007,zwheat2007,zcorn2007)

library(vars)
#Combine into one mv variable

VARselect(y2007,lag.max=90, type="none") #Provides normalized IC measures for VAR(1)-VAR(15)

var_15 <- VAR(y2007, p=15,lag.max=15,K=9, type="none") #Estimate VAR(8)
var_15
summary(var_15) #View Results
causality(var_15, cause=c("zmodis2007","zet2007","zs22007","zs12007","ztmean2007","ztmax2007","zp2007","zcorn2007",
                          "zs32007")) #Test H0: All other variables GC wheat prices

var_8_pred <- predict(var_8, n.ahead=30, ci=.95, dumvar=NULL) #Forecast 30 periods ahead
pred <- var_8_pred$fcst$wheatts[,1]


tmp = VARselect(y2007simple,lag.max=90, type="none")
var_49 <- VAR(y2007simple,p=49,type="none")
summary(var_49)
causality(var_49, cause=c("zp2007","zcorn2007"))



VARselect(y2007diff,lag.max=90, type="none") #Provides normalized IC measures for VAR(1)-VAR(15)

causaltable = data.frame(variable=I(character()),lag=numeric(),fstat=numeric(),df1=numeric(),
                         df2=numeric(),pval=numeric())
#Start over - run it variable by variable
i=1
for (year in 2007:2012){
  for (variable in varlist){
    var1 = get(paste(variable,year,sep=""))
    var1name = paste(variable,year,sep="")
    wheatvar = get(paste("zwheat",year,sep=""))
    outts = ts.union(diff(var1),diff(wheatvar))
    optimallag = VARselect(outts,lag.max=90,type="none")$selection[1]   #Get the optimal VAR lag
    optimalvar = VAR(outts,p=optimallag,type="none")
    causalresult = causality(optimalvar,cause=c("diff.var1."))
    causaltable[i,"variable"]= as.character(paste(variable,year,sep=""))
    causaltable[i,"lag"] = optimallag
    causaltable[i,"fstat"] = causalresult$Granger$statistic[1]
    causaltable[i,"df1"] = causalresult$Granger$parameter[1]
    causaltable[i,"df2"] = causalresult$Granger$parameter[2]
    causaltable[i,"pval"] = causalresult$Granger$p.value[1]
    i = i + 1
  }
}

#Grab variables with p-vals less than 0.05, indicating causality
causesignificant = subset(causaltable,pval<0.05)

#predict ahead 90 days with the significant variables - the time frame of the farmer
for (variable in causesignificant$variable){
  var1 = get(variable)
  year = substr(variable,start=nchar(variable)-3,stop=nchar(variable))
  wheatvar = get(paste("zwheat",year,sep=""))
  outts = ts.union(ts.union(diff(var1),diff(wheatvar)))
  optimallag = VARselect(outts,lag.max=90,type="none")$selection[1]   #Get the optimal VAR lag
  optimalvar = VAR(outts,p=optimallag,type="none")
  varpredict = predict(optimalvar,n.ahead=90,ci=.95,dumvar=NULL)
  forecast = varpredict$fcst$diff.wheatvar.[,1]
  lowerci = varpredict$fcst$diff.wheatvar.[,2]
  upperci = varpredict$fcst$diff.wheatvar.[,3]
  plotfcst = data.frame(forecast,lowerci,upperci)
  plottitlenames = list("zp"="precipitation","zmodis"="modis EVI","ztmean"="Avg Temp",
                        "ztmax"="Max Temp","zet"="ET","zcorn"="Corn Price","zs3"="Soil Moisture at 8 in.")
  plotout = ggplot(plotfcst,aes(x=as.numeric(rownames(plotfcst)),y=forecast))+geom_line()+
    geom_line(aes(x=as.numeric(rownames(plotfcst)),y=lowerci),colour="blue",data=plotfcst)+
    geom_line(aes(x=as.numeric(rownames(plotfcst)),y=upperci),colour="red",data=plotfcst)+
    xlab("Forecast Day")+ylab("Forecasted first differenced wheat price, $/MT")+ theme(legend.position="none")+
    ggtitle(paste("FD Wheat Price Response to",get(gsub(year,"",variable,fixed=TRUE),plottitlenames),"in",year))
  ggsave(filename=paste(variable,".jpg",sep=""),plotout,width=9,height=5)
  
}







# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}