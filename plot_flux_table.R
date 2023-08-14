# This is a script for showing the result from the flux.table both for the daily and diurnal evolution over the selected average windows (box_days + ini_day) 
# date: 2023-04-27 
# author: yi-ying chen 

# load library 
library("lubridate")

# set flux table name
#table.name= c("WuF_2023-03-01 00:00:00_60_flux.table.txt")
#table.name= c("WuF_S1-S9_2019-07-01_2022-12-31_flux.table.txt")
table.name=c("WuF_S10_2023-02-15_2023-05-14_all.table.txt")

# read in the flux.table
flux.table <-  read.csv(table.name) 


flux.table$date.time <- as.POSIXct(flux.table$date.time,format='%Y-%m-%d %H:%M:%S')
flux.table$date.week.hh <- (week(flux.table$date.time)*24 + hour(flux.table$date.time))   
flux.table$flux.ch4[abs(flux.table$flux.ch4) > 0.3] <- 0
flux.table$flux.co2[abs(flux.table$flux.co2) > 50] <- 0
flux.table$flux.le[abs(flux.table$flux.le) > 800] <- 0
flux.table$flux.sh[abs(flux.table$flux.sh) > 800] <- 0
flux.table$flux.ch4[is.na(flux.table$flux.ch4)] <- 0
flux.table$flux.co2[is.na(flux.table$flux.co2)] <- 0


#paddy rice 
#crop_season = c("S2","S3","S8","S9"
crop_season = c("SXX")
#subset the table for the selcted croping season
if (crop_season == "S2") flux.table <- subset(flux.table, ((flux.table$date.time >= "2019-08-15") & (flux.table$date.time <= "2019-11-30" ))) 
if (crop_season == "S3") flux.table <- subset(flux.table, ((flux.table$date.time >= "2020-02-01") & (flux.table$date.time <= "2020-06-30" ))) 

if (crop_season == "S8") flux.table <- subset(flux.table, ((flux.table$date.time >= "2022-02-01") & (flux.table$date.time <= "2022-06-30" ))) 
if (crop_season == "S9") flux.table <- subset(flux.table, ((flux.table$date.time >= "2022-07-01") & (flux.table$date.time <= "2022-10-31" ))) 




 


#set days for box plot 
box_days <- 14
ini_day <-  40

in_date <- flux.table$date.time[1 + ini_day*48 ]  
ed_date <- in_date + days(box_days) 
flux.table$date.hh  <- hour(flux.table$date.time)   
flux.table$date.dd <- yday(flux.table$date.time) 

ld_go <- FALSE

if(ld_go) {
# png(file=paste(substr(in_date,start=1,stop=10),"_",substr(ed_date,start=1,stop=10),"_flux.plot.png",sep=""), width=1024, height=768, res=128)

 par(mfrow=c(4,2),mai=c(0.75,0.75,0.2,0.2), mar=c(2,4.5,1.5,1) )

 #CO2 flux
 ylab.txt <- expression(paste("CO2 flux, (", mu, "mol/", m^2, s,")"))
 plot( x=flux.table$date.time, y=flux.table$flux.co2, type="p",ylab=ylab.txt, xlab="Observation Period [Month date]", xaxt="n",
      col="gray",cex=0.8, pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1) , ylim=c(-3E1,2E1));grid()
 #add gray area 
  rect(xleft=in_date, xright=ed_date, ybottom=par("usr")[3], ytop=par("usr")[4], density=NA, col=adjustcolor("gray", alpha = 0.3)) 
 #rect(xleft=in_date, xright=ed_date, ybottom=par("usr")[3], ytop=par("usr")[4], density=NA, col="gray")
 #set xaex label
  x<- flux.table$date.time ;  at <- seq(min(x), max(x), "month")
  axis(side=1, at=at, labels=format(at, "%b"), cex.axis=1.2) 
  mtext(paste("shade area: ",substr(in_date,start=1,stop=10)," to ",substr(ed_date,start=1,stop=10)),side=3,line=0.5,
             at=par("usr")[1]+0.5*diff(par("usr")[1:2]),
             cex=0.8, col="gray") 
#subset the flux.table 
 plot.table <- subset ( flux.table, (flux.table$date.time >= in_date & flux.table$date.time <= ed_date) )
 boxplot(plot.table$flux.co2 ~ plot.table$date.hh , ylab=ylab.txt, xlab="Local hour [00:00 to 24:00]", col="gray", ylim=c(-3E1,2E1));grid()

 #CH4 Flux
 ylab.txt <- expression(paste("CH4 flux, (", mu, "mol/", m^2, s,")")) 
 plot(x=flux.table$date.time, y=flux.table$flux.ch4, type="p",ylab=ylab.txt, xlab="Observation Period [Month date]",xaxt="n", 
      col="orange",cex=0.8, pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1) , ylim=c(-.5E-1,1.5E-1)) ;grid()
 #add gray area 
 rect(xleft=in_date, xright=ed_date, ybottom=par("usr")[3], ytop=par("usr")[4], density=NA, col=adjustcolor("gray", alpha = 0.3)) 
 axis(side=1, at=at, labels=format(at, "%b"), cex.axis=1.2) 
 boxplot(plot.table$flux.ch4~plot.table$date.hh, ylab=ylab.txt, xlab="Local hour [00:00 to 24:00]", col="orange", ylim=c(-.5E-1,1.5E-1));grid()

 #Latent heat flux
 ylab.txt <- expression(paste("Latent heat, (", "W/", m^2,")"))
 plot(x=flux.table$date.time, y=flux.table$flux.le, type="p",ylab=ylab.txt, xlab="Observation Period [Month date]",xaxt="n", 
      col="blue",cex=0.8, pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1) , ylim=c(-5E1,5E2));grid()
 axis(side=1, at=at, labels=format(at, "%b"), cex.axis=1.2) 
 #add gray area 
 rect(xleft=in_date, xright=ed_date, ybottom=par("usr")[3], ytop=par("usr")[4], density=NA, col=adjustcolor("gray", alpha = 0.3))
 boxplot(plot.table$flux.le~plot.table$date.hh, ylab=ylab.txt, xlab="Local hour [00:00 to 24:00]", col="blue", ylim=c(-5E1,5E2));grid()
 
 #Sensible heat flux
 ylab.txt <- expression(paste("Sensible heat, (", "W/", m^2,")")) 
 plot(x=flux.table$date.time, y=flux.table$flux.sh, type="p",ylab=ylab.txt,xlab="Observation Period [Month date]", xaxt="n",
      col="red",cex=0.8, pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1) , ylim=c(-5E1,5E2));grid()
 axis(side=1, at=at, labels=format(at, "%b"), cex.axis=1.2)  
 #add gray area 
 rect(xleft=in_date, xright=ed_date, ybottom=par("usr")[3], ytop=par("usr")[4], density=NA, col=adjustcolor("gray", alpha = 0.3))
 boxplot(plot.table$flux.sh~plot.table$date.hh, ylab=ylab.txt, xlab="Local hour [00:00 to 24:00]", col="red", ylim=c(-5E1,5E2));grid()
 
# dev.off()
}

#dev.new()

ld_go <- TRUE

if (ld_go) {
  #CO2 flux 
  png(file=paste(substr(table.name,start=5,stop=14),"_CO2_flux.plot.png",sep=""), width=1200, height=850, res=128)
  par(mfrow=c(2,1),mai=c(0.75,0.75,0.2,0.2), mar=c(2.,5.,2.0,1) )
  ylab.txt <- expression(paste("CO")[2]*paste(" flux, (", mu, "mol/", m^2, s,")"))

  # add time series 
  #plot(x=flux.table$date.time, y=flux.table$flux.co2, type="p",ylab=ylab.txt, xlab="Observation Period [Month date]", xaxt="n",
  #    col="black",cex=1.0, ylim=c(-3E1,2E1))
  #set x-axis label
  x<- flux.table$date.dd ;  at <- seq(min(x), max(x), 7)
  #axis(side=1, at=at, labels=format(at, "%b-%d"), cex.axis=1.2) 


  par(cex.axis=1.5); par(cex.lab=1.5) # is for y-axis

  boxplot(main=paste(crop_season,"_cropping season",sep=""), flux.table$flux.co2 ~ flux.table$date.dd, ylab=ylab.txt, ylim=c(-2.5E1,2E1),col="gray", outline=FALSE)
  #x<- flux.table$date.dd ;  at <- seq(min(x), max(x), 7)
  #axis(side=1, at=at, labels=format(at, "%b-%d"), cex.axis=1.2) 
  abline(a=NULL, b=NULL, h=0, v=NULL, col="black")
  grid() 
  means <- tapply(flux.table$flux.co2, flux.table$date.dd, mean)
  
  dd_table <- data.frame(mean=means, cumsum=cumsum(means) )
  points(means, pch=21, cex=0.5, bg="gray")
  #mtext(side=1, "Julian day", line=2.0, cex=1.2) 
 
  boxplot(flux.table$flux.co2 ~ flux.table$date.week.hh, ylab=ylab.txt, ylim=c(-2.5E1,2E1),col="gray",xaxt="n", outline=FALSE)
  abline(a=NULL, b=NULL, h=0, v=NULL, col="black")
 
  x<- flux.table$date.week.hh ; 
  at <- seq(min(x), max(x), 24); xlab<- format(lubridate::ymd( "2022-01-04" ) + lubridate::weeks( at/24 - 1 ),"%b-%d")
  grid()
  par(xpd=T)
  text(x=at-min(x)+12 , y=rep(24,length(at)), paste("[",xlab,")") , cex=1)  
  par(xpd=F)  
  mtext(side=1, "Diurnal-pattern, [weekly average]", line=0.6, cex=1.2) 
  dev.off()

  #CH4 flux
  png(file=paste(substr(table.name,start=5,stop=14),"_CH4_flux.plot.png",sep=""), width=1200, height=850, res=128)
  par(mfrow=c(2,1),mai=c(0.75,0.75,0.2,0.2), mar=c(2.,5.,2.0,1) )
  par(cex.axis=1.5); par(cex.lab=1.5) # is for y-axis

  ylab.txt <- expression(paste("CH")[4]*paste(" flux, (", mu, "mol/", m^2, s,")"))
  boxplot(flux.table$flux.ch4 ~ flux.table$date.dd, ylab=ylab.txt, ylim=c(-.5E-1,1.E-1),col="orange", outline=FALSE)
  means <- tapply(flux.table$flux.ch4, flux.table$date.dd, mean)
  points(means, pch=21, cex=0.5,bg="gray")
  #mtext(side=1, "Julian day", line=2.0, cex=1.2) 
  abline(a=NULL, b=NULL, h=0, v=NULL, col="black")
  grid()
 
  boxplot(flux.table$flux.ch4 ~ flux.table$date.week.hh, ylab=ylab.txt, ylim=c(-.5E-1,1.E-1),col="orange",xaxt="n", outline=FALSE)
  grid()
  x<- flux.table$date.week.hh ; 
  at <- seq(min(x), max(x), 24); xlab<- format(lubridate::ymd( "2022-01-04" ) + lubridate::weeks( at/24 - 1 ),"%b-%d")
  grid()
  par(xpd=T)
  text(x=at-min(x)+12 , y=rep(1.15E-1,length(at)), paste("[",xlab,")") , cex=1)  
  par(xpd=F)  
  mtext(side=1, "Diurnal-pattern, [weekly average]", line=0.6, cex=1.2) 
  dev.off()


  ld_go <- TRUE
  if (ld_go) {
  #Latent heat 
  png(file=paste(substr(table.name,start=5,stop=14),"_Latent_heat_flux.plot.png",sep=""), width=1024, height=768, res=128)
  par(mfrow=c(2,1),mai=c(0.75,0.75,0.2,0.2), mar=c(2,4.5,1.5,1) )

  ylab.txt <- expression(paste("LE")*paste(" flux, (", "W/", m^2,")"))
  boxplot(flux.table$flux.le ~ flux.table$date.dd, ylab=ylab.txt, ylim=c(-5E1,5E2),col="blue", outline=FALSE)
  means <- tapply(flux.table$flux.le, flux.table$date.dd, mean)
  points(means, pch=21, cex=0.5,bg="gray")
  #mtext(side=1, "Julian day", line=2.0, cex=1.2) 
  abline(a=NULL, b=NULL, h=0, v=NULL, col="black")
  grid()
 
  #set x-axis label
  x<- flux.table$date.time ;  at <- seq(min(x), max(x), "week")
  axis(side=1, at=at, labels=format(at, "%b-%d"), cex.axis=1.2) 

  boxplot(flux.table$flux.le ~ flux.table$date.week.hh, ylab=ylab.txt,  ylim=c(-5E1,5E2), col="blue",xaxt="n",outline=FALSE)
  grid()
  mtext(side=1, "Diurnal-pattern, [weekly]", line=0.6, cex=1.2) 
  dev.off()

  #Sensible heat
  png(file=paste(substr(table.name,start=5,stop=14),"_Sensible_heat_flux.plot.png",sep=""), width=1024, height=768, res=128)
  par(mfrow=c(2,1),mai=c(0.75,0.75,0.2,0.2), mar=c(2,4.5,1.5,1) )
 
  ylab.txt <- expression(paste("SH")*paste(" flux, (", "W/", m^2,")"))
  boxplot(flux.table$flux.sh ~ flux.table$date.dd, ylab=ylab.txt, ylim=c(-5E1,5E2),col="red", outline=FALSE)
  means <- tapply(flux.table$flux.sh, flux.table$date.dd, mean)
  points(means, pch=21, cex=0.5,bg="gray")
  #mtext(side=1, "Julian day", line=2.0, cex=1.2) 
  abline(a=NULL, b=NULL, h=0, v=NULL, col="black")
  grid()
  boxplot(flux.table$flux.sh ~ flux.table$date.week.hh, ylab=ylab.txt, ylim=c(-5E1,5E2), col="red",xaxt="n",outline=FALSE)
  grid()
  mtext(side=1, "Diurnal-pattern, [weekly]", line=0.6, cex=1.2) 
  dev.off()
  }


}


