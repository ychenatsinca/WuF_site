# This is a script for showing the result from the flux.table both for the daily and diurnal evolution over the selected average windows (box_days + ini_day) 
# date: 2023-04-27 
# author: yi-ying chen 

# load library 
library("lubridate")

# set flux table name
table.name= c("WuF_2023-01-01 00:00:00_365_flux.table.txt")

# read in the flux.table
flux.table <-  read.csv(table.name) 


flux.table$date.time <- as.POSIXct(flux.table$date.time,format='%Y-%m-%d %H:%M:%S')

#set days for box plot 
box_days <- 14
ini_day <-  95

in_date <- flux.table$date.time[1 + ini_day*48 ]  
ed_date <- in_date + days(box_days) 
flux.table$date.hh <- hour(flux.table$date.time)   

 png(file=paste(substr(in_date,start=1,stop=10),"_",substr(ed_date,start=1,stop=10),"_flux.plot.png",sep=""), width=1024, height=768, res=128)

 par(mfrow=c(4,2),mai=c(0.75,0.75,0.2,0.2), mar=c(2,4.5,1.5,1) )

 #CO2 flux
 ylab.txt <- expression(paste("CO2 flux, (", mu, "mol/", m^2, s,")"))
 plot(x=flux.table$date.time, y=flux.table$flux.co2, type="p",ylab=ylab.txt, xlab="Observation Period [Month date]", xaxt="n",
      col="green",cex=0.8, pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1) , ylim=c(-3E1,2E1));grid()
 #add gray area 
  rect(xleft=in_date, xright=ed_date, ybottom=par("usr")[3], ytop=par("usr")[4], density=NA, col=adjustcolor("gray", alpha = 0.3)) 
 #rect(xleft=in_date, xright=ed_date, ybottom=par("usr")[3], ytop=par("usr")[4], density=NA, col="gray")
 #set xaex label
  x<- flux.table$date.time ;  at <- seq(min(x), max(x), "month")
  axis(side=1, at=at, labels=format(at, "%b"), cex.axis=1.2) 
  mtext(paste("shade area: ",substr(in_date,start=1,stop=10)," to ",substr(ed_date,start=1,stop=10)),side=3,line=0.5,
             at=par("usr")[1]+0.5*diff(par("usr")[1:2]),
             cex=0.8, col="darkgray") 
#subset the flux.table 
 plot.table <- subset ( flux.table, (flux.table$date.time >= in_date & flux.table$date.time <= ed_date) )
 boxplot(plot.table$flux.co2 ~ plot.table$date.hh , ylab=ylab.txt, xlab="Local hour [00:00 to 24:00]", col="green", ylim=c(-3E1,2E1));grid()

 #CH4 Flux
 ylab.txt <- expression(paste("CH4 flux, (", mu, "mol/", m^2, s,")")) 
 plot(x=flux.table$date.time, y=flux.table$flux.ch4, type="p",ylab=ylab.txt, xlab="Observation Period [Month date]",xaxt="n", 
      col="orange",cex=0.8, pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1) , ylim=c(-1E-1,2E-1)) ;grid()
 #add gray area 
 rect(xleft=in_date, xright=ed_date, ybottom=par("usr")[3], ytop=par("usr")[4], density=NA, col=adjustcolor("gray", alpha = 0.3)) 
 axis(side=1, at=at, labels=format(at, "%b"), cex.axis=1.2) 
 boxplot(plot.table$flux.ch4~plot.table$date.hh, ylab=ylab.txt, xlab="Local hour [00:00 to 24:00]", col="orange", ylim=c(-1E-1,2E-1));grid()

 #Latent heat flux
 ylab.txt <- expression(paste("Latent heat, (", "W/", m^2,")"))
 plot(x=flux.table$date.time, y=flux.table$flux.le, type="p",ylab=ylab.txt, xlab="Observation Period [Month date]",xaxt="n", 
      col="lightblue",cex=0.8, pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1) , ylim=c(-5E1,5E2));grid()
 axis(side=1, at=at, labels=format(at, "%b"), cex.axis=1.2) 
 #add gray area 
 rect(xleft=in_date, xright=ed_date, ybottom=par("usr")[3], ytop=par("usr")[4], density=NA, col=adjustcolor("gray", alpha = 0.3))
 boxplot(plot.table$flux.le~plot.table$date.hh, ylab=ylab.txt, xlab="Local hour [00:00 to 24:00]", col="lightblue", ylim=c(-5E1,5E2));grid()
 
 #Sensible heat flux
 ylab.txt <- expression(paste("Sensible heat, (", "W/", m^2,")")) 
 plot(x=flux.table$date.time, y=flux.table$flux.sh, type="p",ylab=ylab.txt,xlab="Observation Period [Month date]", xaxt="n",
      col="red",cex=0.8, pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1) , ylim=c(-5E1,5E2));grid()
 axis(side=1, at=at, labels=format(at, "%b"), cex.axis=1.2)  
 #add gray area 
 rect(xleft=in_date, xright=ed_date, ybottom=par("usr")[3], ytop=par("usr")[4], density=NA, col=adjustcolor("gray", alpha = 0.3))
 boxplot(plot.table$flux.sh~plot.table$date.hh, ylab=ylab.txt, xlab="Local hour [00:00 to 24:00]", col="red", ylim=c(-5E1,5E2));grid()
 
 dev.off()

