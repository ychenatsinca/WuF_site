# ------------------------------------------------------------
#
#     The program calculates the eddy flux from two EC systems
#     open path EC-150 & Young 81000 (H2O, CO2, Ta, u, v, w),
#     closed path: TILDAS(HNO3, N2O, HONO, H2O).
#     Surface fluxes from gradient approach
#     all gd.fluxes reference to the sensible heat flux or N2O flux  
#     Author: Yi-Ying Chen (RCEC, Taiwan) 
#     Contact: yiyingchen@gate.sinica.edu.tw
#     First Date: 2019-07-18
#     Revised: 2020-04-07
#     Version: 1.1
#     Revised: 2023-03-20, 2023-03-30  
#     Version: 1.2 (CO2, H20, CH4, Ta only) without N2O/Tildas 
#     
#     Data Processing Tasks: 
#     1.detrend, 2.despike, 3.lag time fixed, 4.co-spectra check
#     5.coordinate rotation, 6.turbulent test, 7.stationary test
#     8.footprint test, 9.wpl correction(added, 2023-03-30) 
#     10.gradient flux    
# ------------------------------------------------------------

# initial  tables
  flux.table <- data.frame()
  spec.table <- data.frame()
#load source file for the functions used in the program
  source("fun_load_ec150.R")
  #source("fun_load_test.R")
  #source("fun_load_young.R")
  source("fun_despike.R")
  source("fun_ecflux_for.R")
  source("fun_turbulent_test.R")
#load library for the signal processing
  library("signal")
  library("itsmr")
  library("pracma")
  library("Publish")

#--- set the date for analysis 
  site.name = "WuF"
 #  site.name = "MPI"
  
  ini.stamp ="2023-03-16 12:00:00"
  ndays=40
  n30=48*ndays  

  date.id   <- format(seq(as.POSIXct(ini.stamp, tz=""), length.out=n30, by='30 min'),'%Y-%m-%d %H:%M:%S')
# date.id <- format(seq(as.POSIXct("2019-03-14 12:00:00", tz=""), length.out=10, by='15 min'),'%Y_%j %H:%M')
  date.yyyy <- substr(date.id,start=1,stop=4)
  date.mm   <- substr(date.id,start=6,stop=7)
  date.dd   <- substr(date.id,start=9,stop=10)
  date.hhmmss <- substr(date.id,start=12,stop=19)
  time.stamp <- format(paste(substr(date.hhmmss,start=1,stop=2),":",
			   substr(date.hhmmss,start=4,stop=5),":00.0",sep=""))


  #  
  #par(mfrow=c(2,2),mai=c(1,1,1,1) )
  #
  ld.go<-TRUE
  #
  # load library signal for filtering and despike with a crtiria of 3.5-std 
      library("signal")
     # source("fun_despike.R")
     # --- using band pass filter to correct frequency response of closed-path EC measurenment --- 
     #create the band-pass filter
     #bf.h2o <- butter(1, c(0.0,0.9), type="pass")
     #apply the filter to the h2o noisy signal 
     #raw.tildas$H2O.2 <- filtfilt(bf.h2o, as.numeric(raw.tildas$H2O.2)) 
     #create the band-pass filter
     #bf.n2o <- butter(1, c(0.0,0.9), type="pass")
     #apply the filter to the h2o noisy signal 
     #raw.tildas$N2O <- filtfilt(bf.n2o, as.numeric(raw.tildas$N2O)) 
     #create the band-pass filter
     #bf.hno3 <- butter(1, c(0.0,0.9), type="pass")
     #apply the filter to the h2o noisy signal 
     #raw.tildas$HNO3 <- filtfilt(bf.hno3, as.numeric(raw.tildas$HNO3))
     #create the band-pass filter
     #bf.hono <- butter(1, c(0.0,0.9), type="pass")
     #apply the filter to the h2o noisy signal 
     #raw.tildas$HONO <- filtfilt(bf.hono, as.numeric(raw.tildas$HONO)) 

  # load library pracma for detrending
  # library("pracma")
  # load fun_tubulent_test() function(x,ustar,SH,LE,z,d0,Ta,roha,L,flag.itc=TRUE,flag.stat=TRUE)
  # to get nstationary flag and integral turbulent flag 
  #
  # source("fun_turbulent_test.R") 
  # open pdf device 
  # pdf("cssargi_detail.pdf")

  if(ld.go) {
     for (i in 1:length(date.yyyy)) {
     #for (i in 1) {
     # start flux.cal
     #check file exist 
     wrk.yr = date.yyyy[i]
     ec150.dir    <- paste("/lfs/home/ychen/lfs_dir/EC_DATA/",site.name,"/30min/",wrk.yr,sep="")
 
	     file.name <- paste(ec150.dir,"/",site.name,"_",date.yyyy[i],"_",date.mm[i],"_",date.dd[i],"_",
			       substr(date.hhmmss[i],start=1,stop=2),substr(date.hhmmss[i],start=4,stop=5),".dat",sep="")
	     print(file.name)
	     if(file.exists(file.name)==TRUE) {
             raw.data <- fun_load_ec150(
				dir.name="/lfs/home/ychen/lfs_dir/EC_DATA/", site.name=site.name,
                             	dir.subname=date.yyyy[i],
 				yyyy=date.yyyy[i],mm=date.mm[i],dd=date.dd[i],
                                hhmm= paste(substr(date.hhmmss[i],start=1,stop=2),substr(date.hhmmss[i],start=4,stop=5),sep=""),
                             	ld.plot=F,ld.spec=F, ld.na=T, ld.spik=T, dir.plot="./png_plot/")

	     }else{ print("NO EC_DATA file was found!"); next}
    
     # remove un gapfilled NA row in the raw.data
     raw.data <- na.omit(raw.data)
     #create the band-pass filter
     #bf.ch4 <- butter(1, c(0.0,0.025), type="pass")
     #apply the filter to the ch4 noisy signal 
     #raw.data$CH4_mole_fraction <- filtfilt(bf.ch4, as.numeric(raw.data$CH4_mole_fraction)) 

 


     # --- coordinate rotation and covariance calculation via fortran subroutine --- 
    try(
   	obj  <- fun_ecflux_for(raw.data=raw.data)
       )
         #vertical wind after double rotation
     if( exists("obj"))  w_0 <- obj$flux.ec$w
     

    ld.go <-FALSE
     if (ld.go) {
     # --- frequency analysis --- 
     cospec<- spectrum(x=c(as.numeric(raw.data$Uzc),
			   (as.numeric(raw.data$CH4_mole_fraction)-mean(as.numeric(raw.data$CH4_mole_fraction),na.rm=T))),
			   method="pgram",plot=FALSE)
     ch4.spec <- cospec$spec
     #plot(cospec$spec~cospec$freq, log="xy", main="W' HNO3' Co-Spectrum")
     cospec<- spectrum(x=c(as.numeric(raw.data$Uzc),
			   (as.numeric(raw.data$CO2)-mean(as.numeric(raw.data$CO2),na.rm=T))),
		           method="pgram",plot=FALSE)
     co2.spec <- cospec$spec
     #plot(cospec$spec~cospec$freq, log="xy", main="W' HONO' Co-Spectrum")
     cospec<- spectrum(x=c(as.numeric(raw.data$Uzc),
			   (as.numeric(raw.data$H2O)-mean(as.numeric(raw.data$H2O),na.rm=T))),
		           method="pgram",plot=FALSE)
     h2o.spec <- cospec$spec
     #plot(cospec$spec~cospec$freq, log="xy", main="W' N2O' Co-Spectrum")
     cospec<- spectrum(x=c(as.numeric(raw.data$Uzc),
			   (as.numeric(raw.data$Tsc)-mean(as.numeric(raw.data$Tsc),na.rm=T))),
		           method="pgram",plot=FALSE)
     tsc.spec <- cospec$spec
     #plot(cospec$spec~cospec$freq, log="xy", main="W' H2O' Co-Spectrum")
     new.sheet <- data.frame(freq=cospec$freq, ch4.spec= ch4.spec, co2.spec=co2.spec, h2o.spec=h2o.spec, tsc.spec=tsc.spec)
     spec.table <- rbind(spec.table, new.sheet)    
     }


    #add date.time stamp 
    date.time <- as.POSIXct(date.id[i],format='%Y-%m-%d %H:%M:%S') 

    
    #convert charaters to numeric 
    raw.data$Uxc <- as.numeric(raw.data$Uxc)
    raw.data$Uyc <- as.numeric(raw.data$Uyc)
    raw.data$Uzc <- as.numeric(raw.data$Uzc)
    raw.data$Tsc <- as.numeric(raw.data$Tsc)
    raw.data$H2O <- as.numeric(raw.data$H2O)
    raw.data$CO2 <- as.numeric(raw.data$CO2)
    raw.data$cell_press <- as.numeric(raw.data$cell_press)
    #
    if (any(names(raw.data) == "CH4_density")) { 
        raw.data$CH4 <- as.numeric(raw.data$CH4_density)
    }else{
        raw.data$CH4 <- as.numeric(raw.data$CH4)
 
    }

     #adjust the lag time
    library("binhf") 
    #print("Find Lag time ...")
     # --- using acf to find the maximum corelation of two covariance terms ---
     # --- the unit conversion was also applied to all gas specis ---
     # --- H2O ppb(1E-9) to mg/m^3 (1E-6)
     # --- HNO3 ppb(1E-9) ~ g/m3 (1E-9) to nmol(1/63)
     # --- HONO ppb(1E-9)  to nmol (1/47)
     # --- N2O  ppb(1E-9)  to nmol (1/44)
     # copy tildas data to raw data
     #w.id1=1
     #w.id2=100 
     #var1 <- (as.numeric(raw.data$CH4[w.id1:w.id2])-mean(as.numeric(raw.data$CH4[w.id1:w.id2]),na.rm=T)) 
     #var2 <- w_0 
  
     # find  lag time
     #aa <- ccf(x= var1,  y= var2,  lag.max = 100*1.0, plot=FALSE)
     #lag.id = (aa$lag[which.max((aa$acf))])
     #print(paste("CH4 lag time:",lag.id," was found! ",sep=""))
     #raw.data$CH4 <- shift(raw.data$CH4, lag.id, dir="right")


    # covariance calculation 
    # for w_t
    if (any(names(raw.data) == "Tsc") )  { 
        w_t = cov(w_0, raw.data$Tsc)
       }else{
        w_t = NA
        }
    #for w_h2o
    if (any(names(raw.data) == "H2O")  ) { 
        w_h2o = cov(w_0, raw.data$H2O)
       }else{
        w_h2o = NA
        }
    #for w_co2
    if (any(names(raw.data) == "CO2")   ) {  
        w_co2 = cov(w_0, raw.data$CO2)
       }else{
        w_co2 = NA
        }
    #for w_ch4 # unit in [m/s *  mg/m3]
    if ( any(names(raw.data) == "CH4") ) { 
        w_ch4 = cov(w_0, raw.data$CH4*16.)
       }else{
        w_ch4 = NA
        }

    # calculation of the moist air property. air density and air temperature
  #  Rhoa(j)=(Pa*1E5)/(287.05*(273.15+T(j)))        ! unit in [g/m3]      dry air   density   
  #  T(j)=T(j)*(1+0.514*(Rhov(j)/Rhoa(j)))          !-----Change Sonic Temp to Air Temp
    
     #convert sonic temperature to air temperature 
     raw.data$Rhoa <- (((raw.data$cell_press)*1E3)/(287.05*(273.15+raw.data$Tsc)))*1E3   # pressure unit in [Kpa]  temperature unit in [oC] unit in [g/m3]
     raw.data$Ta <- raw.data$Tsc * (1/(1.0 + 0.514*(raw.data$H2O)/raw.data$Rhoa)) 
 
    # NO AIR PRESSURE DATA USE 100.0 Kpa  
     #raw.data$Rhoa <- (100.0*1E6)/(287.05*(273.15+as.numeric(raw.data$Tsc)))
     #raw.data$Ta <- as.numeric(raw.data$Tsc) #* (1- 0.514*(as.numeric(raw.data$H2O)/raw.data$Rhoa))


     air.temp <- mean(raw.data$Ta,na.rm=TRUE)           # [oC]
     air.den  <- mean(raw.data$Rhoa, na.rm=TRUE)        # [g/m3]
     air.press <- mean(raw.data$cell_press, na.rm=TRUE) # [Kpa]
     h2o.den <- mean(raw.data$H2O, na.rm=T)             # [g/m3] 
     co2.den <- mean(raw.data$CO2, na.rm=T)             # [mg/m3]
     ch4.den <- mean(raw.data$CH4, na.rm=T)*16.0        # [mmol/m3]--> [mg/m3]

#get WPL correction factors for CO2 and CH4 flux calculation
Q=(0.*(air.temp**2.)) + (-1.3*1E-7*air.temp) + (3.7*1E-5) 
R=(4.0*1E-8*(air.temp**2.)) + (1.1*1E-5*air.temp) + (2.18*1E-3)
S=(2.0*1E-6*(air.temp**2.)) + (9.8*1E-4*air.temp) + 0.378 
Pe=air.press 
xv=(mean(raw.data$H2O,na.rm=T)/mean(raw.data$Rhoa,na.rm=T))
TkTk= ((-4.0*1E-8*(air.temp**2.) + 1.55*1E-5*air.temp + -7.0*1E-3)*air.press)  + (-4.7*1E-6 * air.temp**2. + 3.0*1E-3*air.temp + 0.927 )  

Big_A <- (Q*Pe**2.) + (R*Pe) + S
Big_B <- 1.0 + (1.0-1.46*xv) * ((-8.2*1E-6*air.temp + 4.3*1E-3)*air.press) + (-1.7*1E-4*air.temp + 0.03)
Big_C <- 1.0 + (1.0-1.*xv)*TkTk  + xv*(Big_B-1.0)

print(paste("Big_A:",formatC(Big_A,digits=4,width=8,format="f",flag="-"),
            "Big_B:",formatC(Big_B,digits=4,width=8,format="f",flag="-"),
            "Big_C:",formatC(Big_C,digits=4,width=8,format="f",flag="-"),sep="  "))


mu=1.61
sigma=(h2o.den/(air.den-h2o.den))

co2.flux = w_co2  + mu*(co2.den/(air.den-h2o.den))*w_h2o + (1.0 + mu*sigma)*co2.den*w_t/(air.temp+273.15)   
ch4.flux = Big_A * ( w_ch4 + Big_B *mu*(ch4.den/(air.den-h2o.den))*w_h2o + (Big_C*(1.0+mu*sigma)*ch4.den*w_t/(air.temp+273.15)   )   )
ch4.flux.0=  w_ch4 
ch4.flux.1=  Big_B *mu*(ch4.den/(air.den-h2o.den))*w_h2o
ch4.flux.2=  Big_C*(1.0+mu*sigma)*ch4.den*w_t/(air.temp+273.15)  
#


# if NAN   
if ( is.na(Big_A) != TRUE) { 
#
   # --- tubulent, stationary, and footprint  test flux 
   flux.test <- fun_turbulent_test(x=raw.data$Uzc, ustar=obj$flux.ec$ustar,uavg=obj$flux.ec$uavg,
				     SH=(air.den/1E3) * 1005. * w_t, 
				     LE=(air.den/1E3) * 2.504 * 1E3 * w_h2o,
				     Ta=air.temp,roha=air.den, zm=0.7,
				     wd=obj$flux.ec$wd)
 }else { 
   # do notthing  
}


    new.row <- data.frame(date.time=date.time, 
                           date.mm=substr(date.time,start=6,stop=7), 
			   date.hh=substr(date.time,start=12,stop=13),
			   u.avg=mean(raw.data$Uxc,na.rm=T),
			   v.avg=mean(raw.data$Uyc,na.rm=T),
                           w.avg=mean(raw.data$Uzc,na.rm=T),
                         co2.avg=co2.den,
                         h2o.avg=h2o.den,
                         ch4.avg=ch4.den,
                         tsc.avg=mean(raw.data$Tsc,na.rm=T),
                         air.press=air.press,
                         air.den=air.den,
                         air.temp= air.temp,
			 flux.ustr=obj$flux.ec$ustar,
                         flux.sh= (air.den/1E3) * 1005. * w_t ,
                         flux.le= (air.den/1E3) * 2.504 * 1E3 * w_h2o,
                         flux.co2= co2.flux/44.*1E3,  
                         flux.ch4= ch4.flux/16.*1E3,
                         flux.ch4.0=ch4.flux.0/16.*1E3,
                         flux.ch4.1=ch4.flux.1/16.*1E3,
                         flux.ch4.2=ch4.flux.2/16.*1E3,
                         flag.xmx=flux.test$xmax,
                         flag.sta=flux.test$flag.stat,
			 flag.itc=flux.test$flag.itc,
			 flag.fpt=flux.test$flag.fpt )
     # show the result
     print(new.row)

     # append the new row to the data.frame
       flux.table <- rbind(flux.table, new.row)   
     
    }#end for 
     #dev.off()
   }#end if

  #library("Publish") 
  # add the units information into the flux.table
  flux.table <- Units(flux.table,list(flux.sh="W/m^2", flux.le="W/m^2", flux.co2="umol(CO2)/m^2-s", flux.ch4="umol(CH4)/m^2-s"))
 
  #write out the flux.table
  write.table(flux.table,file=paste(site.name,"_",ndays,"_flux.table.txt",sep=""), sep=",",row.name=FALSE)

  # calculate mean of the co-spectra 
  #spec.mean <- aggregate( spec.table[,2:5], by=list(spec.table$freq), FUN=mean, na.action = na.omit) 
  #plot(x=as.numeric(spec.mean$Group.1)*10,y=spec.mean$nh2o/max(spec.mean$h2o), log="xy",ylim=c(1e-4,1e0),
  #     main="w'H2O' co-spectra",cex=0.5,ylab="normailzed co-spectra",xlab="frequency",col="blue")
  #plot(x=as.numeric(spec.mean$Group.1)*10,y=spec.mean$n2o/max(spec.mean$n2o), log="xy",ylim=c(1e-4,1e0),
  #     main="w'N2O' co-spectra",cex=0.5,ylab="normailzed co-spectra",xlab="frequency",col="brown")
  #plot(x=as.numeric(spec.mean$Group.1)*10,y=spec.mean$hono/max(spec.mean$hono), log="xy",ylim=c(1e-4,1e0),
  #     main="w'HONO' co-spectra",cex=0.5,ylab="normailzed co-spectra",xlab="frequency",col="black")
  #plot(x=as.numeric(spec.mean$Group.1)*10,y=spec.mean$hno3/max(spec.mean$hno3), log="xy",ylim=c(1e-4,1e0),
  #     main="w'HNO3' co-spectra",cex=0.5,ylab="normailzed co-spectra",xlab="frequency",col="red")

#
# 
#flux.table$flux.ch4  
flux.table$flux.ch4[abs(flux.table$flux.ch4)>=0.2  ] <-  NA


plot.ld <- TRUE
if(plot.ld)
{
# pdf(file=paste("./plot_pdf/flux.table.timeseries.pdf",sep=""))
 par(mfrow=c(4,2),mai=c(0.75,0.75,0.2,0.2) )

 #CO2 flux
 ylab.txt <- expression(paste("CO2 flux, (", mu, "mol/", m^2, s,")"))
 plot(x=flux.table$date.time, y=flux.table$flux.co2, type="p",ylab=ylab.txt, xlab="Observation Period [Month date]", 
      col="gray",cex=0.8, pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1) , ylim=c(-3E1,2E1));grid()
 boxplot(flux.table$flux.co2~flux.table$date.hh, ylab=ylab.txt, xlab="Local hour [00:00 to 24:00]", col="gray", ylim=c(-3E1,2E1));grid()

 #CH4 Flux
 ylab.txt <- expression(paste("CH4 flux, (", mu, "mol/", m^2, s,")")) 
 plot(x=flux.table$date.time, y=flux.table$flux.ch4, type="p",ylab=ylab.txt, xlab="Observation Period [Month date]", 
      col="orange",cex=0.8, pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1) , ylim=c(-1E-1,2E-1));grid()
 boxplot(flux.table$flux.ch4~flux.table$date.hh, ylab=ylab.txt, xlab="Local hour [00:00 to 24:00]", col="orange", ylim=c(-1E-1,2E-1));grid()

 #Latent heat flux
 ylab.txt <- expression(paste("Latent heat flux, (", "W/", m^2,")"))
 plot(x=flux.table$date.time, y=flux.table$flux.le, type="p",ylab=ylab.txt, xlab="Observation Period [Month date]", 
      col="blue",cex=0.8, pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1) , ylim=c(-5E1,5E2));grid()
 boxplot(flux.table$flux.le~flux.table$date.hh, ylab=ylab.txt, xlab="Local hour [00:00 to 24:00]", col="blue", ylim=c(-5E1,5E2));grid()
 
 #Sensible heat flux
 ylab.txt <- expression(paste("Sensible heat flux, (", "W/", m^2,")")) 
 plot(x=flux.table$date.time, y=flux.table$flux.sh, type="p",ylab=ylab.txt,xlab="Observation Period [Month date]", 
      col="red",cex=0.8, pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1) , ylim=c(-5E1,5E2));grid()
 boxplot(flux.table$flux.sh~flux.table$date.hh, ylab=ylab.txt, xlab="Local hour [00:00 to 24:00]", col="red", ylim=c(-5E1,5E2));grid()
 

# dev.off()



} 




