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
#     8.footprint test, 9.wpl correction(fixed, 2023-03-30) 
#     10.gradient flux    
# ------------------------------------------------------------

# initial  tables
  flux.table <- data.frame()
  spec.table <- data.frame()
#load source file for the functions used in the program
  source("fun_load_ec150.R")
  #source("fun_load_young.R")
  source("fun_despike.R")
  source("fun_ecflux_for.R")
  source("fun_turbulent_test.R")
#load library for the signal processing
  library("signal")
  library("itsmr")
  library("pracma")

#--- set the date for analysis 
  site.name = "WuF"
  
  ini.stamp ="2022-10-26 00:00:00"
  ndays=30
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
				dir.name="/lfs/home/ychen/lfs_dir/EC_DATA/", site.name="WuF",
                             	dir.subname=date.yyyy[i],
 				yyyy=date.yyyy[i],mm=date.mm[i],dd=date.dd[i],
                                hhmm= paste(substr(date.hhmmss[i],start=1,stop=2),substr(date.hhmmss[i],start=4,stop=5),sep=""),
                             	ld.plot=F,ld.spec=F, ld.na=TRUE, ld.spik=FALSE, dir.plot="./png_plot/")

	     }else{ print("NO EC-150 file"); next}
    
     # remove un gapfilled NA row in the raw.data
     raw.data <- na.omit(raw.data)
     #create the band-pass filter
     #bf.ch4 <- butter(1, c(0.0,0.025), type="pass")
     #apply the filter to the ch4 noisy signal 
     #raw.data$CH4_mole_fraction <- filtfilt(bf.ch4, as.numeric(raw.data$CH4_mole_fraction)) 

 
     print("Find Lag time ...")
     # --- using acf to find the maximum corelation of two covariance terms ---
     # --- the unit conversion was also applied to all gas specis ---
     # --- H2O ppb(1E-9) to mg/m^3 (1E-6)
     # --- HNO3 ppb(1E-9) ~ g/m3 (1E-9) to nmol(1/63)
     # --- HONO ppb(1E-9)  to nmol (1/47)
     # --- N2O  ppb(1E-9)  to nmol (1/44)
     # copy tildas data to raw data
     #w.id1=1
     #w.id2=100 
     #var1 <- (as.numeric(raw.data$CH4_mole_fraction[w.id1:w.id2])-mean(as.numeric(raw.data$CH4_mole_fraction[w.id1:w.id2]),na.rm=T)) 
     #var2 <- raw.data$Uzc
  
     # find  lag time
    # aa <- ccf(x= var1,  y= var2,  lag.max = 100*1.0, plot=FALSE)
    # lag.id = (aa$lag[which.max((aa$acf))])
    #  print(paste("CH4 lag time:",lag.id," was found! ",sep=""))

    #adjust the lag time
    library("binhf")
    #raw.data$CH4_mole_fraction <- shift(raw.data$CH4_mole_fraction, lag.id, dir="right")



     # --- coordinate rotation and covariance calculation via fortran subroutine --- 
    try(
   	obj  <- fun_ecflux_for(raw.data=raw.data)
       )
         #vertical wind after double rotation
     if( exists("obj"))  w_0 <- obj$flux.ec$w
     


     # --- tubulent, stationary, and footprint  test flux 
  #   flux.test <- fun_turbulent_test(x=raw.data$Uzc, ustar=obj$flux.ec$ustar,uavg=obj$flux.ec$uavg,
#				     SH=obj$air.den * 1005. * obj$flux.ec$f1, 
#				     LE=obj$air.den * 2.504 * 1E3 * obj$flux.ec$f2,
#				     Ta=obj$air.temp,roha=obj$air.den, zm=0.7,
#				     wd=obj$flux.ec$wd)
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
    }


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
    #for w_ch4 
    if ( any(names(raw.data) == "CH4") ) { 
        w_ch4 = cov(w_0, raw.data$CH4)
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

Big_A <- (Q*Pe) + (R*Pe) + S
Big_B <- 1.0 + (1.0-1.46*xv) * ((-8.2*1E-6*air.temp + 4.3*1E-3)*air.press) + (-1.7*1E-4*air.temp + 0.03)
Big_C <- 1.0 + (1.0-1.*xv)*TkTk  + xv*(Big_B-1.0)

mu=1.61
sigma=(h2o.den/(air.den-h2o.den))

co2.flux = w_co2  + mu*(co2.den*1E-3/(air.den-h2o.den))*w_h2o + (1.0 + mu*sigma)*co2.den*w_t/(air.temp+273.15)   
ch4.flux = Big_A * ( w_ch4 + Big_B *mu*(ch4.den*1E3/(air.den-h2o.den))*w_h2o + (Big_C*(1.0+mu*sigma)*ch4.den*w_t/(air.temp+273.15)   )   )

#

    new.row <- data.frame(date.time=date.time, 
                           date.mm=substr(date.time,start=6,stop=7), 
			   date.hh=substr(date.time,start=12,stop=13),
			   u.avg=mean(raw.data$Uxc,na.rm=T),
			   v.avg=mean(raw.data$Uyc,na.rm=T),
                           w.avg=mean(raw.data$Uzc,na.rm=T),
                         co2.avg=mean(raw.data$CO2,na.rm=T),
                         h2o.avg=mean(raw.data$H2O,na.rm=T),
                         ch4.avg=mean(raw.data$CH4,na.rm=T),
                         tsc.avg=mean(raw.data$Tsc,na.rm=T),
                         air.press=air.press,
                         air.den=air.den,
                         air.temp= air.temp,
			 flux.ustr=obj$flux.ec$ustar,
                         flux.sh= (air.den/1E3) * 1005. * w_t ,
                         flux.le= (air.den/1E3) * 2.504 * 1E3 * w_h2o,
                         flux.co2= co2.flux,  
                         flux.ch4= ch4.flux )
                          
                         #flux.db.sh= 1.24* 1005. * cov(obj$flux.ec$w, as.numeric(raw.data$Tsc)),
                         #flux.db.le= 1.24* 2.504 * 1E3* cov(obj$flux.ec$w, as.numeric(raw.data$H2O)),
                         #flux.db.co2= cov(obj$flux.ec$w, as.numeric(raw.data$CO2)) )
                         #flux.db.ch4= cov(obj$flux.ec$w, as.numeric(raw.data$CH4_mole_fraction))  ) 
                         #)
 # to mg 
			#flux.ch4=obj$flux.ec$f4 , # to nmol
 			# air.wdr=obj$flux.ec$wd,
                        # air.utr=obj$flux.ec$ustar, 
                        # air.den=obj$air.den, 
                        # air.pre=obj$air.press,
                        # air.tmp=obj$air.temp,
                        # flag.xmx=flux.test$xmax,
                        # flag.sta=flux.test$flag.stat,
			# flag.itc=flux.test$flag.itc,
			# flag.fpt=flux.test$flag.fpt )
     # show the result
     print(new.row)

     # append the new row to the data.frame
       flux.table <- rbind(flux.table, new.row)   
     
    }#end for 
     #dev.off()
   }#end if
  
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

#write out the flux.table

#
#

plot.ld <- FALSE
if(plot.ld)
{
# pdf(file=paste("./plot_pdf/flux.table.timeseries.pdf",sep=""))
 par(mfrow=c(4,2),mai=c(0.75,0.75,0.2,0.2) )
 plot(x=flux.table$date.time, y=flux.table$lev1.hno3,type="p",xlab="HNO3 flux (nmol m-2s-1)",
      col=ifelse(flux.table$flag.n1n2=="N1","gray","black"),
      pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1) , ylim=c(-1E-3,0.01))

 plot(x=flux.table$date.time, y=flux.table$lev1.hono,type="p",xlab="HONO flux (nmol m-2s-1)",
      col=ifelse(flux.table$flag.n1n2=="N1","gray","black"),
      pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1), ylim=c(-1E-3,0.05))
 
 plot(x=flux.table$date.time, y=flux.table$lev1.n2o,type="p",xlab="N2O flux (nmol m-2s-1)", 
      col=ifelse(flux.table$flag.n1n2=="N1","gray","black"),
      pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1), ylim=c(-5E0,40) )
 
 plot(x=flux.table$date.time, y=flux.table$lev1.le.2,type="p",xlab="Latent heat flux (W m-2; J m-2 s-1)",
      col=ifelse(flux.table$flag.n1n2=="N1","gray","black"),
      pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1),ylim=c(-5E1,600))

 plot(x=flux.table$date.time, y=flux.table$lev1.no.gd,type="p",xlab="NO flux (nmol m-2s-1)", 
      col=ifelse(flux.table$flag.n1n2=="N1","gray","black"),
      pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1),ylim=c(-1E-3,0.1))
 
 plot(x=flux.table$date.time, y=flux.table$lev1.nox.gd,type="p",xlab="(NOx flux (nmolm-2s-1)",
      col=ifelse(flux.table$flag.n1n2=="N1","gray","black"),
      pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1),ylim=c(-1E-3,0.3))
 
 plot(x=flux.table$date.time, y=flux.table$lev1.ustar,type="p",
      col="forestgreen",xlab="Friction velocity (m s-1)",
      pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1),ylim=c(-1E-3,0.4))
 
 plot(x=flux.table$date.time, y=flux.table$lev1.le.1,type="p",xlab="Latent/Sensible heat flux (W m-2)",
      pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1), ylim=c(-5E1,600))
 
 lines(x=flux.table$date.time, y=flux.table$lev1.sh,col="red", 
	pch=ifelse( (flux.table$flag.itc==TRUE)&(flux.table$flag.stat==TRUE),19,1) )

#dev.off()
} 




