# Purpose load the EC 10Hz data
# First Date: 2019-03-06; Revise: 2019-07-12
# Author: Yi-Ying Chen


fun_load_ec150 <- function( dir.name="/lfs/home/ychen/lfs_dir/EC_DATA/", site.name="WuF",
			     dir.subname="2023",yyyy="2023",mm="03",dd="16",hhmm="1300",
		             ld.plot=F,ld.spec=F, dir.plot="./png_plot/", ld.na=TRUE, ld.spik=TRUE){
         
       #  site.name="WuF"; dir.plot="./png_plot/"
#	 dir.name="/lfs/home/ychen/lfs_dir/EC_DATA/"
#         dir.subname="2023"; yyyy="2023";mm="03";dd="16";hhmm="1500";ld.plot="TRUE";ld.spec="TRUE";ld.na=TRUE;ld.spik=TRUE
        
        yyyy <- formatC(yyyy, width=4, flag="0")
        mm <- formatC(mm, width=2, flag="0")
        dd <- formatC(dd, width=2, flag="0")       
        hhmm <- formatC(hhmm, width=4,flag="0")
        
        wrk.time.stamp <- as.POSIXlt(paste(yyyy,"-",mm,"-",dd, " ",  substr(hhmm,start=1,stop=2),":",substr(hhmm,start=3,stop=4),":00",sep=""),tz="")
        ch4.date.stamp <- as.POSIXlt("2022-10-25 23:30:00", tz="")

#	file.name <- paste(dir.name,"/","CSSAGRI_",yyyy,"_",mm,"_",dd,"_",hhmm,".dat",sep="")
       	file.name <- paste(dir.name,"/",site.name,"/30min/",dir.subname,"/",site.name,"_",yyyy,"_",mm,"_",dd,"_",hhmm,".dat",sep="")
        plot.name <- paste(dir.plot,"/",site.name,"_",yyyy,"_",mm,"_",dd,"_",hhmm,sep="") 

       	line.offset <- 2 
        line.skip   <- 1
        print(paste("working on file:", file.name))
     
        #try to load the file
        if (file.exists(file.name) ) { 
 	
	load.data <- read.csv(file=file.name, sep=",", stringsAsFactors=F,
		       skip=line.skip)[-(1:line.offset),]
        # check file length
	  if ( length(load.data$TIMESTAMP)!=18000) {
	     print(paste("please check file length!","Initial lines:",length(load.data$TIMESTAMP),sep=""))
	     load.data<-unique(load.data)
	     print(paste("unique the datastream !", "Final lines:", length(load.data$TIMESTAMP),sep=""))
	   }
         
      
         # Gapfill the data 
         if (ld.spik) {
            # load the source function
            source("./fun_despike.R")   
            #apply the function for despiking the raw data
             if ( wrk.time.stamp > ch4.date.stamp ) { 
             # With LI-7700 data after 2022-10-25 
             load.data$Uzc  <- fun_despike(x=as.numeric(load.data$Uzc),  nstd=c(3.5))
             load.data$Uxc  <- fun_despike(x=as.numeric(load.data$Uxc),  nstd=c(3.5))
             load.data$Uyc  <- fun_despike(x=as.numeric(load.data$Uyc),  nstd=c(3.5))
             load.data$H2O  <- fun_despike(x=as.numeric(load.data$H2O),  nstd=c(3.5))
             load.data$CO2  <- fun_despike(x=as.numeric(load.data$CO2),  nstd=c(3.5))
             load.data$Tsc  <- fun_despike(x=as.numeric(load.data$Tsc),  nstd=c(3.5))
             load.data$CH4_density  <- fun_despike(x=as.numeric(load.data$CH4_density),  nstd=c(3.5)) 
             load.data$CH4_mole_fraction  <- fun_despike(x=as.numeric(load.data$CH4_mole_fraction),  nstd=c(3.5)) 
             # cut off value for CH4
             #load.data$CH4_mole_fraction[ (load.data$CH4_mole_fraction < 0) | (load.data$CH4_mole_fraction > 10)] <- NA
             } else{
              #NO LI-7700 Data before 2022-10-25
             load.data$Uzc  <- fun_despike(x=as.numeric(load.data$Uzc),  nstd=c(3.5))
             load.data$Uxc  <- fun_despike(x=as.numeric(load.data$Uxc),  nstd=c(3.5))
             load.data$Uyc  <- fun_despike(x=as.numeric(load.data$Uyc),  nstd=c(3.5))
             load.data$H2O  <- fun_despike(x=as.numeric(load.data$H2O),  nstd=c(3.5))
             load.data$CO2  <- fun_despike(x=as.numeric(load.data$CO2),  nstd=c(3.5))
             load.data$Tsc  <- fun_despike(x=as.numeric(load.data$Tsc),  nstd=c(3.5))
             } 
          }#ld.spik
 
          # Gapfill the data 
         if (ld.na) {
            # load the source function
            source("./fun_na_fill.R")   
            #apply the function for gapfilling the raw data
            if ( wrk.time.stamp > ch4.date.stamp ) { 
            # With LI-7700 data after 2022-10-25 
            load.data$Uzc <- fun_na_fill(load.data$Uzc,  plot.ld=F, plot.name=c(plot.name))
            load.data$Uxc <- fun_na_fill(load.data$Uxc,  plot.ld=F, plot.name=c(plot.name))
            load.data$Uyc <- fun_na_fill(load.data$Uyc,  plot.ld=F, plot.name=c(plot.name))
            load.data$H2O <- fun_na_fill(load.data$H2O,  plot.ld=F, plot.name=c(plot.name))
            load.data$CO2 <- fun_na_fill(load.data$CO2,  plot.ld=F, plot.name=c(plot.name))
            load.data$Tsc <- fun_na_fill(load.data$Tsc,  plot.ld=F, plot.name=c(plot.name))
            load.data$CH4_density <- fun_na_fill(load.data$CH4_density,  plot.ld=F, plot.name=c(plot.name))
            load.data$CH4_mole_fraction <- fun_na_fill(load.data$CH4_mole_fraction,  plot.ld=F, plot.name=c(plot.name))         
            }else{
            # No LI-7700 before 2022-10-25  
            load.data$Uzc <- fun_na_fill(load.data$Uzc,  plot.ld=F, plot.name=c(plot.name))
            load.data$Uxc <- fun_na_fill(load.data$Uxc,  plot.ld=F, plot.name=c(plot.name))
            load.data$Uyc <- fun_na_fill(load.data$Uyc,  plot.ld=F, plot.name=c(plot.name))
            load.data$H2O <- fun_na_fill(load.data$H2O,  plot.ld=F, plot.name=c(plot.name))
            load.data$CO2 <- fun_na_fill(load.data$CO2,  plot.ld=F, plot.name=c(plot.name))
            load.data$Tsc <- fun_na_fill(load.data$Tsc,  plot.ld=F, plot.name=c(plot.name))
            } 
         }#ld.na
       
      
	#list.out <- list(=load.data)
	if (ld.plot) {
        png(file=paste(plot.name,"_ts",".png",sep=""), width = 6,  height = 6,  units     = "in",  res   = 300,  pointsize = 4)
        par(mfrow=c(3,2),mai = c(0.5, 0.5, 0.01, 0.01), cex.axis=2,cex.lab=2, xaxs="i",yaxs="i")	
        try(plot(na.omit(load.data$Uxc),typ="l",ylim=c(-3,3) ,col="gray" , panel.first = grid()))
        try(plot(na.omit(load.data$H2O),typ="l",col="blue", panel.first = grid()))
        try(plot(na.omit(load.data$Uyc),typ="l",ylim=c(-3,3) ,col="gray" ,panel.first = grid()))
        try(plot(na.omit(load.data$CO2),typ="l" ,col="brown", panel.first = grid() ))
        mtext( " ppm (mmol/mol) =  (mg/m2)×  (24.45/molecular weight) (g/mole) for CO2 convertion factor : 0.5556 ", side =1, line=2 )  
        try(plot(na.omit(load.data$Uzc),typ="l",ylim=c(-1,1) ,col="black", panel.first = grid()))
        plot(na.omit(load.data$CH4_mole_fraction),typ="l",col="red", panel.first = grid())
        legend("right", legend=c("U (m/s)","V (m/s)","W (m/s)", "CO2 (mg/m3)","H2O (g/m3)", "CH4 (umol/mol)  ")) 
 	mtext(paste("Date:",yyyy,"-",mm,"-",dd," ",hhmm,sep=""),side=1,line=2,col="blue")
        dev.off()
        }
       
        #check spectrum 
	if (ld.spec) {
        png(file=paste(plot.name,"_spc",".png",sep=""),  width = 6,  height = 6,  units     = "in",  res   = 300,  pointsize = 4)
        par(mfrow=c(3,2),mai = c(0.5, 0.5, 0.5, 0.5), cex.axis=2,cex.lab=2, xaxs="i",yaxs="i")	
	u.avg <- sqrt(as.numeric(load.data$Uxc)**2. + as.numeric(load.data$Uyc)**2.) 
        asp<-spectrum(na.omit(u.avg),plot=FALSE)
        try(plot(asp$spec~asp$freq, log="xy", main="U_avg Spectrum"))
        x<-asp$freq; y<- x^(-5/3)*0.1
	try(lines(y~x,col="red"))
        grid()

	asp<-spectrum(na.omit(as.numeric(load.data$H2O)),plot=FALSE)
        try(plot(asp$spec~asp$freq, log="xy", main="H2O Spectrum"))
        x<-asp$freq; y<- x^(-5/3)*0.1
	try(lines(y~x,col="red"))
        grid()

        asp<-spectrum(na.omit(as.numeric(load.data$Uzc)),plot=FALSE)
        try(plot(asp$spec~asp$freq, log="xy", main="W Spectrum")); grid()
	x<-asp$freq; y<- x^(-5/3)*0.1
	try(lines(y~x,col="red"))

	asp<-spectrum(na.omit(as.numeric(load.data$CO2)),plot=FALSE)
        try(plot(asp$spec~asp$freq, log="xy", main="CO2 Spectrum")); grid() 
        x<-asp$freq; y<- x^(-5/3)*0.1
	lines(y~x,col="red")
	
        asp<-spectrum( na.omit(u.avg*as.numeric(load.data$Uzc)) ,plot=FALSE)
        try(plot(asp$spec~asp$freq, log="xy", main="Uw Spectrum"));  grid() 
        x<-asp$freq; y<- x^(-5/3)*0.1
	try(lines(y~x,col="red"))	

	asp<-spectrum( na.omit(as.numeric(load.data$CH4_mole_fraction)),plot=FALSE)
        try(plot(asp$spec~asp$freq, log="xy", main="CH4 Spectrum")) ; grid() 
        x<-asp$freq; y<- x^(-5/3)*0.1
	try(lines(y~x,col="red"))
	mtext(paste("Date:",yyyy,"-",mm,"-",dd," ",hhmm,sep=""),side=1,line=2,col="blue")
        dev.off()
        }

	# return the datatable
	return(load.data)
        }else{
         print(paste("Try to load the ", file.name,", but the file is not exist! "), sep="")

        }
}	



