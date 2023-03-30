
# R funtiom flux.cal() 
# 
fun_ecflux_for <- function(raw.data=raw.data ) {
# load the loading funtion load.data()
# Load/make fortran subroutine rot()
# print(system("R CMD SHLIB rot.f90",intern=TRUE))
#
dyn.load("ecflux.so",type="Fortran")
print(is.loaded("ecflux")) 

# load the sharing fortran library for the data processing/ calcualtion 
# the double rotation, the detail can be found in fortran code ROT.f90 
#!-----SITA: Rotation about Z axis (force average V=0)
#!-----PHI:  Rotation about Y axis (force average W=0)
#!-----SPL:  WINDOW SIZE
#!-----U, V, W: Ture velocity (AFTER CORRECTION)
#!-----UAVG, VAVG, WAVG: WINDOW AVERAGE WIND SPEED
#--------------------------------------------------------------------
# Create share library ---> ecflux.so
# by following command 
#> system(R CMD SHLIB ecflux.f90)
# check funtion names in shared library
#> system("nm -D ecflux.so") show function "names" in ecflux.so
#--------------------------------------------------------------------

#        Rhoa(j)=(Pa*1E5)/(287.05*(273.15+T(j)))           ! unit in [g/m3]      dry air     density     Assume T <- sonic=T <- dry <- air
#           T(j)=T(j)*(1-0.514*(Rhov(j)/Rhoa(j)))          !-----Change Sonic Temp to Air Temp

#convert sonic temperature to air temperature 
#raw.data$Rhoa <- (as.numeric(raw.data$cell_press)*1E6)/(287.05*(273.15+as.numeric(raw.data$Tsc)))
#raw.data$Tac <- as.numeric(raw.data$Tsc) * (1- 0.514*(as.numeric(raw.data$H2O)/raw.data$Rhoa))
#raw.data$Rhoa <- (100.0*1E6)/(287.05*(273.15+as.numeric(raw.data$Tsc)))
#raw.data$Tac <- as.numeric(raw.data$Tsc) #* (1- 0.514*(as.numeric(raw.data$H2O)/raw.data$Rhoa))


#air.temp <- mean(raw.data$Tac,na.rm=TRUE)
#air.den  <- mean(raw.data$Rhoa, na.rm=TRUE)
#air.press <- mean(as.numeric(raw.data$cell_press), na.rm=TRUE)


#level for IRGASON(EC150+CSAT3) + LI7700 
wind.rot <- .Fortran("ecflux",
		     spl = as.integer(length(raw.data$TIMESTAMP[!is.na(raw.data$Uxc)])),
	             u = as.double(raw.data$Uxc[!is.na(raw.data$Uxc)]),
	             v = as.double(raw.data$Uyc[!is.na(raw.data$Uxc)]),
	             w = as.double(raw.data$Uzc[!is.na(raw.data$Uxc)]),
	             sita = double(1), phi = double(1),
		     uavg = double(1), wd  = double(1), 
		     ustar= double(1) )

return(list(flux.ec=wind.rot))
}


