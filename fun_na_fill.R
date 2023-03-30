#############################################
# A simple forward gradient approach to gapfill 
# nan values in the ECraw data 
# date: 2023-03-17
# author: yi-ying chen 
##############################################
        fun_na_fill <- function (xx, plot.ld=TRUE, plot.name="gap.plot" ) {
          
          pp <- plot.ld
       # find NAN index
            xx <- as.numeric(xx)
            tmp <- xx

         # id no data exit the function 
        if (length(tmp) > 0) { 

         nan.id <-  which(is.na(xx))
         # any nan values 
         if ( (length(nan.id > 0)) ) {
           #print(nan.id)
          if (nan.id[1] == 1) {
               xx[nan.id[1]] <- mean(xx, na.rm=T)
               xx[nan.id[2]] <- mean(xx, na.rm=T) 
               #fixed the first two id with mean
               nan.id <- nan.id[-1:-2]
             }
             #print(nan.id)
             #using the forward gradient to fill the NAN
             xx[nan.id] <- (xx[nan.id-1]-xx[nan.id-2]) + xx[nan.id-1]
              # plot the gapfillied data
             if (pp) {
                png(file=paste(plot.name,"_na_gf.png",sep=""))
                try(plot(tmp, type="l",col="gray"))
                points(x=which(is.na(tmp)), y=xx[which(is.na(tmp))],col="red")
                dev.off()
             }
            print(paste("NaN was gapfilled !"))
            # return the gapfilled data
            } else{
            print(paste("No NaN was found in the dataset, return to the procedure!",sep=""))
            
            }
            return(xx)
  
         } else { 
           print("No Data!") 
         }
 
      }

  
