

fun_despike <- function(x, nstd=3.5, gap.fill="linear")
{      
      
      #find the index of spikes 
      x <- as.numeric(x)
      #check file length
      #dev.new()
 
      #plot(x,type="l")
      if (length(x) > 0 ) {
        
      n <- length(x)
      # divided into n*sub window 
      sub_n=1

      for (isub in 1:sub_n) {
      
      ix1 <- (n/sub_n)*(isub-1) + 1 
      ix2 <- (n/sub_n)*isub

      #print(ix1)
      #print(ix2)
     
      xx <- x[ix1:ix2]
      xx.mean <- mean(xx,na.rm=T)
      xx.std <- sd(xx,na.rm=T) 

      spike.id <- which ( abs(xx- xx.mean) >= abs(xx.std*nstd) ) 
     
     # print(paste("mean:",xx.mean))
     # print(paste("std:",xx.std))

      #skipe first and the final elenments

      spike.id <- subset(spike.id, (spike.id!=1 & spike.id!=length(xx))) 

      #filled the spike value

      if (is.na(length(spike.id)))  
	  {
           print(paste("No spike found!"))
	   xx <- xx
	  }
          else
	  {	  
           for (i in spike.id ) 
           {
           #print(paste("id:",i,sep=""))
           #print(paste("spik value: ",xx[i],sep=""))
           if (i < n) {
              xx[i] <- ( xx[i-1] +  xx[i+1])/2.
           #print(paste("replace as:",xx[i],sep=""))
           }else{
            xx[i] <- xx.mean
           }
 
      	   }
          print(paste("percentage of spikes:", format(100.0*(length(spike.id)/length(xx)), digits=2,width=5),"%",sep=""))
      #
          }	  
      #return the vector
      
      x[ix1:ix2] <- xx
      #points(x=spike.id,y=x[spike.id],col="red")

      }

      return(as.numeric(x))
     } else {
       print( "No Data")
     }
 
}



