!-----Eddy Covariance Fortran Calculation Module 
!-----Contact: Yi-Ying Chen (yiyingchen@gate.sinica.edu.tw)
!-----Date: 2019-03-17; Revised 2022-03-24
!-----Research Center for Environmental Changes
!-----Coordinate rotation!
!-----SITA: Rotation about Z axis (force average V=0)
!-----PHI:  Rotation about Y axis (force average W=0)
!-----PSI:  New Y and Z Axis Rotate About axis (forcr average vw=0)
!-----SPL:  WINDOW SIZE
!-----U, V, W: Ture velocity (AFTER CORRECTION)
!-----UAVG, VAVG, WAVG: WINDOW AVERAGE WIND SPEED
!-----Commend for R share library > ecflux.so 
!     R CMD SHLIB ecflux.f90
!     or gfortran -c -fPIC ecflux.f90 | gfortran -shared -o ecflux.so ecflux.o 
!-------------------------------------------------------------------
      subroutine ecflux(spl,u,v,w,sita,phi,uavg,wd,ustar)
      implicit none
      integer :: i,j,k,spl
      real*8, external :: cov 
      real*8 :: pi
      real*8 :: vavg,wavg
      real*8 :: vwavg,vvavg,wwavg 
      real*8 :: dum1,dum2,dum3                                        
      real*8, intent(inout) :: u(spl),v(spl),w(spl)
      real*8, intent(out)   :: sita,phi,uavg,wd,ustar
!-----
      pi=4.*atan(float(1))   
!-----
      do i=1,spl
      uavg=u(i)+uavg
      vavg=v(i)+vavg      
      wavg=w(i)+wavg     
!      write(*,*) 'raw data:',ta(i),qa(i),co2(i) 
      end do
!-----
      uavg=uavg/spl
      vavg=vavg/spl
      wavg=wavg/spl
!     write(6,'(3(A5,F8.3))') 'UAVG=',uavg,'VAVG=',vavg,'WAVG=',wavg
!-----      
      if(uavg .eq. 0.)then
            sita=pi/2.
      else 
            sita=atan(abs(vavg/uavg))
      end if 
!-----
      if(sqrt(uavg**2.+vavg**2.) .eq. 0.) then
        phi=pi/2.
        else
        phi=atan(wavg/sqrt(uavg**2.+vavg**2.))
      end if     
!     write(*,'(A9,A6,F8.2,A5,F8.2)') 'ROTATION:','PHI:',(180/PI)*PHI,'SITA:',(180/PI)*SITA
      do i=1,spl
         dum1=u(i) 
         dum2=u(i)*cos(sita)+v(i)*sin(sita)
         u(i)=(u(i)*cos(sita)+v(i)*sin(sita))*cos(phi)+w(i)*sin(phi)  
         v(i)=v(i)*cos(sita)-dum1*sin(sita)
         w(i)=w(i)*cos(phi)-dum2*sin(phi)
      END DO
!-----
      do i=1,spl
      !uavg=u(i)+uavg
      !vavg=v(i)+vavg      
      wavg=w(i)+wavg     
      end do
!-----
      uavg=sqrt(uavg*uavg + vavg*vavg)
      !vavg=vavg/spl
      wavg=wavg/spl
!-----Calculation of ustar 
      ustar=(cov(w,u,spl)**2.+cov(w,v,spl)**2.)**.25 
!-----Calculate surface fluxes 
!-----write output information 
!     write(6,'(3(A5,F8.3))') 'UAVG=',uavg,'VAVG=',vavg,'WAVG=',wavg
!     write(6,*) 'Rotation Angle About Axes:'
!     write(6,*)'Z-Axes=',(180/pi)*sita,'Y-Axes=',(180/pi)*phi
!-----Calculate wind direction
      wd = (180/pi)*sita 
!-----Offset to wind camposs define N=0      
      if (( uavg .ge. 0.) .and. ( vavg .gt. 0. ))  wd = 270. - wd
      if (( uavg .ge. 0.) .and. ( vavg .lt. 0. ))  wd = 270. + abs(wd) 
      if (( uavg .le. 0.) .and. ( vavg .lt. 0. ))  wd = 180. - wd
      if (( uavg .le. 0.) .and. ( vavg .gt. 0. ))  wd = 180. + abs(wd) 
      if ( wd .lt. 0) wd = wd + 360
!     write(6,*)'Wind_dir=', wd
      end subroutine ecflux   
!==== End of ecflux subroutine  ====

!.....Covariance funtion
      real*8 function cov(var1, var2, n)
      implicit none
      real*8, intent(in)    :: var1(n), var2(n)
      integer, intent(in) :: n
      real*8, external :: avg
      real*8 :: dev, acc, avg1, avg2 
      integer :: i
      acc=0.
      dev=0.
      avg1=0.
      avg2=0.
      avg1=avg(var1, n)
      avg2=avg(var2, n)
      do i=1, n
         dev=(var1(i)-avg1)*(var2(i)-avg2)
         acc=dev+acc
      end do
      cov=acc/float(n)      
      return
      end function cov         
!.....Average function    
      real*8 function avg(var, n)
      implicit none
      integer, intent(in) :: n 
      real*8,  intent(in) :: var(n)     
      avg=0.
      avg=sum(var)/float(n)
      return
      end function avg    
!.....Standard deviation funtion
      real*8 function std(var1, n)
      implicit none      
      integer, intent(in) :: n
      real*8,  intent(in) :: var1(n)
      real*8, external :: avg
      real*8           :: acc, dev, avg1       
      integer :: i
      acc=0.
      dev=0.
      avg1=0.
      avg1=avg(var1, n)
      do i=1, n
         dev=(var1(i)-avg1)**2.
         acc=acc+dev
      end do
      acc=acc/float(n)
      std=sqrt(acc)      
      return
      end function std
!=========END OF THE PROGRAM ==========  
