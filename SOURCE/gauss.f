      subroutine gauss(natms,vxx,vyy,vzz)
c     
c*********************************************************************
c     
c     dl_poly subroutine for constructing velocity arrays
c     with a gaussian distribution of unit variance.
c     
c     based on the method described by Allen and Tildesley in
c     Computer Simulation of Liquids , Clarendon Press 1987 P347
c     
c     note - this version uses a universal random number 
c     generator, which generates pseudo-random numbers between
c     0 and 1. it is based on the algorithm of marsaglia, zaman
c     and tsang in: stats and prob. lett. 8 (1990) 35-39.
c     
c     copyright daresbury laboratory 1992
c     author - w. smith         july 1992
c     
c     wl
c     2001/05/30 12:40:07
c     1.2
c     Exp
c     
c*********************************************************************
c     
      
      implicit real*8(a-h,o-z)
      
      dimension vxx(natms),vyy(natms),vzz(natms)
      
      data a1,a3,a5/3.949846138d0,0.252408784d0,0.076542912d0/
      data a7,a9/0.008355968d0,0.029899776d0/
      

#ifdef VAMPIR
      call VTBEGIN(144, ierr)
#endif
c     
c     initialise random number generator

      random = duni()

      do i=1,natms
        
        rrr=(duni()+duni()+duni()+duni()+duni()+duni()
     x    +duni()+duni()+duni()+duni()+duni()+duni()
     x    -6.d0)/4.d0
        rr2=rrr*rrr
        vxx(i)=rrr*(a1+rr2*(a3+rr2*(a5+rr2*(a7+rr2*a9))))
        
        rrr=(duni()+duni()+duni()+duni()+duni()+duni()
     x    +duni()+duni()+duni()+duni()+duni()+duni()
     x    -6.d0)/4.d0
        rr2=rrr*rrr
        vyy(i)=rrr*(a1+rr2*(a3+rr2*(a5+rr2*(a7+rr2*a9))))
        
        rrr=(duni()+duni()+duni()+duni()+duni()+duni()
     x    +duni()+duni()+duni()+duni()+duni()+duni()
     x    -6.d0)/4.d0
        rr2=rrr*rrr
        vzz(i)=rrr*(a1+rr2*(a3+rr2*(a5+rr2*(a7+rr2*a9))))
        
      enddo
      
#ifdef VAMPIR
      call VTEND(144, ierr)
#endif
      return
      end
