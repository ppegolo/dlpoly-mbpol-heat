      subroutine diffsn1
     x  (idnode,natms,ntpatm,mxnode,ltype,numtyp,amsd,
     x  xx0,yy0,zz0,buffer)
      
c***********************************************************************
c     
c     DL_POLY routine for calculating mean squared displacements
c     
c     displacements calculated in diffsn0
c     
c     parallel version - replicated data.
c     
c     copyright daresbury laboratory 1993
c     
c     author - t. forester      june 1993
c     
c     wl
c     2001/05/30 12:40:01
c     1.2
c     Exp
c     
c***********************************************************************
      
#include "dl_params.inc"

      logical newjob

      dimension xx0(mxatms),yy0(mxatms),zz0(mxatms)
      dimension buffer(mxbuff)
      dimension amsd(mxsvdw)
      dimension numtyp(mxsvdw),ltype(mxatms)

      save newjob,iatm1,iatm2

      data newjob/.true./

#ifdef VAMPIR
      call VTBEGIN(133, ierr)
#endif
      if(newjob) then

        newjob=.false.
        iatm1 = (idnode*natms)/mxnode+1
        iatm2 = ((idnode+1)*natms)/mxnode

      endif
c     
c     running sum of squared displacements
      
      do k = 1,ntpatm
        
        amsd(k) = 0.d0
        
      enddo
      
c     
c     calculate square of displacements for each atom type

      do i = iatm1,iatm2
        
        k = ltype(i)
        
        amsd(k)= amsd(k) + xx0(i)**2 + yy0(i)**2 + zz0(i)**2
        
      enddo
      
c     
c     global sum - replicated data strategy
      
      if(mxnode.gt.1) then
        
        do k = 1,ntpatm
          
          buffer(k+ntpatm) = amsd(k)
          
        enddo
        
        call  gdsum(buffer(1+ntpatm),ntpatm,buffer(1))
        
        do k = 1,ntpatm
          
          amsd(k) = buffer(k+ntpatm)
          
        enddo
        
      endif
      
c     
c     mean squared displacement
      
      do k = 1,ntpatm
        
        amsd(k) = amsd(k) / dble(max(numtyp(k),1))
        
      enddo
      
#ifdef VAMPIR
      call VTEND(133, ierr)
#endif
      return
      end


