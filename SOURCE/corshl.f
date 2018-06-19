      subroutine corshl
     x     (idnode,mxnode,ntshl,shlke,listshl,weight,vxx,vyy,vzz,buffer)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating the internal kinetic
c     energy of core-shell units in the shell polarisation model
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith july 1994
c     
c     wl
c     2001/08/31 11:13:43
c     1.5
c     Exp
c
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension listshl(mxshl,3)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension weight(mxatms),buffer(mxbuff)
      
#ifdef VAMPIR
      call VTBEGIN(61, ierr)
#endif
      shlke=0.d0

c
c     block indices

      ishl1 = (idnode*ntshl)/mxnode+1
      ishl2 = ((idnode+1)*ntshl)/mxnode

c     
c     loop over all specified core-shell pairs
      
      m=0

      do k=ishl1,ishl2
        
        m=m+1
        
c     
c     indices of atoms involved
        
        i=listshl(m,2)
        j=listshl(m,3)
c     
c     calculate atom translational kinetic energy
        
        ppp=((weight(i)*vxx(i)+weight(j)*vxx(j))**2
     x      +(weight(i)*vyy(i)+weight(j)*vyy(j))**2
     x      +(weight(i)*vzz(i)+weight(j)*vzz(j))**2)
     x      /(weight(i)+weight(j))
c     
c     calculate individual core and shell kinetic energies
        
        ccc=weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)
        sss=weight(j)*(vxx(j)**2+vyy(j)**2+vzz(j)**2)
        
c     
c     calculate core-shell internal kinetic energy
        
        shlke=shlke+0.5d0*(ccc+sss-ppp)
        
      enddo

c     
c     global average of core-shell internal kinetic energy
        
      if(mxnode.gt.1)call gdsum(shlke,1,buffer)
      
#ifdef VAMPIR
      call VTEND(61, ierr)
#endif
      return
      end
