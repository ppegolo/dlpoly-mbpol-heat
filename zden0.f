      subroutine zden0
     x  (idnode,natms,mxnode,nzden,zlen,ltype,zzz,zdens)
c     
c***********************************************************************
c     
c     dl_poly subroutine for accumulating statistic for density profile
c     zlen = length of cell in z direction
c     
c     double precision accumulators
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1994
c     author    - t. forester    march 1994    
c     
c     wl
c     2000/01/18 14:06:01
c     1.3
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension zzz(mxatms)
      dimension zdens(mxrdf,mxsvdw)
      dimension ltype(mxatms)
#ifdef VAMPIR
      call VTBEGIN(67, ierr)
#endif
c     
c     accumulator

      nzden = nzden+1
c     
c     half of z length

      zleno2 = zlen*0.5d0
c     
c     grid interval for density profiles

      rdelr = dble(mxrdf)/zlen
c     
c     set up atom iatm type

      do iatm=idnode+1,natms,mxnode

        k =ltype(iatm) 

        ll=int((zzz(iatm)+zleno2)*rdelr + 1.0d0)
        ll = max(1,ll)
        ll = min(mxrdf,ll)
c     
c     accumulate statistic

        zdens(ll,k) = zdens(ll,k) + 1.d0

      enddo
      
#ifdef VAMPIR
      call VTEND(67, ierr)
#endif
      return
      end
