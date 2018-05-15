      subroutine zden1
     x  (lpgr,cfgname,unqatm,idnode,mxnode,ntpatm,nzden,
     x  volm,zlen,zdens,buffer)
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating Z density profile
c     from accumulated data.
c     double precision version
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1994
c     author    - t. forester    march 1994    
c     
c     wl
c     2001/05/30 12:40:28
c     1.5
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical lpgr
      
      character*80 cfgname
      character*8 unqatm(mxsite)
      
      dimension zdens(mxrdf,mxsvdw)
      dimension buffer(mxbuff)
      
#ifdef VAMPIR
      call VTBEGIN(168, ierr)
#endif
      if(idnode.eq.0) write(nrite,
     x  "(/,/,12X,'Z DENSITY PROFILES',/,/,
     x  'calculated using ',i10,' configurations')") nzden
      
      if(lpgr) then

c     open Z density file and write headers

        if(idnode.eq.0)then

          open(nzdndt,file='ZDNDAT')

          write(nzdndt,'(a80)')cfgname
          write(nzdndt,'(i10)')mxrdf

        endif
c     
c     grid interval for density profiles
        
        delr = zlen/dble(mxrdf)
c     
c     volume of z strip

        dvolz=(volm/zlen)*delr
c     
c     normalisation factor 
        
        nzden = max(nzden,1)
        factor = 1.d0/(dble(nzden)*dvolz)

        do k = 1,ntpatm
          
          if(idnode.eq.0) then

             write(nrite,
     x      "(/,'rho(r)  :',a8,/,/,8x,'r',6x,'rho',9x,'n(r)',/)")
     x      unqatm(k)
             write(nzdndt,'(a8)')unqatm(k)

          endif
c     
c     global sum of data on all nodes
          
          if(mxnode.gt.1) call gdsum(zdens(1,k),mxrdf,buffer)
c     
c     running integration of z-density
          
          sum = 0.d0
c     
c     loop over distances
          
          do j = 1,mxrdf
            
            rrr = (dble(j)-0.5d0)*delr - zlen*0.5d0
            rho = zdens(j,k)*factor
            sum = sum + rho*dvolz
c     
c     print out information
            
            if(idnode.eq.0) then

              write(nrite,"(f10.4,1p,2e14.6)") rrr,rho,sum
              write(nzdndt,"(1p,2e14.6)") rrr,rho

            endif
            
          enddo
          
        enddo
        
        if(idnode.eq.0)close (nzdndt)

      endif
      
#ifdef VAMPIR
      call VTEND(168, ierr)
#endif
      return
      end
