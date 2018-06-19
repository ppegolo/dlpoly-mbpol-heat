      subroutine ewald4cp
     x   (iatm,ilst,engcpe,vircpe,alpha,epsq,nexatm2,
     x   lexatm2,chge,xdf,ydf,zdf,fxx,fyy,fzz,stress,
     x   dipx,dipy,dipz,rcut,efieldkx,efieldky,efieldkz)
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces in a
c     periodic system using ewald's method. This is the subrotine
c     for correction of the excluded pairs.
c     
c     parallel replicated data version (part 3)
c     
c     copyright - daresbury laboratory 1992
c                 voth group
c     author    - w. smith dec 1992.
c                 c. j. burnham  sept 2003.
c                 t. yan dec 2003
c     
c     stress stensor added t.forester may 1994
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension nexatm2(msatms),lexatm2(msatms,mxexcl)
      dimension chge(mxatms)
      dimension efieldkx(mxatms),efieldky(mxatms),efieldkz(mxatms)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension stress(9)

#ifdef VAMPIR
      call VTBEGIN(103, ierr)
#endif

c     
c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0
c     
c     start of primary loop for forces evaluation
      
      do m=1,nexatm2(ilst)
c     
c     atomic index and charge product
        
        jatm=lexatm2(ilst,m)
c     
c     calculate interatomic distance
        
        rsq=xdf(m)**2+ydf(m)**2+zdf(m)**2
        rrr = dsqrt(rsq)

ccc        if (rcut.gt.rrr) then

           r3i = 1.d0/(rrr**3)
           r5i = r3i/rsq
c     exclude intramolecular charge field

           eij_ki = -r4pie0*chge(jatm)*r3i/epsq
           eij_kj = -r4pie0*chge(iatm)*r3i/epsq

           efieldkx(iatm) = efieldkx(iatm) + eij_ki*xdf(m)
           efieldky(iatm) = efieldky(iatm) + eij_ki*ydf(m)
           efieldkz(iatm) = efieldkz(iatm) + eij_ki*zdf(m)

           efieldkx(jatm) = efieldkx(jatm) - eij_kj*xdf(m)
           efieldky(jatm) = efieldky(jatm) - eij_kj*ydf(m)
           efieldkz(jatm) = efieldkz(jatm) - eij_kj*zdf(m)
c
c
c     exclude charge-dipole interactions

           rdmui = dipx(iatm)*xdf(m)+dipy(iatm)*ydf(m)+
     x             dipz(iatm)*zdf(m)
           rdmuj = dipx(jatm)*xdf(m)+dipy(jatm)*ydf(m)+
     x             dipz(jatm)*zdf(m)

           facmu = (chge(jatm)*rdmui-chge(iatm)*rdmuj)
     x             *3.d0*r5i*r4pie0/epsq

           dfacmui = -chge(jatm)*r3i*r4pie0/epsq
           dfacmuj =  chge(iatm)*r3i*r4pie0/epsq

           dfacx = dfacmui*dipx(iatm)+dfacmuj*dipx(jatm)
           dfacy = dfacmui*dipy(iatm)+dfacmuj*dipy(jatm)
           dfacz = dfacmui*dipz(iatm)+dfacmuj*dipz(jatm)
c
c     increment forces

           ffac = facmu

           virdum = dfacx*xdf(m)+dfacy*ydf(m)+dfacz*zdf(m)
           vircpe = vircpe - ffac*rsq - virdum

           fx = ffac*xdf(m)+dfacx
           fy = ffac*ydf(m)+dfacy
           fz = ffac*zdf(m)+dfacz


           fxx(iatm) = fxx(iatm)+fx
           fyy(iatm) = fyy(iatm)+fy
           fzz(iatm) = fzz(iatm)+fz

           fxx(jatm) = fxx(jatm)-fx
           fyy(jatm) = fyy(jatm)-fy
           fzz(jatm) = fzz(jatm)-fz
        
#ifdef STRESS
c     
c     calculate stress tensor
        
           stress(1) = stress(1)+ xdf(m)*fx
           stress(2) = stress(2)+ xdf(m)*fy
           stress(3) = stress(3)+ xdf(m)*fz
        
           stress(5) = stress(5)+ ydf(m)*fy
           stress(6) = stress(6)+ ydf(m)*fz
        
           stress(9) = stress(9)+ zdf(m)*fz

#endif        

ccc        endif

      enddo
#ifdef STRESS
c     
c     complete stress tensor
      
      stress(4) = stress(2)
      stress(7) = stress(3)
      stress(8) = stress(6)
#endif              
      
#ifdef VAMPIR
      call VTEND(103, ierr)
#endif
      return
      end



