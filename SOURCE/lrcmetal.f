      subroutine lrcmetal
     x  (idnode,imcon,mxnode,natms,ntpatm,engunit,rcut,volm,
     x  lstvdw,ltpvdw,ltype,numtyp,prmvdw,dens,elrcm,vlrcm)
      
c*************************************************************************
c     
c     DL_POLY subroutine to evaluate long-range corrections to
c     pressure and energy in a periodic metal system.
c     
c     copyright daresbury laboratory 1995
c     
c     author -       w. smith   june 1995
c     
c     wl
c     2002/08/01 14:46:26
c     1.4
c     Exp
c     
c***************************************************************************
      
#include "dl_params.inc"
      
      logical newjob
      dimension numtyp(mxsvdw),ltype(mxatms)
      dimension prmvdw(mxvdw,mxpvdw),ltpvdw(mxvdw),lstvdw(mxvdw)
      dimension dens(mxsvdw),elrcm(0:mxsvdw),vlrcm(0:mxsvdw)
      
      save newjob
      data newjob/.true./
c     
c     sutton-chen potentials
      
      vv1(r,a,b,c)=(a*b**3*(b/r)**(c-3.d0))/(c-3.d0)
      gg1(r,a,b,c)=(a*c*b**3*(b/r)**(c-3.d0))/(c-3.d0)
      
#ifdef VAMPIR
      call VTBEGIN(147, ierr)
#endif
      twopi = 2.0d0*pi
      forpi = 4.0d0*pi
c     
c     initalise counter arrays
      
      do i = 1,mxsvdw
        
        numtyp(i) = 0
        
      enddo
c     
c     evaluate species populations in system
      
      do i = 1, natms
        
        ka = ltype(i)
        numtyp(ka) = numtyp(ka)+1
        
      enddo
      
c     
c     number densities
      
      do i = 1,ntpatm
        
        dens(i) = dble(numtyp(i))/volm
        
      enddo
      
c     
c     long range corrections to energy and pressure
      
      do i=0,mxsvdw

        elrcm(i) = 0.d0
        vlrcm(i) = 0.d0

      enddo
      
      if(imcon.ne.0.and.imcon.ne.6) then
        
        kvdw = 0
        elrcsum=0.d0
        
        do i = 1, ntpatm
          
          do j = 1,i
            
            elrc0=0.d0
            elrc1=0.d0
            elrc2=0.d0
            vlrc0=0.d0
            vlrc1=0.d0
            vlrc2=0.d0
            
            kvdw = kvdw + 1
            k0 = lstvdw(kvdw)
            
            if(ltpvdw(k0).eq.100) then
              
              elrc0 = vv1(rcut,prmvdw(k0,1),prmvdw(k0,2),prmvdw(k0,3))
              vlrc0 = gg1(rcut,prmvdw(k0,1),prmvdw(k0,2),prmvdw(k0,3))
              
              if(i.eq.j) then
                
                elrc1 = vv1(rcut,1.d0,prmvdw(k0,2),
     x            prmvdw(k0,4))*(prmvdw(k0,1)*prmvdw(k0,5))**2
                elrcm(i)=elrcm(i)+forpi*dens(i)*elrc1
                elrcsum=elrcsum+twopi*volm*dens(i)**2*elrc1
                vlrc1 = gg1(rcut,1.d0,prmvdw(k0,2),
     x            prmvdw(k0,4))*(prmvdw(k0,1)*prmvdw(k0,5))**2
                vlrcm(i)=vlrcm(i)+twopi*dens(i)*vlrc1

              else

                k1=lstvdw((i*(i+1))/2)
                k2=lstvdw((j*(j+1))/2)
                elrc1 = vv1(rcut,1.d0,prmvdw(k0,2),
     x            prmvdw(k0,4))*(prmvdw(k1,1)*prmvdw(k1,5))**2
                elrc2 = vv1(rcut,1.d0,prmvdw(k0,2),
     x            prmvdw(k0,4))*(prmvdw(k2,1)*prmvdw(k2,5))**2
                elrcm(i)=elrcm(i)+forpi*dens(j)*elrc1
                elrcm(j)=elrcm(j)+forpi*dens(i)*elrc2
                elrcsum=elrcsum+twopi*volm*dens(i)*dens(j)*(elrc1+elrc2)
                vlrc1 = gg1(rcut,1.d0,prmvdw(k0,2),
     x            prmvdw(k0,4))*(prmvdw(k1,1)*prmvdw(k1,5))**2
                vlrc2 = gg1(rcut,1.d0,prmvdw(k0,2),
     x            prmvdw(k0,4))*(prmvdw(k2,1)*prmvdw(k2,5))**2
                vlrcm(i)=vlrcm(i)+twopi*dens(j)*vlrc1
                vlrcm(j)=vlrcm(j)+twopi*dens(i)*vlrc2

              endif
              
            endif
            
            if(i.ne.j) then
              elrc0 = elrc0*2.d0
              vlrc0 = vlrc0*2.d0
            endif
            
            elrcm(0) = elrcm(0) + twopi*volm*dens(i)*dens(j)*elrc0
            vlrcm(0) = vlrcm(0) - twopi*volm*dens(i)*dens(j)*vlrc0
            
          enddo
          
        enddo
        
      endif

      if(newjob) then
        newjob =.false.
        if(idnode.eq.0) write(nrite,
     x    "(/,/,1x,1p,
     x    'long range correction to metal energy    ',e15.6,/,
     x    1x,'lr correction for metal atom density     ',e15.6,/,
     x    1x,'1st partial lr correction to metal virial',e15.6)")
     x    elrcm(0)/engunit,elrcsum/engunit**2,
     x    vlrcm(0)/engunit
        
      endif
      
#ifdef VAMPIR
      call VTEND(147, ierr)
#endif
      return
      end


