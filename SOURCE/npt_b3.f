      subroutine npt_b3
     x  (safe,lshmov,idnode,imcon,mxnode,natms,ntpatm,nscons,
     x  ntcons,elrc,engke,virlrc,press,taup,taut,sigma,tolnce,
     x  tstep,vircon,volm,lashap,lishap,listcon,listme,
     x  lstfrz,listot,buffer,cell,dens,dxt,dxx,dyt,dyy,dzt,dzz,
     x  fxx,fyy,fzz,prmcon,txx,tyy,tzz,uxx,uyy,uzz,vxx,vyy,vzz,
     x  weight,xdf,xxt,xxx,ydf,yyt,yyy,zdf,zzt,zzz,vx1,vy1,vz1,
     x  stress,eta,dens0)
      
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Berendsen
c     thermostat and anisotropic pressure control
c     isothermal compressibility (beta) set to that of liquid water
c     = 0.007372 dlpoly units
c     
c     parallel replicated data version
c     
c     for systems using bond CONSTRAINTS. Frozen atoms feb 1994
c     
c     copyright - daresbury laboratory 1993
c     author    -    t. forester december 1993
c     
c     wl
c     2000/01/18 14:05:43
c     1.7
c     $State Exp $
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical safe,lshmov,newjob
      
      dimension listme(mxatms),lstfrz(mxatms)
      dimension lishap(mxlshp)
      dimension lashap(mxproc),listot(mxatms)
      dimension listcon(mxcons,3),buffer(mxbuff)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension uxx(mxatms),uyy(mxatms),uzz(mxatms)
      dimension dxx(mxcons),dyy(mxcons),dzz(mxcons)
      dimension dxt(mxcons),dyt(mxcons),dzt(mxcons)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf)
      dimension weight(mxatms),prmcon(mxtcon),cell(9)
      dimension dens(mxsvdw),dens0(mxsvdw),celn(9)
      dimension vx1(msatms),vy1(msatms),vz1(msatms)
      dimension stress(9),stres1(9),stres2(9),eta(9)

      save newjob,volm0,elrc0,virlrc0

      data beta/7.3728d-3/
      data newjob/.true./
#ifndef STRESS
c     
c     this must be compiled with the stress flag on

      call error(idnode,434)
#endif

#ifdef VAMPIR
      call VTBEGIN(35, ierr)
#endif
      safe=.false.
      nbuff=mxbuff
c     
c     store initial values of volume, long range corrections etc

      if(newjob) then

        volm0 = volm
        elrc0 = elrc
        virlrc0 = virlrc
        do i = 1,ntpatm
          dens0(i) = dens(i)
        enddo
        newjob = .false.

      endif
c     
c     constraint virial

      vircon=0.d0
c     
c     temporary stress tensor accumulators
      do i = 1,9
        stres1(i) = 0.d0
        stres2(i) = 0.d0
      enddo
c     
c     set up block indices
      iatm0 = (idnode*natms)/mxnode + 1
      iatm1 = ((idnode+1)*natms)/mxnode
      
c     reciprocal of timestep
      rstep = 1.d0/tstep
      rtsq = rstep**2 
c     
c     store initial values of position and velocity

      j=0
      do i=iatm0,iatm1
        
        j=j+1
        xdf(j)=xxx(i)
        ydf(j)=yyy(i)
        zdf(j)=zzz(i)

        vx1(j) = vxx(i)
        vy1(j) = vyy(i)
        vz1(j) = vzz(i)

      enddo
c     
c     construct current bond vectors
      
      do k=1,nscons
c     
c     indices of atoms in bond
        
        i=listcon(k,2)
        j=listcon(k,3)
c     
c     calculate current bond vector
        
        dxx(k)=xxx(i)-xxx(j)
        dyy(k)=yyy(i)-yyy(j)
        dzz(k)=zzz(i)-zzz(j)
        
      enddo
c     
c     periodic boundary condition for bond vectors
      
      call images(imcon,0,1,nscons,cell,dxx,dyy,dzz)
c     
c     estimate kinetic energy and stress tensor at current timestep

      engke = 0.d0
      do i = iatm0,iatm1

        vxt = vxx(i) + 0.5d0*tstep*(fxx(i)/weight(i))
        vyt = vyy(i) + 0.5d0*tstep*(fyy(i)/weight(i))
        vzt = vzz(i) + 0.5d0*tstep*(fzz(i)/weight(i))
c     
c     2* kinetic energy at current timestep

        engke = engke + weight(i)*(vxt**2+vyt**2+vzt**2)
c     
c     kinetic contribution to stress tensor

        stres1(1) = stres1(1) + weight(i)*vxt*vxt
        stres1(2) = stres1(2) + weight(i)*vxt*vyt
        stres1(3) = stres1(3) + weight(i)*vxt*vzt
        stres1(5) = stres1(5) + weight(i)*vyt*vyt
        stres1(6) = stres1(6) + weight(i)*vyt*vzt
        stres1(9) = stres1(9) + weight(i)*vzt*vzt
c     
c     first estimate of new velocities - no thermostat

        vxx(i) = vxx(i) + tstep*(fxx(i)/weight(i))
        vyy(i) = vyy(i) + tstep*(fyy(i)/weight(i))
        vzz(i) = vzz(i) + tstep*(fzz(i)/weight(i))
c     
c     first estimate of new positions

        xxx(i) = xxx(i) + tstep*vxx(i)
        yyy(i) = yyy(i) + tstep*vyy(i)
        zzz(i) = zzz(i) + tstep*vzz(i)

      enddo
c     
c     kinetic energy

      engke =engke*0.5d0

      if(mxnode.gt.1) then 
        call gdsum(engke,1,buffer)
        call gdsum(stres1,9,buffer)
      endif
c     
c     full stress tensor

      stres1(4) = stres1(2)
      stres1(7) = stres1(3)
      stres1(8) = stres1(6)

      do i = 1,9
        stres1(i) = stres1(i) + stress(i)
      enddo
c     
c     find eta

      eta1 = beta*tstep/taup*(stres1(1)/volm - press)
      eta5 = beta*tstep/taup*(stres1(5)/volm - press)
      eta9 = beta*tstep/taup*(stres1(9)/volm - press)
      eta2 = beta*tstep/taup*(stres1(2)/volm)
      eta3 = beta*tstep/taup*(stres1(3)/volm)
      eta6 = beta*tstep/taup*(stres1(6)/volm)
c     
c     temperature scaling  coefficient - taut is the decay constant
      
      chit0 =  sqrt(1.d0 + tstep/taut*(sigma/engke -1.d0))
c     
c     begin iterations !!-----------------------------------------------

      maxit = 5
      if(ntcons.eq.0) maxit=maxit -1
      
      do iter = 1,maxit
c     
c     unconstrained new positions
        
        j = 0
        do i = iatm0,iatm1
          j = j+1
c     
c     advance velocity using leapfrog

          vxx(i)=(vx1(j) + tstep*(fxx(i)/weight(i)))*chit0
          vyy(i)=(vy1(j) + tstep*(fyy(i)/weight(i)))*chit0
          vzz(i)=(vz1(j) + tstep*(fzz(i)/weight(i)))*chit0
c     
c     update positions : 
          
          xxx(i)= tstep*vxx(i) + (1.d0+eta1)*xdf(j) + eta2*ydf(j) +
     x      eta3*zdf(j)
          yyy(i)= tstep*vyy(i) + eta2*xdf(j) + (1.d0+eta5)*ydf(j) +
     x      eta6*zdf(j)
          zzz(i)= tstep*vzz(i) + eta3*xdf(j) + eta6*ydf(j) +
     x      (1.d0+eta9)*zdf(j)

        enddo

        if(ntcons.eq.0) safe =.true.
        if(ntcons.gt.0) then
c     
c     store integrated positions

          j= 0
          do i=iatm0,iatm1
            j=j+1
            uxx(j) = xxx(i)
            uyy(j) = yyy(i)
            uzz(j) = zzz(i)
          enddo
c     
c     apply bond constraint procedures
c     
c     global exchange of configuration data
          
          if(mxnode.gt.1) then
            
            call merge(idnode,mxnode,natms,nbuff,xxx,yyy,zzz,buffer)
            
          endif
c
c     estimate new cell parameters
      
          a1 = eta1*cell(1) + eta2*cell(4) + eta3*cell(7)
          a2 = eta1*cell(2) + eta2*cell(5) + eta3*cell(8)
          a3 = eta1*cell(3) + eta2*cell(6) + eta3*cell(9)
          a4 = eta2*cell(1) + eta5*cell(4) + eta6*cell(7)
          a5 = eta2*cell(2) + eta5*cell(5) + eta6*cell(8)
          a6 = eta2*cell(3) + eta5*cell(6) + eta6*cell(9)
          a7 = eta3*cell(1) + eta6*cell(4) + eta9*cell(7)
          a8 = eta3*cell(2) + eta6*cell(5) + eta9*cell(8)
          a9 = eta3*cell(3) + eta6*cell(6) + eta9*cell(9)

          celn(1) = cell(1) + a1
          celn(2) = cell(2) + a2
          celn(3) = cell(3) + a3
          celn(4) = cell(4) + a4
          celn(5) = cell(5) + a5
          celn(6) = cell(6) + a6
          celn(7) = cell(7) + a7
          celn(8) = cell(8) + a8
          celn(9) = cell(9) + a9
c     
c     apply constraint correction
          
          call rdshake_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,viracc,lashap,lishap,listcon,listme,
     x      listot,lstfrz,buffer,celn,dxt,dxx,dyt,dyy,dzt,dzz,
     x      prmcon,txx,tyy,tzz,weight,xxt,xxx,yyt,yyy,zzt,zzz,
     x      stres2)
c     
c     contribution to constraint virial 
          
          vircon = vircon + viracc
c     
c     calculate force and velocity correction
          
          j=0
          do i=iatm0,iatm1

            j=j+1
            fxx(i) = fxx(i) + (xxx(i)-uxx(j))*weight(i)*rtsq
            fyy(i) = fyy(i) + (yyy(i)-uyy(j))*weight(i)*rtsq
            fzz(i) = fzz(i) + (zzz(i)-uzz(j))*weight(i)*rtsq

            vxx(i) = vxx(i) + (xxx(i)-uxx(j))*rstep
            vyy(i) = vyy(i) + (yyy(i)-uyy(j))*rstep
            vzz(i) = vzz(i) + (zzz(i)-uzz(j))*rstep

          enddo
c     
c     end of shake corrections
          
        endif
c     
c     calculate new kinetic energy and stress tensor at current timestep
        
        engke=0.d0
        do i = 1,9
          stres1(i) = 0.d0
        enddo

        j = 0
        do i = iatm0,iatm1

          j=j+1
          vxt = 0.5d0*(vxx(i) + vx1(j))
          vyt = 0.5d0*(vyy(i) + vy1(j))
          vzt = 0.5d0*(vzz(i) + vz1(j))
c     
c     2*kinetic energy at current timestep

          engke = engke + weight(i)*(vxt**2+vyt**2+vzt**2)
c     
c     kinetic contribution stress tensor

          stres1(1) = stres1(1) + weight(i)*vxt*vxt
          stres1(2) = stres1(2) + weight(i)*vxt*vyt
          stres1(3) = stres1(3) + weight(i)*vxt*vzt
          stres1(5) = stres1(5) + weight(i)*vyt*vyt
          stres1(6) = stres1(6) + weight(i)*vyt*vzt
          stres1(9) = stres1(9) + weight(i)*vzt*vzt

        enddo

        engke = engke*0.5d0
        if(mxnode.gt.1) call gdsum(engke,1,buffer)
c     
c     complete stress tensor - add constraint contributions

        do i = 1,9
          stres1(i) = stres1(i) + stres2(i)
        enddo

        stres1(4) = stres1(2)
        stres1(7) = stres1(3)
        stres1(8) = stres1(6)

        if(mxnode.gt.1.) call gdsum(stres1,9,buffer)
c     
c     full stress tensor

        do i = 1,9
          stres1(i) = stres1(i) + stress(i)
        enddo
c     
c     improved prediction of eta and chit 

        eta1 = beta*tstep/taup*(stres1(1)/volm - press)
        eta5 = beta*tstep/taup*(stres1(5)/volm - press)
        eta9 = beta*tstep/taup*(stres1(9)/volm - press)
        eta2 = beta*tstep/taup*(stres1(2)/volm)
        eta3 = beta*tstep/taup*(stres1(3)/volm)
        eta6 = beta*tstep/taup*(stres1(6)/volm)
        chit0 =  sqrt(1.d0 + tstep/taut*(sigma/engke -1.d0))
c     
c     end of thermostat iterations

      enddo
c     
c     update stress tensor

      do i = 1,9
        stress(i) = stres1(i)
      enddo
c     
c     update volume
      
      volm = volm*(1.d0+eta1)*(1.d0+eta5)*(1.d0+eta9)
c     
c     adjust cell vectors - anisotropic

      a1 = eta1*cell(1) + eta2*cell(4) + eta3*cell(7)
      a2 = eta1*cell(2) + eta2*cell(5) + eta3*cell(8)
      a3 = eta1*cell(3) + eta2*cell(6) + eta3*cell(9)
      a4 = eta2*cell(1) + eta5*cell(4) + eta6*cell(7)
      a5 = eta2*cell(2) + eta5*cell(5) + eta6*cell(8)
      a6 = eta2*cell(3) + eta5*cell(6) + eta6*cell(9)
      a7 = eta3*cell(1) + eta6*cell(4) + eta9*cell(7)
      a8 = eta3*cell(2) + eta6*cell(5) + eta9*cell(8)
      a9 = eta3*cell(3) + eta6*cell(6) + eta9*cell(9)

      cell(1) = cell(1) + a1
      cell(2) = cell(2) + a2
      cell(3) = cell(3) + a3
      cell(4) = cell(4) + a4
      cell(5) = cell(5) + a5
      cell(6) = cell(6) + a6
      cell(7) = cell(7) + a7
      cell(8) = cell(8) + a8
      cell(9) = cell(9) + a9
c     
c     adjust long range corrections and number density
      
      elrc = elrc0*(volm0/volm)
      virlrc = virlrc0*(volm0/volm)
      
      do kk = 1,ntpatm
        
        dens(kk) = dens0(kk)*(volm0/volm)
        
      enddo
c
c     construct scaling tensor (for later!)

      eta(1) = eta1 + 1.d0
      eta(2) = eta2
      eta(3) = eta3
      eta(4) = eta2
      eta(5) = eta5 + 1.d0
      eta(6) = eta6
      eta(7) = eta3
      eta(8) = eta6
      eta(9) = eta9 + 1.d0
c     
c     periodic boundary condition
      
      call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
c     
c     global exchange of configuration data
      
      if(mxnode.gt.1)then
        
        nbuff=mxbuff
        call merge(idnode,mxnode,natms,nbuff,xxx,yyy,zzz,buffer)
        call merge(idnode,mxnode,natms,nbuff,vxx,vyy,vzz,buffer)
        
        if(ntcons.gt.0) then
          call merge(idnode,mxnode,natms,nbuff,fxx,fyy,fzz,buffer)
        endif
        
      endif
#ifdef VAMPIR
      call VTEND(35, ierr)
#endif
      return
      end


