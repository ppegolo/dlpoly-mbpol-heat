      subroutine npt_h0
     x  (idnode,imcon,mxnode,natms,ntpatm,chip,chit,conint,consv,elrc,
     x  engke,virlrc,press,taup,taut,sigma,temp,tstep,virtot,volm,
     x  buffer,cell,dens,fxx,fyy,fzz,vxx,vyy,vzz,weight,xdf,
     x  xxx,ydf,yyy,zdf,zzz,vx1,vy1,vz1,stress,eta,dens0)
      
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Hoover 
c     thermostat + piston.
c     
c     reference: Melchionna, Ciccotti and Holian,
c     Mol Phys 1993, 78, p533
c     
c     parallel replicated data version
c     
c     copyright daresbury laboratory 1995/1999
c     author    -    s. melchionna   april 1995
c     and       -    t. forester     april 1995
c     adapted   -    w.smith for atomic systems oct 1999
c     
c     wl
c     2003/05/08 08:45:11
c     1.3
c     $Sate: Exp $
c     
c***********************************************************************
      
#include "dl_params.inc"
      
      logical newjob
      
      dimension buffer(mxbuff)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf)
      dimension weight(mxatms),cell(9),cell0(9)
      dimension dens(mxsvdw),dens0(mxsvdw),celn(9)
      dimension vx1(msatms),vy1(msatms),vz1(msatms)
      dimension stress(9),stres1(9),eta(9)
      
      save newjob,volm0,elrc0,virlrc0,cell0

      data newjob/.true./

#ifdef VAMPIR
      call VTBEGIN(36, ierr)
#endif
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
        do i = 1,9
          cell0(i) = cell(i)
        enddo
        newjob = .false.

      endif
#ifdef STRESS
c     
c     temporary stress tensor accumulators
      do i = 1,9
        stres1(i) = 0.d0
      enddo
#endif
c     
c     set up block indices
      iatm0 = (idnode*natms)/mxnode + 1
      iatm1 = ((idnode+1)*natms)/mxnode
c      
c     reciprocal of timestep
      rstep = 1.d0/tstep
      rtsq = rstep**2 
c     
c     inertia parameter for Nose-Hoover thermostat
      
      pmass = natms*boltz*temp*taup**2
c     
c     store initial values of position and velocity
      
      totmas = 0.d0
      j=0
      do i=iatm0,iatm1
        
        j=j+1
        xdf(j)=xxx(i)
        ydf(j)=yyy(i)
        zdf(j)=zzz(i)
        totmas=totmas+weight(i)

        vx1(j) = vxx(i)
        vy1(j) = vyy(i)
        vz1(j) = vzz(i)

      enddo

      if(mxnode.gt.1) call gdsum(totmas,1,buffer)
c     
c     calculate centre of mass

      xcmo = 0.d0
      ycmo = 0.d0
      zcmo = 0.d0
      
      j=0
      do i=iatm0,iatm1
        
        j=j+1
        xcmo = xcmo + weight(i)*xxx(i)/totmas
        ycmo = ycmo + weight(i)*yyy(i)/totmas
        zcmo = zcmo + weight(i)*zzz(i)/totmas

      enddo

      if(mxnode.gt.1) then

        buffer(4) = xcmo
        buffer(5) = ycmo
        buffer(6) = zcmo

        call gdsum(buffer(4),3,buffer(1))

        xcmo = buffer(4)
        ycmo = buffer(5)
        zcmo = buffer(6)

      endif
c     
c     estimate kinetic energy at current timestep

      engke = 0.d0
      do i = iatm0,iatm1

        vxt = vxx(i) + 0.5d0*tstep*(fxx(i)/weight(i))
        vyt = vyy(i) + 0.5d0*tstep*(fyy(i)/weight(i))
        vzt = vzz(i) + 0.5d0*tstep*(fzz(i)/weight(i))
c     
c     kinetic energy at current timestep

        engke = engke + weight(i)*(vxt**2+vyt**2+vzt**2)
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
      if(mxnode.gt.1) call gdsum(engke,1,buffer)
c     
c     propagate chip

      chipp = ((2.d0*engke-virtot)/3.d0 - press*volm)/pmass
      chipnew = chip + tstep*chipp
      chip0 = 0.5d0*(chip+chipnew)
c     
c     propagate chit

      chitp = (engke/sigma - 1.d0)/taut**2
      chitnew = chit + tstep*chitp
      chit0 = 0.5d0*(chit+chitnew)
c     
c     begin iterations !!-----------------------------------------------

      maxit = 4
      
      do iter = 1,maxit
c     
c     unconstrained new positions
        
        j = 0
        do i = iatm0,iatm1
          j = j+1
c     
c     advance velocity using leapfrog

          vxx(i) = vx1(j) + tstep*(fxx(i)/weight(i) -(chit0+chip0)*
     x      0.5d0*(vxx(i)+vx1(j)))
          vyy(i) = vy1(j) + tstep*(fyy(i)/weight(i) -(chit0+chip0)*
     x      0.5d0*(vyy(i)+vy1(j)))
          vzz(i) = vz1(j) + tstep*(fzz(i)/weight(i) -(chit0+chip0)*
     x      0.5d0*(vzz(i)+vz1(j)))
c     
c     advance positions using leapfrog

          xxx(i) = xdf(j) + tstep*(vxx(i) +
     x      chipnew*((xxx(i)+xdf(j))*0.5d0 - xcmo))
          yyy(i) = ydf(j) + tstep*(vyy(i) +
     x      chipnew*((yyy(i)+ydf(j))*0.5d0 - ycmo))
          zzz(i) = zdf(j) + tstep*(vzz(i) +
     x      chipnew*((zzz(i)+zdf(j))*0.5d0 - zcmo))

        enddo
c     
c     calculate new kinetic energy at current timestep
        
        engke=0.d0

        j = 0
        do i = iatm0,iatm1

          j=j+1
          vxt = 0.5d0*(vxx(i) + vx1(j))
          vyt = 0.5d0*(vyy(i) + vy1(j))
          vzt = 0.5d0*(vzz(i) + vz1(j))
c     
c     kinetic energy at current timestep

          engke = engke + weight(i)*(vxt**2+vyt**2+vzt**2)

        enddo

        engke = engke*0.5d0
        if(mxnode.gt.1) call gdsum(engke,1,buffer)
c     
c     improved prediction of chip and chit 

        chipp = ((2.d0*engke-virtot)/3.d0 - press*volm)/pmass
        chipnew = chip + tstep*chipp
        chip0 = 0.5d0*(chip+chipnew)

        chitp = (engke/sigma - 1.d0)/taut**2
        chitnew = chit + tstep*chitp
        chit0 = 0.5d0*(chit+chitnew)
c     
c     end of thermostat iterations

      enddo
#ifdef STRESS
c     
c     kinetic contribution to stress tensor
      j = 0
      do i = iatm0,iatm1

        j = j+1
        vxt = 0.5d0*(vxx(i)+vx1(j))
        vyt = 0.5d0*(vyy(i)+vy1(j))
        vzt = 0.5d0*(vzz(i)+vz1(j))

        stres1(1) = stres1(1) + weight(i)*vxt*vxt
        stres1(2) = stres1(2) + weight(i)*vxt*vyt
        stres1(3) = stres1(3) + weight(i)*vxt*vzt
        stres1(5) = stres1(5) + weight(i)*vyt*vyt
        stres1(6) = stres1(6) + weight(i)*vyt*vzt
        stres1(9) = stres1(9) + weight(i)*vzt*vzt

      enddo
#endif
c     
c     update volume
      
      vold = volm
      volm = volm*exp(3.d0*tstep*chipnew)
c     
c     scale cell vectors - isotropic

      scale= (volm/volm0)**(1.d0/3.d0)
      do i = 1,9
        cell(i) = cell0(i)*scale
      enddo
c     
c     construct scaling tensor (for later!)

      do i = 2,8
        eta(i) = 0.d0
      enddo
      eta(1) = chipnew
      eta(5) = chipnew
      eta(9) = chipnew
c     
c     adjust long range corrections and number density
      
      elrc = elrc0*(volm0/volm)
      virlrc = virlrc0*(volm0/volm)
      
      do kk = 1,ntpatm
        
        dens(kk) = dens0(kk)*(volm0/volm)
        
      enddo
c     
c     update thermostat and barostat variables

      chit = chitnew
      chip = chipnew
c     
c     conserved quantity less kinetic and potential energy terms
      
      conint = conint + 2.d0*sigma*tstep*chit0
      cons0 = conint
      cons1 = sigma*(taut*chit0)**2
      cons2 = press*vold
      cons3 = 1.5d0*pmass*(chip0)**2
      consv = cons0 + cons1 + cons2 + cons3
c     
c     ensure total momentum is zero

      vxm = 0.d0
      vym = 0.d0
      vzm = 0.d0
      do i = iatm0,iatm1
        vxm  = vxm + vxx(i)*weight(i)
        vym  = vym + vyy(i)*weight(i)
        vzm  = vzm + vzz(i)*weight(i)
      enddo

      if(mxnode.gt.1) then

        buffer(4) = vxm
        buffer(5) = vym
        buffer(6) = vzm
        call gdsum(buffer(4),3,buffer(1))
        vxm  =  buffer(4)
        vym  =  buffer(5)
        vzm  =  buffer(6)

      endif
c     
c     correction to velocities

      vxm = vxm/totmas
      vym = vym/totmas
      vzm = vzm/totmas

      do i = iatm0,iatm1
        vxx(i) = vxx(i) - vxm
        vyy(i) = vyy(i) - vym
        vzz(i) = vzz(i) - vzm
      enddo
c     
c     periodic boundary condition
      
      call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
c     
c     global exchange of configuration data
      
      if(mxnode.gt.1)then
        
        nbuff=mxbuff
        call merge(idnode,mxnode,natms,nbuff,xxx,yyy,zzz,buffer)
        call merge(idnode,mxnode,natms,nbuff,vxx,vyy,vzz,buffer)
        
      endif
#ifdef STRESS
c     
c     complete stress tensor

      if(mxnode.gt.1) call gdsum(stres1,9,buffer)

      stress(1) = stress(1) + stres1(1)
      stress(2) = stress(2) + stres1(2)
      stress(3) = stress(3) + stres1(3)
      stress(4) = stress(4) + stres1(2)
      stress(5) = stress(5) + stres1(5)
      stress(6) = stress(6) + stres1(6)
      stress(7) = stress(7) + stres1(3)
      stress(8) = stress(8) + stres1(6)
      stress(9) = stress(9) + stres1(9)
#endif
#ifdef VAMPIR
      call VTEND(36, ierr)
#endif
      return
      end


