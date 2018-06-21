      subroutine nst_h0
     x  (idnode,imcon,mxnode,natms,ntpatm,chit,conint,consv,
     x  elrc,engke,virlrc,press,taup,taut,sigma,temp,tstep,volm,
     x  lstfrz,buffer,cell,dens,fxx,fyy,fzz,vxx,vyy,vzz,weight,
     x  xdf,xxx,ydf,yyy,zdf,zzz,vx1,vy1,vz1,stress,eta,dens0)
      
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Hoover 
c     thermostat + piston.
c     
c     Parrinello- Rahman type : changing cell shape.
c     
c     reference: Melchionna, Ciccotti and Holian,
c     Mol Phys 1993, 78, p533
c     
c     parallel replicated data version
c     
c     copyright daresbury laboratory 1995/1999
c     author    - t. forester     june  1995
c     adapted   - w. smith atomic systems nov 1999
c     
c     wl
c     2001/05/30 12:40:17
c     1.3
c     $Sate: Exp $
c     
c***********************************************************************
      
#include "dl_params.inc"
      
      logical newjob
      
      dimension lstfrz(mxatms),buffer(mxbuff)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf)
      dimension dens(mxsvdw),dens0(mxsvdw),celn(9)
      dimension vx1(msatms),vy1(msatms),vz1(msatms)
      dimension stress(9),stres1(9),stres2(9)
      dimension eta(9),etanew(9),eta0(9)
      dimension weight(mxatms),cell(9)

      save newjob,volm0,elrc0,virlrc0

      data newjob/.true./
#ifndef STRESS
c     
c     this must be compiled with the stress flag on

      call error(idnode,434)
#endif

#ifdef VAMPIR
      call VTBEGIN(167, ierr)
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
        newjob = .false.

      endif
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
c      
c     reciprocal of timestep
      rstep = 1.d0/tstep
      rtsq = rstep**2 
c     
c     inertia parameter for Nose-Hoover thermostat
      
c     qmass = 2.0d0*sigma*taut**2
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
c     propagate eta

      do i = 1,9

        etap = 0.d0
        if(i.eq.1) then
          etap = press*volm
        elseif(i.eq.5) then
          etap = press*volm
        elseif(i.eq.9) then
          etap = press*volm
        endif

        etanew(i) = eta(i) + tstep*(stres1(i) - etap)/pmass
        eta0(i) = 0.5d0*(etanew(i)+eta(i)) 

      enddo
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

          vxt = 0.5d0*(vxx(i)+vx1(j))
          vyt = 0.5d0*(vyy(i)+vy1(j))
          vzt = 0.5d0*(vzz(i)+vz1(j))

          vxx(i) = vx1(j) + tstep*(fxx(i)/weight(i) -
     x      (chit0+eta0(1))*vxt - eta0(2)*vyt - eta0(3)*vzt)
          vyy(i) = vy1(j) + tstep*(fyy(i)/weight(i) -
     x      eta0(2)*vxt - (eta0(5)+chit0)*vyt - eta0(6)*vzt)
          vzz(i) = vz1(j) + tstep*(fzz(i)/weight(i) -
     x      eta0(3)*vxt - eta0(6)*vyt - (eta0(9)+chit0)*vzt)
c     
c     advance positions using leapfrog

          xxa = (xxx(i)+xdf(j))*0.5d0 - xcmo
          yya = (yyy(i)+ydf(j))*0.5d0 - ycmo
          zza = (zzz(i)+zdf(j))*0.5d0 - zcmo

          xxx(i) = xdf(j) + tstep*(vxx(i) +
     x      etanew(1)*xxa + etanew(2)*yya + etanew(3)*zza)
          yyy(i) = ydf(j) + tstep*(vyy(i) +
     x      etanew(2)*xxa + etanew(5)*yya + etanew(6)*zza)
          zzz(i) = zdf(j) + tstep*(vzz(i) +
     x      etanew(3)*xxa + etanew(6)*yya + etanew(9)*zza)

        enddo
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

        do i = 1,9

          etap = 0.d0
          if(i.eq.1) then
            etap = press*volm
          elseif(i.eq.5) then
            etap = press*volm
          elseif(i.eq.9) then
            etap = press*volm
          endif

          etanew(i) = eta(i) + tstep*(stres1(i) - etap)/pmass
          eta0(i) = 0.5d0*(etanew(i)+eta(i)) 

        enddo

        chitp = (engke/sigma - 1.d0)/taut**2
        chitnew = chit + tstep*chitp
        chit0 = 0.5d0*(chit+chitnew)
c     
c     end of thermostat iterations

      enddo
c     
c     update stress tensor

      do i = 1,9
        stress(i) = stres1(i)
      enddo
c     
c     update thermostat and barostat variables

      chit = chitnew
      do i = 1,9
        eta(i) = etanew(i)
      enddo
c     
c     update volume
      
      chip = eta(1)+eta(5)+eta(9)
      vold = volm
      volm = volm*exp(tstep*chip)
c     
c     adjust cell vectors - anisotropic

      a1 = (tstep*eta(1))
      a2 = (tstep*eta(2))
      a3 = (tstep*eta(3))
      a5 = (tstep*eta(5))
      a6 = (tstep*eta(6))
      a9 = (tstep*eta(9))

      b1 = (a1*a1 + a2*a2 + a3*a3)*0.5d0 + a1
      b2 = (a1*a2 + a2*a5 + a3*a6)*0.5d0 + a2
      b3 = (a1*a3 + a2*a6 + a3*a9)*0.5d0 + a3
      b5 = (a2*a2 + a5*a5 + a6*a6)*0.5d0 + a5
      b6 = (a2*a3 + a5*a6 + a6*a9)*0.5d0 + a6
      b9 = (a3*a3 + a6*a6 + a9*a9)*0.5d0 + a9

      c1 = b1*cell(1) + b2*cell(4) + b3*cell(7)
      c2 = b1*cell(2) + b2*cell(5) + b3*cell(8)
      c3 = b1*cell(3) + b2*cell(6) + b3*cell(9)
      c4 = b2*cell(1) + b5*cell(4) + b6*cell(7)
      c5 = b2*cell(2) + b5*cell(5) + b6*cell(8)
      c6 = b2*cell(3) + b5*cell(6) + b6*cell(9)
      c7 = b3*cell(1) + b6*cell(4) + b9*cell(7)
      c8 = b3*cell(2) + b6*cell(5) + b9*cell(8)
      c9 = b3*cell(3) + b6*cell(6) + b9*cell(9)

      cell(1) = c1 + cell(1)
      cell(2) = c2 + cell(2)
      cell(3) = c3 + cell(3)
      cell(4) = c4 + cell(4)
      cell(5) = c5 + cell(5)
      cell(6) = c6 + cell(6)
      cell(7) = c7 + cell(7)
      cell(8) = c8 + cell(8)
      cell(9) = c9 + cell(9)
c     
c     adjust long range corrections and number density
      
      elrc = elrc0*(volm0/volm)
      virlrc = virlrc0*(volm0/volm)
      
      do kk = 1,ntpatm
        
        dens(kk) = dens0(kk)*(volm0/volm)
        
      enddo
c     
c     conserved quantity less kinetic and potential energy
      
      conint = conint + 2.d0*sigma*tstep*chit0
      cons0 = conint
      cons1 = sigma*(taut*chit0)**2
      cons2 = press*vold
      eta0(4) = eta0(2)
      eta0(7) = eta0(3)
      eta0(8) = eta0(6)
      tr = 0.d0
      do i = 1,9
        tr = tr + eta0(i)**2
      enddo
      cons3 = 0.5d0*pmass*tr
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
#ifdef VAMPIR
      call VTEND(167, ierr)
#endif
      return
      end


