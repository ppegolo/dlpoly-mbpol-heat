      subroutine npt_b0
     x  (idnode,imcon,mxnode,natms,ntpatm,elrc,engke,virlrc,press,
     x  taup,taut,sigma,tstep,virtot,volm,lstfrz,buffer,cell,dens,
     x  fxx,fyy,fzz,vxx,vyy,vzz,weight,xdf,xxx,ydf,yyy,zdf,zzz,
     x  vx1,vy1,vz1,stress,eta,dens0)
      
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Berendsen
c     thermostat amd isotropic pressure control
c     isothermal compressibility (beta) set to that of liquid water
c     = 0.007372 dlpoly units
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1993/1999
c     author    -    t. forester december 1993
c     adapted - atomic systems w.smith oct 1999
c     
c     wl
c     2001/05/30 12:40:14
c     1.3
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical newjob
      
      dimension lstfrz(mxatms),buffer(mxbuff)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf)
      dimension weight(mxatms),cell(9)
      dimension dens(mxsvdw),dens0(mxsvdw),celn(9)
      dimension vx1(msatms),vy1(msatms),vz1(msatms)
      dimension stress(9),stres1(9),eta(9)

      save newjob,volm0,elrc0,virlrc0

      data newjob/.true./
      data beta/7.3728d-3/
      
#ifdef VAMPIR
      call VTBEGIN(165, ierr)
#endif
      nbuff=mxbuff
c
c     store initial values of volume and long range corrections

      if(newjob) then

        volm0 = volm
        elrc0 = elrc
        virlrc0 = virlrc
        do i = 1,ntpatm
          dens0(i) = dens(i)
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
c     estimate kinetic energy and momentum at current timestep

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
c     pressure control variable

      pr=(2.d0*engke-virtot)/(3.d0*volm)
      
      chip0 = 1.d0 + beta*tstep*(pr - press)/taup
      scale = chip0**(1.d0/3.d0)
c     
c     temperature scaling  coefficient - taut is the decay constant
      
      chit0 =  sqrt(1.d0 + tstep/taut*(sigma/engke -1.d0))
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

          vxx(i)=(vx1(j) + tstep*(fxx(i)/weight(i)))*chit0
          vyy(i)=(vy1(j) + tstep*(fyy(i)/weight(i)))*chit0
          vzz(i)=(vz1(j) + tstep*(fzz(i)/weight(i)))*chit0
c     
c     update positions : 
          
          xxx(i)= tstep*vxx(i) + scale*xdf(j)
          yyy(i)= tstep*vyy(i) + scale*ydf(j)
          zzz(i)= tstep*vzz(i) + scale*zdf(j)

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

        pr=(2.d0*engke-virtot)/(3.d0*volm)
        chip0 = 1.d0 + beta*tstep*(pr - press)/taup
        scale = chip0**(1.d0/3.d0)

        chit0 =  sqrt(1.d0 + tstep/taut*(sigma/engke -1.d0))
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

      volm = volm*chip0
c     
c     scale cell vectors - isotropic

      do i = 1,9
        cell(i) = cell(i)*scale
      enddo
c
c     construct scaling tensor (for later!)

      do i = 2,8
        eta(i) = 0.d0
      enddo
      eta(1) = scale
      eta(5) = scale
      eta(9) = scale
c     
c     adjust long range corrections and number density
      
      elrc = elrc0*(volm0/volm)
      virlrc = virlrc0*(volm0/volm)
      
      do kk = 1,ntpatm
        
        dens(kk) = dens0(kk)*(volm0/volm)
        
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
      call VTEND(165, ierr)
#endif
      return
      end


