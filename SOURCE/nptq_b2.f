      subroutine nptq_b2
     x  (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,ntpatm,elrc,engke,engrot,virlrc,press,
     x  quattol,sigma,taup,taut,tolnce,tstep,vircom,vircon,
     x  virtot,volm,lstme,lishap,lashap,lstfrz,listcon,listme,
     x  lstfre,lstcsit,lstbod,lstgtp,lstrgd,listot,numgsit,cell,
     x  gcmx1,gcmy1,gcmz1,dens,dxx,dyy,dzz,dxt,dyt,dzt,uxx,uyy,
     x  uzz,prmcon,omx,omy,omz,opx,opy,opz,oqx,oqy,oqz,q0,q1,q2,
     x  q3,qn0,qn1,qn2,qn3,gxx,gyy,gzz,xxx,yyy,zzz,vxx,vyy,vzz,
     x  gvxx,gvyy,gvzz,fxx,fyy,fzz,tqx,tqy,tqz,weight,rotinx,
     x  rotiny,rotinz,gcmx,gcmy,gcmz,gmass,buffer,xxt,yyt,zzt,
     x  txx,tyy,tzz,vx1,vy1,vz1,stress,gvx1,gvy1,gvz1,omx1,omy1,
     x  omz1,redmass,eta,dens0,esig1)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using implicit leapfrog quaternion algorithm
c     plus bond constraints - rigid body sites and constraint sites
c     may coincide.
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Berendsen
c     thermostat and barostat.
c     
c     parallel replicated data version : block data
c
c     tqx,tqy,tqz = torque in lab fixed frame (input)
c     omx,omy,omz = angular velocity in body fixed frame (principle axis)
c     rotinx,y,z  = rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory 1995
c     author      t.forester june 1995
c     amended     w.smith sep 1999 : euler equation
c     
c     wl
c     2001/05/30 12:40:14
c     1.7
c     $Sate: Exp $
c     
c**********************************************************************
      
#include "dl_params.inc"
      
      logical safe,safeq,lshmov,newstep,newjob
      
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension rot(9),cell(9),cell0(9),celn(9)
      dimension buffer(mxbuff),weight(mxatms)
      dimension tqx(msgrp),tqy(msgrp),tqz(msgrp)
      dimension omx(mxgrp),omy(mxgrp),omz(mxgrp)
      dimension opx(msgrp),opy(msgrp),opz(msgrp)
      dimension oqx(msgrp),oqy(msgrp),oqz(msgrp)
      dimension gcmx(mxgrp),gcmy(mxgrp),gcmz(mxgrp)
      dimension gcmx1(msgrp),gcmy1(msgrp),gcmz1(msgrp)
      dimension gvxx(mxgrp),gvyy(mxgrp),gvzz(mxgrp),gmass(mxungp)
      dimension q0(mxgrp),q1(mxgrp),q2(mxgrp),q3(mxgrp)
      dimension qn0(msgrp),qn1(msgrp),qn2(msgrp),qn3(msgrp)
      dimension gxx(mxungp,mxngp),gyy(mxungp,mxngp),gzz(mxungp,mxngp)
      dimension rotinx(mxungp,2),rotiny(mxungp,2),rotinz(mxungp,2)
      dimension lstfre(mxatms),lstrgd(mxgatm),numgsit(mxungp)
      dimension lstgtp(mxgrp),lstme(mxatms),lstfrz(mxatms)
      dimension dxx(mxcons),dyy(mxcons),dzz(mxcons)
      dimension uxx(mxatms),uyy(mxatms),uzz(mxatms)
      dimension listcon(mxcons,3)
      dimension listme(mxatms)
      dimension lishap(mxlshp)
      dimension lashap(mxproc),listot(mxatms)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension dxt(mxcons),dyt(mxcons),dzt(mxcons)
      dimension prmcon(mxtcon),esig1(mxcons)
      dimension vx1(msatms),vy1(msatms),vz1(msatms)
      dimension stress(9),stres1(9),stres2(9),eta(9)
      dimension lstcsit(2*mxcons)
      dimension gvx1(msgrp),gvy1(msgrp),gvz1(msgrp)
      dimension lstbod(mxatms)
      dimension omx1(msgrp),omy1(msgrp),omz1(msgrp)
      dimension redmass(mxcons),dens(mxsvdw),dens0(mxsvdw)
      
      save newjob,volm0,elrc0,virlrc0,cell0,chit0,chip0
      
      data newjob/.true./
      data pt5/0.5d0/
      data beta/7.3728d-3/
#ifdef VAMPIR
      call VTBEGIN(39, ierr)
#endif
c     
c     store initial values of volume and long range corrections

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

        chit0 = 1.d0
        chip0 = 1.d0

      endif
      
      nbuff=mxbuff
c     
c     group block indices
      
      igrp1 = (idnode*ngrp)/mxnode + 1
      igrp2 = ((idnode+1)*ngrp)/mxnode
c
c     check work arrays are large enough

      safe =  (igrp2-igrp1+1.le.msgrp) 
      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) then 
        igrp = igrp2-igrp1+1
        call gimax(igrp,1,idum)
        if(idnode.eq.0) write(nrite,*) ' make msgrp >= ',igrp
        call  error(idnode,506)
      endif
      safe=.false.
c     
c     free atom block indices
      
      ifre1 = (idnode*ntfree)/mxnode + 1
      ifre2 = ((idnode+1)*ntfree)/mxnode
c     
c     store initial values of position and velocity
      
      do i = 1,natms
        uxx(i) = xxx(i)
        uyy(i) = yyy(i)
        uzz(i) = zzz(i)
      enddo
      
      j = 0
      do ifre=ifre1,ifre2
        
        i=lstfre(ifre)
        j = j+1
        
        vx1(j) = vxx(i)
        vy1(j) = vyy(i)
        vz1(j) = vzz(i)
        
      enddo
      
      jg = 0
      do ig = igrp1,igrp2
        
        jg = jg+1
        
        gcmx1(jg) = gcmx(ig)
        gcmy1(jg) = gcmy(ig)
        gcmz1(jg) = gcmz(ig)
        
        gvx1(jg) = gvxx(ig)
        gvy1(jg) = gvyy(ig)
        gvz1(jg) = gvzz(ig)
        
        omx1(jg) = omx(ig)
        omy1(jg) = omy(ig)
        omz1(jg) = omz(ig)
        
      enddo
      
c     
c     construct current bond vectors - required by shake
      
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
c     length scaling parameter

      scale = chip0**(1.d0/3.d0)
c     
c     shake iterations and thermostat iterations here!

      do i = 1,9
        stres2(i) = 0.d0
      enddo
      vircon = 0.d0
      mxshak1 = mxshak
      if(ntcons.eq.0) mxshak1 = 3
      do icyc = 1,mxshak1

#ifdef STRESS
c
c     temporary stress tensor accumulators
      do i = 1,9
        stres1(i) = 0.d0
      enddo
#endif
c     
c    integrate 'free' particles
      
      j = 0
      do ifre = ifre1,ifre2
        
        j = j+1
        i = lstfre(ifre)
c     
c     advance velocity by leapfrog

        vxx(i) = (vx1(j) + tstep*(fxx(i)/weight(i)))*chit0
        vyy(i) = (vy1(j) + tstep*(fyy(i)/weight(i)))*chit0
        vzz(i) = (vz1(j) + tstep*(fzz(i)/weight(i)))*chit0
c     
c     advance position by leapfrog

        xxx(i) = uxx(i)*scale + tstep*vxx(i)
        yyy(i) = uyy(i)*scale + tstep*vyy(i)
        zzz(i) = uzz(i)*scale + tstep*vzz(i)
        
      enddo

c     
c     calculate new kinetic energy at current timestep
        
        engke=0.d0
        j = 0
        do ifre = ifre1,ifre2
          i = lstfre(ifre)
          j=j+1
          vxt = 0.5d0*(vxx(i) + vx1(j))
          vyt = 0.5d0*(vyy(i) + vy1(j))
          vzt = 0.5d0*(vzz(i) + vz1(j))
c     
c     kinetic energy at current timestep

          engke = engke + weight(i)*(vxt**2+vyt**2+vzt**2)*pt5
#ifdef STRESS
c
c     kinetic contribution to stress tensor

          stres1(1) = stres1(1) + weight(i)*vxt*vxt
          stres1(2) = stres1(2) + weight(i)*vxt*vyt
          stres1(3) = stres1(3) + weight(i)*vxt*vzt
          stres1(5) = stres1(5) + weight(i)*vyt*vyt
          stres1(6) = stres1(6) + weight(i)*vyt*vzt
          stres1(9) = stres1(9) + weight(i)*vzt*vzt
#endif

        enddo
c
c     *************  Rigid body motion ****************************

c     
c     translational kinetic energy
      
      engtrn = 0.d0
      jg = 0
      jr=0
      do ig = igrp1,igrp2
        
        jg = jg+1
c     
c     at this point qn0-2 are work arrays holding force on c.o.m.
        
        qn0(jg) = 0.d0
        qn1(jg) = 0.d0
        qn2(jg) = 0.d0
        
        id =lstgtp(ig)
        do j = 1,numgsit(id)
          
          jr = jr+1
          i = lstrgd(jr)
c     
c     forces on com
          
          qn0(jg) = qn0(jg) + fxx(i)
          qn1(jg) = qn1(jg) + fyy(i)
          qn2(jg) = qn2(jg) + fzz(i)
          
        enddo
c
c     centre of mass velocities at half-step
        
        tmp = pt5*tstep/gmass(id)
        vxt = gvx1(jg) + tmp*qn0(jg)
        vyt = gvy1(jg) + tmp*qn1(jg)
        vzt = gvz1(jg) + tmp*qn2(jg)
c     
c     translational kinetic energy 
        
        engtrn = engtrn+gmass(id)*(vxt*vxt + vyt*vyt + vzt*vzt)*pt5
#ifdef STRESS
c
c     kinetic contribution to stress tensor

        stres1(1) = stres1(1) + gmass(id)*vxt*vxt
        stres1(2) = stres1(2) + gmass(id)*vxt*vyt
        stres1(3) = stres1(3) + gmass(id)*vxt*vzt
        stres1(5) = stres1(5) + gmass(id)*vyt*vyt
        stres1(6) = stres1(6) + gmass(id)*vyt*vzt
        stres1(9) = stres1(9) + gmass(id)*vzt*vzt
#endif
c     
c     advance velocity by leapfrog
        
        tmp = tstep/gmass(id)
        gvxx(ig) = (gvx1(jg) + tmp*qn0(jg))*chit0
        gvyy(ig) = (gvy1(jg) + tmp*qn1(jg))*chit0
        gvzz(ig) = (gvz1(jg) + tmp*qn2(jg))*chit0
c     
c     advance position by leapfrog

        gcmx(ig) = gcmx1(jg)*scale + tstep*gvxx(ig)
        gcmy(ig) = gcmy1(jg)*scale + tstep*gvyy(ig)
        gcmz(ig) = gcmz1(jg)*scale + tstep*gvzz(ig)
        
      enddo
c     
c     rotational motion .....
      
      engrot = 0.d0
      jg = 0
      jr = 0
      do ig = igrp1,igrp2
        
        id = lstgtp(ig)
        jg =jg+1

        do j = 1,numgsit(id)
          
          jr = jr +1
          i = lstrgd(jr)
          
          xxt(jr) = uxx(i) - gcmx1(jg)
          yyt(jr) = uyy(i) - gcmy1(jg)
          zzt(jr) = uzz(i) - gcmz1(jg)
          
        enddo
        
      enddo
c     
c     minimum images
      
      call images(imcon,0,1,jr,cell,xxt,yyt,zzt)
c     
c     convert atomic virial to molecular one
c     note convention: virial(atom-atom) = -sum(Ri.Fi)
c     : virial(com-com)   = -sum(Rcom.Fcom)
c     so  virial(com-com) = virial(atom-atom) + sum((Ri-Rcom).Fi)
      
      
      vircom = 0.d0
      jr=0
      do ig=igrp1,igrp2
        
        id = lstgtp(ig)
        
        do j = 1,numgsit(id)
          
          jr = jr+1
          i = lstrgd(jr)
          
          vircom=vircom+(xxt(jr)*fxx(i)+yyt(jr)*fyy(i)+zzt(jr)*fzz(i))
#ifdef STRESS
c
c     stress tensor : rigid body contributions

          stres1(1) = stres1(1) - xxt(jr)*fxx(i)
          stres1(2) = stres1(2) - xxt(jr)*fyy(i)
          stres1(3) = stres1(3) - xxt(jr)*fzz(i)
          stres1(4) = stres1(4) - yyt(jr)*fxx(i)
          stres1(5) = stres1(5) - yyt(jr)*fyy(i)
          stres1(6) = stres1(6) - yyt(jr)*fzz(i)
          stres1(7) = stres1(7) - zzt(jr)*fxx(i)
          stres1(8) = stres1(8) - zzt(jr)*fyy(i)
          stres1(9) = stres1(9) - zzt(jr)*fzz(i)
#endif
        enddo

      enddo
c     
c     torques in lab frame
      
      jr=0
      jg=0
      do ig = igrp1,igrp2
        
        jg = jg+1
        id = lstgtp(ig)
        
        tqx(jg) = 0.d0
        tqy(jg) = 0.d0
        tqz(jg) = 0.d0
        
        do j = 1,numgsit(id)
          
          jr = jr +1
          i = lstrgd(jr)
          
          tqx(jg) = tqx(jg) + yyt(jr)*fzz(i) - zzt(jr)*fyy(i)
          tqy(jg) = tqy(jg) + zzt(jr)*fxx(i) - xxt(jr)*fzz(i)
          tqz(jg) = tqz(jg) + xxt(jr)*fyy(i) - yyt(jr)*fxx(i)
          
        enddo
        
c     
c     current rotational matrix 
        
        rot(1) = q0(ig)**2+q1(ig)**2-q2(ig)**2-q3(ig)**2
        rot(2) = 2.d0*(q1(ig)*q2(ig) - q0(ig)*q3(ig))
        rot(3) = 2.d0*(q1(ig)*q3(ig) + q0(ig)*q2(ig))
        rot(4) = 2.d0*(q1(ig)*q2(ig) + q0(ig)*q3(ig))
        rot(5) = q0(ig)**2-q1(ig)**2+q2(ig)**2-q3(ig)**2
        rot(6) = 2.d0*(q2(ig)*q3(ig) - q0(ig)*q1(ig))
        rot(7) = 2.d0*(q1(ig)*q3(ig) - q0(ig)*q2(ig))
        rot(8) = 2.d0*(q2(ig)*q3(ig) + q0(ig)*q1(ig))
        rot(9) = q0(ig)**2-q1(ig)**2-q2(ig)**2+q3(ig)**2
c     
c     angular velocity at time step n (first guess)

        opx(jg) = omx1(jg)
        opy(jg) = omy1(jg)
        opz(jg) = omz1(jg)
c     
c     iterate angular velocity for time step n (e. yezdimer)
        
        do i=1,5
          
          trx=(tqx(jg)*rot(1)+tqy(jg)*rot(4)+tqz(jg)*rot(7))*
     x      rotinx(id,2)
     x      +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*
     x      rotinx(id,2)
          try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x      rotiny(id,2)
     x      +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*
     x      rotiny(id,2)
          trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x      rotinz(id,2)
     x      +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*
     x      rotinz(id,2)

          delx = tstep*trx
          dely = tstep*try
          delz = tstep*trz
c     
c     improved angular velocity at time step n
          
          opx(jg)=(omx1(jg)+pt5*delx)
          opy(jg)=(omy1(jg)+pt5*dely)
          opz(jg)=(omz1(jg)+pt5*delz)
          
        enddo
c
c     scaled angular velocity at time step n

        opx(jg)=opx(jg)*chit0
        opy(jg)=opy(jg)*chit0
        opz(jg)=opz(jg)*chit0
c     
c     angular velocity at time step n+1  (needed for quat algorithm)
          
        oqx(jg) = (omx1(jg) + delx*1.5d0)*chit0
        oqy(jg) = (omy1(jg) + dely*1.5d0)*chit0
        oqz(jg) = (omz1(jg) + delz*1.5d0)*chit0
c     
c     angular velocity at time step n+1/2

        omx(ig) = (omx1(jg) + delx)*chit0
        omy(ig) = (omy1(jg) + dely)*chit0
        omz(ig) = (omz1(jg) + delz)*chit0
c     
c     rotational kinetic energy
        
        engrot = engrot+ pt5*(rotinx(id,1)*opx(jg)**2
     x    +rotiny(id,1)*opy(jg)**2
     x    +rotinz(id,1)*opz(jg)**2)
        
      enddo
c     
c     assign new quaternions

      safeq = .true.
      jg = 0
      jr=0
      do ig = igrp1,igrp2

        id = lstgtp(ig)
        jg = jg+1
c     
c     first iteration of new quaternions (lab fixed)
        
        qn0(jg)=q0(ig)+(-q1(ig)*opx(jg)-q2(ig)*opy(jg)-q3(ig)*opz(jg))
     x    *tstep*pt5
        qn1(jg)=q1(ig)+( q0(ig)*opx(jg)-q3(ig)*opy(jg)+q2(ig)*opz(jg))
     x    *tstep*pt5
        qn2(jg)=q2(ig)+( q3(ig)*opx(jg)+q0(ig)*opy(jg)-q1(ig)*opz(jg))
     x    *tstep*pt5
        qn3(jg)=q3(ig)+(-q2(ig)*opx(jg)+q1(ig)*opy(jg)+q0(ig)*opz(jg))
     x    *tstep*pt5
        
        qn0b =0.d0
        qn1b =0.d0
        qn2b =0.d0
        qn3b =0.d0
        
        itq= 0
  100   itq=itq+1
        
        qn0a = pt5*( -q1(ig)*opx(jg)- q2(ig)*opy(jg)- q3(ig)*opz(jg))
     x    + pt5*(-qn1(jg)*oqx(jg)-qn2(jg)*oqy(jg)-qn3(jg)*oqz(jg))
        qn1a = pt5*(  q0(ig)*opx(jg)- q3(ig)*opy(jg)+ q2(ig)*opz(jg))
     x    +    pt5*( qn0(jg)*oqx(jg)-qn3(jg)*oqy(jg)+qn2(jg)*oqz(jg))
        qn2a = pt5*(  q3(ig)*opx(jg)+ q0(ig)*opy(jg)- q1(ig)*opz(jg))
     x    +    pt5*( qn3(jg)*oqx(jg)+qn0(jg)*oqy(jg)-qn1(jg)*oqz(jg))
        qn3a = pt5*( -q2(ig)*opx(jg)+ q1(ig)*opy(jg)+ q0(ig)*opz(jg))
     x    +    pt5*(-qn2(jg)*oqx(jg)+qn1(jg)*oqy(jg)+qn0(jg)*oqz(jg))
        
        qn0(jg) = q0(ig) + pt5*qn0a*tstep
        qn1(jg) = q1(ig) + pt5*qn1a*tstep
        qn2(jg) = q2(ig) + pt5*qn2a*tstep
        qn3(jg) = q3(ig) + pt5*qn3a*tstep
        
        rnorm = 1.d0/sqrt(qn0(jg)**2+qn1(jg)**2+qn2(jg)**2+qn3(jg)**2)
        qn0(jg) = qn0(jg)*rnorm
        qn1(jg) = qn1(jg)*rnorm
        qn2(jg) = qn2(jg)*rnorm
        qn3(jg) = qn3(jg)*rnorm
c     
c     convergence test 
        
        eps = sqrt(((qn0a-qn0b)**2+(qn1a-qn1b)**2+(qn2a-qn2b)**2
     x    +(qn3a-qn3b)**2)*tstep**2)
        
        qn0b = qn0a
        qn1b = qn1a
        qn2b = qn2a
        qn3b = qn3a
        
        if((itq.lt.mxquat).and.(eps.gt.quattol)) goto 100
        if(itq.ge.mxquat) safeq= .false.
        
      enddo
c     
c     minimum images of group positions and particle positions
      
      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
      call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)
c     
c     new atomic positions for atoms in rigid bodies
      
      jg =0
      jr=0
      do ig = igrp1,igrp2

        jg =jg+1
c     
c     group type
        
        id = lstgtp(ig)
c     
c     new rotational matrix
        
        rot(1) = qn0(jg)**2+qn1(jg)**2-qn2(jg)**2-qn3(jg)**2
        rot(2) = 2.d0*(qn1(jg)*qn2(jg) - qn0(jg)*qn3(jg))
        rot(3) = 2.d0*(qn1(jg)*qn3(jg) + qn0(jg)*qn2(jg))
        rot(4) = 2.d0*(qn1(jg)*qn2(jg) + qn0(jg)*qn3(jg))
        rot(5) = qn0(jg)**2-qn1(jg)**2+qn2(jg)**2-qn3(jg)**2
        rot(6) = 2.d0*(qn2(jg)*qn3(jg) - qn0(jg)*qn1(jg))
        rot(7) = 2.d0*(qn1(jg)*qn3(jg) - qn0(jg)*qn2(jg))
        rot(8) = 2.d0*(qn2(jg)*qn3(jg) + qn0(jg)*qn1(jg))
        rot(9) = qn0(jg)**2-qn1(jg)**2-qn2(jg)**2+qn3(jg)**2
        
        do j = 1,numgsit(id)
          
          jr = jr +1
          i = lstrgd(jr)
          
          xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x      +gcmx(ig)
          yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x      +gcmy(ig)
          zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x      +gcmz(ig)
c
c     new atomic velocites in body frame
          
          vaa = omy(ig)*gzz(id,j) - omz(ig)*gyy(id,j)
          vbb = omz(ig)*gxx(id,j) - omx(ig)*gzz(id,j)
          vcc = omx(ig)*gyy(id,j) - omy(ig)*gxx(id,j)
c     
c     new atomic velocites in lab frame
          
          vxx(i) = rot(1)*vaa+rot(2)*vbb+rot(3)*vcc + gvxx(ig)
          vyy(i) = rot(4)*vaa+rot(5)*vbb+rot(6)*vcc + gvyy(ig)
          vzz(i) = rot(7)*vaa+rot(8)*vbb+rot(9)*vcc + gvzz(ig)
          
        enddo
        
      enddo

      if(mxnode.gt.1) then
c     
c     merge new atomic coordinates 
      
        call merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
      
      endif
c
c     globally sum kinetic energies

      if(mxnode.gt.1) then

        buffer(5) = engke
        buffer(6) = engtrn
        buffer(7) = engrot
        buffer(8) = vircom
        
        call  gdsum(buffer(5),4,buffer(1))
        
        engke =  buffer(5)
        engtrn = buffer(6)
        engrot = buffer(7)
        vircom = buffer(8)
        
      endif
c     
c     pressure control variable

      engke = engke + engtrn
      pr=(2.d0*engke-virtot-vircon-vircom)/(3.d0*volm)
      chip0 = 1.d0 + beta*tstep*(pr - press)/taup
      scale = chip0**(1.d0/3.d0)
c
c     new estimate of chit
      
      engtot = engke + engrot 
      chit0 =  sqrt(1.d0 + tstep/taut*(sigma/engtot -1.d0))

      if(ntcons.gt.0) then
c     
c     apply constraint correction
          
        newstep = .false.
        if(icyc.eq.1) newstep = .true.
c
c     guess new cell parameters 

        do i = 1,9
          celn(i) = scale*cell(i)
        enddo

        call qshake
     x    (newstep,safe,lshmov,idnode,imcon,mxnode,natms,
     x    nscons,tolnce,tstep,viracc,lashap,lishap,listcon,
     x    listme,listot,lstfrz,lstbod,lstgtp,lstcsit,
     x    buffer,celn,dxt,dxx,dyt,dyy,dzt,dzz,prmcon,txx,
     x    tyy,tzz,weight,xxt,xxx,yyt,yyy,zzt,zzz,fxx,fyy,
     x    fzz,stres2,q0,q1,q2,q3,gxx,gyy,gzz,gmass,rotinx,
     x    rotiny,rotinz,redmass,esig1)

        if(viracc.eq.0.d0.and.icyc.gt.3) goto 110
        vircon = vircon + viracc
c     
c     end of shake corrections
          
      endif

      enddo
      
  110  continue
c
c     update quaternions

      jg = 0
      do ig = igrp1,igrp2
        
        jg = jg+1
        q0(ig) = qn0(jg)
        q1(ig) = qn1(jg)
        q2(ig) = qn2(jg)
        q3(ig) = qn3(jg)

      enddo
c     
c     scale cell vectors

      scale = ((chip0*volm)/volm0)**(1.d0/3.d0)

      do i = 1,9
        cell(i) = scale*cell0(i)
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

      if(mxnode.gt.1) then
c
c     merge new group coordinates and velocities

        nbuff = mxbuff
        call merge(idnode,mxnode,ngrp,nbuff,gcmx,gcmy,gcmz,buffer)
        call merge(idnode,mxnode,ngrp,nbuff,gvxx,gvyy,gvzz,buffer)
c     
c     merge new atomic velocities
      
        call merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
c
c     merge new quaternions

        call merge4(idnode,mxnode,ngrp,nbuff,q0,q1,q2,q3,buffer)

      endif
c
c     ensure all atoms are within cell boundaries

      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)

#ifdef STRESS
c     
c     complete stress tensor

      do i = 1,9
        stres1(i) = stres1(i) + stres2(i)
      enddo

      if(mxnode.gt.1) call gdsum(stres1,9,buffer)

      do i = 1,9
        stress(i) = stress(i) + stres1(i)
      enddo
c
c     symmetrise stress tensor

      stress(2) = 0.5d0*(stress(2)+stress(4))
      stress(4) = stress(2)
      stress(3) = 0.5d0*(stress(3)+stress(7))
      stress(7) = stress(3)
      stress(6) = 0.5d0*(stress(6)+stress(8))
      stress(8) = stress(6)
#endif
#ifdef VAMPIR
      call VTEND(39, ierr)
#endif
      return
      end
