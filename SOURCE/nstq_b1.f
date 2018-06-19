      subroutine nstq_b1
     x  (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,ntpatm,elrc,engke,engrot,virlrc,press,
     x  quattol,sigma,taup,taut,tolnce,tstep,vircom,vircon,volm,
     x  lstme,lishap,lashap,lstfrz,listcon,listme,lstfre,lstgtp,
     x  lstrgd,listot,numgsit,cell,dens,gcmx1,gcmy1,gcmz1,dxx,
     x  dyy,dzz,dxt,dyt,dzt,gvx1,gvy1,gvz1,uxx,uyy,uzz,prmcon,
     x  xdf,ydf,zdf,omx,omy,omz,omx1,omy1,omz1,opx,opy,opz,oqx,
     x  oqy,oqz,q0,q1,q2,q3,qn0,qn1,qn2,qn3,gxx,gyy,gzz,xxx,yyy,
     x  zzz,vxx,vyy,vzz,gvxx,gvyy,gvzz,fxx,fyy,fzz,tqx,tqy,
     x  tqz,weight,rotinx,rotiny,rotinz,gcmx,gcmy,gcmz,gmass,
     x  buffer,xxt,yyt,zzt,txx,tyy,tzz,vx1,vy1,vz1,stress,eta,
     x  dens0)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using implicit leapfrog quaternion algorithm
c     plus bond constraints - PROVIDED rigid body sites
c     and constraint sites do not coincide.
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Berendsen
c     thermostat and barostat. (cell may change shape).
c     isothermal compressibility (beta) set to that of liquid water
c     = 0.007372 dlpoly units
c     
c     parallel replicated data version : block data
c     
c     for systems using bond constraints
c     
c     parallel replicated data version : block data
c     
c     tqx,tqy,tqz = torque in lab fixed frame (input)
c     omx,omy,omz = angular velocity in body fixed frame (principl axis)
c     rotinx,y,z  = rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory 1993
c     author      t.forester october 1993
c     amended     t.forester dec 1994 : block data
c     amended     w.smith sep 1999 : euler equation
c     
c     wl
c     2003/03/28 16:21:14
c     1.1
c     $Sate: Exp $
c     
c**********************************************************************
      
#include "dl_params.inc"
      
      logical safe,lshmov,newjob,safeq
      
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension rot(9),cell(9),celn(9)
      dimension buffer(mxbuff),weight(mxatms)
      dimension tqx(msgrp),tqy(msgrp),tqz(msgrp)
      dimension omx(mxgrp),omy(mxgrp),omz(mxgrp)
      dimension omx1(msgrp),omy1(msgrp),omz1(msgrp)
      dimension opx(msgrp),opy(msgrp),opz(msgrp)
      dimension oqx(msgrp),oqy(msgrp),oqz(msgrp)
      dimension gcmx(mxgrp),gcmy(mxgrp),gcmz(mxgrp)
      dimension gcmx1(msgrp),gcmy1(msgrp),gcmz1(msgrp)
      dimension gvxx(mxgrp),gvyy(mxgrp),gvzz(mxgrp),gmass(mxungp)
      dimension gvx1(msgrp),gvy1(msgrp),gvz1(msgrp)
      dimension q0(mxgrp),q1(mxgrp),q2(mxgrp),q3(mxgrp)
      dimension qn0(msgrp),qn1(msgrp),qn2(msgrp),qn3(msgrp)
      dimension gxx(mxungp,mxngp),gyy(mxungp,mxngp),gzz(mxungp,mxngp)
      dimension rotinx(mxungp,2),rotiny(mxungp,2),rotinz(mxungp,2)
      dimension lstfre(mxatms),lstrgd(mxgatm),numgsit(mxungp)
      dimension lstgtp(mxgrp),lstme(mxatms),lstfrz(mxatms)
      dimension dxx(mxcons),dyy(mxcons),dzz(mxcons)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf)
      dimension uxx(mxatms),uyy(mxatms),uzz(mxatms)
      dimension listcon(mxcons,3)
      dimension listme(mxatms)
      dimension lishap(mxlshp)
      dimension lashap(mxproc),listot(mxatms)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension dxt(mxcons),dyt(mxcons),dzt(mxcons)
      dimension prmcon(mxtcon),dens(mxsvdw),dens0(mxsvdw)
      dimension vx1(msatms),vy1(msatms),vz1(msatms)
      dimension stress(9),stres1(9),stres2(9),eta(9)

      save newjob,volm0,elrc0,virlrc0

      data newjob/.true./
      data pt5/0.5d0/
      data beta/7.3728d-3/
#ifndef STRESS
c     
c     this must be compiled with the stress flag on

      call error(idnode,434)
#endif

#ifdef VAMPIR
      call VTBEGIN(40, ierr)
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
c     constraint virial

      vircon=0.d0
c     
c     temporary stress tensor accumulators and new cell
      do i = 1,9
        stres1(i) = 0.d0
        stres2(i) = 0.d0
        celn(i) = cell(i)
      enddo
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
c     reciprocal of timestep
      rstep = 1.d0/tstep
      rtsq = rstep**2
c     
c     store initial values of position and velocity

      j = 0
      do ifre=ifre1,ifre2
        
        i=lstfre(ifre)
        j = j+1
        
        xdf(j)=xxx(i)
        ydf(j)=yyy(i)
        zdf(j)=zzz(i)

        vx1(j) = vxx(i)
        vy1(j) = vyy(i)
        vz1(j) = vzz(i)
        
      enddo
      
      jg = 0
      do ig = igrp1,igrp2
        
        jg = jg+1
        omx1(jg) = omx(ig)
        omy1(jg) = omy(ig)
        omz1(jg) = omz(ig)
        
        gcmx1(jg) = gcmx(ig)
        gcmy1(jg) = gcmy(ig)
        gcmz1(jg) = gcmz(ig)
        
        gvx1(jg) = gvxx(ig)
        gvy1(jg) = gvyy(ig)
        gvz1(jg) = gvzz(ig)

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
c     estimate velocity and temperature at half-time step
      
      engke = 0.d0
      j =0
      do ifre = ifre1,ifre2
        
        i = lstfre(ifre)
        j = j+1
c     
c     estimate new velocity
        
        tmp = tstep/weight(i)*pt5
        vxt = vxx(i) + tmp*fxx(i)
        vyt = vyy(i) + tmp*fyy(i)
        vzt = vzz(i) + tmp*fzz(i)
c     
c     kinetic energy * 2
        
        engke=engke+weight(i)*(vxt*vxt + vyt*vyt + vzt*vzt)
c
c     stress tensor

        stres1(1) = stres1(1) + weight(i)*vxt*vxt
        stres1(2) = stres1(2) + weight(i)*vxt*vyt
        stres1(3) = stres1(3) + weight(i)*vxt*vzt
        stres1(5) = stres1(5) + weight(i)*vyt*vyt
        stres1(6) = stres1(6) + weight(i)*vyt*vzt
        stres1(9) = stres1(9) + weight(i)*vzt*vzt
        
      enddo
c
c     kinetic energy

      engke = engke*pt5
c     
c     estimate translational kinetic energy
      
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
        vxt = gvxx(ig) + tmp*qn0(jg)
        vyt = gvyy(ig) + tmp*qn1(jg)
        vzt = gvzz(ig) + tmp*qn2(jg)
c     
c     translational kinetic energy * 2
        
        engtrn = engtrn+gmass(id)*(vxt*vxt + vyt*vyt + vzt*vzt)
c
c     stress tensor

        stres1(1) = stres1(1) + gmass(id)*vxt*vxt
        stres1(2) = stres1(2) + gmass(id)*vxt*vyt
        stres1(3) = stres1(3) + gmass(id)*vxt*vzt
        stres1(5) = stres1(5) + gmass(id)*vyt*vyt
        stres1(6) = stres1(6) + gmass(id)*vyt*vzt
        stres1(9) = stres1(9) + gmass(id)*vzt*vzt
        
      enddo
c
c     translation kinetic energy of rigid bodies

      engtrn = engtrn*pt5
c     
c     estimate rotational kinetic energy
      
      engrot = 0.d0
      jg = 0
      jr = 0
      do ig = igrp1,igrp2
        
        jg = jg+1 
        id = lstgtp(ig)
        do j = 1,numgsit(id)
          
          jr = jr +1
          i = lstrgd(jr)
          
          xxt(jr) = xxx(i) - gcmx1(jg)
          yyt(jr) = yyy(i) - gcmy1(jg)
          zzt(jr) = zzz(i) - gcmz1(jg)
          
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
c
c     stress tensor : rigid body contributions

          stres2(1) = stres2(1) - xxt(jr)*fxx(i)
          stres2(2) = stres2(2) - xxt(jr)*fyy(i)
          stres2(3) = stres2(3) - xxt(jr)*fzz(i)
          stres2(4) = stres2(4) - yyt(jr)*fxx(i)
          stres2(5) = stres2(5) - yyt(jr)*fyy(i)
          stres2(6) = stres2(6) - yyt(jr)*fzz(i)
          stres2(7) = stres2(7) - zzt(jr)*fxx(i)
          stres2(8) = stres2(8) - zzt(jr)*fyy(i)
          stres2(9) = stres2(9) - zzt(jr)*fzz(i)

        enddo

      enddo
c
c     make stres2 symmetric

      stres2(2) = 0.5d0*(stres2(2)+stres2(4))
      stres2(4) = stres2(2)
      stres2(3) = 0.5d0*(stres2(3)+stres2(7))
      stres2(7) = stres2(3)
      stres2(6) = 0.5d0*(stres2(6)+stres2(8))
      stres2(8) = stres2(6)
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

        opx(jg) = omx(ig)
        opy(jg) = omy(ig)
        opz(jg) = omz(ig)
c
c       iterate angular velocity for time step n (e. yezdimer)

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
c     
c     improved angular velocity at time step n

          opx(jg) = omx(ig) + pt5*tstep*trx
          opy(jg) = omy(ig) + pt5*tstep*try
          opz(jg) = omz(ig) + pt5*tstep*trz

        enddo
c     
c     rotational kinetic energy
        
        engrot = engrot+ pt5*(rotinx(id,1)*opx(jg)**2
     x    +rotiny(id,1)*opy(jg)**2
     x    +rotinz(id,1)*opz(jg)**2)
        
      enddo
c
c     complete stress tensor
      do i = 1,9
        stres1(i) = stres1(i) + stres2(i)
      enddo
      
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

        call gdsum(stres1,9,buffer)
        
      endif

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
      
      engke = engke + engtrn
      engtot = engke + engrot
      chit0 =  sqrt(1.d0 + tstep/taut*(sigma/engtot -1.d0))
c     
c     begin iterations !!-----------------------------------------------

      mxiter = 5
      if(ntcons.eq.0) mxiter = mxiter -1
      
      do iter = 1,mxiter
c     
c     unconstrained new positions
        
        j = 0
        do ifre = ifre1,ifre2

          i = lstfre(ifre)
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
c
c     estimate new cell tensor

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

        if(ntcons.eq.0) safe =.true.
        if(ntcons.gt.0) then
c     
c     store integrated positions

          j= 0
          do ifre=ifre1,ifre2
            i = lstfre(ifre)
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
            
             call merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
            
          endif

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
c     calculate force correction
          
          j=0
          do ifre = ifre1,ifre2
            i = lstfre(ifre)
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
        do ifre = ifre1,ifre2
          i = lstfre(ifre)
          j=j+1
          vxt = 0.5d0*(vxx(i) + vx1(j))
          vyt = 0.5d0*(vyy(i) + vy1(j))
          vzt = 0.5d0*(vzz(i) + vz1(j))
c     
c     kinetic energy at current timestep

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

c     
c     ********: rigid body motion - thermostated  :*****************************
        
c     
c     ***** step 1 : integrate centre of mass motion *********
        
        jg  = 0
        do ig = igrp1,igrp2

          jg = jg+1
c     
c     calculate thermostated velocities
          
          id =lstgtp(ig)
          
          rmass = 1.d0/gmass(id)

          gvxx(ig) = (gvx1(jg) + tstep*(qn0(jg)*rmass))*chit0
          gvyy(ig) = (gvy1(jg) + tstep*(qn1(jg)*rmass))*chit0
          gvzz(ig) = (gvz1(jg) + tstep*(qn2(jg)*rmass))*chit0
c     
c     update positions : 

          gcmx(ig)= tstep*gvxx(ig)+(1.d0+eta1)*gcmx1(jg)+eta2*gcmy1(jg)+
     x      eta3*gcmz1(jg)
          gcmy(ig)= tstep*gvyy(ig)+eta2*gcmx1(jg)+(1.d0+eta5)*gcmy1(jg)+
     x      eta6*gcmz1(jg)
          gcmz(ig)= tstep*gvzz(ig)+eta3*gcmx1(jg)+eta6*gcmy1(jg)+
     x      (1.d0+eta9)*gcmz1(jg)
          
        enddo
c     
c     calculate kinetic energy
        
        engtrn = 0.d0
        jg = 0
        do ig = igrp1,igrp2
          
          jg = jg+1
          id = lstgtp(ig)
c     
c     velocity at half step
          
          vxt = pt5*(gvx1(jg)+gvxx(ig))
          vyt = pt5*(gvy1(jg)+gvyy(ig))
          vzt = pt5*(gvz1(jg)+gvzz(ig))
          
          engtrn = engtrn+pt5*gmass(id)*(vxt*vxt+vyt*vyt+vzt*vzt)
c     
c     kinetic contribution stress tensor

          stres1(1) = stres1(1) + gmass(id)*vxt*vxt
          stres1(2) = stres1(2) + gmass(id)*vxt*vyt
          stres1(3) = stres1(3) + gmass(id)*vxt*vzt
          stres1(5) = stres1(5) + gmass(id)*vyt*vyt
          stres1(6) = stres1(6) + gmass(id)*vyt*vzt
          stres1(9) = stres1(9) + gmass(id)*vzt*vzt

        enddo
c     
c     ****** step 2 : integrate rotational motion **********
        
        safeq = .true.
        engrot = 0.d0
        
        jg=0
        do ig = igrp1,igrp2
          
          jg = jg+1
          id = lstgtp(ig)
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
c     
c     scaled angular velocity at time step n
            
          opx(jg)=(omx1(jg)+pt5*tstep*trx)*chit0
          opy(jg)=(omy1(jg)+pt5*tstep*try)*chit0
          opz(jg)=(omz1(jg)+pt5*tstep*trz)*chit0
c     
c     angular velocity at time step n+1/2
          
          omx(ig)=(omx1(jg)+tstep*trx)*chit0
          omy(ig)=(omy1(jg)+tstep*try)*chit0
          omz(ig)=(omz1(jg)+tstep*trz)*chit0
c     
c     angular velocity at time step n+1 (needed for quat algorithm)
          
          oqx(jg)=(omx1(jg)+1.5d0*tstep*trx)*chit0
          oqy(jg)=(omy1(jg)+1.5d0*tstep*try)*chit0
          oqz(jg)=(omz1(jg)+1.5d0*tstep*trz)*chit0
c     
c     rotational kinetic energy
          
          engrot = engrot+ pt5*(rotinx(id,1)*opx(jg)**2
     x      +rotiny(id,1)*opy(jg)**2
     x      +rotinz(id,1)*opz(jg)**2)
          
        enddo
        
        if(mxnode.gt.1) then
          
          buffer(4) = engke
          buffer(5) = engtrn
          buffer(6) = engrot
          
          call  gdsum(buffer(4),3,buffer(1))
          
          engke = buffer(4)
          engtrn = buffer(5)
          engrot = buffer(6)
          
        endif
c     
c     complete stress tensor - add constraint and c.o.m. contributions

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

        engke = engke+engtrn
        engtot = engke+engrot
        chit0 =  sqrt(1.d0 + tstep/taut*(sigma/engtot -1.d0))
c     
c     end of thermostat/barostat iterations

      enddo
c     
c     update stress tensor and cell

      do i = 1,9
        stress(i) = stres1(i)
        cell(i) = celn(i)
      enddo
c     
c     update volume
      
      volm = volm*(1.d0+eta1)*(1.d0+eta5)*(1.d0+eta9)
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
c     adjust long range corrections and number density
      
      elrc = elrc0*(volm0/volm)
      virlrc = virlrc0*(volm0/volm)
      
      do kk = 1,ntpatm
        
        dens(kk) = dens0(kk)*(volm0/volm)
        
      enddo
c     
c     assign new quaternions

      jg = 0
      jr=0
      do ig = igrp1,igrp2

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
     x       + pt5*(-qn1(jg)*oqx(jg)-qn2(jg)*oqy(jg)-qn3(jg)*oqz(jg))
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
        
        eps = ((qn0a-qn0b)**2+(qn1a-qn1b)**2+(qn2a-qn2b)**2
     x    +(qn3a-qn3b)**2)*tstep**2
        
        qn0b = qn0a
        qn1b = qn1a
        qn2b = qn2a
        qn3b = qn3a
        
        if((itq.lt.mxquat).and.(eps.gt.quattol)) goto 100
        if(itq.ge.mxquat) safeq= .false.
c     
c     store new quaternions
        
        q0(ig) = qn0(jg)
        q1(ig) = qn1(jg)
        q2(ig) = qn2(jg)
        q3(ig) = qn3(jg)
        
      enddo
c     
c     minimum images 
      
      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
      call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)
c     
c     new atomic positions for atoms in rigid bodies
      
      jr=0
      do ig = igrp1,igrp2
c     
c     group type
        
        id = lstgtp(ig)
c     
c     new rotational matrix
        
        rot(1) = q0(ig)**2+q1(ig)**2-q2(ig)**2-q3(ig)**2
        rot(2) = 2.d0*(q1(ig)*q2(ig) - q0(ig)*q3(ig))
        rot(3) = 2.d0*(q1(ig)*q3(ig) + q0(ig)*q2(ig))
        rot(4) = 2.d0*(q1(ig)*q2(ig) + q0(ig)*q3(ig))
        rot(5) = q0(ig)**2-q1(ig)**2+q2(ig)**2-q3(ig)**2
        rot(6) = 2.d0*(q2(ig)*q3(ig) - q0(ig)*q1(ig))
        rot(7) = 2.d0*(q1(ig)*q3(ig) - q0(ig)*q2(ig))
        rot(8) = 2.d0*(q2(ig)*q3(ig) + q0(ig)*q1(ig))
        rot(9) = q0(ig)**2-q1(ig)**2-q2(ig)**2+q3(ig)**2
        
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
c     merge new group coordinates and velocities
      
        nbuff = mxbuff
        call merge(idnode,mxnode,ngrp,nbuff,gcmx,gcmy,gcmz,buffer)
        call merge(idnode,mxnode,ngrp,nbuff,gvxx,gvyy,gvzz,buffer)
      
c     
c     merge new atomic coordinates and velocities
      
        call merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
      
      endif
c
c     ensure all atoms are within cell boundaries

      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)

#ifdef VAMPIR
      call VTEND(40, ierr)
#endif
      return
      end


