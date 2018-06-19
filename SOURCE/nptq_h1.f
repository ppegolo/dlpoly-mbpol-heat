      subroutine nptq_h1
     x (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,ntpatm,chip,chit,consv,conint,elrc,engke,
     x  engrot,virlrc,press,quattol,sigma,taup,taut,temp,tolnce,
     x  tstep,virtot,vircom,vircon,volm,lstme,lishap,lashap,
     x  lstfrz,listcon,listme,lstfre,lstgtp,lstrgd,listot,
     x  numgsit,cell,dens,dxx,dyy,dzz,dxt,dyt,dzt,gcmx1,gcmy1,
     x  gcmz1,gvx1,gvy1,gvz1,uxx,uyy,uzz,prmcon,xdf,ydf,zdf,omx,
     x  omy,omz,omx1,omy1,omz1,opx,opy,opz,oqx,oqy,oqz,q0,q1,q2,
     x  q3,qn0,qn1,qn2,qn3,gxx,gyy,gzz,xxx,yyy,zzz,vxx,vyy,vzz,
     x  gvxx,gvyy,gvzz,fxx,fyy,fzz,tqx,tqy,tqz,weight,rotinx,
     x  rotiny,rotinz,gcmx,gcmy,gcmz,gmass,buffer,xxt,yyt,zzt,
     x  txx,tyy,tzz,vx1,vy1,vz1,stress,eta,dens0)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using implicit leapfrog quaternion algorithm
c     plus bond constraints - PROVIDED rigid body sites
c     and constraint sites do not coincide.
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Hoover
c     thermostat and barostat (Melchionna et al variant)
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
c     2001/05/30 12:40:15
c     1.8
c     $Sate: Exp $
c     
c**********************************************************************
      
#include "dl_params.inc"
      
      logical safe,safeq,lshmov,newjob
      
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension rot(9),cell(9),cell0(9),celn(9)
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
      dimension stress(9),stres1(9),eta(9)

      save newjob,volm0,elrc0,virlrc0,cell0
      
      data newjob/.true./
      data pt5/0.5d0/
#ifdef VAMPIR
      call VTBEGIN(42, ierr)
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

      endif

      nbuff=mxbuff
      vircon = 0.d0
#ifdef STRESS
c
c     temporary stress tensor accumulators
      do i = 1,9
        stres1(i) = 0.d0
      enddo
#endif
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
c     inertia parameter for Nose-Hoover thermostat
      
c     qmass = 2.0d0*sigma*taut**2
      pmass = natms*boltz*temp*taup**2

c     
c     store initial values of position and velocity

      totmas = 0.d0
      j = 0
      do ifre=ifre1,ifre2
        
        i=lstfre(ifre)
        j = j+1
        
        xdf(j)=xxx(i)
        ydf(j)=yyy(i)
        zdf(j)=zzz(i)
        totmas = totmas+weight(i)

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
        totmas = totmas+gmass(lstgtp(ig))

        gcmx1(jg) = gcmx(ig)
        gcmy1(jg) = gcmy(ig)
        gcmz1(jg) = gcmz(ig)
        
        gvx1(jg) = gvxx(ig)
        gvy1(jg) = gvyy(ig)
        gvz1(jg) = gvzz(ig)

      enddo

      if(mxnode.gt.1) call gdsum(totmas,1,buffer)
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
c     calculate centre of mass

      xcmo = 0.d0
      ycmo = 0.d0
      zcmo = 0.d0
      
      j=0
      iatm0 = (idnode*natms)/mxnode + 1
      iatm1 = ((idnode+1)*natms)/mxnode
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
c     first estimate of new velocities - no thermostat
        
        tmp = tstep/gmass(id)
        gvxx(ig) = gvxx(ig) + tmp*qn0(jg)
        gvyy(ig) = gvyy(ig) + tmp*qn1(jg)
        gvzz(ig) = gvzz(ig) + tmp*qn2(jg)
c     
c     first estimate of new positions

        gcmx(ig) = gcmx(ig) + tstep*gvxx(ig)
        gcmy(ig) = gcmy(ig) + tstep*gvyy(ig)
        gcmz(ig) = gcmz(ig) + tstep*gvzz(ig)
        
      enddo
      engtrn = engtrn*pt5
c     
c     estimate rotational kinetic energy
      
      engrot = 0.d0
      jg = 0
      jr = 0
      do ig = igrp1,igrp2
        
        id = lstgtp(ig)
        jg =jg+1

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
c     
c     first estimate of new angular velocities - no thermostat

        omx(ig) = omx(ig) + tstep*trx
        omy(ig) = omy(ig) + tstep*try
        omz(ig) = omz(ig) + tstep*trz

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
        
      endif
c     
c     propagate chip

      engke = engke + engtrn
      chipp=((2.d0*engke-virtot-vircon-vircom)/3.d0-press*volm)/pmass
      chipnew = chip + tstep*chipp
      chip0 = 0.5d0*(chip+chipnew)
c     
c     propagate chit

      engtot = engke+engrot
      chitp = (engtot/sigma - 1.d0)/taut**2
      chitnew = chit + tstep*chitp
      chit0 = 0.5d0*(chit+chitnew)
c     
c     begin iterations !!-----------------------------------------------

      mxiter = 5
      if(ntcons.eq.0) mxiter=mxiter -1
      
      do iter = 1,mxiter
c     
c     unconstrained new positions
        
        j = 0
        do ifre = ifre1,ifre2

          i = lstfre(ifre)
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
c     estimate new cell tensor

          volnew = volm*exp(3.d0*tstep*chipnew)
          scale= (volnew/volm0)**(1.d0/3.d0)
          do i = 1,9
            celn(i) = cell0(i)*scale
          enddo
c     
c     apply constraint correction
          
          call rdshake_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,viracc,lashap,lishap,listcon,listme,
     x      listot,lstfrz,buffer,celn,dxt,dxx,dyt,dyy,dzt,dzz,
     x      prmcon,txx,tyy,tzz,weight,xxt,xxx,yyt,yyy,zzt,zzz,
     x      stres1)
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

          enddo
c     
c     end of shake corrections
          
        endif
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

          engke = engke + weight(i)*(vxt**2+vyt**2+vzt**2)

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
c     
c     advance velocity using leapfrog

          gvxx(ig) = gvx1(jg) + tstep*(qn0(jg)*rmass -(chit0+chip0)*
     x      0.5d0*(gvxx(ig)+gvx1(jg)))
          gvyy(ig) = gvy1(jg) + tstep*(qn1(jg)*rmass -(chit0+chip0)*
     x      0.5d0*(gvyy(ig)+gvy1(jg)))
          gvzz(ig) = gvz1(jg) + tstep*(qn2(jg)*rmass -(chit0+chip0)*
     x      0.5d0*(gvzz(ig)+gvz1(jg)))
c     
c     advance positions using leapfrog

          gcmx(ig) = gcmx1(jg) + tstep*(gvxx(ig) +
     x      chipnew*((gcmx(ig)+gcmx1(jg))*0.5d0 - xcmo))
          gcmy(ig) = gcmy1(jg) + tstep*(gvyy(ig) +
     x      chipnew*((gcmy(ig)+gcmy1(jg))*0.5d0 - ycmo))
          gcmz(ig) = gcmz1(jg) + tstep*(gvzz(ig) +
     x      chipnew*((gcmz(ig)+gcmz1(jg))*0.5d0 - zcmo))
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
          
        enddo
c     
c     ****** step 2 : integrate rotational motion **********
        
        safeq = .true.
        engrot = 0.d0
        
        jg = 0
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
          
          delx = tstep*(trx - chit0*pt5*(omx(ig)+omx1(jg)))
          dely = tstep*(try - chit0*pt5*(omy(ig)+omy1(jg)))
          delz = tstep*(trz - chit0*pt5*(omz(ig)+omz1(jg)))
c     
c     angular velocity at time step n
          
          opx(jg) = omx1(jg) + delx*pt5
          opy(jg) = omy1(jg) + dely*pt5
          opz(jg) = omz1(jg) + delz*pt5
c     
c     angular velocity at time step n+1/2
          
          omx(ig) = omx1(jg) + delx
          omy(ig) = omy1(jg) + dely
          omz(ig) = omz1(jg) + delz
c     
c     angular velocity at time step n+1 (needed for quat algorithm)
          
          oqx(jg) = omx1(jg) + delx*1.5d0
          oqy(jg) = omy1(jg) + dely*1.5d0
          oqz(jg) = omz1(jg) + delz*1.5d0
c     
c     rotational kinetic energy
          
          engrot = engrot+ pt5*(rotinx(id,1)*opx(jg)**2
     x      +rotiny(id,1)*opy(jg)**2
     x      +rotinz(id,1)*opz(jg)**2)
          
        enddo
c
c     global sum of kinetic energy terms

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
c     improved prediction of chip and chit 

        engke = engke+engtrn
        chipp=((2.d0*engke-virtot-vircom-vircon)/3.d0-press*volm)/pmass
        chipnew = chip + tstep*chipp
        chip0 = 0.5d0*(chip+chipnew)

        engtot = engke+engrot
        chitp = (engtot/sigma - 1.d0)/taut**2
        chitnew = chit + tstep*chitp
        chit0 = 0.5d0*(chit+chitnew)

c     
c     end of thermostat iterations

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
c     ensure total momentum is zero

      vxm = 0.d0
      vym = 0.d0
      vzm = 0.d0
      do ifre = ifre1,ifre2
        i = lstfre(ifre)
        vxm  = vxm + vxx(i)*weight(i)
        vym  = vym + vyy(i)*weight(i)
        vzm  = vzm + vzz(i)*weight(i)
      enddo
      do ig = igrp1,igrp2

        id = lstgtp(ig)
        vxm = vxm + gvxx(ig)*gmass(id)
        vym = vym + gvyy(ig)*gmass(id)
        vzm = vzm + gvzz(ig)*gmass(id)
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

      do ifre = ifre1,ifre2
        i = lstfre(ifre)
        vxx(i) = vxx(i) - vxm
        vyy(i) = vyy(i) - vym
        vzz(i) = vzz(i) - vzm
      enddo
      do ig = igrp1,igrp2
        gvxx(ig) = gvxx(ig) - vxm
        gvyy(ig) = gvyy(ig) - vym
        gvzz(ig) = gvzz(ig) - vzm
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
c     xs
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

#ifdef STRESS
c
c     kinetic contribution to stress tensor
      j = 0
      do ifre = ifre1,ifre2
        i = lstfre(ifre)
        j=j+1
        vxt = 0.5d0*(vxx(i) + vx1(j))
        vyt = 0.5d0*(vyy(i) + vy1(j))
        vzt = 0.5d0*(vzz(i) + vz1(j))

        stres1(1) = stres1(1) + weight(i)*vxt*vxt
        stres1(2) = stres1(2) + weight(i)*vxt*vyt
        stres1(3) = stres1(3) + weight(i)*vxt*vzt
        stres1(5) = stres1(5) + weight(i)*vyt*vyt
        stres1(6) = stres1(6) + weight(i)*vyt*vzt
        stres1(9) = stres1(9) + weight(i)*vzt*vzt
      enddo

      jg = 0
      do ig = igrp1,igrp2
          
        jg = jg+1
        id = lstgtp(ig)
          
        vxt = pt5*(gvx1(jg)+gvxx(ig))
        vyt = pt5*(gvy1(jg)+gvyy(ig))
        vzt = pt5*(gvz1(jg)+gvzz(ig))
          
        stres1(1) = stres1(1) + gmass(id)*vxt*vxt
        stres1(2) = stres1(2) + gmass(id)*vxt*vyt
        stres1(3) = stres1(3) + gmass(id)*vxt*vzt
        stres1(5) = stres1(5) + gmass(id)*vyt*vyt
        stres1(6) = stres1(6) + gmass(id)*vyt*vzt
        stres1(9) = stres1(9) + gmass(id)*vzt*vzt
      enddo
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
c
c     symmetrise  stress tensor

      stress(2) = 0.5d0*(stress(2)+stress(4))
      stress(4) = stress(2)
      stress(3) = 0.5d0*(stress(3)+stress(7))
      stress(7) = stress(3)
      stress(6) = 0.5d0*(stress(6)+stress(8))
      stress(8) = stress(6)
#endif
#ifdef VAMPIR
      call VTEND(42, ierr)
#endif
      return
      end
