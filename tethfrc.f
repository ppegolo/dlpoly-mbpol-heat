      subroutine tethfrc
     x  (idnode,mxnode,imcon,natms,nstep,ntteth,keytet,listtet,
     x  engtet,virtet,buffer,cell,fxx,fyy,fzz,prmtet,xxx,yyy,
     x  zzz,xxs,yys,zzs,xdab,ydab,zdab,stress)

c***********************************************************************
c     
c     dl_poly routine to tether atoms to initial positions
c     includes stress tensor
c     
c     replicated data version : block data
c     
c     copyright daresbury laboratory 1994
c     author     t.forester feb 1994
c     amended    t.forester dec 1994 : block data
c     
c     wl
c     2000/01/18 14:05:59
c     1.4
c     Exp
c     
c***********************************************************************

#include "dl_params.inc"

      logical safe
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension xxs(mxatms),yys(mxatms),zzs(mxatms)
      dimension xdab(msbad),ydab(msbad),zdab(msbad)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension prmtet(mxteth,mxpbnd),cell(9)
      dimension listtet(msteth,2),keytet(mxteth)
      dimension buffer(mxbuff)
      dimension stress(9)

      data safe/.true./
c     
c     define tethering functions using parameters in prmtet
c     
c     harmonic function

      vharm(r,k) = 0.5d0*prmtet(k,1)*r*r
      dharm(r,k) = prmtet(k,1)*r
c
c     restrained harmonic: 

      vt2(r,k) = 0.5d0*prmtet(k,1)*(min(r,prmtet(k,2)))**2
     x  + prmtet(k,1)*prmtet(k,2)*
     x  (sign(max(r-prmtet(k,2),0.d0),r))
      gt2(r,k) = prmtet(k,1)*(sign(min(r,prmtet(k,2)),r))

c
c     quartic potential

      vt3(r,k)=0.5d0*prmtet(k,1)*r**2 +
     x  1.d0/3.d0*prmtet(k,2)*r**3 + 
     x  0.25d0*prmtet(k,3)*r**4
      gt3(r,k)=prmtet(k,1)*r +
     x  prmtet(k,2)*r**2 +
     x  prmtet(k,3)*r**3 

#ifdef VAMPIR
      call VTBEGIN(25, ierr)
#endif
c     
c     set up reference positions at start of job

      if(nstep.le.1) then

        do i = 1,natms
          xxs(i) = xxx(i)
          yys(i) = yyy(i)
          zzs(i) = zzz(i)
        enddo

      endif
c
c     check size of work arrays

      if((ntteth-mxnode+1)/mxnode.gt.msbad) call error(idnode,420)
c     
c     block indices

      itet1 = (idnode*ntteth)/mxnode + 1
      itet2 = ((idnode+1)*ntteth)/mxnode
      
      ii = 0
      do i = itet1,itet2

        ii=ii+1
c     
c     atomic indice

        ia= listtet(ii,2)
c     
c     tether vector

        xdab(ii) = xxx(ia)-xxs(ia)
        ydab(ii) = yyy(ia)-yys(ia)
        zdab(ii) = zzz(ia)-zzs(ia)

      enddo

c     
c     ignore  periodic boundary condition
      
      call images(imcon,0,1,ii,cell,xdab,ydab,zdab)
c     
c     zero tether energy and virial accumulators
      
      engtet=0.d0
      virtet=0.d0
c     
c     loop over all specified tethered atoms

      ii=0
      do i=itet1,itet2
        
        ii=ii+1
c     
c     define components of bond vector
        
        rab=sqrt(xdab(ii)**2+ydab(ii)**2+zdab(ii)**2)

c     
c     check for possible zero length vector

        if(rab.lt.1.d-10) then

          rrab =0.d0

        else

          rrab = 1.d0/rab

        endif
c     
c     index of potential function parameters

        kk=listtet(ii,1)
c     
c     calculate scalar constant terms

        if(keytet(kk).eq.1)then

          omega=vharm(rab,kk)
          gamma=dharm(rab,kk)*rrab

        elseif(keytet(kk).eq.2)then

          omega=vt2(rab,kk)
          gamma=gt2(rab,kk)*rrab

        elseif(keytet(kk).eq.3)then

          omega=vt3(rab,kk)
          gamma=gt3(rab,kk)*rrab

        else
          safe = .false.
          omega=0.d0
          gamma=0.d0
        endif
        
        gamma = -gamma
c     
c     calculate tether energy and virial

        engtet=engtet+omega
        virtet=virtet-gamma*rab*rab
        
c     
c     indice of atom
        
        ia=listtet(ii,2)

c     
c     calculate atomic forces
        
        fxx(ia)=fxx(ia)+gamma*xdab(ii)
        fyy(ia)=fyy(ia)+gamma*ydab(ii)
        fzz(ia)=fzz(ia)+gamma*zdab(ii)
#ifdef STRESS
c     
c     stress tensor 

        stress(1) = stress(1) + xdab(ii)*gamma*xdab(ii)
        stress(2) = stress(2) + xdab(ii)*gamma*ydab(ii)
        stress(3) = stress(3) + xdab(ii)*gamma*zdab(ii)
        stress(4) = stress(4) + ydab(ii)*gamma*xdab(ii)
        stress(5) = stress(5) + ydab(ii)*gamma*ydab(ii)
        stress(6) = stress(6) + ydab(ii)*gamma*zdab(ii)
        stress(7) = stress(7) + zdab(ii)*gamma*xdab(ii)
        stress(8) = stress(8) + zdab(ii)*gamma*ydab(ii)
        stress(9) = stress(9) + zdab(ii)*gamma*zdab(ii)
#endif

      enddo

c
c     check for undefined potentials

      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,450)
c     
c     sum contributions to potential and virial

      if(mxnode.gt.1) then

        buffer(3)=engtet
        buffer(4)=virtet

        call gdsum(buffer(3),2,buffer(1))

        engtet=buffer(3)
        virtet=buffer(4)

      endif

#ifdef VAMPIR
      call VTEND(25, ierr)
#endif
      return
      end
