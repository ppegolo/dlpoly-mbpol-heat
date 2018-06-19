! 24 NOV 06 - IUCHI - COMMENT OUT VX,Y,Z0_AT_T
! 12 NOV 06 - IUCHI - INTRODUCE VX,Y,Z0_AT_T FOR WITHOUT SHAKE
! 14 SEP 06 - IUCHI - INTRODUCE VX,Y,Z_AT_T FOR VELOCITY BELONGED TO NT
! 
      subroutine nve_1
     x  (safe,lshmov,idnode,imcon,mxnode,natms,nscons,ntcons,
     x  engke,tolnce,tstep,vircon,lashap,lishap,listcon,
     x  listme,listot,lstfrz,buffer,cell,dxt,dxx,dyt,
     x  dyy,dzt,dzz,fxx,fyy,fzz,prmcon,txx,tyy,tzz,uxx,
     x  uyy,uzz,vxx,vyy,vzz,weight,xdf,xxt,xxx,ydf,yyt,
     x  yyy,zdf,zzt,zzz,stress,lttm)
c     
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics. Verlet leapfrog With RD-SHAKE
c     
c     parallel replicated data version : block data
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith august 1992.
c     amended   - t.forester sept 1994
c     amended   - t.forester dec  1994 : block data
c     
c     wl
c     2000/01/18 14:05:47
c     1.5
c     Exp
!
!     Last updated: 24 Nov 2006
c     
c***********************************************************************
c     
!  from MODULE
      use ttm_forces, only: vx_at_t, vy_at_t, vz_at_t
c$$$     $     ,vx0_at_t, vy0_at_t, vz0_at_t

#include "dl_params.inc"
      
      logical safe,lshmov

!  for TTM2 case
      logical lttm

      dimension listme(mxatms),lishap(mxlshp),lstfrz(mxatms)
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
      dimension stress(9),stres1(9)
      
#ifdef VAMPIR
      call VTBEGIN(33, ierr)
#endif
      safe=.false.
      nbuff=mxbuff
c     
c     block indices

      iatm0 = (idnode*natms)/mxnode + 1
      iatm1 = ((idnode+1)*natms)/mxnode
c     
c     constraint virial

      vircon=0.d0
#ifdef STRESS
c
c     temporary stress tensor accumulators
      do i = 1,9
        stres1(i) = 0.d0
      enddo
#endif
c     
c     store initial values of position
      
      j=0
      do i=iatm0,iatm1
        
        j=j+1
        xdf(j)=xxx(i)
        ydf(j)=yyy(i)
        zdf(j)=zzz(i)
        
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
c     move atoms by leapfrog algorithm
      
      j = 0
      do i=iatm0,iatm1
        j = j+1
c     
c     update velocities
        
        uxx(i)=vxx(i) + tstep*fxx(i)/weight(i)
        uyy(i)=vyy(i) + tstep*fyy(i)/weight(i)
        uzz(i)=vzz(i) + tstep*fzz(i)/weight(i)
c     
c     update positions
        
        xxx(i)=xdf(j)+tstep*uxx(i)
        yyy(i)=ydf(j)+tstep*uyy(i)
        zzz(i)=zdf(j)+tstep*uzz(i)
        
      enddo

c$$$! current velocity without SHAKE
c$$$      if( lttm ) then
c$$$         do i=1,natms
c$$$            vx0_at_t(i) = 0.5d0 * ( uxx(i) + vxx(i) )
c$$$            vy0_at_t(i) = 0.5d0 * ( uyy(i) + vyy(i) )
c$$$            vz0_at_t(i) = 0.5d0 * ( uzz(i) + vzz(i) )
c$$$         enddo
c$$$      endif
      
      if(ntcons.eq.0) safe =.true.
      if(ntcons.gt.0) then
c     
c     RDSHAKE procedure 
c     
c     global exchange of configuration data
        
        if(mxnode.gt.1) then
          
          call merge(idnode,mxnode,natms,nbuff,xxx,yyy,zzz,buffer)
          
        endif
c     
c     apply constraint correction
        
        call rdshake_1
     x    (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x    tolnce,tstep,vircon,lashap,lishap,listcon,listme,
     x    listot,lstfrz,buffer,cell,dxt,dxx,dyt,dyy,dzt,dzz,
     x    prmcon,txx,tyy,tzz,weight,xxt,xxx,yyt,yyy,zzt,zzz,
     x    stres1)
c     
c     calculate velocity correction
        
        rstep = 1.d0/tstep
        j=0
        do i=iatm0,iatm1
          
c     
c     update corrected velocity
          
          j=j+1
          uxx(i)=(xxx(i)-xdf(j))*rstep
          uyy(i)=(yyy(i)-ydf(j))*rstep
          uzz(i)=(zzz(i)-zdf(j))*rstep
          
c     
c     calculate the corrected forces
          
          fxx(i)=(uxx(i)-vxx(i))*weight(i)*rstep
          fyy(i)=(uyy(i)-vyy(i))*weight(i)*rstep
          fzz(i)=(uzz(i)-vzz(i))*weight(i)*rstep

        enddo
        
      endif
c     
c     updated velocity and kinetic energy
      
      engke =0.d0
      do i = iatm0,iatm1

        vxt = 0.5d0*(vxx(i)+uxx(i))
        vyt = 0.5d0*(vyy(i)+uyy(i))
        vzt = 0.5d0*(vzz(i)+uzz(i))

        engke=engke+0.5d0*weight(i)*(vxt*vxt+vyt*vyt+vzt*vzt)

        vxx(i) = uxx(i)
        vyy(i) = uyy(i)
        vzz(i) = uzz(i)
! 
!      current velocities

        if( lttm ) then
           
           vx_at_t(i) = vxt
           vy_at_t(i) = vyt
           vz_at_t(i) = vzt

        endif

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
c     periodic boundary condition
      
      call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
      
c     
c     global exchange of configuration data
      
      if(mxnode.gt.1)then

        nbuff=mxbuff
        call gdsum(engke,1,buffer)
        call merge(idnode,mxnode,natms,nbuff,xxx,yyy,zzz,buffer)
        call merge(idnode,mxnode,natms,nbuff,vxx,vyy,vzz,buffer)
        
        if(ntcons.gt.0) then
          
          call merge(idnode,mxnode,natms,nbuff,fxx,fyy,fzz,buffer)
          
        endif
        
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
      call VTEND(33, ierr)
#endif
      return
      end



