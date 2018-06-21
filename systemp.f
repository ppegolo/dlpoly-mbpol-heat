! 02 AUG 06 - IUCHI - REPLACE KEYRES.GT.1 BY KEYRES.EQ.2
!
      subroutine systemp
     x  (idnode,imcon,keyens,keyres,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,ntshl,levcfg,degfre,degshl,degrot,temp,
     x  tolnce,lashap,lishap,listcon,listme,listot,lstfrz,lstgtp,
     x  lstrgd,lstme,numgsit,listshl,buffer,cell,dxt,dyt,dzt,
     x  fxx,fyy,fzz,gcmx,gcmy,gcmz,gmass,gvxx,gvyy,gvzz,
     x  gxx,gyy,gzz,q0,q1,q2,q3,rotinx,rotiny,rotinz,weight,
     x  uxx,uyy,uzz,vxx,vyy,vzz,xxt,xxx,yyt,yyy,zzt,zzz,
     x  omx,omy,omz,lttm,nfict)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for setting the initial system temperature
c     
c     copyright - daresbury laboratory 1997
c     author    - w. smith         july 1997
c     
c     wl
c     2001/08/31 11:13:53
c     1.6
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"

      logical lttm
      dimension listshl(mxshl,3),lstme(mxatms)
      dimension listcon(mxcons,3),lstfrz(mxatms)
      dimension listme(mxatms),lashap(mxproc)
      dimension lishap(mxlshp),listot(mxatms)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension uxx(mxatms),uyy(mxatms),uzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension dxt(mxcons),dyt(mxcons),dzt(mxcons)
      dimension omx(mxgrp),omy(mxgrp),omz(mxgrp)
      dimension gcmx(mxgrp),gcmy(mxgrp),gcmz(mxgrp)
      dimension gvxx(mxgrp),gvyy(mxgrp),gvzz(mxgrp),gmass(mxungp)
      dimension q0(mxgrp),q1(mxgrp),q2(mxgrp),q3(mxgrp)
      dimension gxx(mxungp,mxngp),gyy(mxungp,mxngp),gzz(mxungp,mxngp)
      dimension rotinx(mxungp,2),rotiny(mxungp,2),rotinz(mxungp,2)
      dimension cell(9),weight(mxatms),buffer(mxbuff)
      dimension lstrgd(mxgatm),numgsit(mxungp),lstgtp(mxgrp)
#ifdef VAMPIR
      call VTBEGIN(6, ierr)
#endif
c     
c     number of degrees of freedom 
c     3 for com translation
c     3 for angular momentum about origin (non-periodic systems only)
      
      degfre = dble(3*(ntfree-ntshl)-3-ntcons)+degfre
      if(imcon.eq.0.or.imcon.eq.6) degfre = degfre - 3.0d0
      if(imcon.eq.0.or.imcon.eq.6) degrot = max(0.d0,degrot - 3.0d0)
      degshl= dble(3*ntshl)
c
c     remove m-site degree of freedom for ttm2 water

      if(lttm)degfre=degfre-3*nfict
c     
c     lose one degree of freedom if temperature constrained
c     gaussian constraints
      
c     if(keyens.eq.1) degfre = degfre - 1.d0
      
      if(idnode.eq.0) 
     x  write(nrite,"(/,/,' total degrees of freedom       ',f20.0,/,
     x  ' rotational degrees of freedom  ',f20.0,/,
     x  ' shell pseudo degrees of freedom',f20.0)")
     x  degfre,degrot,degshl
      if(degfre.lt.1.d0) call error(idnode,350)
      
c     
c     generate starting velocities
      
      sigma=temp*boltz*degfre*0.5d0
      
      if(keyres.eq.0)then
        
        call gauss(natms,vxx,vyy,vzz)
        
        do i=1,natms
          
          if(weight(i).eq.0.d0) then
            
            rsq = 0.d0
            
          else
            
            rsq = 1.d0/sqrt(weight(i))
            
          endif
          
          vxx(i)=vxx(i)*rsq
          vyy(i)=vyy(i)*rsq
          vzz(i)=vzz(i)*rsq
          
        enddo
        
        if(ntcons.gt.0)call quench
     x    (imcon,idnode,mxnode,natms,nscons,tolnce,
     x    lashap,lishap,listcon,listme,listot,
     x    buffer,cell,dxt,dyt,dzt,uxx,uyy,uzz,vxx,vyy,vzz,
     x    weight,xxt,xxx,yyt,yyy,zzt,zzz)
        
        if(ngrp.gt.0) call quatqnch
     x    (idnode,imcon,mxnode,natms,ngrp,lstgtp,lstrgd,numgsit,
     x    lstme,gxx,gyy,gzz,buffer,cell,xxt,yyt,zzt,
     x    gcmx,gcmy,gcmz,gmass,gvxx,gvyy,gvzz,q0,q1,q2,q3,omx,omy,
     x    omz,rotinx,rotiny,rotinz,vxx,vyy,vzz,weight,xxx,yyy,zzz)
        
        if(ntshl.gt.0) then
          
          do k=1,4
            
            call vscaleg
     x        (idnode,mxnode,imcon,natms,ngrp,sigma,lstfrz,lstgtp,
     x        lstrgd,numgsit,lstme,buffer,cell,gcmx,gcmy,gcmz,gmass,
     x        gvxx,gvyy,gvzz,gxx,gyy,gzz,q0,q1,q2,q3,omx,omy,omz,weight,
     x        vxx,vyy,vzz,xxx,yyy,zzz,xxt,yyt,zzt,rotinx,rotiny,rotinz)
            
            call shlqnch
     x        (idnode,mxnode,ntshl,temp,listshl,weight,vxx,vyy,vzz,
     x        buffer)
            
          enddo
          
        else
          
          call vscaleg
     x      (idnode,mxnode,imcon,natms,ngrp,sigma,lstfrz,lstgtp,
     x      lstrgd,numgsit,lstme,buffer,cell,gcmx,gcmy,gcmz,gmass,
     x      gvxx,gvyy,gvzz,gxx,gyy,gzz,q0,q1,q2,q3,omx,omy,omz,weight,
     x      vxx,vyy,vzz,xxx,yyy,zzz,xxt,yyt,zzt,rotinx,rotiny,rotinz)

        endif
        
      elseif(keyres.eq.1) then 
        
        if(ngrp.gt.0) call quatqnch
     x    (idnode,imcon,mxnode,natms,ngrp,lstgtp,lstrgd,numgsit,
     x    lstme,gxx,gyy,gzz,buffer,cell,xxt,yyt,zzt,
     x    gcmx,gcmy,gcmz,gmass,gvxx,gvyy,gvzz,q0,q1,q2,q3,omx,omy,
     x    omz,rotinx,rotiny,rotinz,vxx,vyy,vzz,weight,xxx,yyy,zzz)
        
      else if(keyres.gt.1) then
        
        if(ngrp.gt.0) then 
          
          call vscaleg
     x      (idnode,mxnode,imcon,natms,ngrp,sigma,lstfrz,lstgtp,
     x      lstrgd,numgsit,lstme,buffer,cell,gcmx,gcmy,gcmz,gmass,
     x      gvxx,gvyy,gvzz,gxx,gyy,gzz,q0,q1,q2,q3,omx,omy,omz,weight,
     x      vxx,vyy,vzz,xxx,yyy,zzz,xxt,yyt,zzt,rotinx,rotiny,rotinz)
          
        elseif(ntshl.gt.0) then
          
          do k=1,4
            
            call vscaleg
     x        (idnode,mxnode,imcon,natms,ngrp,sigma,lstfrz,lstgtp,
     x        lstrgd,numgsit,lstme,buffer,cell,gcmx,gcmy,gcmz,gmass,
     x        gvxx,gvyy,gvzz,gxx,gyy,gzz,q0,q1,q2,q3,omx,omy,omz,weight,
     x        vxx,vyy,vzz,xxx,yyy,zzz,xxt,yyt,zzt,rotinx,rotiny,rotinz)

            call shlqnch
     x        (idnode,mxnode,ntshl,temp,listshl,weight,vxx,vyy,vzz,
     x        buffer)
            
          enddo
          
        else
 
          call vscaleg
     x      (idnode,mxnode,imcon,natms,ngrp,sigma,lstfrz,lstgtp,
     x      lstrgd,numgsit,lstme,buffer,cell,gcmx,gcmy,gcmz,gmass,
     x      gvxx,gvyy,gvzz,gxx,gyy,gzz,q0,q1,q2,q3,omx,omy,omz,weight,
     x      vxx,vyy,vzz,xxx,yyy,zzz,xxt,yyt,zzt,rotinx,rotiny,rotinz)

        endif

      endif
      
c     
c     print out sample of initial configuration 
      
      if(idnode.eq.0) write(nrite,
     x  "(/,/,1x,'sample of starting configuration',/)")
      
      io=(natms+19)/20
      if((levcfg.le.1).and.(idnode.eq.0)) 
     x  write(nrite,"(6x,'i',7x,'x(i)',8x,'y(i)',8x,'z(i)',
     x  7x,'vx(i)',7x,'vy(i)',7x,'vz(i)',/,/)")
      if((levcfg.eq.2).and.(idnode.eq.0)) 
     x  write(nrite,"(6x,'i',7x,'x(i)',8x,'y(i)',8x,'z(i)',
     x  7x,'vx(i)',7x,'vy(i)',7x,'vz(i)',
     x  7x,'fx(i)',7x,'fy(i)',7x,'fz(i)',/,/)")

      do i=1,natms,io
        
        if(levcfg.le.1)then
          
          if(idnode.eq.0) write(nrite,
     x      "(1x,i6,1p,3e12.4,3e12.4,3e12.4)")
     x      i,xxx(i),yyy(i),zzz(i),vxx(i),vyy(i),vzz(i)
          
        else if(levcfg.eq.2)then
          
          if(idnode.eq.0) write(nrite,
     x      "(1x,i6,1p,3e12.4,3e12.4,3e12.4)")
     x      i,xxx(i),yyy(i),zzz(i),
     x      vxx(i),vyy(i),vzz(i),fxx(i),fyy(i),fzz(i)
          
        endif
        
      enddo
#ifdef VAMPIR
      call VTEND(6, ierr)
#endif
      return
      end
