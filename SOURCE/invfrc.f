      subroutine invfrc
     x  (idnode,imcon,mxnode,ntinv,enginv,virinv,keyinv,listinv,
     x  buffer,cell,fxx,fyy,fzz,prminv,xxx,yyy,zzz,xdab,ydab,zdab,
     x  xdac,ydac,zdac,xdad,ydad,zdad,stress)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating inversion energy and force 
c     terms in molecular dynamics.
c     
c     copyright - daresbury laboratory 1996
c     author    - w. smith       may   1996
c     
c     wl
c     2002/05/31 14:01:12
c     1.6
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical safe
      dimension keyinv(mxtinv),listinv(mxinv,5),prminv(mxtinv,mxpinv)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms),cell(9)
      dimension xdab(msbad),ydab(msbad),zdab(msbad)
      dimension xdac(msbad),ydac(msbad),zdac(msbad)
      dimension xdad(msbad),ydad(msbad),zdad(msbad)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension buffer(mxbuff),stress(9)
      
c     
c     define angular potential function and derivative
c     using the parameters in array prminv
      
c     
c     key=1 for harmonic inversion potential
      
      vharm(phi,k)=0.5d0*prminv(k,1)*(phi-prminv(k,2))**2
      dharm(phi,k)=prminv(k,1)*(phi-prminv(k,2))
c     
c     key=2 for harmonic cosine inversion potential
      
      vcos(cphi,k)=0.5d0*prminv(k,1)*(cphi-prminv(k,2))**2
      dcoz(cphi,k)=prminv(k,1)*(cphi-prminv(k,2))
c     
c     key=3 for planar inversion potentials
      
      vplan(cphi,k)=prminv(k,1)*(1-cphi)
      dplan(k)=prminv(k,1)
      
#ifdef VAMPIR
      call VTBEGIN(24, ierr)
#endif
c     
c     check size of work arrays
      
      if((ntinv-mxnode+1)/mxnode.gt.msbad) call error(idnode,427)
c     
c     block indices
      
      inv1 = (idnode*ntinv)/mxnode+1
      inv2 = ((idnode+1)*ntinv)/mxnode
      
      twopi = 2.d0*pi
      rtwopi=1.d0/twopi
      safe = .true.
      
#ifdef STRESS
c     
c     initialise stress tensor accumulators
      
      strs1 = 0.d0
      strs2 = 0.d0
      strs3 = 0.d0
      strs5 = 0.d0
      strs6 = 0.d0
      strs9 = 0.d0
#endif
c     
c     calculate bond vectors
      
      ii=0
      do i=inv1,inv2
        
        ii=ii+1
c     
c     indices of atoms involved
        
        ia=listinv(ii,2)
        ib=listinv(ii,3)
        ic=listinv(ii,4)
        id=listinv(ii,5)
        
c     
c     define components of bond vectors
        
        xdab(ii)=xxx(ib)-xxx(ia)
        ydab(ii)=yyy(ib)-yyy(ia)
        zdab(ii)=zzz(ib)-zzz(ia)
        
        xdac(ii)=xxx(ic)-xxx(ia)
        ydac(ii)=yyy(ic)-yyy(ia)
        zdac(ii)=zzz(ic)-zzz(ia)
        
        xdad(ii)=xxx(id)-xxx(ia)
        ydad(ii)=yyy(id)-yyy(ia)
        zdad(ii)=zzz(id)-zzz(ia)
        
      enddo
      
c     
c     periodic boundary condition
      
      call images(imcon,0,1,ii,cell,xdab,ydab,zdab)
      call images(imcon,0,1,ii,cell,xdac,ydac,zdac)
      call images(imcon,0,1,ii,cell,xdad,ydad,zdad)
      
c     
c     zero inversion energy accumulator
      
      enginv=0.d0
      virinv=0.d0
      
c     
c     loop over all specified inversions
      
      ii=0
      do i=inv1,inv2
        
c     
c     define components of bond vectors
        
        ii=ii+1
        
        xab=xdab(ii)
        yab=ydab(ii)
        zab=zdab(ii)
        rab2=xab*xab+yab*yab+zab*zab
        rrab=1.d0/sqrt(rab2)
        
        xac=xdac(ii)
        yac=ydac(ii)
        zac=zdac(ii)
        rac2=xac*xac+yac*yac+zac*zac
        rrac=1.d0/sqrt(rac2)
        
        xad=xdad(ii)
        yad=ydad(ii)
        zad=zdad(ii)
        rad2=xad*xad+yad*yad+zad*zad
        rrad=1.d0/sqrt(rad2)
        
        rbc=xab*xac+yab*yac+zab*zac
        rcd=xac*xad+yac*yad+zac*zad
        rdb=xad*xab+yad*yab+zad*zab
        
c     
c     calculate bond-angle-plane vectors
        
        ubx=xac*rrac+xad*rrad
        uby=yac*rrac+yad*rrad
        ubz=zac*rrac+zad*rrad
        ubn=1.d0/sqrt(ubx**2+uby**2+ubz**2)
        ubx=ubn*ubx
        uby=ubn*uby
        ubz=ubn*ubz
        rub=xab*ubx+yab*uby+zab*ubz
        
        vbx=xac*rrac-xad*rrad
        vby=yac*rrac-yad*rrad
        vbz=zac*rrac-zad*rrad
        vbn=1.d0/sqrt(vbx**2+vby**2+vbz**2)
        vbx=vbn*vbx
        vby=vbn*vby
        vbz=vbn*vbz
        rvb=xab*vbx+yab*vby+zab*vbz
        wwb=sqrt(rub**2+rvb**2)
        
        ucx=xad*rrad+xab*rrab
        ucy=yad*rrad+yab*rrab
        ucz=zad*rrad+zab*rrab
        ucn=1.d0/sqrt(ucx**2+ucy**2+ucz**2)
        ucx=ucn*ucx
        ucy=ucn*ucy
        ucz=ucn*ucz
        ruc=xac*ucx+yac*ucy+zac*ucz
        
        vcx=xad*rrad-xab*rrab
        vcy=yad*rrad-yab*rrab
        vcz=zad*rrad-zab*rrab
        vcn=1.d0/sqrt(vcx**2+vcy**2+vcz**2)
        vcx=vcn*vcx
        vcy=vcn*vcy
        vcz=vcn*vcz
        rvc=xac*vcx+yac*vcy+zac*vcz
        wwc=sqrt(ruc**2+rvc**2)
        
        udx=xab*rrab+xac*rrac
        udy=yab*rrab+yac*rrac
        udz=zab*rrab+zac*rrac
        udn=1.d0/sqrt(udx**2+udy**2+udz**2)
        udx=udn*udx
        udy=udn*udy
        udz=udn*udz
        rud=xad*udx+yad*udy+zad*udz
        
        vdx=xab*rrab-xac*rrac
        vdy=yab*rrab-yac*rrac
        vdz=zab*rrab-zac*rrac
        vdn=1.d0/sqrt(vdx**2+vdy**2+vdz**2)
        vdx=vdn*vdx
        vdy=vdn*vdy
        vdz=vdn*vdz
        rvd=xad*vdx+yad*vdy+zad*vdz
        wwd=sqrt(rud**2+rvd**2)
        
c     
c     calculate inversion angle cosines
        
        cosb=wwb*rrab
        cosc=wwc*rrac
        cosd=wwd*rrad
        
c     
c     select potential energy function type
        
        kk=listinv(ii,1)
        
c     
c     calculate potential energy and scalar force term
        
        if(keyinv(kk).eq.1)then
          
          thb=acos(cosb)
          thc=acos(cosc)
          thd=acos(cosd)
          enginv=enginv+
     x      (vharm(thb,kk)+vharm(thc,kk)+vharm(thd,kk))/3.d0
          gamb=0.d0
          if(abs(thb).gt.1.d-12)gamb=dharm(thb,kk)/(3.d0*sin(thb))
          gamc=0.d0
          if(abs(thc).gt.1.d-12)gamc=dharm(thc,kk)/(3.d0*sin(thc))
          gamd=0.d0
          if(abs(thd).gt.1.d-12)gamd=dharm(thd,kk)/(3.d0*sin(thd))
          
        else if(keyinv(kk).eq.2)then
          
          enginv=enginv+
     x      (vcos(cosb,kk)+vcos(cosc,kk)+vcos(cosd,kk))/3.d0
          gamb=-dcoz(cosb,kk)/3.d0
          gamc=-dcoz(cosc,kk)/3.d0
          gamd=-dcoz(cosd,kk)/3.d0
          
        else if(keyinv(kk).eq.3)then
          
          enginv=enginv+
     x      (vplan(cosb,kk)+vplan(cosc,kk)+vplan(cosd,kk))/3.d0
          gamb=dplan(kk)/3.d0
          gamc=dplan(kk)/3.d0
          gamd=dplan(kk)/3.d0
          
        else
c     
c     undefined potential
          
          safe = .false.
          gamb=0.d0
          gamc=0.d0
          gamd=0.d0
          
        endif
        
c     
c     indices of atoms involved
        
        ia=listinv(ii,2)
        ib=listinv(ii,3)
        ic=listinv(ii,4)
        id=listinv(ii,5)
        
c
c     calculate bond and u,v scalar products

        rubc=xab*ucx+yab*ucy+zab*ucz
        rubd=xab*udx+yab*udy+zab*udz
        rucd=xac*udx+yac*udy+zac*udz
        rucb=xac*ubx+yac*uby+zac*ubz
        rudb=xad*ubx+yad*uby+zad*ubz
        rudc=xad*ucx+yad*ucy+zad*ucz

        rvbc=xab*vcx+yab*vcy+zab*vcz
        rvbd=xab*vdx+yab*vdy+zab*vdz
        rvcd=xac*vdx+yac*vdy+zac*vdz
        rvcb=xac*vbx+yac*vby+zac*vbz
        rvdb=xad*vbx+yad*vby+zad*vbz
        rvdc=xad*vcx+yad*vcy+zad*vcz

c     
c     calculate atomic forces
        
        fbx = gamb*(-cosb*xab*rrab**2+rrab*(rub*ubx+rvb*vbx)/wwb)
     x    +(ruc*ucn*rrab*(xac-ruc*ucx-(rbc-ruc*rubc)*xab*rrab**2)
     x    - rvc*vcn*rrab*(xac-rvc*vcx-(rbc-rvc*rvbc)*xab*rrab**2))
     x    * gamc*rrac/wwc
     x    +(rud*udn*rrab*(xad-rud*udx-(rdb-rud*rubd)*xab*rrab**2)
     x    + rvd*vdn*rrab*(xad-rvd*vdx-(rdb-rvd*rvbd)*xab*rrab**2))
     x    * gamd*rrad/wwd
        
        fby = gamb*(-cosb*yab*rrab**2+rrab*(rub*uby+rvb*vby)/wwb)
     x    +(ruc*ucn*rrab*(yac-ruc*ucy-(rbc-ruc*rubc)*yab*rrab**2)
     x    - rvc*vcn*rrab*(yac-rvc*vcy-(rbc-rvc*rvbc)*yab*rrab**2))
     x    * gamc*rrac/wwc
     x    +(rud*udn*rrab*(yad-rud*udy-(rdb-rud*rubd)*yab*rrab**2)
     x    + rvd*vdn*rrab*(yad-rvd*vdy-(rdb-rvd*rvbd)*yab*rrab**2))
     x    * gamd*rrad/wwd
        
        fbz = gamb*(-cosb*zab*rrab**2+rrab*(rub*ubz+rvb*vbz)/wwb)
     x    +(ruc*ucn*rrab*(zac-ruc*ucz-(rbc-ruc*rubc)*zab*rrab**2)
     x    - rvc*vcn*rrab*(zac-rvc*vcz-(rbc-rvc*rvbc)*zab*rrab**2))
     x    * gamc*rrac/wwc
     x    +(rud*udn*rrab*(zad-rud*udz-(rdb-rud*rubd)*zab*rrab**2)
     x    + rvd*vdn*rrab*(zad-rvd*vdz-(rdb-rvd*rvbd)*zab*rrab**2))
     x    * gamd*rrad/wwd
        
        fcx = gamc*(-cosc*xac*rrac**2+rrac*(ruc*ucx+rvc*vcx)/wwc)
     x    +(rud*udn*rrac*(xad-rud*udx-(rcd-rud*rucd)*xac*rrac**2)
     x    - rvd*vdn*rrac*(xad-rvd*vdx-(rcd-rvd*rvcd)*xac*rrac**2))
     x    * gamd*rrad/wwd
     x    +(rub*ubn*rrac*(xab-rub*ubx-(rbc-rub*rucb)*xac*rrac**2)
     x    + rvb*vbn*rrac*(xab-rvb*vbx-(rbc-rvb*rvcb)*xac*rrac**2))
     x    * gamb*rrab/wwb
        
        fcy = gamc*(-cosc*yac*rrac**2+rrac*(ruc*ucy+rvc*vcy)/wwc)
     x    +(rud*udn*rrac*(yad-rud*udy-(rcd-rud*rucd)*yac*rrac**2)
     x    - rvd*vdn*rrac*(yad-rvd*vdy-(rcd-rvd*rvcd)*yac*rrac**2))
     x    * gamd*rrad/wwd
     x    +(rub*ubn*rrac*(yab-rub*uby-(rbc-rub*rucb)*yac*rrac**2)
     x    + rvb*vbn*rrac*(yab-rvb*vby-(rbc-rvb*rvcb)*yac*rrac**2))
     x    * gamb*rrab/wwb
        
        fcz = gamc*(-cosc*zac*rrac**2+rrac*(ruc*ucz+rvc*vcz)/wwc)
     x    +(rud*udn*rrac*(zad-rud*udz-(rcd-rud*rucd)*zac*rrac**2)
     x    - rvd*vdn*rrac*(zad-rvd*vdz-(rcd-rvd*rvcd)*zac*rrac**2))
     x    * gamd*rrad/wwd
     x    +(rub*ubn*rrac*(zab-rub*ubz-(rbc-rub*rucb)*zac*rrac**2)
     x    + rvb*vbn*rrac*(zab-rvb*vbz-(rbc-rvb*rvcb)*zac*rrac**2))
     x    * gamb*rrab/wwb
        
        fdx = gamd*(-cosd*xad*rrad**2+rrad*(rud*udx+rvd*vdx)/wwd)
     x    +(rub*ubn*rrad*(xab-rub*ubx-(rdb-rub*rudb)*xad*rrad**2)
     x    - rvb*vbn*rrad*(xab-rvb*vbx-(rdb-rvb*rvdb)*xad*rrad**2))
     x    * gamb*rrab/wwb
     x    +(ruc*ucn*rrad*(xac-ruc*ucx-(rcd-ruc*rudc)*xad*rrad**2)
     x    + rvc*vcn*rrad*(xac-rvc*vcx-(rcd-rvc*rvdc)*xad*rrad**2))
     x    * gamc*rrac/wwc
        
        fdy = gamd*(-cosd*yad*rrad**2+rrad*(rud*udy+rvd*vdy)/wwd)
     x    +(rub*ubn*rrad*(yab-rub*uby-(rdb-rub*rudb)*yad*rrad**2)
     x    - rvb*vbn*rrad*(yab-rvb*vby-(rdb-rvb*rvdb)*yad*rrad**2))
     x    * gamb*rrab/wwb
     x    +(ruc*ucn*rrad*(yac-ruc*ucy-(rcd-ruc*rudc)*yad*rrad**2)
     x    + rvc*vcn*rrad*(yac-rvc*vcy-(rcd-rvc*rvdc)*yad*rrad**2))
     x    * gamc*rrac/wwc
        
        fdz = gamd*(-cosd*zad*rrad**2+rrad*(rud*udz+rvd*vdz)/wwd)
     x    +(rub*ubn*rrad*(zab-rub*ubz-(rdb-rub*rudb)*zad*rrad**2)
     x    - rvb*vbn*rrad*(zab-rvb*vbz-(rdb-rvb*rvdb)*zad*rrad**2))
     x    * gamb*rrab/wwb
     x    +(ruc*ucn*rrad*(zac-ruc*ucz-(rcd-ruc*rudc)*zad*rrad**2)
     x    + rvc*vcn*rrad*(zac-rvc*vcz-(rcd-rvc*rvdc)*zad*rrad**2))
     x    * gamc*rrac/wwc
        
        fax = -(fbx+fcx+fdx)
        fay = -(fby+fcy+fdy)
        faz = -(fbz+fcz+fdz)
        
        fxx(ia)=fxx(ia)+fax
        fyy(ia)=fyy(ia)+fay
        fzz(ia)=fzz(ia)+faz
        
        fxx(ib)=fxx(ib)+fbx
        fyy(ib)=fyy(ib)+fby
        fzz(ib)=fzz(ib)+fbz
        
        fxx(ic)=fxx(ic)+fcx
        fyy(ic)=fyy(ic)+fcy
        fzz(ic)=fzz(ic)+fcz
        
        fxx(id)=fxx(id)+fdx
        fyy(id)=fyy(id)+fdy
        fzz(id)=fzz(id)+fdz

#ifdef STRESS
c     stress tensor calculation for inversion terms
        
        strs1 = strs1 + xab*fbx + xac*fcx + xad*fdx 
        strs2 = strs2 + yab*fbx + yac*fcx + yad*fdx 
        strs3 = strs3 + zab*fbx + zac*fcx + zad*fdx 
        strs5 = strs5 + yab*fby + yac*fcy + yad*fdy 
        strs6 = strs6 + yab*fbz + yac*fcz + yad*fdz 
        strs9 = strs9 + zab*fbz + zac*fcz + zad*fdz 
#endif

      enddo
      
#ifdef STRESS
c     
c     complete stress tensor
      
      stress(1) = stress(1) + strs1
      stress(2) = stress(2) + strs2
      stress(3) = stress(3) + strs3
      stress(4) = stress(4) + strs2
      stress(5) = stress(5) + strs5
      stress(6) = stress(6) + strs6
      stress(7) = stress(7) + strs3
      stress(8) = stress(8) + strs6
      stress(9) = stress(9) + strs9
#endif
c     
c     check for undefined potentials
      
      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,449)
      
c     
c     sum contributions over all nodes
      
      if(mxnode.gt.1) then
        
        buffer(2) = enginv

        call gdsum(buffer(2),1,buffer(1))
        
        enginv = buffer(2)
        
      endif

#ifdef VAMPIR
      call VTEND(24, ierr)
#endif
      return
      end




