      subroutine fbpfrc
     x  (idnode,mxnode,natms,imcon,rcutfb,engfbp,virfbp,
     x  latinx,ltype,lst,lct,link,lstfbp,ltpfbp,prmfbp,cell,
     x  xxx,yyy,zzz,fxx,fyy,fzz,rcut4b,stress,buffer)
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating four body inversion forces
c     arising from the inversion angle between three atoms around a
c     nominated central atom
c
c     Note: the subroutine converts coordinates to reduced units
c     to avoid a call to images.f. The link cell algorithm used
c     here necessitates a parallelepiped cell geometry
c     
c     copyright - daresbury laboratory 1996
c     author   - w.smith july 1996
c     
c     wl
c     2002/05/31 13:59:31
c     1.7
c     Exp
c     
c***********************************************************************
c     

#include "dl_params.inc"
      
      logical safe
      dimension lst(mxcell),lct(mxcell),link(mxatms)
      dimension nix(27),niy(27),niz(27)
      dimension latinx(mxatms),ltype(mxatms)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension lstfbp(mxfbp),ltpfbp(mxfbp),prmfbp(mxfbp,mxpfbp)
      dimension cell(9),rcut4b(mxfbp),stress(9),buffer(mxbuff)
      dimension rcell(9),cprp(10)
      
      data nix/ 0,-1,-1,-1, 0, 0,-1, 1,-1, 0, 1,-1, 0, 1,
     x  1, 1, 1, 0, 0, 1,-1, 1, 0,-1, 1, 0,-1/
      data niy/ 0, 0,-1, 1, 1, 0, 0, 0,-1,-1,-1, 1, 1, 1,
     x  0, 1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1/
      data niz/ 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     x  0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1/
      
c     
c     key=1 for harmonic inversion potential
      
      vharm(phi,k)=0.5d0*prmfbp(k,1)*(phi-prmfbp(k,2))**2
      dharm(phi,k)=prmfbp(k,1)*(phi-prmfbp(k,2))
c     
c     key=2 for harmonic cosine inversion potential
      
      vcos(cphi,k)=0.5d0*prmfbp(k,1)*(cphi-prmfbp(k,2))**2
      dcoz(cphi,k)=prmfbp(k,1)*(cphi-prmfbp(k,2))
c     
c     key=3 for planar inversion potentials
      
      vplan(cphi,k)=prmfbp(k,1)*(1-cphi)
      dplan(k)=prmfbp(k,1)
      
#ifdef VAMPIR
      call VTBEGIN(20, ierr)
#endif
c
c     flag for undefined potentials

      safe = .true.

c     
c     initialise potential energy and virial

      engfbp=0.d0
      virfbp=0.d0
c
c     create mock cell vectors for non-periodic system
      if(imcon.eq.0) then
        xm = 0.d0
        ym = 0.d0
        zm = 0.d0
        do i = 1,natms
          xm = max(xm,abs(xxx(i)))
          ym = max(ym,abs(yyy(i)))
          zm = max(zm,abs(zzz(i)))
        enddo

        cell(1) = 2.d0*xm+rcutfb
        cell(2) = 0.d0
        cell(3) = 0.d0
        cell(4) = 0.d0
        cell(5) = 2.d0*ym+rcutfb
        cell(6) = 0.d0
        cell(7) = 0.d0
        cell(8) = 0.d0
        cell(9) = 2.d0*zm+rcutfb
      
      endif
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
c     check for appropriate boundary conditions
      
      if(imcon.gt.3)call error(idnode,79)
      call invert(cell,rcell,det)
      call dcell(cell,cprp)
c     
c     calculate link cell numbers
      
      nbx=int(cprp(7)/(rcutfb+1.d-6))
      nby=int(cprp(8)/(rcutfb+1.d-6))
      nbz=int(cprp(9)/(rcutfb+1.d-6))
      ncells=nbx*nby*nbz
      if(ncells.gt.mxcell) then

        if(idnode.eq.0) write(nrite,'(a,i6)')
     x    'number of required link cells in routine fbpfrc is ',ncells
          write(nrite,'(a,i6)')
     x    'number of default link cells in routine fbpfrc is ',mxcell
        call error(idnode,87)

      endif

c     
c     transform atomic coordinates and construct link cells
      
      do l=1,ncells
        
        lct(l)=0
        lst(l)=0
        
      enddo
      
      xdc=dble(nbx)
      ydc=dble(nby)
      zdc=dble(nbz)

      do i=1,natms

        sxx=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
        syy=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
        szz=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)
        
        xxx(i)=sxx
        yyy(i)=syy
        zzz(i)=szz

        ix=int(xdc*(sxx+0.5d0))
        iy=int(ydc*(syy+0.5d0))
        iz=int(zdc*(szz+0.5d0))
        k=1+ix+nbx*(iy+nby*iz)
        lst(k)=lst(k)+1
        link(i)=lct(k)
        lct(k)=i
        
      enddo

c     
c     loop over central atoms of inversion
      
      ix=0
      iy=1
      iz=1
      do icell=1,ncells

        ix=ix+1
        if(ix.gt.nbx)then
          ix=1
          iy=iy+1
          if(iy.gt.nby)then
            iy=1
            iz=iz+1
          endif
        endif
        
        k=0
        do kk=1,27

          jx=ix+nix(kk)
          jy=iy+niy(kk)
          jz=iz+niz(kk)
          
          if(jx.gt.nbx)jx=1
          if(jy.gt.nby)jy=1
          if(jz.gt.nbz)jz=1
          if(jx.lt.1)jx=jx+nbx
          if(jy.lt.1)jy=jy+nby
          if(jz.lt.1)jz=jz+nbz
          
          jcell=jx+nbx*(jy-1+nby*(jz-1))
          j=lct(jcell)
          
          do ii=1,lst(jcell)
            
            k=k+1
            latinx(k)=j
            j=link(j)
            
          enddo
          
        enddo

        limit=k
        
        do ii=1,lst(icell)
          
          ia=latinx(ii)
          ifbp=mx3fbp*(ltype(ia)-1)
          if(mod(ia,mxnode).eq.idnode.and.ltpfbp(ifbp+1).ge.0)then
            
          do jj=1,limit-2

          ib=latinx(jj)

          do kk=jj+1,limit-1

          ic=latinx(kk)

          do ll=kk+1,limit

          id=latinx(ll)
                  
          jfbp=max(ltype(ib),ltype(ic),ltype(id))
          lfbp=min(ltype(ib),ltype(ic),ltype(id))
          kfbp=ltype(ib)+ltype(ic)+ltype(id)-jfbp-lfbp
          jklbd=ifbp+lfbp+(kfbp*(kfbp-1))/2+(jfbp*(jfbp**2-1))/6
          kkfbp=lstfbp(jklbd)
          if(kkfbp.gt.0)then
                    
          sxab = xxx(ib)-xxx(ia)
          sxab = sxab-nint(sxab)
          syab = yyy(ib)-yyy(ia)
          syab = syab-nint(syab)
          szab = zzz(ib)-zzz(ia)
          szab = szab-nint(szab)

          xab=cell(1)*sxab+cell(4)*syab+cell(7)*szab
          if(abs(xab).lt.rcutfb)then

          yab=cell(2)*sxab+cell(5)*syab+cell(8)*szab
          if(abs(yab).lt.rcutfb)then

          zab=cell(3)*sxab+cell(6)*syab+cell(9)*szab
          if(abs(zab).lt.rcutfb)then
                          
          rab2=xab*xab+yab*yab+zab*zab

          sxac = xxx(ic)-xxx(ia)
          sxac = sxac-nint(sxac)
          syac = yyy(ic)-yyy(ia)
          syac = syac-nint(syac)
          szac = zzz(ic)-zzz(ia)
          szac = szac-nint(szac)

          xac=cell(1)*sxac+cell(4)*syac+cell(7)*szac
          if(abs(xac).lt.rcutfb)then
                            
          yac=cell(2)*sxac+cell(5)*syac+cell(8)*szac
          if(abs(yac).lt.rcutfb)then

          zac=cell(3)*sxac+cell(6)*syac+cell(9)*szac
          if(abs(zac).lt.rcutfb)then

          rac2=xac*xac+yac*yac+zac*zac

          sxad = xxx(id)-xxx(ia)
          sxad = sxad-nint(sxad)
          syad = yyy(id)-yyy(ia)
          syad = syad-nint(syad)
          szad = zzz(id)-zzz(ia)
          szad = szad-nint(szad)

          xad=cell(1)*sxad+cell(4)*syad+cell(7)*szad
          if(abs(xad).lt.rcutfb)then
                            
          yad=cell(2)*sxad+cell(5)*syad+cell(8)*szad
          if(abs(yad).lt.rcutfb)then

          zad=cell(3)*sxad+cell(6)*syad+cell(9)*szad
          if(abs(zad).lt.rcutfb)then

          rad2=xad*xad+yad*yad+zad*zad

          if(rcut4b(kkfbp)**2.ge.max(rab2,rac2,rad2))then

          rrab=1.d0/sqrt(rab2)
          rrac=1.d0/sqrt(rac2)
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
          if(abs(cosb).gt.1.d0)cosb=sign(1.d0,cosb)
          if(abs(cosc).gt.1.d0)cosc=sign(1.d0,cosc)
          if(abs(cosd).gt.1.d0)cosd=sign(1.d0,cosd)
          
c     
c     select potential energy function type
          
          ktyp=ltpfbp(kkfbp)

c     
c     calculate potential energy and scalar force term
          
          if(ktyp.eq.1)then
            
            thb=acos(cosb)
            thc=acos(cosc)
            thd=acos(cosd)
            engfbp=engfbp+
     x        (vharm(thb,kkfbp)+vharm(thc,kkfbp)+vharm(thd,kkfbp))/3.d0
            gamb=0.d0
            if(abs(thb).gt.1.d-12)gamb=dharm(thb,kkfbp)/(3.d0*sin(thb))
            gamc=0.d0
            if(abs(thc).gt.1.d-12)gamc=dharm(thc,kkfbp)/(3.d0*sin(thc))
            gamd=0.d0
            if(abs(thd).gt.1.d-12)gamd=dharm(thd,kkfbp)/(3.d0*sin(thd))
            
          else if(ktyp.eq.2)then
            
            engfbp=engfbp+
     x        (vcos(cosb,kkfbp)+vcos(cosc,kkfbp)+vcos(cosd,kkfbp))
     x        /3.d0
            gamb=-dcoz(cosb,kkfbp)/3.d0
            gamc=-dcoz(cosc,kkfbp)/3.d0
            gamd=-dcoz(cosd,kkfbp)/3.d0

          else if(ktyp.eq.3)then
            
            engfbp=engfbp+
     x        (vplan(cosb,kkfbp)+vplan(cosc,kkfbp)+vplan(cosd,kkfbp))
     x        /3.d0
            gamb=-dplan(kkfbp)/3.d0
            gamc=-dplan(kkfbp)/3.d0
            gamd=-dplan(kkfbp)/3.d0
            
          else
c     
c     undefined potential
            
            safe = .false.
            gamb=0.d0
            gamc=0.d0
            gamd=0.d0
            
          endif
          
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
     x      +(ruc*ucn*rrab*(xac-ruc*ucx-(rbc-ruc*rubc)*xab*rrab**2)
     x      - rvc*vcn*rrab*(xac-rvc*vcx-(rbc-rvc*rvbc)*xab*rrab**2))
     x      * gamc*rrac/wwc
     x      +(rud*udn*rrab*(xad-rud*udx-(rdb-rud*rubd)*xab*rrab**2)
     x      + rvd*vdn*rrab*(xad-rvd*vdx-(rdb-rvd*rvbd)*xab*rrab**2))
     x      * gamd*rrad/wwd
          
          fby = gamb*(-cosb*yab*rrab**2+rrab*(rub*uby+rvb*vby)/wwb)
     x      +(ruc*ucn*rrab*(yac-ruc*ucy-(rbc-ruc*rubc)*yab*rrab**2)
     x      - rvc*vcn*rrab*(yac-rvc*vcy-(rbc-rvc*rvbc)*yab*rrab**2))
     x      * gamc*rrac/wwc
     x      +(rud*udn*rrab*(yad-rud*udy-(rdb-rud*rubd)*yab*rrab**2)
     x      + rvd*vdn*rrab*(yad-rvd*vdy-(rdb-rvd*rvbd)*yab*rrab**2))
     x      * gamd*rrad/wwd
          
          fbz = gamb*(-cosb*zab*rrab**2+rrab*(rub*ubz+rvb*vbz)/wwb)
     x      +(ruc*ucn*rrab*(zac-ruc*ucz-(rbc-ruc*rubc)*zab*rrab**2)
     x      - rvc*vcn*rrab*(zac-rvc*vcz-(rbc-rvc*rvbc)*zab*rrab**2))
     x      * gamc*rrac/wwc
     x      +(rud*udn*rrab*(zad-rud*udz-(rdb-rud*rubd)*zab*rrab**2)
     x      + rvd*vdn*rrab*(zad-rvd*vdz-(rdb-rvd*rvbd)*zab*rrab**2))
     x      * gamd*rrad/wwd
          
          fcx = gamc*(-cosc*xac*rrac**2+rrac*(ruc*ucx+rvc*vcx)/wwc)
     x      +(rud*udn*rrac*(xad-rud*udx-(rcd-rud*rucd)*xac*rrac**2)
     x      - rvd*vdn*rrac*(xad-rvd*vdx-(rcd-rvd*rvcd)*xac*rrac**2))
     x      * gamd*rrad/wwd
     x      +(rub*ubn*rrac*(xab-rub*ubx-(rbc-rub*rucb)*xac*rrac**2)
     x      + rvb*vbn*rrac*(xab-rvb*vbx-(rbc-rvb*rvcb)*xac*rrac**2))
     x      * gamb*rrab/wwb
          
          fcy = gamc*(-cosc*yac*rrac**2+rrac*(ruc*ucy+rvc*vcy)/wwc)
     x      +(rud*udn*rrac*(yad-rud*udy-(rcd-rud*rucd)*yac*rrac**2)
     x      - rvd*vdn*rrac*(yad-rvd*vdy-(rcd-rvd*rvcd)*yac*rrac**2))
     x      * gamd*rrad/wwd
     x      +(rub*ubn*rrac*(yab-rub*uby-(rbc-rub*rucb)*yac*rrac**2)
     x      + rvb*vbn*rrac*(yab-rvb*vby-(rbc-rvb*rvcb)*yac*rrac**2))
     x      * gamb*rrab/wwb
          
          fcz = gamc*(-cosc*zac*rrac**2+rrac*(ruc*ucz+rvc*vcz)/wwc)
     x      +(rud*udn*rrac*(zad-rud*udz-(rcd-rud*rucd)*zac*rrac**2)
     x      - rvd*vdn*rrac*(zad-rvd*vdz-(rcd-rvd*rvcd)*zac*rrac**2))
     x      * gamd*rrad/wwd
     x      +(rub*ubn*rrac*(zab-rub*ubz-(rbc-rub*rucb)*zac*rrac**2)
     x      + rvb*vbn*rrac*(zab-rvb*vbz-(rbc-rvb*rvcb)*zac*rrac**2))
     x      * gamb*rrab/wwb
          
          fdx = gamd*(-cosd*xad*rrad**2+rrad*(rud*udx+rvd*vdx)/wwd)
     x      +(rub*ubn*rrad*(xab-rub*ubx-(rdb-rub*rudb)*xad*rrad**2)
     x      - rvb*vbn*rrad*(xab-rvb*vbx-(rdb-rvb*rvdb)*xad*rrad**2))
     x      * gamb*rrab/wwb
     x      +(ruc*ucn*rrad*(xac-ruc*ucx-(rcd-ruc*rudc)*xad*rrad**2)
     x      + rvc*vcn*rrad*(xac-rvc*vcx-(rcd-rvc*rvdc)*xad*rrad**2))
     x      * gamc*rrac/wwc
          
          fdy = gamd*(-cosd*yad*rrad**2+rrad*(rud*udy+rvd*vdy)/wwd)
     x      +(rub*ubn*rrad*(yab-rub*uby-(rdb-rub*rudb)*yad*rrad**2)
     x      - rvb*vbn*rrad*(yab-rvb*vby-(rdb-rvb*rvdb)*yad*rrad**2))
     x      * gamb*rrab/wwb
     x      +(ruc*ucn*rrad*(yac-ruc*ucy-(rcd-ruc*rudc)*yad*rrad**2)
     x      + rvc*vcn*rrad*(yac-rvc*vcy-(rcd-rvc*rvdc)*yad*rrad**2))
     x      * gamc*rrac/wwc
          
          fdz = gamd*(-cosd*zad*rrad**2+rrad*(rud*udz+rvd*vdz)/wwd)
     x      +(rub*ubn*rrad*(zab-rub*ubz-(rdb-rub*rudb)*zad*rrad**2)
     x      - rvb*vbn*rrad*(zab-rvb*vbz-(rdb-rvb*rvdb)*zad*rrad**2))
     x      * gamb*rrab/wwb
     x      +(ruc*ucn*rrad*(zac-ruc*ucz-(rcd-ruc*rudc)*zad*rrad**2)
     x      + rvc*vcn*rrad*(zac-rvc*vcz-(rcd-rvc*rvdc)*zad*rrad**2))
     x      * gamc*rrac/wwc
          
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

          endif
          endif
          endif
          endif
          endif
          endif
          endif
          endif
          endif
          endif
          endif

          enddo
          enddo
          enddo
          
          endif

        enddo

      enddo
c
c     check for undefined potentials

      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,443)

c     
c     global sum of four body potential and virial

      buffer(1)=engfbp
      buffer(2)=virfbp
      call gdsum(buffer(1),2,buffer(3))
      engfbp=buffer(1)
      virfbp=buffer(2)

c     
c     restore coordinate array to original representation

      do i=1,natms

        sxx=xxx(i)
        syy=yyy(i)
        szz=zzz(i)

        xxx(i)=cell(1)*sxx+cell(4)*syy+cell(7)*szz
        yyy(i)=cell(2)*sxx+cell(5)*syy+cell(8)*szz
        zzz(i)=cell(3)*sxx+cell(6)*syy+cell(9)*szz

      enddo
c
c     restore cell vector
      if(imcon.eq.0) then
        cell(1) = 0.d0
        cell(5) = 0.d0
        cell(9) = 0.d0
      endif
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

#ifdef VAMPIR
      call VTEND(20, ierr)
#endif
      return

      end
