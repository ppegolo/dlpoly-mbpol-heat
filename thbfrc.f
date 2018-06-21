      subroutine thbfrc
     x  (idnode,mxnode,natms,imcon,rcuttb,engtbp,virtbp,
     x  latinx,ltype,lst,lct,link,lsttbp,ltptbp,prmtbp,cell,
     x  xxx,yyy,zzz,fxx,fyy,fzz,rcut3b,stress,buffer)
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating three body forces arising 
c     from the included angle between three atoms
c
c     Note: the subroutine converts coordinates to reduced units
c     to avoid a call to images.f. The link cell algorithm used
c     here necessitates a parallelepiped cell geometry
c     
c     copyright - daresbury laboratory 1994
c     author   - w.smith march 1994 
c     
c     wl
c     2002/05/31 14:46:41
c     1.11
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
      dimension lsttbp(mxtbp),ltptbp(mxtbp),prmtbp(mxtbp,mxptbp)
      dimension cell(9),rcut3b(mxtbp),stress(9),buffer(mxbuff)
      dimension rcell(9),cprp(10)
      
      data nix/ 0,-1,-1,-1, 0, 0,-1, 1,-1, 0, 1,-1, 0, 1,
     x  1, 1, 1, 0, 0, 1,-1, 1, 0,-1, 1, 0,-1/
      data niy/ 0, 0,-1, 1, 1, 0, 0, 0,-1,-1,-1, 1, 1, 1,
     x  0, 1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1/
      data niz/ 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     x  0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1/
      
c     
c     truncation and screening functions

      sw1(x,y,a)=exp(-(x**8+y**8)/a**8)
      sw2(x,y,a,b)=exp(-(x/a+y/b))

c     
c     truncated harmonic valence angle potential - use with sw1
      
      vt1(t,a,b)=0.5d0*a*(t-b)**2
      gt1(t,a,b)=a*(t-b)
c     
c     screened harmonic valence angle potential - use with sw2
      
      vt2(t,a,b)=0.5d0*a*(t-b)**2
      gt2(t,a,b)=a*(t-b)
c     
c     screened vessal potential type 1 - use with sw2
      
      vt3(t,a,b)=a/(8.d0*(b-pi)**2)*(((b-pi)**2-(t-pi)**2)**2)
      gt3(t,a,b)=a/(2.d0*(b-pi)**2)*((b-pi)**2-(t-pi)**2)*(t-pi)
c     
c     truncated vessal potential type 2 - use with sw1
      
      vt4(t,a,b,c)=a*(t**c*(t-b)**2*(t+b-2.d0*pi)**2-
     x  0.5d0*c*pi**(c-1.d0)*(t-b)**2*(pi-b)**3)
      gt4(t,a,b,c)=a*(t**(c-1.d0)*(t-b)*(t+b-2.d0*pi)*
     x  ((c+4.d0)*t**2-2.d0*pi*(c+2.d0)*t+c*b*(2.d0*pi-b))-
     x  c*pi**(c-1.d0)*(t-b)*(pi-b)**3)
c
c     dreiding/charmm hydrogen bond

      shb(r,a)=(5.d0*(a/r)**2-6.d0)*(a/r)**10
      vt5(a,b)=a*b**4
      gt5(a,b)=4.d0*a*b**3

#ifdef VAMPIR
      call VTBEGIN(19, ierr)
#endif
c
c     flag for undefined potentials

      safe = .true.
c     
c     initialise potential energy and virial

      engtbp=0.d0
      virtbp=0.d0
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

        cell(1) = 2.d0*xm+rcuttb
        cell(2) = 0.d0
        cell(3) = 0.d0
        cell(4) = 0.d0
        cell(5) = 2.d0*ym+rcuttb
        cell(6) = 0.d0
        cell(7) = 0.d0
        cell(8) = 0.d0
        cell(9) = 2.d0*zm+rcuttb
      
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
      
      if(imcon.gt.3)call error(idnode,67)
      call invert(cell,rcell,det)
      call dcell(cell,cprp)
c     
c     calculate link cell numbers
      
      nbx=int(cprp(7)/(rcuttb+1.d-6))
      nby=int(cprp(8)/(rcuttb+1.d-6))
      nbz=int(cprp(9)/(rcuttb+1.d-6))
      ncells=nbx*nby*nbz
      if(nbx.lt.3.or.nby.lt.3.or.nbz.lt.3)then

        call error(idnode,68)

      endif
      if(ncells.gt.mxcell) then
        
        if(idnode.eq.0) then

          write(nrite,'(a,i6)')
     x    'number of required link cells in routine thbfrc is ',ncells
          write(nrite,'(a,i6)')
     x    'number of default link cells in routine thbfrc is ',mxcell

        endif

        call error(idnode,69)

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
c     loop over central atoms of angles
      
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
        
c     construct mini-list of neighbour cell contents

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
          
          i=latinx(ii)
          itbp=mx2tbp*(ltype(i)-1)
          if(mod(i,mxnode).eq.idnode.and.lsttbp(itbp+1).ge.0)then

          last=limit

          do kk=1,limit/2
              
          if(kk.gt.(limit-1)/2)last=limit/2

          do jj=1,last

          j=latinx(jj)
          jk=jj+kk
          if(jk.gt.limit)jk=jk-limit
          k=latinx(jk)
          if(i.ne.j.and.i.ne.k)then
                  
          jtbp=max(ltype(j),ltype(k))
          ktbp=min(ltype(j),ltype(k))
          jktbp=itbp+(jtbp*(jtbp-1))/2+ktbp
          kktbp=lsttbp(jktbp)
          if(kktbp.gt.0)then
c     
c     make labels etc consistent with angfrc.f

          ia = j
          ib = i
          ic = k
                    
          sxab = xxx(ia)-xxx(ib)
          sxab = sxab-nint(sxab)
          syab = yyy(ia)-yyy(ib)
          syab = syab-nint(syab)
          szab = zzz(ia)-zzz(ib)
          szab = szab-nint(szab)

          xab=cell(1)*sxab+cell(4)*syab+cell(7)*szab
          if(abs(xab).lt.rcuttb)then

          yab=cell(2)*sxab+cell(5)*syab+cell(8)*szab
          if(abs(yab).lt.rcuttb)then

          zab=cell(3)*sxab+cell(6)*syab+cell(9)*szab
          if(abs(zab).lt.rcuttb)then
                          
          sxbc = xxx(ic)-xxx(ib)
          sxbc = sxbc-nint(sxbc)
          sybc = yyy(ic)-yyy(ib)
          sybc = sybc-nint(sybc)
          szbc = zzz(ic)-zzz(ib)
          szbc = szbc-nint(szbc)

          xbc=cell(1)*sxbc+cell(4)*sybc+cell(7)*szbc
          if(abs(xbc).lt.rcuttb)then
                            
          ybc=cell(2)*sxbc+cell(5)*sybc+cell(8)*szbc
          if(abs(ybc).lt.rcuttb)then

          zbc=cell(3)*sxbc+cell(6)*sybc+cell(9)*szbc
          if(abs(zbc).lt.rcuttb)then

          ktyp=ltptbp(kktbp)
          rab=sqrt(xab*xab+yab*yab+zab*zab)
          rbc=sqrt(xbc*xbc+ybc*ybc+zbc*zbc)

          if(rcut3b(kktbp).ge.max(rab,rbc))then

          xac = xab - xbc
          yac = yab - ybc
          zac = zab - zbc
          rac=sqrt(xac*xac+yac*yac+zac*zac)

          rrab = 1.d0/rab
          rrbc = 1.d0/rbc
          rrac = 1.d0/rac
c     
c     normalise direction vectors

          xab = xab*rrab
          yab = yab*rrab
          zab = zab*rrab

          xbc = xbc*rrbc
          ybc = ybc*rrbc
          zbc = zbc*rrbc

          xac = xac*rrac
          yac = yac*rrac
          zac = zac*rrac

          cost=(xab*xbc+yab*ybc+zab*zbc)
          if(abs(cost).gt.1.d0)cost=sign(1.d0,cost)
          if(ktyp.ne.5)then

            sint=max(1.d-8,sqrt(1.d0-cost*cost))
            theta=acos(cost)

          endif
                                  
          if(ktyp.eq.0)then

          pterm=vt1(theta,prmtbp(kktbp,1),prmtbp(kktbp,2))
          gamma=gt1(theta,prmtbp(kktbp,1),prmtbp(kktbp,2))/sint
          vterm=0.d0
          gamsa=0.d0
          gamsc=0.d0
          gamsb=0.d0

          elseif(ktyp.eq.1)then

          scrn=sw1(rab,rbc,prmtbp(kktbp,3))
          pterm=scrn*vt1(theta,prmtbp(kktbp,1),prmtbp(kktbp,2))
          vterm=-8.d0*pterm*(rab**8+rbc**8)/prmtbp(kktbp,3)**8
          gamma=scrn*gt1(theta,prmtbp(kktbp,1),prmtbp(kktbp,2))/sint
          gamsa=(8.d0*pterm/prmtbp(kktbp,3)**8)*rab**7
          gamsc=(8.d0*pterm/prmtbp(kktbp,3)**8)*rbc**7
          gamsb=0.d0

          elseif(ktyp.eq.2)then

          scrn=sw2(rab,rbc,prmtbp(kktbp,3),prmtbp(kktbp,4))
          pterm=scrn*vt2(theta,prmtbp(kktbp,1),prmtbp(kktbp,2))
          vterm=-pterm*(rab/prmtbp(kktbp,3)+rbc/prmtbp(kktbp,4))
          gamma=scrn*gt2(theta,prmtbp(kktbp,1),prmtbp(kktbp,2))/sint
          gamsa=(pterm/prmtbp(kktbp,3))
          gamsc=(pterm/prmtbp(kktbp,4))
          gamsb=0.d0

          elseif(ktyp.eq.3)then

          scrn=sw2(rab,rbc,prmtbp(kktbp,3),prmtbp(kktbp,4))
          pterm=scrn*vt3(theta,prmtbp(kktbp,1),prmtbp(kktbp,2))
          vterm=-pterm*(rab/prmtbp(kktbp,3)+rbc/prmtbp(kktbp,4))
          gamma=scrn*gt3(theta,prmtbp(kktbp,1),prmtbp(kktbp,2))/sint
          gamsa=(pterm/prmtbp(kktbp,3))
          gamsc=(pterm/prmtbp(kktbp,4))
          gamsb=0.d0

          elseif(ktyp.eq.4)then
                                    
          scrn=sw1(rab,rbc,prmtbp(kktbp,4))
          pterm=scrn*vt4(theta,prmtbp(kktbp,1),prmtbp(kktbp,2),
     x      prmtbp(kktbp,3))
          vterm=-8.d0*pterm*(rab**8+rbc**8)/prmtbp(kktbp,4)**8
          gamma=scrn*gt4(theta,prmtbp(kktbp,1),prmtbp(kktbp,2),
     x      prmtbp(kktbp,3))/sint
          gamsa=(8.d0*pterm/prmtbp(kktbp,4)**8)*rab**7
          gamsc=(8.d0*pterm/prmtbp(kktbp,4)**8)*rbc**7
          gamsb=0.d0

          elseif(ktyp.eq.5)then

          if(min(rab,rbc).lt.1.5d0.and.rac.le.rcut3b(kktbp))then

          scrn=shb(rac,prmtbp(kktbp,2))
          tterm=vt5(prmtbp(kktbp,1),cost)
          pterm=tterm*scrn
          uterm=60.d0*((prmtbp(kktbp,2)/rac)**2-1.d0)*
     x      (prmtbp(kktbp,2)/rac)**10
          vterm=-tterm*uterm
          gamma=-scrn*gt5(prmtbp(kktbp,1),cost)
          gamsb=tterm*uterm/rac
          gamsa=0.d0
          gamsc=0.d0

          else

          pterm=0.d0
          vterm=0.d0
          gamma=0.d0
          gamsa=0.d0
          gamsb=0.d0
          gamsc=0.d0

          endif

          else

          safe = .false.
          pterm=0.d0
          vterm=0.d0
          gamma=0.d0
          gamsa=0.d0
          gamsb=0.d0
          gamsc=0.d0
                     
          endif

          engtbp=engtbp+pterm
          virtbp=virtbp+vterm

c     
c     calculate atomic forces

          fxa = gamma*(xbc-xab*cost)*rrab+gamsa*xab+gamsb*xac
          fya = gamma*(ybc-yab*cost)*rrab+gamsa*yab+gamsb*yac
          fza = gamma*(zbc-zab*cost)*rrab+gamsa*zab+gamsb*zac
                                
          fxc = gamma*(xab-xbc*cost)*rrbc+gamsc*xbc-gamsb*xac
          fyc = gamma*(yab-ybc*cost)*rrbc+gamsc*ybc-gamsb*yac
          fzc = gamma*(zab-zbc*cost)*rrbc+gamsc*zbc-gamsb*zac
                                  
          fxx(ia)=fxx(ia)+fxa
          fyy(ia)=fyy(ia)+fya
          fzz(ia)=fzz(ia)+fza
                                  
          fxx(ib)=fxx(ib)-fxa-fxc
          fyy(ib)=fyy(ib)-fya-fyc
          fzz(ib)=fzz(ib)-fza-fzc
                                  
          fxx(ic)=fxx(ic)+fxc
          fyy(ic)=fyy(ic)+fyc
          fzz(ic)=fzz(ic)+fzc
#ifdef STRESS
c     
c     calculate stress tensor
              
          strs1 = strs1 + rab*xab*fxa + rbc*xbc*fxc
          strs2 = strs2 + rab*xab*fya + rbc*xbc*fyc
          strs3 = strs3 + rab*xab*fza + rbc*xbc*fzc
          strs5 = strs5 + rab*yab*fya + rbc*ybc*fyc
          strs6 = strs6 + rab*yab*fza + rbc*ybc*fzc
          strs9 = strs9 + rab*zab*fza + rbc*zbc*fzc
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

          enddo
          enddo
          
          endif

        enddo

      enddo
c
c     check for undefined potentials

      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,442)

c     
c     global sum of three body potential and virial

      if(mxnode.gt.1)then

        buffer(1)=engtbp
        buffer(2)=virtbp
        call gdsum(buffer(1),2,buffer(3))
        engtbp=buffer(1)
        virtbp=buffer(2)

      endif

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
      call VTEND(19, ierr)
#endif
      return

      end
