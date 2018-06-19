      subroutine dihfrc
     x  (idnode,imcon,mxnode,ntdihd,keyfce,dlrpot,epsq,engcpe,
     x  engdih,engsrp,rcut,rvdw,vircpe,virdih,virsrp,keydih,listdih,
     x  ltype,lstvdw,buffer,cell,chge,fxx,fyy,fzz,prmdih,xxx,yyy,
     x  zzz,xdab,ydab,zdab,xdbc,ydbc,zdbc,xdcd,ydcd,zdcd,vvv,ggg,
     x  stress)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating dihedral energy and force 
c     terms in molecular dynamics.
c     
c     version 3: scale factors for reduces electrostatic and vdw
c     1-4 interactions.
c     
c     NOTE: assumes 1-4 interactions are in the exclude list
c     
c     block as opposed to stride version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith       march 1992
c     
c     version 3  december 1993
c     author    - t. forester.
c     
c     stress tensor added t.forester june 1995
c
c     ryckaert-bellemans potential added a.smondyrev may 2000
c
c     fluorinated ryckaert-bellemans potential added a.smondyrev may2000
c     
c     wl
c     2003/05/08 08:45:10
c     1.10
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical safe,csafe
      dimension keydih(mxtdih),listdih(mxdihd,5),prmdih(mxtdih,mxpdih)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms),cell(9)
      dimension xdab(msbad),ydab(msbad),zdab(msbad)
      dimension xdbc(msbad),ydbc(msbad),zdbc(msbad)
      dimension xdcd(msbad),ydcd(msbad),zdcd(msbad)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension chge(mxatms)
      dimension buffer(mxbuff)
      dimension vvv(mxgrid,mxvdw),ggg(mxgrid,mxvdw)
      dimension ltype(mxatms),lstvdw(mxvdw)
      dimension stress(9)
      
c     
c     define angular potential function and derivative
c     using the parameters in array prmdih

c
c     key=1 for torsion dihedral potential

      vdihd(phi,k)=prmdih(k,1)*(1.d0+cos(prmdih(k,3)*phi-prmdih(k,2)))
      ddihd(phi,k)=-prmdih(k,1)*prmdih(k,3)*sin(prmdih(k,3)*phi-
     x  prmdih(k,2))
c
c     key=2 for harmonic improper dihedral

      vharm(phi,k)=0.5d0*prmdih(k,1)*(phi*phi)
      dharm(phi,k)=prmdih(k,1)*(phi)
c
c     key=3 for harmonic cosine dihedral (note sint cancelled in dhcos)

      vhcos(phi,k)=0.5d0*prmdih(k,1)*(cos(phi)-cos(prmdih(k,2)))**2
      dhcos(phi,k)=-prmdih(k,1)*(cos(phi)-cos(prmdih(k,2)))
c
c     key=4 for 3-term cosine dihedral

      vcos3(phi,k)=0.5d0*(prmdih(k,1)*(1.d0+cos(phi))+prmdih(k,2)*
     x  (1.d0-cos(2.d0*phi))+prmdih(k,3)*(1.d0+cos(3.d0*phi)))
      dcos3(phi,k)=-0.5d0*(prmdih(k,1)*sin(phi)-2.d0*prmdih(k,2)*
     x  sin(2.d0*phi)+3.d0*prmdih(k,3)*sin(3.d0*phi))

c     key=5 for ryckaert-bellemans potential
c           chem.phys.lett., vol.30, p.123, 1975.
c           ATTENTION !!! Modified to have trans configuration
c           correspond to theta=180 rather than
c           theta=0 as in original form.

      vryck(phi,k)=prmdih(k,1)*(9.28d0-12.16d0*cos(phi)-
     x  13.12d0*(cos(phi))**2+3.06d0*(cos(phi))**3+
     x  26.24d0*(cos(phi))**4+31.5d0*(cos(phi))**5)
      dryck(phi,k)=prmdih(k,1)*(12.16d0+
     x  26.24d0*(cos(phi))-9.18d0*(cos(phi))**2-
     x  104.96d0*(cos(phi))**3-157.5d0*(cos(phi))**4)


c     dryck(phi,k)=prmdih(1,k)*(1.462d0+3.156d0*cos(phi)-
c    x  1.104d0*(cos(phi))**2-12.624d0*(cos(phi))**3-
c    x  18.94d0*(cos(phi))**4)

c     key=6 for fluorinated ryckaert-bellemans potential
c           Rice at al., JCP 104, 2101, (1996).

      vrbf(phi,k)=prmdih(k,1)*(3.55d0-2.78d0*cos(phi)-
     x  3.56d0*(cos(phi))**2-1.64d0*(cos(phi))**3+
     x  7.13d0*(cos(phi))**4+12.84d0*(cos(phi))**5+
     x  9.67d0*exp(-56.d0*(phi-pi)**2))
      drbf(phi,k)=prmdih(k,1)*(2.78d0+7.12d0*cos(phi)+
     x  4.92d0*(cos(phi))**2-28.52d0*(cos(phi))**3-
     x  64.2d0*(cos(phi))**4)
      erbf(phi,k)=-1083.04d0*(phi-pi)*exp(-56.0*(phi-pi)**2)

c     key=7 for opls cosine dihedral

      vopls(phi,k)=prmdih(k,1)+
     x  0.5d0*(prmdih(k,2)*(1.d0+cos(phi))+prmdih(k,3)*
     x  (1.d0-cos(2.d0*phi))+prmdih(k,6)*(1.d0+cos(3.d0*phi)))
      dopls(phi,k)=-0.5d0*(prmdih(k,2)*sin(phi)-2.d0*prmdih(k,3)*
     x  sin(2.d0*phi)+3.d0*prmdih(k,6)*sin(3.d0*phi))

#ifdef VAMPIR
      call VTBEGIN(23, ierr)
#endif
c
c     check size of work arrays

      if((ntdihd-mxnode+1)/mxnode.gt.msbad) call error(idnode,421)
c
c     block indices

      idih1 = (idnode*ntdihd)/mxnode+1
      idih2 = ((idnode+1)*ntdihd)/mxnode

      twopi = 2.d0*pi
      rtwopi=1.d0/twopi
      safe = .true.
      csafe = .true.

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
      do i=idih1,idih2

        ii=ii+1
c
c     indices of atoms involved

        ia=listdih(ii,2)
        ib=listdih(ii,3)
        ic=listdih(ii,4)
        id=listdih(ii,5)

c
c     define components of bond vectors

        xdab(ii)=xxx(ia)-xxx(ib)
        ydab(ii)=yyy(ia)-yyy(ib)
        zdab(ii)=zzz(ia)-zzz(ib)

        xdbc(ii)=xxx(ib)-xxx(ic)
        ydbc(ii)=yyy(ib)-yyy(ic)
        zdbc(ii)=zzz(ib)-zzz(ic)

        xdcd(ii)=xxx(ic)-xxx(id)
        ydcd(ii)=yyy(ic)-yyy(id)
        zdcd(ii)=zzz(ic)-zzz(id)

      enddo

c
c     periodic boundary condition

      call images(imcon,0,1,ii,cell,xdab,ydab,zdab)
      call images(imcon,0,1,ii,cell,xdbc,ydbc,zdbc)
      call images(imcon,0,1,ii,cell,xdcd,ydcd,zdcd)

c
c     zero dihedral energy accumulator

      engdih=0.d0
      virdih=0.d0

c
c     zero scaled 1-4 electrostatic and short range potential accumulators

      engc14 = 0.d0
      virc14 = 0.d0
      engs14 = 0.d0
      virs14 = 0.d0
c
c     loop over all specified dihedrals

      ii=0
      do i=idih1,idih2

c
c     define components of bond vectors

        ii=ii+1

        xab=xdab(ii)
        yab=ydab(ii)
        zab=zdab(ii)

        xbc=xdbc(ii)
        ybc=ydbc(ii)
        zbc=zdbc(ii)
        rrbc=1.d0/sqrt(xbc*xbc+ybc*ybc+zbc*zbc)

        xcd=xdcd(ii)
        ycd=ydcd(ii)
        zcd=zdcd(ii)

        xac=xab+xbc
        yac=yab+ybc
        zac=zab+zbc
c
c     construct first dihedral vector

        pbx=yab*zbc-zab*ybc
        pby=zab*xbc-xab*zbc
        pbz=xab*ybc-yab*xbc
        pb2=pbx*pbx+pby*pby+pbz*pbz
        rpb1=1.d0/sqrt(pb2)
        rpb2 =rpb1*rpb1
c
c     construct second dihedral vector

        pcx=ybc*zcd-zbc*ycd
        pcy=zbc*xcd-xbc*zcd
        pcz=xbc*ycd-ybc*xcd
        pc2=pcx*pcx+pcy*pcy+pcz*pcz
        rpc1=1.d0/sqrt(pc2)
        rpc2 = rpc1*rpc1
c
c     determine dihedral angle

        pbpc=pbx*pcx+pby*pcy+pbz*pcz
        cost=pbpc*rpb1*rpc1
        sint=(xbc*(pcy*pbz-pcz*pby)+ybc*(pbx*pcz-pbz*pcx)+
     x    zbc*(pcx*pby-pcy*pbx))*(rpb1*rpc1*rrbc)

        theta=atan2(sint,cost)
c
c     avoid singularity in sint

        sint = sign(max(1.d-8,abs(sint)),sint)
        rsint = 1.d0/sint
c
c     selection of potential energy function type

        kk=listdih(ii,1)

c
c     calculate potential energy and scalar force term

        if(keydih(kk).eq.1)then

          engdih=engdih+vdihd(theta,kk)
          gamma=ddihd(theta,kk)*rpb1*rpc1*rsint

        else if(keydih(kk).eq.2)then
c
c     find deviation from equilibrium angle

          theta = theta - prmdih(kk,2)
          theta = theta - nint(theta*rtwopi)*twopi

          engdih=engdih+vharm(theta,kk)
          gamma=dharm(theta,kk)*rpb1*rpc1*rsint

        else if(keydih(kk).eq.3)then

          engdih=engdih+vhcos(theta,kk)
          gamma=dhcos(theta,kk)*rpb1*rpc1

        else if(keydih(kk).eq.4)then

          engdih=engdih+vcos3(theta,kk)
          gamma=dcos3(theta,kk)*rpb1*rpc1*rsint

        else if(keydih(kk).eq.5)then

          engdih=engdih+vryck(theta,kk)
          gamma=dryck(theta,kk)*rpb1*rpc1

        else if(keydih(kk).eq.6)then

          engdih=engdih+vrbf(theta,kk)
          gamma=(drbf(theta,kk)+erbf(theta,kk)*rsint)*rpb1*rpc1

        else if(keydih(kk).eq.7)then

          theta = theta - prmdih(kk,7)
          engdih=engdih+vopls(theta,kk)
          gamma=dopls(theta,kk)*rpb1*rpc1*rsint

        else
c
c     undefined potential

          safe = .false.
          gamma = 0.d0

        endif

c
c     indices of atoms involved

        ia=listdih(ii,2)
        ib=listdih(ii,3)
        ic=listdih(ii,4)
        id=listdih(ii,5)

c
c     calculate atomic forces

        fax = gamma*((-pcy*zbc+pcz*ybc)-pbpc*rpb2*(-pby*zbc+pbz*ybc))
        fay = gamma*(( pcx*zbc-pcz*xbc)-pbpc*rpb2*( pbx*zbc-pbz*xbc))
        faz = gamma*((-pcx*ybc+pcy*xbc)-pbpc*rpb2*(-pbx*ybc+pby*xbc))

        fcx = gamma*((-pcy*zab+pcz*yab)-pbpc*rpb2*(-pby*zab+pbz*yab))
        fcy = gamma*(( pcx*zab-pcz*xab)-pbpc*rpb2*( pbx*zab-pbz*xab))
        fcz = gamma*((-pcx*yab+pcy*xab)-pbpc*rpb2*(-pbx*yab+pby*xab))

        fb1x= gamma*((-pby*zcd+pbz*ycd)-pbpc*rpc2*(-pcy*zcd+pcz*ycd))
        fb1y= gamma*(( pbx*zcd-pbz*xcd)-pbpc*rpc2*( pcx*zcd-pcz*xcd))
        fb1z= gamma*((-pbx*ycd+pby*xcd)-pbpc*rpc2*(-pcx*ycd+pcy*xcd))

        fd1x= gamma*((-pby*zbc+pbz*ybc)-pbpc*rpc2*(-pcy*zbc+pcz*ybc))
        fd1y= gamma*(( pbx*zbc-pbz*xbc)-pbpc*rpc2*( pcx*zbc-pcz*xbc))
        fd1z= gamma*((-pbx*ybc+pby*xbc)-pbpc*rpc2*(-pcx*ybc+pcy*xbc))

        fxx(ia)=fxx(ia)+fax
        fyy(ia)=fyy(ia)+fay
        fzz(ia)=fzz(ia)+faz

        fxx(ib)=fxx(ib)-fax-fcx+fb1x
        fyy(ib)=fyy(ib)-fay-fcy+fb1y
        fzz(ib)=fzz(ib)-faz-fcz+fb1z

        fxx(ic)=fxx(ic)+fcx-fb1x-fd1x
        fyy(ic)=fyy(ic)+fcy-fb1y-fd1y
        fzz(ic)=fzz(ic)+fcz-fb1z-fd1z

        fxx(id)=fxx(id)+fd1x
        fyy(id)=fyy(id)+fd1y
        fzz(id)=fzz(id)+fd1z
#ifdef STRESS
c     stress tensor calculation for dihedral terms
        
        strs1 = strs1 + xab*fax + xbc*(fb1x-fcx) - xcd*fd1x 
        strs2 = strs2 + yab*fax + ybc*(fb1x-fcx) - ycd*fd1x 
        strs3 = strs3 + zab*fax + zbc*(fb1x-fcx) - zcd*fd1x 
        strs5 = strs5 + yab*fay + ybc*(fb1y-fcy) - ycd*fd1y 
        strs6 = strs6 + yab*faz + ybc*(fb1z-fcz) - ycd*fd1z 
        strs9 = strs9 + zab*faz + zbc*(fb1z-fcz) - zcd*fd1z 
#endif
c     
c     1-4 electrostatics : adjust by weighting factor
c     assumes 1-4 interactions are in the exclude list
        
        kk = listdih(ii,1)
        scale = prmdih(kk,4)
        
        xad = xac+xcd
        yad = yac+ycd
        zad = zac+zcd
        
        rad = sqrt(xad**2 + yad**2 + zad**2)
c     
c     scaled charge product*dielectric
        
        chgprd = scale*chge(ia)*chge(id)*r4pie0
        coul = 0.d0
        fcoul = 0.d0
        
c     
c     truncation of potential for all schemes except ewald sum
        
        if(chgprd.ne.0.d0.and.keyfce.gt.0) then
          
c     
c     Electrostatics by ewald sum
          
          if(keyfce/2.eq.1.or.keyfce/2.eq.6.or.keyfce/2.eq.7) then
            
            coul = chgprd/(epsq*rad)
            fcoul = coul/(rad**2)
c     
c     distance dependent dielectric
            
          elseif(rcut.gt.rad) then
            
            if(keyfce/2.eq.2) then
              
              coul = chgprd/(epsq*rad**2)
              fcoul = 2.0d0*coul/(rad**2)
c     
c     coulombic
              
            else if(keyfce/2.eq.3) then
              
              coul = chgprd/(epsq*rad)
              fcoul = coul/(rad**2)
              
c     
c     truncated and shifted coulombic
              
            else if(keyfce/2.eq.4) then
              
              coul = chgprd*(rcut-rad)/(epsq*rad*rcut)
              fcoul = chgprd/(epsq*rad**3)
c
c     reaction field (not zero shifted)
              
            else if(keyfce/2.eq.5) then
              
              b0 = 2.d0*(epsq - 1.d0)/(2.d0*epsq + 1.d0)
              rfld0 = b0/rcut**3
              rfld2 = rfld0*0.5d0
              coul  = chgprd*(1.d0/rad + rfld2*rad*rad)
              fcoul = chgprd*(1.d0/rad**3 - rfld0)
              
            elseif (keyfce/2.eq.0) then

              coul = 0.d0
              fcoul = 0.d0
              
            else
              
              call error(idnode,446)
              
            endif
            
          endif
c
c     correction to electrostatic energy and virial
          
          engc14 = engc14 + coul
          virc14 = virc14 - fcoul*rad**2
          
          fx = fcoul*xad
          fy = fcoul*yad
          fz = fcoul*zad
          
          fxx(ia) = fxx(ia) + fx
          fyy(ia) = fyy(ia) + fy
          fzz(ia) = fzz(ia) + fz
          
          fxx(id) = fxx(id) - fx
          fyy(id) = fyy(id) - fy
          fzz(id) = fzz(id) - fz
#ifdef STRESS
c     
c     calculate stress tensor
          
          strs1 = strs1 + xad*fx
          strs2 = strs2 + xad*fy
          strs3 = strs3 + xad*fz
          strs5 = strs5 + yad*fy
          strs6 = strs6 + yad*fz
          strs9 = strs9 + zad*fz
#endif
        endif
c     
c     1--4 short-range interactions scaled
        
        scale = prmdih(kk,5)
        gamma = 0.d0
        
        if(mod(keyfce,2).eq.1) then
c     
c     atomic and potential function indices
          
          ka=max(ltype(ia),ltype(id))
          kb=min(ltype(ia),ltype(id))
          k=lstvdw((ka*(ka-1))/2+kb)
          
          if(scale*vvv(1,k).ne.0.d0) then 
c     
c     apply truncation of potential

            if(rvdw.gt.rad)then
              
c     
c     determine interpolation panel for force arrays
              
              l=int(rad/dlrpot)
              ppp=rad/dlrpot-dble(l)
              
              
c     calculate interaction energy using 3-point interpolation
              
              vk = vvv(l,k)
              vk1 = vvv(l+1,k)
              vk2 = vvv(l+2,k)
              
              t1 = vk + (vk1-vk)*ppp
              t2 = vk1 +(vk2 - vk1)*(ppp - 1.0d0)
              
              engs14 = (t1 + (t2-t1)*ppp*0.5d0)*scale + engs14
c     
c     calculate forces using 3-point interpolation
              
              gk = ggg(l,k)
              gk1 = ggg(l+1,k)
              gk2 = ggg(l+2,k)
              
              t1 = gk + (gk1-gk)*ppp
              t2 = gk1 + (gk2-gk1)*(ppp - 1.0d0)
              
              gamma = scale*(t1 +(t2-t1)*ppp*0.5d0)/(rad**2)
              
              fx = gamma*xad
              fy = gamma*yad
              fz = gamma*zad
              
              fxx(ia)=fxx(ia)+fx
              fyy(ia)=fyy(ia)+fy
              fzz(ia)=fzz(ia)+fz
              
              fxx(id)=fxx(id)-fx
              fyy(id)=fyy(id)-fy
              fzz(id)=fzz(id)-fz
#ifdef STRESS
c     
c     calculate stress tensor

              strs1 = strs1 + xad*fx
              strs2 = strs2 + xad*fy
              strs3 = strs3 + xad*fz
              strs5 = strs5 + yad*fy
              strs6 = strs6 + yad*fz
              strs9 = strs9 + zad*fz
#endif
c     
c     calculate scaled 1-4 short-range potential virial
              
              virs14=virs14-gamma*rad**2
              
            else

              csafe=.false.

            endif
            
          endif

        endif
        
      enddo
      
c     
c     sum contributions to potentials
      
      if(mxnode.gt.1) then
        
        buffer(6) = engdih
        buffer(7) = engc14
        buffer(8) = virc14
        buffer(9) = engs14
        buffer(10) = virs14

        call gdsum(buffer(6),5,buffer(1))
        
        engdih = buffer(6)
        engc14 = buffer(7)
        virc14 = buffer(8)
        engs14 = buffer(9)
        virs14 = buffer(10)
        
      endif
      
      engcpe = engcpe + engc14
      vircpe = vircpe + virc14
      engsrp = engsrp + engs14
      virsrp = virsrp + virs14
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
c     check for undefined potentials and 1-4 range errors
      
      if(mxnode.gt.1)then

        call gstate(safe)
        call gstate(csafe)

      endif
      if(.not.safe) call error(idnode,448)
      if(.not.csafe) call error(idnode,447)
      
#ifdef VAMPIR
      call VTEND(23, ierr)
#endif
      return
      end




