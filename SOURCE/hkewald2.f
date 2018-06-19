      subroutine hkewald2
     x  (idnode,mxnode,nhko,nlatt,imcon,natms,engcpe,
     x  vircpe,drewd,rcut,epsq,cell,chge,ahk,zzn,zzd,
     x  hon,dhn,xxx,yyy,zzz,fxx,fyy,fzz,stress,xdf,ydf,
     x  zdf,sss,rsqdf) 
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating real-space contributions to
c     the hautman-klein-ewald electrostatic method
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 2000
c     author    - w. smith  may 2000
c     
c     wl
c     2001/08/31 11:13:46
c     1.2
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension ahk(0:mxhko),zzn(mxxdf),zzd(mxxdf)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension chge(mxatms),stress(9),cell(9)
      dimension hon(mxgrid,0:mxhko),sss(mxxdf)
      dimension dhn(mxgrid,0:mxhko),rcell(9)
      dimension nix(25),niy(25),rsqdf(mxxdf)
      data nix/ 0, 1, 1, 0,-1,-1,-1, 0, 1, 2, 2,
     x  2, 1, 0,-1,-2,-2,-2,-2,-2,-1, 0, 1, 2, 2/
      data niy/ 0, 0, 1, 1, 1, 0,-1,-1,-1, 0, 1,
     x  2, 2, 2, 2, 2, 1, 0,-1,-2,-2,-2,-2,-2,-1/
CDIR$ CACHE_ALIGN fi
#ifdef VAMPIR
      call VTBEGIN(173, ierr)
#endif
c     
c     check boundary condition

      if(imcon.ne.6)call error(idnode,66)
c
c     number of neighbouring real space cells

      if(nlatt.gt.2)call error(idnode,488)
      step=dble(2*nlatt+1)
      nboxes=(2*nlatt+1)**2
c     
c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0
c     
c     set cutoff condition for pair forces
      
      rcsq=(step*rcut)**2
c     
c     reciprocal of interpolation interval
      
      rdrewd = 1.d0/drewd
c     
c     reciprocal cell

      call invert(cell,rcell,det)
      do i=1,9
        rcell(i)=rcell(i)/step
      enddo
#ifdef STRESS
c     
c     initialise stress tensor accumulators
      strs3 = 0.d0
      strs6 = 0.d0
      strs9 = 0.d0
      strs1 = 0.d0
      strs2 = 0.d0
      strs5 = 0.d0
#endif
c     
c     loop over image cells, starting with central cell

      ma=1
      mpm2=natms/2
      npm2=(natms-1)/2

      do k=1,nboxes

        last=natms
        dcx=dble(nix(k))
        dcy=dble(niy(k))
        udx=cell(1)*dcx+cell(4)*dcy
        udy=cell(2)*dcx+cell(5)*dcy
c     
c     outer loop over atoms

        do m=ma,mpm2

          fac=r4pie0/epsq
          if(m.eq.0)fac=fac*0.5d0
          if(m.gt.npm2)last=mpm2
c     
c     set initial array values

          ii=0
          do i=idnode+1,last,mxnode
            
            ii=ii+1
            chgea = fac*chge(i)
            
            if(chgea.ne.0.d0) then
              
              j=i+m
              if(j.gt.natms)j=j-natms
              
              chgprd=chgea*chge(j)
              
              if(chgprd.ne.0.d0) then

                zzn(ii)=1.d0
                zzd(ii)=0.d0
c     
c     calculate interatomic separation
                
                ddx=xxx(i)-xxx(j)+udx
                ddy=yyy(i)-yyy(j)+udy
                ssx=rcell(1)*ddx+rcell(4)*ddy
                ssy=rcell(2)*ddx+rcell(5)*ddy
                ssx=ssx-nint(ssx)
                ssy=ssy-nint(ssy)
                xdf(ii)=step*(ssx*cell(1)+ssy*cell(4))
                ydf(ii)=step*(ssx*cell(2)+ssy*cell(5))
                zdf(ii)=zzz(i)-zzz(j)
                rsqdf(ii)=xdf(ii)**2+ydf(ii)**2+zdf(ii)**2
                
              endif

            endif

          enddo
c     
c     loop over HK orders

          do n=0,nhko
c     
c     inner loop over atoms

            ii=0
            do i=idnode+1,last,mxnode

              ii=ii+1
              chgea = fac*chge(i)

              if(chgea.ne.0.d0) then
                
                j=i+m
                if(j.gt.natms)j=j-natms

                chgprd=chgea*chge(j)
                
                if(chgprd.ne.0.d0) then
c     
c     apply truncation of potential
                  
                  ssq=rsqdf(ii)-zdf(ii)*zdf(ii)

                  if(rcsq.gt.ssq)then
c     
c     calculate potential energy and virial
                    
                    coul=0.d0
                    fcoul=0.d0
                    rrr = sqrt(rsqdf(ii))
                    sss(ii)=sqrt(ssq)

                    if(n.eq.0)then

                      coul = chgprd/rrr
                      fcoul = coul/rsqdf(ii)

                    endif
c     
c     interpolation parameters

                    l0=int(sss(ii)*rdrewd)
                    ppp=sss(ii)*rdrewd-dble(l0)
                    l0=l0+1
                    l1=l0+1
                    l2=l0+2
c     
c     calculate interaction energy using 3-point interpolation
                    
                    vk0 = hon(l0,n)
                    vk1 = hon(l1,n)
                    vk2 = hon(l2,n)
                    t1 = vk0 + (vk1 - vk0)*ppp
                    t2 = vk1 + (vk2 - vk1)*(ppp - 1.0d0)
                    
                    eterm=(t1+(t2-t1)*ppp*0.5d0)*ahk(n)*chgprd
                    engcpe=engcpe+coul-eterm*zzn(ii)
c     
c     calculate forces using 3pt interpolation
                    
                    vk0 = dhn(l0,n)
                    vk1 = dhn(l1,n)
                    vk2 = dhn(l2,n)
                    
                    t1 = vk0 + (vk1 - vk0)*ppp
                    t2 = vk1 + (vk2 - vk1)*(ppp - 1.0d0)
c     
c     calculate in-plane forces
                    
                    egamma=fcoul+
     x                (t1+(t2-t1)*ppp*0.5d0)*chgprd*zzn(ii)*ahk(n)
                    fx=egamma*xdf(ii)
                    fy=egamma*ydf(ii)
c     
c     calculate perpendicular forces
                    
                    fz=fcoul*zdf(ii)+2.d0*dble(n)*eterm*zzd(ii)
c
c     add to force accumulators

                    fxx(i)=fxx(i)+fx
                    fyy(i)=fyy(i)+fy
                    fzz(i)=fzz(i)+fz

                    fxx(j)=fxx(j)-fx
                    fyy(j)=fyy(j)-fy
                    fzz(j)=fzz(j)-fz
c     
c     reset zzn array for next order of convergence function
              
                    zzd(ii)=zzn(ii)*zdf(ii)
                    zzn(ii)=zzd(ii)*zdf(ii)
              
#ifdef STRESS
c     
c     calculate stress tensor
                    
                    strs1 = strs1 + xdf(ii)*fx
                    strs2 = strs2 + xdf(ii)*fy
                    strs3 = strs3 + xdf(ii)*fz
                    
                    strs5 = strs5 + ydf(ii)*fy
                    strs6 = strs6 + ydf(ii)*fz
                    
                    strs9 = strs9 + zdf(ii)*fz
#endif
                  endif
                  
                endif
                
              endif

            enddo

          enddo

        enddo

        ma=0

      enddo
c
c     calculate virial

      vircpe=-(strs1+strs5+strs9)

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
      call VTEND(173, ierr)
#endif
      return
      end
