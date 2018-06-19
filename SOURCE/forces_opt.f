      subroutine forces_opt
     x  (loglnk,lgofr,lzeql,idnode,imcon,keyfce,kmax1,
     x  kmax2,kmax3,nhko,nlatt,mxnode,ntpvdw,natms,nstbgr,nstep,
     x  nsteql,numrdf,nospl,fplan,bplan,alpha,dlrpot,drewd,engcpe,
     x  engsrp,epsq,rcut,rvdw,vircpe,virsrp,volm,ilist,jlist,
     x  lentry,lexatm,list,lstvdw,ltpvdw,ltype,nexatm,nauxfft,
     x  lexatm2,nexatm2,key1,key2,key3,buffer,cell,chge,polr,polr2,
     x  toler,ckc,cks,clm,elc,els,emc,ckr,skr,ercp,
     x  ems,enc,ens,erc,fer,fxx,fyy,fzz,ggg,rdf,rsqdf,slm,vvv,
     x  xdf,xxx,ydf,yyy,zdf,zzz,stress,rho,elrcm,vlrcm,ewlbuf,
     x  csp,qqc,txx,tyy,tzz,bspx,bspy,bspz,bsdx,bsdy,bsdz,qqq,
     x  bscx,bscy,bscz,ffttable,ww1,ww2,ww3,ahk,zzn,zzd,sss,hon,
     x  dhn,pp,znp,zgc,zgs,crn,lcp,weight,engunit,
     x  lpolar,lthole,athole,dipx,dipy,dipz,emux,emuy,emuz,
     x  efieldkx,efieldky,efieldkz,efdcrecx,efdcrecy,efdcrecz,
     x  efddmurecx,efddmurecy,efddmurecz,potcc,
     x  uuu)
c     
      
#include "dl_params.inc"
      
      logical lgofr,lzeql,loglnk,lewald,lspme,lhke
      logical lpolar,lcp
      complex*16 qqq,ww1,ww2,ww3,bscx,bscy,bscz

      dimension ahk(0:mxhko),zzn(mxxdf),zzd(mxxdf),sss(mxxdf)
      dimension hon(mxegrd,0:mxhko),dhn(mxegrd,0:mxhko)
      dimension pp(2*mxhko),znp(mxhke,0:2*mxhko)
      dimension zgc(0:2*mxhko),zgs(0:2*mxhko),crn(0:mxhko,0:mxhko)
      dimension ilist(mxxdf),jlist(mxxdf),ltpvdw(mxvdw),nauxfft(4)
      dimension key1(kmaxd),key2(kmaxe),key3(kmaxf)
      dimension lentry(msatms),list(msatms,mxlist)
      dimension nexatm(msatms),lexatm(msatms,mxexcl)
      dimension nexatm2(msatms),lexatm2(msatms,mxexcl)
      dimension lstvdw(mxvdw),ltype(mxatms)
      dimension ckc(mxewld),cks(mxewld),clm(mxewld),slm(mxewld)
      dimension ckr(mxewld),skr(mxewld)
      dimension elc(mxewld,0:1),els(mxewld,0:1),ewlbuf(mxebuf)
      dimension emc(mxewld,0:kmaxb),ems(mxewld,0:kmaxb)
      dimension enc(mxewld,0:kmaxc),ens(mxewld,0:kmaxc)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension vvv(mxgrid,mxvdw),ggg(mxgrid,mxvdw)
      dimension cell(9),chge(mxatms),buffer(mxbuff)
      dimension erc(mxegrd),fer(mxegrd),rho(mxatms)
      dimension ercp(mxegrd,0:3),weight(mxatms)
      dimension rdf(mxrdf,mxvdw),stress(9),elrcm(2),vlrcm(2)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension csp(mxspl),ffttable(mxftab)
      dimension ww1(kmaxd),ww2(kmaxe),ww3(kmaxf)
      dimension bspx(mxspme,mxspl),bspy(mxspme,mxspl)
      dimension bspz(mxspme,mxspl),bsdx(mxspme,mxspl)
      dimension bsdy(mxspme,mxspl),bsdz(mxspme,mxspl)
      dimension qqc(kmaxd,kmaxe,kmaxf)
      dimension qqq(kmaxd,kmaxe,kmaxf)
      dimension bscx(kmaxd),bscy(kmaxe),bscz(kmaxf)
      dimension polr(mxatms),polr2(mxatms),potcc(mxatms)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension emux(mxatms),emuy(mxatms),emuz(mxatms)
      dimension efieldkx(mxatms),efieldky(mxatms),efieldkz(mxatms)
      dimension efdcrecx(mxatms),efdcrecy(mxatms),efdcrecz(mxatms)
      dimension efddmurecx(mxatms),efddmurecy(mxatms),
     x          efddmurecz(mxatms)

#ifdef FFTW
      FFTW_PLAN_TYPE fplan,bplan
#else
      integer fplan, bplan
#endif

#ifdef VAMPIR
      call VTBEGIN(15, ierr)
#endif
      lhke=(keyfce/2.eq.7)
      lspme=(keyfce/2.eq.6)
      lewald=(keyfce/2.eq.1)
c     
c     initialise force arrays
      
      do i=1,natms
        fxx(i)=0.d0
        fyy(i)=0.d0
        fzz(i)=0.d0
      enddo

c     
c     initialise energy and virial accumulators
      
      engcpe=0.d0
      engsrp=0.d0
      vircpe=0.d0
      virsrp=0.d0
c
c     fourier contribution to coulombic forces in Ewald sum

      if(lewald .and. .not.lpolar)then
        
        call ewald1
     x    (idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3,
     x    engac1,viracc,alpha,volm,epsq,cell,chge,xxx,yyy,zzz,
     x    fxx,fyy,fzz,elc,emc,enc,els,ems,ens,ckc,cks,clm,slm,
     x    stress,buffer,ewlbuf)
        
        engcpe=engcpe+engac1
        vircpe=vircpe+viracc

      endif
c
c     hautman-klein-ewald method

      if(lhke)then
c
c     fourier terms of hk-ewald

        call hkewald1
     x    (idnode,mxnode,natms,imcon,nhko,kmax1,kmax2,
     x    engacc,viracc,alpha,epsq,cell,ahk,chge,xxx,
     x    yyy,zzz,fxx,fyy,fzz,elc,emc,els,ems,ckc,cks,
     x    stress,crn,pp,znp,zgc,zgs,buffer)
        
        engcpe=engcpe+engacc
        vircpe=vircpe+viracc
c
c     real space terms of hk-ewald

        call hkewald2
     x    (idnode,mxnode,nhko,nlatt,imcon,natms,engacc,
     x    viracc,drewd,rcut,epsq,cell,chge,ahk,zzn,zzd,
     x    hon,dhn,xxx,yyy,zzz,fxx,fyy,fzz,stress,xdf,ydf,
     x    zdf,sss,rsqdf)

        engcpe=engcpe+engacc
        vircpe=vircpe+viracc
        
      endif
c
c     smooth particle mesh ewald

      if(lspme)then
        
        call ewald_spme
     x    (idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3,nospl,fplan,
     x     bplan,engac1,viracc,alpha,volm,epsq,nauxfft,key1,key2,key3,
     x     cell,chge,xxx,yyy,zzz,txx,tyy,tzz,fxx,fyy,fzz,stress,buffer,
     x     csp,ww1,ww2,ww3,bspx,bspy,bspz,bsdx,bsdy,bsdz,qqc,
     x     qqq,bscx,bscy,bscz,ffttable)

        engcpe=engcpe+engac1
        vircpe=vircpe+viracc

      endif
c
c     calculate induced dipole iteratively

      if (lewald .and. lpolar .and. (.not.lcp)) then
c
c     clear the fields

ccc         toler=4.8d-15

         iloop = 1
         iflag = 0
c
c     start iteration

  101    continue

         do i=1,natms

            if (iloop.eq.1) then

cc               potcc(i)=0.d0

               efieldkx(i)=0.d0
               efieldky(i)=0.d0
               efieldkz(i)=0.d0

               efdcrecx(i)=0.d0
               efdcrecy(i)=0.d0
               efdcrecz(i)=0.d0

            endif

            efddmurecx(i)=0.d0
            efddmurecy(i)=0.d0
            efddmurecz(i)=0.d0

            emux(i)=0.d0
            emuy(i)=0.d0
            emuz(i)=0.d0

         enddo
c
c     electric field contributed by r-space

         eps = epsq
         if(loglnk) eps = eps*2.0d0
         eps2 = epsq*2.0d0

         call rfield
     x    (idnode,imcon,mxnode,natms,lentry,list,ilist,rcut,eps,
     x    nexatm,lexatm,jlist,chge,polr,polr2,xxx,yyy,zzz,alpha,cell,
     x    buffer,efieldkx,efieldky,efieldkz,ercp,drewd,eps2,
     x    dipx,dipy,dipz,emux,emuy,emuz,iloop,lexatm2,nexatm2,
     x    lthole,athole)
c
c     electric field contributed by k-space

         call ewald1p
     x    (lpolar,idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3,
     x    engac1,viracc,alpha,volm,epsq,cell,chge,xxx,yyy,zzz,
     x    fxx,fyy,fzz,elc,emc,enc,els,ems,ens,ckc,cks,clm,slm,
     x    stress,buffer,ewlbuf,efdcrecx,efdcrecy,efdcrecz,
     x    efddmurecx,efddmurecy,efddmurecz,dipx,dipy,dipz,polr,polr2,
     x    ckr,skr,iloop,iflag)
c
c     calculate the converged induced-dipole iteratively

         call iter
     x    (iflag,idnode,mxnode,natms,imcon,
     x     dipx,dipy,dipz,efieldkx,efieldky,efieldkz,
     x     efdcrecx,efdcrecy,efdcrecz,emux,emuy,emuz,
     x     efddmurecx,efddmurecy,efddmurecz,polr,polr2,toler)


        iloop=iloop+1
c
c     check self-consistence

        if (mod(iloop-2,10).eq.0 .and. idnode.eq.0) then

           open(28)
           write(28,*)'iloop = ',iloop-1
           do ii=1,natms
              if (polr2(ii).gt.1.d-6) write(28,'(I6,3X,3f28.20)')
     x           ii,dipx(ii),dipy(ii),dipz(ii)
           enddo
           close(28)

        endif

        if (iloop.gt. 2581) then

           write(nrite,*)'dipole diverged!'
           do ii=1,natms
              if (polr2(ii).gt.1.d-6) write(nrite,'(I6,3X,3f18.6)')
     x            ii,dipx(ii),dipy(ii),dipz(ii)
           enddo

           stop

        endif

        if (iflag.eq.0) goto 101
c
c     calculate coulombic energy in a polar system for
c     charge-dipole and dipole-dipole interactions

        call  coulomb_polar
     x    (idnode,imcon,mxnode,natms,engacp,viracp,
     x    efieldkx,efieldky,efieldkz,
     x    efdcrecx,efdcrecy,efdcrecz,
     x    efddmurecx,efddmurecy,efddmurecz,
     x    dipx,dipy,dipz,emux,emuy,emuz,polr,polr2)

        engcpe=engcpe+engacp

        call ewald1p
     x    (lpolar,idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3,
     x    engac1,viracc,alpha,volm,epsq,cell,chge,xxx,yyy,zzz,
     x    fxx,fyy,fzz,elc,emc,enc,els,ems,ens,ckc,cks,clm,slm,
     x    stress,buffer,ewlbuf,efdcrecx,efdcrecy,efdcrecz,
     x    efddmurecx,efddmurecy,efddmurecz,dipx,dipy,dipz,polr,polr2,
     x    ckr,skr,iloop,iflag)

        engcpe=engcpe+engac1
        vircpe=vircpe+viracc

      endif
c
c     initialize electric field and dipole force arrays
c     if polarizability is used with car-parrinello

      if (lewald .and. lpolar .and. lcp) then

         engcpe=0.d0
         vircpe=0.d0

         do i=1,natms

            efieldkx(i)=0.d0
            efieldky(i)=0.d0
            efieldkz(i)=0.d0

            efdcrecx(i)=0.d0
            efdcrecy(i)=0.d0
            efdcrecz(i)=0.d0

            efddmurecx(i)=0.d0
            efddmurecy(i)=0.d0
            efddmurecz(i)=0.d0

            emux(i)=0.d0
            emuy(i)=0.d0
            emuz(i)=0.d0

         enddo

        call ewald1cp
     x    (lpolar,idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3,
     x    engac1,viracc,alpha,volm,epsq,cell,chge,xxx,yyy,zzz,
     x    fxx,fyy,fzz,elc,emc,enc,els,ems,ens,ckc,cks,clm,slm,
     x    stress,buffer,ewlbuf,efdcrecx,efdcrecy,efdcrecz,
     x    efddmurecx,efddmurecy,efddmurecz,dipx,dipy,dipz,polr,polr2,
     x    ckr,skr)

        engcpe=engcpe+engac1
        vircpe=vircpe+viracc

      endif
c     
c     outer loop over atoms
      
      ii=0

      do i=idnode+1,natms,mxnode
        
        ii=ii+1
c     
c     calculate interatomic distances
        
        do k=1,lentry(ii)
          
          j=list(ii,k)
          ilist(k) = j
          
          xdf(k)=xxx(i)-xxx(j)
          ydf(k)=yyy(i)-yyy(j)
          zdf(k)=zzz(i)-zzz(j)

        enddo
c     
c     periodic boundary conditions

        call images(imcon,0,1,lentry(ii),cell,xdf,ydf,zdf)
c     
c     square of distances

        do k=1,lentry(ii)
          
          rsqdf(k) = xdf(k)**2+ydf(k)**2+zdf(k)**2
          
        enddo
c
c     calculate short range force and potential terms
        
        if(mod(keyfce,2).eq.1) then
          
          call srfrce
     x      (i,lentry(ii),engacc,viracc,rvdw,dlrpot,ilist,ltype,lstvdw,
     x      ltpvdw,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,vvv,ggg,stress)
          
          engsrp=engsrp+engacc
          virsrp=virsrp+viracc

        endif
c     
c     calculate coulombic force and potential terms
c     (real space contributions to ewald sum)
        
        if ((lewald.or.lspme) .and. (.not.lpolar))then
          
          call ewald2
     x      (i,lentry(ii),engacc,viracc,drewd,rcut,epsq,ilist,
     x      chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,erc,fer,stress)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc

        elseif (lewald .and. lpolar .and. (.not.lcp)) then

          call ewald2p
     x         (i,lentry(ii),engacc,viracc,rcut,epsq,ilist,
     x         chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress,
     x         dipx,dipy,dipz,alpha,ercp,drewd)

          engcpe=engcpe+engacc
          vircpe=vircpe+viracc

        elseif (lewald .and. lpolar .and. lcp) then

          call ewald2cp
     x         (i,lentry(ii),engacc,viracc,rcut,epsq,ilist,
     x         chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress,
     x         dipx,dipy,dipz,alpha,ercp,drewd,
     x         efieldkx,efieldky,efieldkz,emux,emuy,emuz)

          engcpe=engcpe+engacc
          vircpe=vircpe+viracc

        elseif (keyfce/2.eq.2) then
          
          call coul2(i,lentry(ii),engacc,viracc,rcut,epsq,
     x      ilist,chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        elseif (keyfce/2.eq.3) then
          
          call coul0(i,lentry(ii),engacc,viracc,rcut,epsq,
     x      ilist,chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        elseif(keyfce/2.eq.4) then
          
          call coul4(i,lentry(ii),engacc,viracc,rcut,epsq,
     x      ilist,chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        elseif(keyfce/2.eq.5) then
          
          call coul3 (i,lentry(ii),engacc,viracc,rcut,epsq,
     x      ilist,chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress)  
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc

        endif
c     
c     accumulate radial distribution functions
        
        if(lgofr.and.((.not.lzeql).or.(nstep.gt.nsteql))) then

          if(mod(nstep,nstbgr).eq.0) then
            
            call rdf0(i,lentry(ii),rcut,ilist,ltype,lstvdw,rdf,rsqdf)
            
          endif

        endif
        
      enddo
c     
c     calculate corrections for intramolecular coulomb terms in
c     Ewald sum
      
      if(lewald.or.lspme.or.lhke)then

        eps = epsq
        if(loglnk) eps = eps*2.0d0
        eps2 = epsq*2.0d0
c     
c     outer loop over atoms
        
        ii=0
        
        do i=idnode+1,natms,mxnode
          
          ii=ii+1
          
c     
c     calculate interatomic distances
          
          do k=1,nexatm(ii)
            
            j=lexatm(ii,k)
            jlist(k)=j
            
            xdf(k)=xxx(i)-xxx(j)
            ydf(k)=yyy(i)-yyy(j)
            zdf(k)=zzz(i)-zzz(j)
            
          enddo
          
c     
c     periodic boundary condition
          
          call images(imcon,0,1,nexatm(ii),cell,xdf,ydf,zdf)
          
c     
c     calculate correction terms
          
          if(lhke)then

            call hkewald3
     x        (i,ii,engacc,viracc,eps,nexatm,lexatm,
     x        chge,xdf,ydf,zdf,fxx,fyy,fzz,stress)
            
          elseif ((lewald.or.lspme) .and. (.not.lpolar))then

            call ewald3
     x        (i,ii,engacc,viracc,alpha,eps,nexatm,
     x        lexatm,chge,xdf,ydf,zdf,fxx,fyy,fzz,stress)

          elseif (lpolar .and. (.not.lcp)) then

            call ewald3p
     x        (i,ii,engacc,viracc,alpha,eps,nexatm,
     x        lexatm,chge,xdf,ydf,zdf,fxx,fyy,fzz,stress,
     x        dipx,dipy,dipz,polr,polr2,athole,lthole)

          elseif (lpolar .and. lcp) then

            call ewald3cp
     x        (i,ii,engacc,viracc,alpha,eps,nexatm,
     x        lexatm,chge,xdf,ydf,zdf,fxx,fyy,fzz,stress,
     x        dipx,dipy,dipz,polr,polr2,athole,lthole,
     x        efieldkx,efieldky,efieldkz,emux,emuy,emuz)

          endif

          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        enddo
c
c     exclude intramolecular charge-dipole interactions

        if (lpolar) then

c
c     outer loop over atoms

          ii=0

          do i=idnode+1,natms,mxnode

            ii=ii+1

c
c     calculate interatomic distances

            do k=1,nexatm2(ii)

              j=lexatm2(ii,k)
              jlist(k)=j

              xdf(k)=xxx(i)-xxx(j)
              ydf(k)=yyy(i)-yyy(j)
              zdf(k)=zzz(i)-zzz(j)

            enddo

c
c     periodic boundary condition

          call images(imcon,0,1,nexatm2(ii),cell,xdf,ydf,zdf)

            if (.not.lcp) then

              call ewald4p
     x          (i,ii,engacc,viracc,alpha,eps2,nexatm2,
     x          lexatm2,chge,xdf,ydf,zdf,fxx,fyy,fzz,stress,
     x          dipx,dipy,dipz,rcut)

            elseif (lcp) then

              call ewald4cp
     x         (i,ii,engacc,viracc,alpha,eps2,nexatm2,
     x         lexatm2,chge,xdf,ydf,zdf,fxx,fyy,fzz,stress,
     x         dipx,dipy,dipz,rcut,efieldkx,efieldky,efieldkz)

            endif

            vircpe=vircpe+viracc

          enddo

        endif
        
      endif

      if(lpolar .and. lcp) then
c
c     global sum of electirc fields for polarizable model
c     with car-parrinello extended lagrangian

        if(mxnode.gt.1) then

          j=0
          do i=1,natms

            buffer(j+1)=efieldkx(i)
            buffer(j+2)=efieldky(i)
            buffer(j+3)=efieldkz(i)
            buffer(j+4)=emux(i)
            buffer(j+5)=emuy(i)
            buffer(j+6)=emuz(i)
            j=j+6

          enddo

          call gdsum(buffer(1),6*natms,buffer(6*natms+1))

          j=0
          do i=1,natms

            efieldkx(i)=buffer(j+1)
            efieldky(i)=buffer(j+2)
            efieldkz(i)=buffer(j+3)
            emux(i)=buffer(j+4)
            emuy(i)=buffer(j+5)
            emuz(i)=buffer(j+6)
            j=j+6

          enddo

        endif

        call  coulomb_polar
     x    (idnode,imcon,mxnode,natms,engacp,viracp,
     x    efieldkx,efieldky,efieldkz,
     x    efdcrecx,efdcrecy,efdcrecz,
     x    efddmurecx,efddmurecy,efddmurecz,
     x    dipx,dipy,dipz,emux,emuy,emuz,polr,polr2)

        engcpe=engcpe+engacp

      endif
c     
c     counter for rdf statistics outside loop structure
      
      if(lgofr.and.((.not.lzeql).or.(nstep.gt.nsteql))) then

        if(mod(nstep,nstbgr).eq.0) then

          numrdf = numrdf + 1

        endif

      endif

c     
c     sum up contributions to short range and coulombic potential
      
      if(mxnode.gt.1) then
        
        buffer(5)=engsrp
        buffer(6)=virsrp
        buffer(7)=engcpe
        buffer(8)=vircpe
        
        call gdsum(buffer(5),4,buffer(1))
        
        engsrp=buffer(5)
        virsrp=buffer(6)
        engcpe=buffer(7)
        vircpe=buffer(8)
        
      endif

      uuu = engsrp + engcpe

      write(*,*)'in forces_opt:'
      write(*,*)'engsrp + engcpe'
      write(*,'(10f8.3)')engsrp/engunit,engcpe/engunit
      write(*,*)'in forces_opt',uuu/engunit

#ifdef VAMPIR
      call VTEND(15, ierr)
#endif
      return
      end
