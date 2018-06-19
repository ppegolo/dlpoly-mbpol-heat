#define USE_PDSUM 1
! 11 JUL 06 - IUCHI - INTRODUCING SRFRCE2 FOR SAFETY
! 10 NOV 05 - IUCHI - ADD LTTM2 IN ARGUMENTS OF SRFRCE
! 07 NOV 05 - IUCHI - POLARIZATION TURNS OFF WHEN NECESSARY (LPOLOFF)
! 07 NOV 05 - IUCHI - ADD USE MODULE TTM_FORCES
! 07 NOV 05 - IUCHI - FIRST ITERATION IS VALIDIATED UNLESS KEYRES=1
! 20 OCT 05 - IUCHI - ADD EATD FROM MODULE
!
      subroutine forces
     x  (lmetal,loglnk,lgofr,lzeql,idnode,imcon,keyfce,kmax1,
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
     x  dhn,pp,znp,zgc,zgs,crn,lcp,lttm,nttm2,listttm2,lads,lacs,
     x  lpolar,lthole,athole,athole12,athole13,
     x  athole_ion,ithole,n_ions,athole_ionwat,
     x  dipx,dipy,dipz,
     x  emux,emuy,emuz,nthole,ascc,ascd,efieldkx,efieldky,efieldkz,
     x  efdcrecx,efdcrecy,efdcrecz,
     x  efddmurecx,efddmurecy,efddmurecz,potcc,keyres,
     x  n_water,gammattm,restartd,laspc,sor_omega,aspc_iter_nstep,
     x  this_step)
c
c***********************************************************************
c
c     dl_poly subroutine for calculating interatomic forces
c     using the verlet neighbour list
c
c     parallel replicated data version
c
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     amended  - t. forester sept 1994
c     amended  - w. smith june 1995 for metal potentialsc
c     amended  - t. yan dec 2003 for polarizable model
c
c     key:
c     keyfce = odd  ------ short range potentials calculated : srfrce
c     = 0,1  ------ no electrostatics
c     = 2,3  ------ Ewald sum                         : ewald1,2,3
c     = 4,5  ------ distance dependent dielectric     : coul2
c     = 6,7  ------ coulombic                         : coul0
c     = 8,9  ------ truncated and shifted coulombic   : coul4
c     = 10,11 ----- reaction field                    : coul3
c     = 12,13 ----- Smooth Particle Mesh Ewald        : ewald[_spme,2,3]
c     = 14,15 ----- Hautman-Klein-Ewald               : hkewald1,2,3
c
c     wl
c     2001/08/31 11:13:45
c     1.10
c     Exp
!
!     Last updated: 11 July 2006 by S. Iuchi
c
c***********************************************************************
c
! from module
      use ttm_forces,      only: lpoloff
      use unit_parameters, only: eatd
      use ps_type_dms,     only: vesp, vesp_k, vesp_r, vesp_s, vesp_c,
     x                           vesp_dc_k, vesp_dc_r, vesp_dc_c, ldms
          ! PP_:
      use heatcurrent, only: check_forces

#include "dl_params.inc"

      logical lgofr,lzeql,loglnk,lmetal,lewald,lspme,lhke
      logical lpolar,lcp,lttm,lthole,lacs,lads,laspc
      logical lfirst
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
      dimension ercp(mxegrd,0:3)
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
      dimension listttm2(mxatms)

      logical restartd

      real(8) :: sor_omega
      integer :: aspc_iter_nstep
      integer :: this_step

      data lfirst/.true./
      save lfirst


      ! PP_:
      real(8) :: fxx_tmp(1024)=0.d0, fxx_tmp_old(1024)

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

c     set up atoms numbers for nodes
      iatm0 = (idnode*natms)/mxnode + 1
      iatm1 = ((idnode+1)*natms)/mxnode

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

#if 0
        ii=0

        do i=1+idnode,natms,mxnode
          ii=ii+1

          write(50+idnode,*) 'list for i=',i,nexatm(ii)
          do k=1,nexatm(ii)
            j=lexatm(ii,k)
            write(50+idnode,*) i,j
          enddo
        enddo
        call gsync()
        stop
#endif

c
c     calculate local density in metals

      if(lmetal) then

        call scdens
     x    (idnode,imcon,mxnode,ntpvdw,natms,rvdw,dlrpot,engacc,virden,
     x    ilist,lentry,list,lstvdw,ltpvdw,ltype,buffer,cell,rho,vvv,
     x    rsqdf,xdf,xxx,ydf,yyy,zdf,zzz,elrcm,vlrcm)

        engsrp=engsrp+engacc
        virsrp=virsrp+rsrden*vlrcm(2)/dble(mxnode)
#ifdef STRESS
        stress(1) = stress(1)-rsrden*vlrcm(2)/(3.d0*dble(mxnode))
        stress(5) = stress(5)-rsrden*vlrcm(2)/(3.d0*dble(mxnode))
        stress(9) = stress(9)-rsrden*vlrcm(2)/(3.d0*dble(mxnode))
#endif
      endif
c
c     fourier contribution to coulombic forces in Ewald sum

      if(lewald .and. .not.lpolar)then

        !         PP_:
              /* fxx_tmp_old = fxx */
        call ewald1
     x    (idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3,
     x    engac1,viracc,alpha,volm,epsq,cell,chge,xxx,yyy,zzz,
     x    fxx,fyy,fzz,elc,emc,enc,els,ems,ens,ckc,cks,clm,slm,
     x    stress,buffer,ewlbuf)

            ! PP_:
          /* fxx_tmp = fxx_tmp + fxx - fxx_tmp_old */

        engcpe=engcpe+engac1
        vircpe=vircpe+viracc

        !PP_:
        write(1111,*) engac1

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

      if( lpoloff ) goto 999    ! skip when pol. evaluation is unnecessary

!     calculate induced dipole iteratively
      if (((lewald .and. lpolar .and. (.not.lcp)) .or.
     x     (lewald .and. lpolar .and. lfirst .and. keyres.ne.1))
     x     .and..not.restartd) then

!     write whether first iteration is necessary or not in case of lcp
      if (idnode.eq.0.and.lcp) write(nrite,*)
     x ' ** NOTE ** : iter. is used for 1st CP/dipoles step'

       iloop = 1
       iflag = 0

!      clear the fields

       efieldkx(1:natms)=0.d0
       efieldky(1:natms)=0.d0
       efieldkz(1:natms)=0.d0

       efdcrecx(iatm0:iatm1)=0.d0
       efdcrecy(iatm0:iatm1)=0.d0
       efdcrecz(iatm0:iatm1)=0.d0

!     start iteration

  101    continue

         efddmurecx(iatm0:iatm1)=0.d0
         efddmurecy(iatm0:iatm1)=0.d0
         efddmurecz(iatm0:iatm1)=0.d0

         emux(1:natms)=0.d0
         emuy(1:natms)=0.d0
         emuz(1:natms)=0.d0

!     electric field contributed by r-space

         eps = epsq
         if(loglnk) eps = eps*2.0d0
         eps2 = epsq*2.0d0

         call rfield
     x    (idnode,imcon,mxnode,natms,lentry,list,ilist,rcut,eps,
     x    nexatm,lexatm,jlist,chge,polr,polr2,xxx,yyy,zzz,alpha,cell,
     x    buffer,efieldkx,efieldky,efieldkz,ercp,drewd,eps2,
     x    dipx,dipy,dipz,emux,emuy,emuz,iloop,lexatm2,nexatm2,
     x    lthole,athole,athole12,athole13,
     x    athole_ion,ithole,n_ions,athole_ionwat,
     x    nthole,lttm,nttm2,listttm2,
     x    lads,ascd,n_water)

!     electric field contributed by k-space

         call ewald1p
     x    (lpolar,idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3,
     x    engac1,viracc,alpha,volm,epsq,cell,chge,xxx,yyy,zzz,
     x    fxx,fyy,fzz,elc,emc,enc,els,ems,ens,ckc,cks,clm,slm,
     x    stress,buffer,ewlbuf,efdcrecx,efdcrecy,efdcrecz,
     x    efddmurecx,efddmurecy,efddmurecz,dipx,dipy,dipz,polr,polr2,
     x    ckr,skr,iloop,iflag)

!     calculate the converged induced-dipole iteratively

         ! when doing iterations, use sor_omega of 0.3
         this_sor_omega=sor_omega
         if(laspc) then
             if(this_step .lt. 5) then
                 this_sor_omega=0.3
             endif
         endif

         call iter
     x    (iflag,idnode,mxnode,natms,imcon,
     x     dipx,dipy,dipz,efieldkx,efieldky,efieldkz,
     x     efdcrecx,efdcrecy,efdcrecz,emux,emuy,emuz,
     x     efddmurecx,efddmurecy,efddmurecz,polr2,toler,this_sor_omega,
     x     buffer)

        iloop=iloop+1

!     check self-consistence

        if (iloop.gt. 22581) then
           write(nrite,*)'dipole diverged!'
           do ii=iatm0,iatm1
              if (polr2(ii).gt.1.d-6) write(nrite,'(I6,3X,3f18.6)')
     x           ii,dipx(ii),dipy(ii),dipz(ii)
           enddo

           stop
        endif

        if(laspc) then
            if(this_step .le. aspc_iter_nstep) then
!                write(89,*) this_step
                if ((iflag.eq.0)) goto 101 ! no convergence, loop is over
            else
                iflag = 1 ! mark as converged after one iteration
            endif
        else if(iflag.eq.0)then
             goto 101 ! no convergence, loop is over
         endif

        lfirst = .false.
        if(lcp) goto 102

!     dipoles converged!
!     calculate coulombic energy in a polar system for
!     charge-dipole and dipole-dipole interactions

        call  coulomb_polar
     x    (idnode,imcon,mxnode,natms,engacp,viracp,
     x    efieldkx,efieldky,efieldkz,
     x    efdcrecx,efdcrecy,efdcrecz,
     x    efddmurecx,efddmurecy,efddmurecz,
     x    dipx,dipy,dipz,emux,emuy,emuz,polr2)

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

  102 continue

!     initialize electric field and dipole force arrays
!     if polarizability is used with car-parrinello

      if (lewald .and. lpolar .and. lcp) then

         engcpe=0.d0
         vircpe=0.d0

         efdcrecx(iatm0:iatm1)=0.d0
         efdcrecy(iatm0:iatm1)=0.d0
         efdcrecz(iatm0:iatm1)=0.d0

         efddmurecx(iatm0:iatm1)=0.d0
         efddmurecy(iatm0:iatm1)=0.d0
         efddmurecz(iatm0:iatm1)=0.d0

         efieldkx(1:natms)=0.d0
         efieldky(1:natms)=0.d0
         efieldkz(1:natms)=0.d0

         emux(1:natms)=0.d0
         emuy(1:natms)=0.d0
         emuz(1:natms)=0.d0

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
!
!     when lpoloff

 999  continue

      if (lewald .and. lpolar .and. lpoloff ) then

         efdcrecx(iatm0:iatm1) = 0.d0
         efdcrecy(iatm0:iatm1) = 0.d0
         efdcrecz(iatm0:iatm1) = 0.d0

         efddmurecx(iatm0:iatm1)=0.d0
         efddmurecy(iatm0:iatm1)=0.d0
         efddmurecz(iatm0:iatm1)=0.d0

         dipx(iatm0:iatm1) = 0.0d0
         dipy(iatm0:iatm1) = 0.0d0
         dipz(iatm0:iatm1) = 0.0d0

         if( .not. lcp ) then

            iloop=0
            iflag=1 ! equivalent to converged case in iteration method

            call ewald1p
     x           (lpolar,idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3,
     x           engac1,viracc,alpha,volm,epsq,cell,chge,xxx,yyy,zzz,
     x           fxx,fyy,fzz,elc,emc,enc,els,ems,ens,ckc,cks,clm,slm,
     x           stress,buffer,ewlbuf,efdcrecx,efdcrecy,efdcrecz,
     x           efddmurecx,efddmurecy,efddmurecz,dipx,dipy,dipz,polr,
     x           polr2,ckr,skr,iloop,iflag)

#ifndef DISABLE_ASSERT
            do i=1,natms

               if (abs(efdcrecx(i)).gt.1.0d-10)
     x            call afailed('forces.f',__LINE__)
               if (abs(efdcrecy(i)).gt.1.0d-10)
     x            call afailed('forces.f',__LINE__)
               if (abs(efdcrecz(i)).gt.1.0d-10)
     x            call afailed('forces.f',__LINE__)

               if (abs(efddmurecx(i)).gt.1.0d-10)
     x            call afailed('forces.f',__LINE__)
               if (abs(efddmurecy(i)).gt.1.0d-10)
     x            call afailed('forces.f',__LINE__)
               if (abs(efddmurecz(i)).gt.1.0d-10)
     x            call afailed('forces.f',__LINE__)

            enddo
#endif /* DISABLE_ASSERT */

            engcpe=engcpe+engac1
            vircpe=vircpe+viracc

         end if

!!!! add lcp case here in the future

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
c     calculate metal forces and potential

        if(lmetal) then

          call suttchen
     x      (i,lentry(ii),ntpvdw,engacc,viracc,rvdw,dlrpot,ilist,ltype,
     x      lstvdw,ltpvdw,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,vvv,ggg,stress,
     x      rho)

          engsrp=engsrp+engacc
          virsrp=virsrp+viracc

        endif
c
c     calculate short range force and potential terms

        if(mod(keyfce,2).eq.1) then

           if( lttm ) then
              call srfrce2
     $             (i,lentry(ii),engacc,viracc,rvdw,dlrpot,ilist,ltype,
     $             lstvdw,ltpvdw,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,vvv,
     $             ggg,stress,lttm)
           else


              call srfrce
     x             (i,lentry(ii),engacc,viracc,rvdw,dlrpot,ilist,ltype,
     $             lstvdw,ltpvdw,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,vvv,
     $             ggg,stress)


           endif

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

           if( lpoloff ) then

              dipx(iatm0:iatm1) = 0.0d0
              dipy(iatm0:iatm1) = 0.0d0
              dipz(iatm0:iatm1) = 0.0d0

           endif

           call ewald2p
     x          (i,lentry(ii),engacc,viracc,rcut,epsq,ilist,
     x          chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress,
     x          polr,polr2,athole,athole_ion,ithole,n_ions,
     x          athole_ionwat,
     x          nthole,lthole,lacs,lads,
     x          ascc,ascd,dipx,dipy,dipz,alpha,ercp,drewd,
     x          lttm,nttm2,listttm2)

           engcpe=engcpe+engacc
           vircpe=vircpe+viracc

        elseif (lewald .and. lpolar .and. lcp) then

          call ewald2cp
     x         (i,lentry(ii),engacc,viracc,rcut,epsq,ilist,
     x         chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress,
     x         dipx,dipy,dipz,alpha,ercp,drewd,polr,polr2,lacs,
     x         efieldkx,efieldky,efieldkz,emux,emuy,emuz,
     x         lthole,athole,nthole,lttm,nttm2,listttm2,
     x         ascc,ascd,lads)

          engcpe=engcpe+engacc
          vircpe=vircpe+viracc

        elseif (keyfce/2.eq.2) then

          call coul2(i,lentry(ii),engacc,viracc,rcut,epsq,
     x      ilist,chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress)

          engcpe=engcpe+engacc
          vircpe=vircpe+viracc

        elseif (keyfce/2.eq.3) then

c        write(nrite,*)'i,lentry(ii)',i,lentry(ii)
c        do m=1,lentry(ii)
c        write(nrite,*)i,ilist(m)
c        enddo

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

             if( lpoloff ) then

                dipx(iatm0:iatm1) = 0.0d0
                dipy(iatm0:iatm1) = 0.0d0
                dipz(iatm0:iatm1) = 0.0d0

             end if

             call ewald3p
     x            (i,ii,engacc,viracc,alpha,eps,nexatm,
     x            lexatm,chge,xdf,ydf,zdf,fxx,fyy,fzz,stress,
     x            dipx,dipy,dipz,polr,polr2,athole12,athole13,
     x            nthole,lthole,lttm,nttm2,listttm2,n_water,natms)

          elseif (lpolar .and. lcp) then

            call ewald3cp
     x        (i,ii,engacc,viracc,alpha,eps,nexatm,
     x        lexatm,chge,xdf,ydf,zdf,fxx,fyy,fzz,stress,
     x        dipx,dipy,dipz,polr,polr2,athole12,athole13,
     x        nthole,lthole,lttm,nttm2,listttm2,
     x        efieldkx,efieldky,efieldkz,emux,emuy,emuz)

          endif

          engcpe=engcpe+engacc
          vircpe=vircpe+viracc

        enddo

      endif ! lewald.or.lspme.or.lhke
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

      endif ! lpolar

      if(lpolar .and. lcp) then
c
c     global sum of electirc fields for polarizable model
c     with car-parrinello extended lagrangian

        if(mxnode.gt.1) then
#if USE_PDSUM
           call pdsum6(idnode,mxnode,natms,
     x         efieldkx,efieldky,efieldkz,emux,emuy,emuz,mxbuff,buffer)
#else
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
#endif /* partial reduce */
        endif

        call  coulomb_polar
     x    (idnode,imcon,mxnode,natms,engacp,viracp,
     x    efieldkx,efieldky,efieldkz,
     x    efdcrecx,efdcrecy,efdcrecz,
     x    efddmurecx,efddmurecy,efddmurecz,
     x    dipx,dipy,dipz,emux,emuy,emuz,polr2)

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
c     contributions of charge derivatives to forces
      if (ldms) then

         do i = 1, natms
            vesp(i) = vesp_k(i) + vesp_r(i) + vesp_s(i) + vesp_c(i)
     x              + vesp_dc_k(i) + vesp_dc_r(i) + vesp_dc_c(i)
!            write(101,'(8f15.6)') vesp(i),vesp_k(i), vesp_r(i),
!     x                            vesp_s(i), vesp_c(i), vesp_dc_k(i),
!     x                            vesp_dc_r(i), vesp_dc_c(i)
         end do

         if (mxnode.gt.1) call gdsum(vesp,natms,buffer)

         if (lttm) then
            virqdf=0.d0
            call qdforce(listttm2,n_water,nttm2,natms,imcon,cell,
     x           gammattm,xxx,yyy,zzz,fxx,fyy,fzz,stress,virqdf)
            vircpe=vircpe+virqdf
         end if ! lttm

      end if ! ldms
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

      ! PP_:
      /* call check_forces(0.d0,fxx_tmp) */

#ifdef VAMPIR
      call VTEND(15, ierr)
#endif
      return
      end
