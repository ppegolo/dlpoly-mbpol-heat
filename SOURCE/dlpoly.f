! 13 APR 07 - IUCHI - ADD SITE_DIPOLE FILE FOR IND. DIPOLE (NON TTM WATER)
! 23 FEB 07 - IUCHI-  ADD if( nstep > 0 ) flag_com = .false. BEFORE MD
! 27 NOV 06 - IUCHI - DEBUG: ADD LOGICAL LRETURN
! 24 NOV 06 - IUCHI - COMMENT OUT VX,Y,Z0_AT_T
! 22 NOV 06 - IUCHI - SKIP DIPOLE GUESS AT FIRST FOUR STEPS
! 12 NOV 06 - IUCHI - INTRODUCE VX,Y,Z0_AT_T FOR VELOCITY WITHOUT SHAKE
! 14 SEP 06 - IUCHI - INTRODUCE VX,Y,Z_AT_T FOR VELOCITY BELONGED TO NT
! 13 SEP 06 - IUCHI - CHANGE PLACE OF OUTPUT_TTM2_FORCE
! 28 AUG 06 - IUCHI - COMMENT OUT SET_GUESS_DIP
! 05 AUG 06 - IUCHI - INTRODUCE DIPX, Y, Z_OLD FOR LCP CASE
! 07 JUL 06 - IUCHI - GOOD INITIAL GUESS: dip*xyz, set_guess_dip
! 06 JUL 06 - IUCHI - CHANGE PLACES OF OUTPUT_TTM* ROUNES
! 28 APR 06 - IUCHI - ADD STRESS AND VIRIAL IN QDFORCE
! 06 MAR 06 - IUCHI - FLAG LUNFORMAT TO CONTROL FORMAT OF HISTORY FILE
! 23 FEB 06 - IUCHI - REMOVE UNNECESSARY ARGUMENTS (OUTPUT_TTM2_DIP, QDFORCE) 
!                     SEE MAKEFILE $DEBUG1 -warn
! 09 FEB 06 - IUCHI - MODIFY ARGUMENTS IN OUTPUT_TTM2_DIP
! 03 FEB 06 - IUCHI - XXOLD, YYOLD, ZZOLD FOR OUTPUT_TTM2_COM 
! 14 DEC 05 - IUCHI - MODIFY CALL OUTPUT_TTM2_FORCE FOR DISTRIBUTE FORCES
! 12 DEC 05 - IUCHI - XXOLD, YYOLD, ZZOLD FOR HISTORY AND OUTPUT_TTM2_FORCE
! 06 DEC 05 - IUCHI - MODIFY CALL OUTPUT_TTM2_FORCE
! 01 DEC 05 - IUCHI - INTRODUCING STEEPEST-DESCENT LIKE OPTIMIZATION
! 28 NOV 05 - IUCHI - CALL QDFORCE IF LDMS
! 18 NOV 05 - IUCHI - CALL INIT AND END_INTRA_GRAD_DMS IF LDMS
! 17 NOV 05 - IUCHI - CALL QDFORCE IF LDMS
! 17 NOV 05 - IUCHI - ADD MODULE PS_TYPE_DMS
! 16 NOV 05 - IUCHI - CALL QDETERMINE
! 16 NOV 05 - IUCHI - CALL INIT_ AND END_DIPOLE_MOMENTS
! 16 NOV 05 - IUCHI - ADD USE MODULE DIPOLE_MOMENTS
! 14 NOV 05 - IUCHI - INTRODUCING FLAG LDECOMP
! 10 NOV 05 - IUCHI - ADD LTTM2 IN ARGUMENTS OF BNDFRC AND ANGFRC
! 09 NOV 05 - IUCHI - CALL INIT_ AND END_TTM_FORCES
! 09 NOV 05 - IUCHI - INTRODUCING FLAG LTIP4P_GEO
! 07 NOV 05 - IUCHI - ADD USE MODULE TTM_FORCES
! 02 NOV 05 - IUCHI - CALL CHECKING COM OF THE SYSTEM
! 01 NOV 05 - IUCHI - OUTPUT INTRA-POTENTIAL FOR TTM MODEL
! 28 OCT 05 - IUCHI - ADD P-S TYPE PES FOR TTM MODEL
! 25 OCT 05 - IUCHI - M SITES INFORMATION FOR REVCON 
! 25 OCT 05 - IUCHI - ADD FINITE DIFFERENCE METHOD: FORDEBUG.F90 
!        05 - IUCHI - LHEAD OPTION
! 17 OCT 05 - IUCHI - CALL OUTPUT_TTM2_DIPOLE, _COM IF LTTM2
!
      program dlpoly
c     
c***********************************************************************
c     
c     dl_poly is an cclrc/ccp5 program package for the dynamical 
c     simulation of molecular systems.
c     
c     dl_poly is the property of the cclrc daresbury laboratory, 
c     daresbury, warrington wa4 4ad. no part of the package may
c     be redistributed to third parties without the consent of
c     daresbury laboratory.
c     
c     dl_poly is available free of charge to academic institutions
c     engaged in non-commercial research only. potential users not
c     in this category must consult the ccp5 program librarian at 
c     daresbury to negotiate terms of use.
c     
c     neither the cclrc, daresbury laboratory, ccp5 nor the authors
c     of this package claim that it is free from errors and do not
c     accept liability for any loss or damage that may arise from
c     its use. it is the users responsibility to verify that the 
c     package dl_poly is fit for the purpose the user intends for
c     it.
c     
c     users of this package are recommended to consult the dl_poly
c     user and reference manuals for the full terms and conditions
c     of its use.
c     
c     authors: w.smith and t.r.forester 1995
c     copyright daresbury laboratory 1995
c     
c     release 2.0 
c     
c     wl
c     2001/08/31 11:13:44
c     1.18
c     Exp
!
!     Last updated: Apr 7, 2007 by S. Iuchi
c     
c***********************************************************************
c     
! from module
      use dipole_moments,  only: init_dipole_moments,
     $     end_dipole_moments, dip1x, dip1y, dip1z, dip2x, dip2y, dip2z,
     $     dip3x, dip3y, dip3z
      use format_history,  only: lunformat
      use param_type_ttm,  only: gammattm
      use ps_type_dms,     only: ldms, init_intra_grad_dms,
     $                           end_intra_grad_dms
      use ttm_forces,      only: lpoloff, ltip4p_geo, ldecomp,
     $                           fx_ttm_per, fy_ttm_per, fz_ttm_per, 
     $                           fx_ttm_ind, fy_ttm_ind, fz_ttm_ind, 
     $                           init_ttm_forces, end_ttm_forces,
     $                           vx_at_t, vy_at_t, vz_at_t
c$$$     $                          ,vx0_at_t, vy0_at_t, vz0_at_t
      use unit_parameters, only: eatd
      use optimize,        only: lgeoopt

      use aspc_module, only:
     x        aspc_init, aspc_predict, aspc_set, aspc_fini
      

      use multibead, only: mb_init, mb_finalize, bead_rank, bead_size

      parameter (mxmc=4,mxmi=36,mxmr=89)

#include "dl_params.inc"
      
      logical lhead
      logical ltscal,lzeql,loptim,lopt2,ltraj,lgofr,lpgr,lfcap,safe,
     x        safeq
      logical newjob,newlst,lneut,loglnk,lnsq,lzden,lshmov,lcnb,lmetal
      logical safep,lpolar,lthole,lacs,lads,lcp,ldpts,lttm,lttm3
      logical lfd, llfirst, lfdx, lfdy, lfdz, lreturn  ! for lfd
      logical flag_com          ! for COM velocity check
      
      data llfirst/.true./
      data lfdx, lfdy, lfdz/.true.,.true.,.true./
      save llfirs
!c
!c     1 e.Angstrom = 4.803239 Debye
!      data eatd/4.803239d0/
      
      character*80  sysname
      character*80  cfgname

      dimension memc(mxmc),memi(mxmi),memr(mxmr),nauxfft(4)
      dimension cell(9),celprp(10),elrcm(2),eta(9),stresl(9)
      dimension stress(9),vlrcm(2),npmf(2)

      logical, allocatable :: lms(:)

      character*40, allocatable :: molnam(:)
      character*8, allocatable :: atmnam(:)
      character*8, allocatable :: sitnam(:)
      character*8, allocatable :: unqatm(:)

#ifdef SHMEM
      pointer (buf_ptr,buffer)
      pointer (tbf_ptr,tbuffer)
      real, dimension(1) :: buffer
      real, dimension(1) :: tbuffer
#elif SGISHMEM
      pointer (buf_ptr,buffer)
      pointer (tbf_ptr,tbuffer)
      real*8, dimension(1) :: buffer
      real*8, dimension(1) :: tbuffer
#else
      real*8, allocatable :: buffer(:)
      real*8, allocatable :: tbuffer(:)
#endif
      real*8, allocatable :: ahk(:),hon(:,:),dhn(:,:)
      real*8, allocatable :: fon(:,:),pp(:)
      real*8, allocatable :: accum(:),rotmin(:),gaxs(:,:)
      real*8, allocatable :: amsd(:),dens0(:),esig1(:)
      real*8, allocatable :: chge(:),rcut3b(:),ewlbuf(:)
      real*8, allocatable :: ckc(:),cks(:),clm(:),slm(:)
      real*8, allocatable :: dxp(:),dyp(:),dzp(:)
      real*8, allocatable :: dxt(:),dyt(:),dzt(:)
      real*8, allocatable :: dxx(:),dyy(:),dzz(:)
      real*8, allocatable :: elc(:,:),els(:,:)
      real*8, allocatable :: emc(:,:),ems(:,:)
      real*8, allocatable :: enc(:,:),ens(:,:)
      real*8, allocatable :: erc(:),fer(:),dens(:)
      real*8, allocatable :: ffttable(:)
      real*8, allocatable :: flx(:),fly(:),flz(:)
      real*8, allocatable :: fpx(:),fpy(:),fpz(:)
      real*8, allocatable :: fxx(:),fyy(:),fzz(:)
      real*8, allocatable :: gcmx(:),gcmy(:),gcmz(:)
      real*8, allocatable :: gcmx1(:),gcmy1(:),gcmz1(:)
      real*8, allocatable :: gvx1(:),gvy1(:),gvz1(:)
      real*8, allocatable :: gvxx(:),gvyy(:),gvzz(:),gmass(:)
      real*8, allocatable :: gxx(:,:),gyy(:,:),gzz(:,:)
      real*8, allocatable :: omx(:),omy(:),omz(:)
      real*8, allocatable :: omx1(:),omy1(:),omz1(:)
      real*8, allocatable :: opx(:),opy(:),opz(:)
      real*8, allocatable :: oqx(:),oqy(:),oqz(:)
      real*8, allocatable :: pmfwght(:),redmass(:),dsq(:)
      real*8, allocatable :: prmang(:,:),prmfld(:)
      real*8, allocatable :: prmcon(:),prmbnd(:,:),rcut4b(:)
      real*8, allocatable :: prmdih(:,:),prmvdw(:,:)
      real*8, allocatable :: prmtet(:,:),prmshl(:),prminv(:,:)
      real*8, allocatable :: q0(:),q1(:),q2(:),q3(:)
      real*8, allocatable :: qn0(:),qn1(:),qn2(:),qn3(:)
      real*8, allocatable :: rdf(:,:),prmtbp(:,:),prmfbp(:,:)
      real*8, allocatable :: rotinx(:,:),rotiny(:,:),rotinz(:,:)
      real*8, allocatable :: ssqval(:),stkval(:,:)
      real*8, allocatable :: stpval(:),sumval(:)
      real*8, allocatable :: tqx(:),tqy(:),tqz(:)
      real*8, allocatable :: xold(:),yold(:),zold(:)
      real*8, allocatable :: txx(:),tyy(:),tzz(:)
      real*8, allocatable :: uxx(:),uyy(:),uzz(:)
      real*8, allocatable :: vvv(:,:),ggg(:,:)
      real*8, allocatable :: vx1(:),vy1(:),vz1(:)
      real*8, allocatable :: vxx(:),vyy(:),vzz(:)
      real*8, allocatable :: weight(:),chgsit(:),wgtsit(:)
      real*8, allocatable :: polarsit(:),polr(:),potcc(:)
      real*8, allocatable :: polarsit2(:),polr2(:)
      real*8, allocatable :: dpms(:)
      real*8, allocatable :: dipx(:),dipy(:),dipz(:)
      real*8, allocatable :: dipx_old(:), dipy_old(:), dipz_old(:) ! lcp case
      real*8, allocatable :: vdxx(:),vdyy(:),vdzz(:)
      real*8, allocatable :: fdxx(:),fdyy(:),fdzz(:)
      real*8, allocatable :: udxx(:),udyy(:),udzz(:)
      real*8, allocatable :: emux(:),emuy(:),emuz(:)
      real*8, allocatable :: efieldkx(:),efieldky(:),efieldkz(:)
      real*8, allocatable :: efdcrecx(:),efdcrecy(:),efdcrecz(:)
      real*8, allocatable :: efddmurecx(:),efddmurecy(:),
     x                       efddmurecz(:)
      real*8, allocatable :: ckr(:),skr(:)
      real*8, allocatable :: p(:),xi(:),g(:)
      real*8, allocatable :: h(:),pcom(:),xicom(:)
      real*8, allocatable :: xt(:),ooo(:,:)
      real*8, allocatable :: ercp(:,:)
      real*8, allocatable :: xa(:,:),ya(:,:),za(:,:)
      real*8, allocatable :: xdab(:),ydab(:),zdab(:)
      real*8, allocatable :: xdab2(:),ydab2(:),zdab2(:)
      real*8, allocatable :: listin2(:),listbnd2(:),listbnd3(:)
      real*8, allocatable :: xdac(:),ydac(:),zdac(:)
      real*8, allocatable :: xdad(:),ydad(:),zdad(:)
      real*8, allocatable :: xdbc(:),ydbc(:),zdbc(:)
      real*8, allocatable :: xdcd(:),ydcd(:),zdcd(:)
      real*8, allocatable :: xdf(:),ydf(:),zdf(:),rsqdf(:)
      real*8, allocatable :: xx0(:),yy0(:),zz0(:)
      real*8, allocatable :: xxs(:),yys(:),zzs(:)
      real*8, allocatable :: xxt(:),yyt(:),zzt(:)
      real*8, allocatable :: xxx(:),yyy(:),zzz(:)
      real*8, allocatable :: xxold(:),yyold(:),zzold(:) ! for HISTORY etc.
      real*8, allocatable :: xxx1(:),yyy1(:),zzz1(:)  
      real*8, allocatable :: xcm(:),ycm(:),zcm(:)
      real*8, allocatable :: vxcm(:),vycm(:),vzcm(:)
      real*8, allocatable :: xmsd(:),ymsd(:),zmsd(:)
      real*8, allocatable :: zdens(:,:),rho(:)
      real*8, allocatable :: zumval(:),ravval(:)
      real*8, allocatable :: csp(:),qqc(:,:,:)
      real*8, allocatable :: bspx(:,:), bspy(:,:), bspz(:,:)
      real*8, allocatable :: bsdx(:,:), bsdy(:,:), bsdz(:,:)
      real*8, allocatable :: crn(:,:),znp(:,:),zgc(:),zgs(:)
      real*8, allocatable :: zzn(:),zzd(:),sss(:)
      real*8, allocatable :: fanaly(:,:) ! for lfd

      integer, allocatable :: ilist(:),jlist(:),ind(:,:)
      integer, allocatable :: keybnd(:),keyang(:),keydih(:)
      integer, allocatable :: lentry(:),list(:,:),ltype(:)
      integer, allocatable :: lexatm(:,:),nexatm(:),noxatm(:)
      integer, allocatable :: lexatm2(:,:),nexatm2(:),noxatm2(:)
      integer, allocatable :: lexsit(:,:),nexsit(:)
      integer, allocatable :: link(:),lct(:),lst(:)
      integer, allocatable :: lishap(:),listot(:),numtyp(:)
      integer, allocatable :: listme(:),lashap(:),listin(:)
      integer, allocatable :: listpm(:),lstpmt(:),numfrz(:)
      integer, allocatable :: listtet(:,:),keytet(:),lsttet(:)
      integer, allocatable :: lstang(:,:),listang(:,:)
      integer, allocatable :: lstbnd(:,:),listbnd(:,:),keyinv(:)
      integer, allocatable :: lstcon(:,:),listcon(:,:)
      integer, allocatable :: lstcsit(:),lstpmf(:,:)
      integer, allocatable :: lstdih(:,:),listdih(:,:),lfzsit(:)
      integer, allocatable :: lstfbp(:),ltpfbp(:)
      integer, allocatable :: lstfre(:),lstrgd(:),numgsit(:)
      integer, allocatable :: lstfrz(:),lstme(:)
      integer, allocatable :: lstgtp(:),lstgst(:,:)
      integer, allocatable :: lstinv(:,:),listinv(:,:),numinv(:)
      integer, allocatable :: lstneu(:),neulst(:),nugrp(:)
      integer, allocatable :: lstout(:),lstbod(:)
      integer, allocatable :: lstshl(:,:),listshl(:,:)
      integer, allocatable :: lsttbp(:),ltptbp(:)
      integer, allocatable :: lstvdw(:),ltpvdw(:),ltpsit(:)
      integer, allocatable :: numang(:),numdih(:),numcon(:)
      integer, allocatable :: numbonds(:),nummols(:),numsit(:)
      integer, allocatable :: numgrp(:),listyp(:)
      integer, allocatable :: numpmf(:),indpmf(:)
      integer, allocatable :: numteth(:),numshl(:)
      integer, allocatable :: itest(:),index(:),kscons(:)
      integer, allocatable :: msite(:),mconst(:)
      integer, allocatable :: key1(:),key2(:),key3(:)
      integer, allocatable :: listttm2(:)
      complex*16, allocatable :: ww1(:), ww2(:), ww3(:)
      complex*16, allocatable :: qqq(:,:,:)
      complex*16, allocatable :: bscx(:), bscy(:), bscz(:)
#ifdef FFTW
      FFTW_PLAN_TYPE fplan, bplan 
#else
      integer fplan, bplan
#endif
      data memc/mxmc*0/,memi/mxmi*0/,memr/mxmr*0/,ifail/0/
      data npage,lines/8,0/,newjob/.true./,safe/.true./safeq/.true./
      data safep/.true./

c     ASPC GRM
      integer :: aspc_k_param
      integer :: aspc_iter_nstep
      real*8 :: sor_omega
      logical :: laspc

#ifdef VAMPIR
#include "VT.inc"
#endif
c
c     set up the communications

      call mb_init(1)
!      call gsync()
#ifdef VAMPIR
c
c     set up VAMPIR

      call VTSETUP()
      call VTTRACEON(ierr)
      call VTBEGIN(99, ierr)
#endif
c     
c     set idnode/mxnode used across the code

      idnode = bead_rank
      mxnode = bead_size
!
! determine lhead
      if( idnode == 0 ) then
         lhead=.true.
      else
         lhead=.false.
      endif
c     
c     open main printing file
      
      if(idnode.eq.0)open(nrite,file='OUTPUT')
      if(idnode.eq.0) write (nrite,
     x  "(/,20x,'DL_POLY Version 2.0',
     x  /,/,30x,'Running on ',i4,' nodes',/,/)") mxnode
c
c     define all major array sizes

#ifdef SHMEM
      call shpalloc(tbf_ptr,10,memr(1),ifail)
#elif SGISHMEM
      call shpalloc(tbf_ptr,20,memr(1),ifail)
#else
      allocate (tbuffer(10),stat=memr(1))
#endif
      call parset(idnode,mxnode,tbuffer)
c
c     allocate arrays

      allocate (lms(mxneut),stat=meml)
      if(meml.gt.0)call error(idnode,35)

      allocate (molnam(mxtmls),stat=memc(1))
      allocate (atmnam(mxatms),stat=memc(2))
      allocate (sitnam(mxsite),stat=memc(3))
      allocate (unqatm(mxsite),stat=memc(4))
c
c     check character memory allocation

      do i=1,mxmc

        if(memc(i).ne.0)safe=.false.

      enddo
      if(.not.safe)then

        if(idnode.eq.0)write(nrite,'(10i5)')memc
        call error(idnode,34)

      endif

#ifdef SHMEM
      call shpalloc(buf_ptr,mxbuff,memr(2),ifail)
#elif SGISHMEM
      call shpalloc(buf_ptr,mxbuff*2,memr(2),ifail)
#else
      allocate (buffer(mxbuff),stat=memr(2))
#endif
      allocate (accum(mxungp),rotmin(mxungp),gaxs(mxungp,9),
     x         stat=memr(3))
      allocate (amsd(mxsvdw),dens0(mxsvdw),esig1(mxcons),stat=memr(4))
      allocate (chge(mxatms),rcut3b(mxtbp),ewlbuf(mxebuf),stat=memr(5))
      allocate (ckc(mxewld),cks(mxewld),clm(mxewld),slm(mxewld),
     x         stat=memr(6))
      allocate (dxp(mspmf),dyp(mspmf),dzp(mspmf),stat=memr(7))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=memr(8))
      allocate (dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=memr(9))
      allocate (elc(mxewld,0:1),els(mxewld,0:1),stat=memr(10))
      allocate (emc(mxewld,0:kmaxb),ems(mxewld,0:kmaxb),stat=memr(11))
      allocate (enc(mxewld,0:kmaxc),ens(mxewld,0:kmaxc),stat=memr(12))
      allocate (erc(mxegrd),fer(mxegrd),dens(mxsvdw),stat=memr(13))
      allocate (flx(mxatms),fly(mxatms),flz(mxatms),stat=memr(14))
      allocate (fpx(mxatms),fpy(mxatms),fpz(mxatms),stat=memr(15))
      allocate (fxx(mxatms),fyy(mxatms),fzz(mxatms),stat=memr(16))
      allocate (gcmx(mxgrp),gcmy(mxgrp),gcmz(mxgrp),stat=memr(17))
      allocate (gcmx1(msgrp),gcmy1(msgrp),gcmz1(msgrp),stat=memr(18))
      allocate (gvx1(msgrp),gvy1(msgrp),gvz1(msgrp),stat=memr(19))
      allocate (gvxx(mxgrp),gvyy(mxgrp),gvzz(mxgrp),gmass(mxungp),
     x         stat=memr(20))
      allocate (gxx(mxungp,mxngp),gyy(mxungp,mxngp),gzz(mxungp,mxngp),
     x         stat=memr(21))
      allocate (omx(mxgrp),omy(mxgrp),omz(mxgrp),stat=memr(22))
      allocate (omx1(msgrp),omy1(msgrp),omz1(msgrp),stat=memr(23))
      allocate (opx(msgrp),opy(msgrp),opz(msgrp),stat=memr(24))
      allocate (oqx(msgrp),oqy(msgrp),oqz(msgrp),stat=memr(25))
      allocate (pmfwght(mxspmf),redmass(mxcons),dsq(mspmf),
     x         stat=memr(26))
      allocate (prmang(mxtang,mxpang),prmfld(mxfld),stat=memr(27))
      allocate (prmcon(mxtcon),prmbnd(mxtbnd,mxpbnd),rcut4b(mxfbp),
     x         stat=memr(28))
      allocate (prmdih(mxtdih,mxpdih),prmvdw(mxvdw,mxpvdw),
     x         stat=memr(29))
      allocate (prmtet(mxteth,mxpbnd),prmshl(mxtshl),stat=memr(30))
      allocate (q0(mxgrp),q1(mxgrp),q2(mxgrp),q3(mxgrp),stat=memr(31))
      allocate (qn0(msgrp),qn1(msgrp),qn2(msgrp),qn3(msgrp),
     x         stat=memr(32))
      allocate (rdf(mxrdf,mxvdw),prminv(mxtinv,mxpinv),prmtbp(mxtbp,
     x          mxptbp),stat=memr(33))
      allocate (rotinx(mxungp,2),rotiny(mxungp,2),rotinz(mxungp,2),
     x         stat=memr(34))
      allocate (ssqval(mxnstk),stkval(mxstak,mxnstk),stat=memr(35))
      allocate (stpval(mxnstk),sumval(mxnstk),stat=memr(36))
      allocate (tqx(msgrp),tqy(msgrp),tqz(msgrp),stat=memr(37))
      allocate (xold(msatms),yold(msatms),zold(msatms),stat=memr(38))
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=memr(39))
      allocate (uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=memr(40))
      allocate (vvv(mxgrid,mxvdw),ggg(mxgrid,mxvdw),stat=memr(41))
      allocate (vx1(msatms),vy1(msatms),vz1(msatms),stat=memr(42))
      allocate (vxx(mxatms),vyy(mxatms),vzz(mxatms),stat=memr(43))
      allocate (weight(mxatms),chgsit(mxsite),wgtsit(mxsite),
     x         stat=memr(44))
      allocate (polarsit(mxsite),polarsit2(mxsite),
     x          polr(mxatms),polr2(mxatms),potcc(mxatms),
     x         stat=memr(70))
      allocate (dipx(mxatms),dipy(mxatms),dipz(mxatms),stat=memr(71))
      allocate (dipx_old(mxatms),dipy_old(mxatms),dipz_old(mxatms)) ! for lcp
      allocate (vdxx(mxatms),vdyy(mxatms),vdzz(mxatms),stat=memr(78))
      allocate (udxx(mxatms),udyy(mxatms),udzz(mxatms),stat=memr(79))
      allocate (fdxx(mxatms),fdyy(mxatms),fdzz(mxatms),stat=memr(80))
      allocate (dpms(mxatms),stat=memr(81)) 
      allocate (emux(mxatms),emuy(mxatms),emuz(mxatms),stat=memr(72))
      allocate (efieldkx(mxatms),efieldky(mxatms),efieldkz(mxatms),
     x         stat=memr(75))
      allocate (efdcrecx(mxatms),efdcrecy(mxatms),efdcrecz(mxatms),
     x         stat=memr(73))
      allocate (efddmurecx(mxatms),efddmurecy(mxatms),
     x          efddmurecz(mxatms),stat=memr(74))
      allocate (ckr(mxewld),skr(mxewld),stat=memr(76))
      allocate (p(6*mxatms),xi(6*mxatms),g(6*mxatms),stat=memr(85))
      allocate (h(6*mxatms),pcom(6*mxatms),xicom(6*mxatms),
     x          stat=memr(86))
      allocate (xt(6*mxatms),ooo(22,6*mxatms),stat=memr(87))
      allocate (ercp(mxegrd,0:3),stat=memr(77))
      allocate (xa(2,mspmf),ya(2,mspmf),za(2,mspmf),stat=memr(45))
      allocate (xdab(msbad),ydab(msbad),zdab(msbad),stat=memr(46))
      allocate (xdab2(mxatms),ydab2(mxatms),zdab2(mxatms),stat=memr(84))
      allocate (xdac(msbad),ydac(msbad),zdac(msbad),stat=memr(47))
      allocate (xdad(msbad),ydad(msbad),zdad(msbad),stat=memr(48))
      allocate (xdbc(msbad),ydbc(msbad),zdbc(msbad),stat=memr(49))
      allocate (xdcd(msbad),ydcd(msbad),zdcd(msbad),stat=memr(50))
      allocate (xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf),
     x         stat=memr(51))
!      allocate (xx0(mxatms),yy0(mxatms),zz0(mxatms),stat=memr(52))
      allocate (xx0(mxatms),yy0(mxatms),zz0(mxatms),
     $     xxold(mxatms),yyold(mxatms),zzold(mxatms), stat=memr(52))
      allocate (xxs(mxatms),yys(mxatms),zzs(mxatms),stat=memr(53))
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=memr(54))
      allocate (xxx(mxatms),yyy(mxatms),zzz(mxatms),stat=memr(55))
      allocate (xxx1(mxatms),yyy1(mxatms),zzz1(mxatms),stat=memr(81))
      allocate (xcm(mxatms),ycm(mxatms),zcm(mxatms),stat=memr(82))
      allocate (vxcm(mxatms),vycm(mxatms),vzcm(mxatms),stat=memr(83))
      allocate (xmsd(mxatms),ymsd(mxatms),zmsd(mxatms),stat=memr(88))
      allocate (zdens(mxrdf,mxsvdw),rho(mxatms),prmfbp(mxfbp,mxpfbp),
     x          stat=memr(56))
      allocate (zumval(mxnstk),ravval(mxnstk),csp(mxspl),stat=memr(57))
      allocate (bspx(mxspme,mxspl),bspy(mxspme,mxspl),stat=memr(58))
      allocate (bspz(mxspme,mxspl),bsdx(mxspme,mxspl),stat=memr(59))
      allocate (bsdy(mxspme,mxspl),bsdz(mxspme,mxspl),stat=memr(60))
      allocate (ww1(kmaxd),ww2(kmaxe),ww3(kmaxf),stat=memr(61))
      allocate (qqc(kmaxd,kmaxe,kmaxf),qqq(kmaxd,kmaxe,kmaxf),
     x          stat=memr(62))
      allocate (bscx(kmaxd),bscy(kmaxe),bscz(kmaxf),stat=memr(63))
      allocate (ffttable(mxftab),ahk(0:mxhko),stat=memr(64))
      allocate (hon(mxegrd,0:mxhko),dhn(mxegrd,0:mxhko),stat=memr(65))
      allocate (fon(mxegrd,0:7),pp(2*mxhko),stat=memr(66))
      allocate (crn(0:mxhko,0:mxhko),stat=memr(67))
      allocate (znp(mxhke,0:2*mxhko),zgc(0:2*mxhko),zgs(0:2*mxhko),
     x          stat=memr(68))
      allocate (zzn(mxxdf),zzd(mxxdf),sss(mxxdf),stat=memr(69))

      nauxfft = (/ 3,0,0,0 /)

c
c     check real memory allocation

      do i=1,mxmr

        if(memr(i).ne.0)safe=.false.

      enddo
      if(.not.safe)then

        if(idnode.eq.0)write(nrite,'(10i5)')memr
        call error(idnode,33)

      endif

      allocate (ilist(mxxdf),jlist(mxxdf),ind(mxgrp,3),stat=memi(1))
      allocate (keybnd(mxtbnd),keyang(mxtang),keydih(mxtdih),
     x         stat=memi(2))
      allocate (lentry(msatms),list(msatms,mxlist),ltype(mxatms),
     x         stat=memi(3))
      allocate (lexatm(msatms,mxexcl),nexatm(msatms),noxatm(msatms),
     x         stat=memi(4))
      allocate (lexatm2(msatms,mxexcl),nexatm2(msatms),noxatm2(msatms),
     x         stat=memi(34))
      allocate (lexsit(mxsite,mxexcl),nexsit(mxsite),stat=memi(5))
      allocate (link(mxatms),lct(mxcell),lst(mxcell),stat=memi(6))
      allocate (lishap(mxlshp),listot(mxatms),numtyp(mxsvdw),
     x         stat=memi(7))
      allocate (listme(mxatms),lashap(mxproc),listin(mxatms),
     x         stat=memi(8))
      allocate (listpm(mxpmf),lstpmt(mxpmf),numfrz(mxsvdw),
     x         stat=memi(9))
      allocate (listtet(msteth,2),keytet(mxteth),lsttet(mxteth),
     x          stat=memi(10))
      allocate (lstang(mxtang,3),listang(mxangl,4),stat=memi(11))
      allocate (lstbnd(mxtbnd,2),listbnd(mxbond,3),keyinv(mxtinv),
     x         stat=memi(12))
      allocate (listin2(4*mxatms),listbnd2(mxatms),listbnd3(mxatms),
     x         stat=memi(35))
      allocate (listttm2(mxatms),stat=memi(36))
      allocate (lstcon(mxtcon,2),listcon(mxcons,3),stat=memi(13))
      allocate (lstcsit(2*mxcons),lstpmf(mxspmf,mspmf),stat=memi(14))
      allocate (lstdih(mxtdih,4),listdih(mxdihd,5),lfzsit(mxsite),
     x         stat=memi(15))
      allocate (lstfbp(mxfbp),ltpfbp(mxfbp),stat=memi(16))
      allocate (lstfre(mxatms),lstrgd(mxgatm),numgsit(mxungp),
     x         stat=memi(17))
      allocate (lstfrz(mxatms),lstme(mxatms),stat=memi(18))
      allocate (lstgtp(mxgrp),lstgst(mxungp,mxngp),stat=memi(19))
      allocate (lstinv(mxtinv,4),listinv(mxinv,5),numinv(mxtmls),
     x         stat=memi(20))
      allocate (lstneu(mxatms),neulst(mxneut),nugrp(mxsite),
     x         stat=memi(21))
      allocate (lstout(mxatms),lstbod(mxatms),stat=memi(22))
      allocate (lstshl(mxtshl,2),listshl(mxshl,3),stat=memi(23))
      allocate (lsttbp(mxtbp),ltptbp(mxtbp),stat=memi(24))
      allocate (lstvdw(mxvdw),ltpvdw(mxvdw),ltpsit(mxsite),
     x         stat=memi(25))
      allocate (numang(mxtmls),numdih(mxtmls),numcon(mxtmls),
     x         stat=memi(26))
      allocate (numbonds(mxtmls),nummols(mxtmls),numsit(mxtmls),
     x          stat=memi(27))
      allocate (numgrp(mxtmls),listyp(mxungp),stat=memi(28))
      allocate (numpmf(mxtmls),indpmf(mxspmf),stat=memi(29))
      allocate (numteth(mxtmls),numshl(mxtmls),stat=memi(30))
      allocate (itest(mxtmls),index(mxtmls),kscons(0:mxproc-1),
     x         stat=memi(31))
      allocate (msite(mxtmls),mconst(mxtmls),stat=memi(32))
      allocate (key1(kmaxd),key2(kmaxe),key3(kmaxf),stat=memi(33))
c
c     check integer memory allocation

      do i=1,mxmi

        if(memi(i).ne.0)safe=.false.

      enddo
      if(.not.safe)then

        if(idnode.eq.0)write(nrite,'(10i5)')memi
        call error(idnode,32)

      endif
c     
c     start clock

      call timchk(0,tzer0)
c     
c     input the control parameters defining the simulation

      call simdef
     x  (lfcap,lgofr,lnsq,loptim,lpgr,ltraj,ltscal,lzeql,lzden,lpolar,
     x  lthole,sysname,idnode,mxnode,intsta,istraj,keyens,keyfce,lopt2,
     x  lfd,keyres,keytrj,kmax1,kmax2,kmax3,multt,nstack,nstbgr,nstbpo,
     x  nhko,nlatt,nstbts,nsteql,nstraj,nstrun,nospl,fplan,bplan,ftol,
     x  alpha,delr,epsq,fmax,press,quattol,rcut,rprim,rvdw,taup,taut,
     x  athole,athole12,athole13,temp,timcls,timjob,tolnce,tstep,
     x  ffttable,
     x  lcp,ldpts,lads,lacs,lttm,lttm3,nthole,dipmas,diptmp,tautd,toler,
     x  npartem,temgap,
     x  delta,ascc,ascd,gammattm,laspc,sor_omega,aspc_k_param,
     x  aspc_iter_nstep)
c     
c     input the system force field

      call sysdef
     x  (lneut,lmetal,lnsq,molnam,sitnam,unqatm,idnode,mxnode,keyfce,
     x  keyfld,natms,ngrp,ntpatm,ntpmls,ntpvdw,ntptbp,ntpfbp,nshels,
     x  nhko,nlatt,alpha,dlrpot,drewd,engunit,prmpmf,rcut,rvdw,rcuttb,
     x  rcutfb,indpmf,keyang,keybnd,keydih,keyinv,keytet,lfzsit,listyp,
     x  lstang,lstbnd,lstcon,lstdih,lstinv,lstgst,lstvdw,lsttbp,lstfbp,
     x  lstshl,ltpsit,lsttet,ltpvdw,ltptbp,ltpfbp,npmf,nugrp,numang,
     x  numbonds,numcon,numdih,numinv,numgrp,numgsit,nummols,numsit,
     x  numpmf,numteth,numshl,chgsit,lpolar,polarsit,polarsit2,
     x  erc,fer,ercp,pmfwght,ggg,prmang,prmbnd,
     x  prmcon,prmdih,prminv,prmfld,prmtet,prmvdw,prmtbp,prmfbp,prmshl,
     x  vvv,wgtsit,rcut3b,rcut4b,buffer,ahk,hon,dhn,fon)

      if(lmetal.and.multt.gt.1)call error(idnode,153)
c     
c     construct initial configuration of system

      call sysgen
     x  (lhead,loglnk,lneut,lopt2,atmnam,cfgname,sitnam,idnode,imcon,
     x  keyens,keyfce,keyres,levcfg,multt,mxnode,ntpmls,delr,rcut,
     x  volm,lfzsit,lstfrz,ltpsit,ltype,lpolar,nummols,numsit,nugrp,
     x  lstneu,buffer,cell,chge,chgsit,polr,polr2,polarsit,polarsit2,
     x  fxx,fyy,fzz,vxx,vyy,vzz,
     x  weight,wgtsit,xxx,yyy,zzz)
c
c     construct initial bookkeeping arrays

      call sysbook
     x  (loglnk,lneut,lshmov,lcnb,idnode,imcon,mxnode,natms,
     x  nneut,ngrp,nscons,ntangl,ntbond,ntcons,ntdihd,ntinv,
     x  ntpmls,ntpmf,nspmf,ntfree,ntteth,ntshl,degfre,degrot,
     x  keybnd,keyang,lashap,lexsit,lexatm,lexatm2,lishap,listang,
     x  listbnd,listcon,listdih,listinv,listin,listme,listot,
     x  listyp,lstang,lstbnd,lstcon,lstdih,lstinv,lstfre,
     x  lstfrz,lstgtp,lstgst,lstrgd,lstcsit,nexatm,nexatm2,nexsit,
     x  numang,lstme,lstbod,lstout,lstshl,numbonds,numcon,
     x  numdih,numinv,numgsit,indpmf,lstpmf,numpmf,npmf,ind,
     x  lstpmt,listpm,numgrp,nummols,numsit,numteth,numshl,
     x  itest,index,kscons,msite,mconst,lsttet,listtet,
     x  listshl,neulst,lstneu,buffer,cell,gcmx,gcmy,gcmz,
     x  gmass,gxx,gyy,gzz,prmdih,q0,q1,q2,q3,rotinx,rotiny,
     x  rotinz,txx,tyy,tzz,weight,xxt,xxx,yyt,yyy,zzt,zzz,
     x  accum,gaxs,rotmin)

!
!     check whether polarization is turned off in the case of TTM2 model

      if( lpoloff ) then

         if( .not.lttm ) then
            
            write(nrite,*) 'POLOFF is only valid for TTM2 model'
            stop
            
         else if( loglnk ) then
            
            write(nrite,*) 'POLOFF does not support link cell algorithm'
            stop
            
         else if( lcp ) then
            
            write(nrite,*) 'POLOFF does not support CP'
            stop

         endif
         
      endif

      nfict=0
      if (lttm) then

        call ttm2list
     x  (nttm2,ntpmls,nummols,numsit,listttm2,weight)

        call exclude_ttm
     x  (idnode,mxnode,natms,nexatm,lexatm,nttm2,listttm2)

c
c     save real and fictitious degree of freedom

        natms_f=natms
        natms_r=natms-nummols(ntpmls)
        nfict=nummols(ntpmls)
c
c     real degree of freedom

        natms=natms_r
! 
!     for good guess
        if( .not. lcp ) then
           allocate( dip1x(mxatms), dip1y(mxatms), dip1z(mxatms) )
           allocate( dip2x(mxatms), dip2y(mxatms), dip2z(mxatms) )
           allocate( dip3x(mxatms), dip3y(mxatms), dip3z(mxatms) )
        end if
!
!     for force decomposition

        call init_ttm_forces(natms_f)
!
!     molecular dipole moments

        n_water = nummols(ntpmls)
        call init_dipole_moments(n_water)
!
!     gradients with respect to intra. coordinates 

        if( ldms ) call init_intra_grad_dms(n_water)

      endif ! lttm
c
c     set initial system temperature

      if(.not.lopt2) call systemp
     x  (idnode,imcon,keyens,keyres,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,ntshl,levcfg,degfre,degshl,degrot,temp,
     x  tolnce,lashap,lishap,listcon,listme,listot,lstfrz,lstgtp,
     x  lstrgd,lstme,numgsit,listshl,buffer,cell,dxt,dyt,dzt,
     x  fxx,fyy,fzz,gcmx,gcmy,gcmz,gmass,gvxx,gvyy,gvzz,
     x  gxx,gyy,gzz,q0,q1,q2,q3,rotinx,rotiny,rotinz,weight,
     x  uxx,uyy,uzz,vxx,vyy,vzz,xxt,xxx,yyt,yyy,zzt,zzz,
     x  omx,omy,omz,lttm,nfict)
c
c     read thermodynamic and structural data from restart file

      call sysinit
     x  (lgofr,lzden,lmetal,idnode,imcon,keyfce,keyres,mxnode,
     x  natms,nstep,numacc,numrdf,ntpatm,nzden,chip,chit,chitd,
     x  conint,elrc,engunit,virlrc,rvdw,volm,lstvdw,ltpvdw,lstfrz,
     x  ltype,numtyp,numfrz,buffer,cell,dens,prmvdw,ravval,rdf,
     x  ssqval,stkval,stpval,sumval,xx0,yy0,zz0,zumval,zdens,
     x  xxs,yys,zzs,elrcm,vlrcm,eta,dipx,dipy,dipz,
     x  vdxx,vdyy,vdzz,xmsd,ymsd,zmsd,lpolar,lcp,ldpts,conintd)

      if(laspc) then
        !check for errors in CONTROL file
        if(laspc.and.lcp) call error(idnode, 817)
        if(laspc.and.(.not.lpolar)) call error(idnode, 818)

        call aspc_init(aspc_k_param, aspc_iter_nstep)
      endif

      call timchk(1,tzero)
!
!    for optimization

      if( loptim ) then
         
         call input_opt
         
         if( lhead .and. lgeoopt ) call output_opt
         
      end if

c***********************************************************************
c     start of molecular dynamics calculations
c***********************************************************************
!
!     for check COM velocity
      flag_com = .true.       
      if( nstep > 0 ) flag_com = .false.

c     zero long range component of stress
      do i = 1,9
        stresl(i) = 0.d0
      enddo
c
c     geometry optimization

      if (lopt2) then

        call engforce
     x (idnode,imcon,mxnode,natms,ntpatm,engunit,rvdw,volm,
     x  lstvdw,ltpvdw,ltype,numtyp,prmvdw,dens,elrcm,vlrcm,
     x  lmetal,keyens,cell,celprp,newlst,delr,tstep,vxx,vyy,
     x  vzz,xold,yold,zold,engcp3,epsq,rcut,vircp3,lexatm,
     x  nexatm,noxatm,lentry,list,lstfrz,xxx,yyy,zzz,xdf,ydf,
     x  zdf,flx,fly,flz,chge,stresl,lct,link,uxx,uyy,uzz,
     x  buffer,lneut,lms,nneut,neulst,lstneu,lnsq,lgofr,lzeql,
     x  keyfce,multt,nstep,nstbgr,nsteql,numrdf,dlrpot,engcpe,
     x  engsrp,rprim,ilist,fpx,fpy,fpz,ggg,
     x  rdf,rsqdf,vvv,loglnk,kmax1,kmax2,kmax3,nhko,nlatt,
     x  ntpvdw,nospl,fplan,bplan,alpha,drewd,jlist,nauxfft,
     x  lexatm2,nexatm2,key1,key2,key3,polr,polr2,rho,ewlbuf,toler,
     x  ckc,cks,clm,elc,els,emc,ckr,skr,ercp,ntptbp,rcut3b,
     x  ems,enc,ens,erc,fer,fxx,fyy,fzz,slm,ntbond,ntpfbp,
     x  csp,qqc,txx,tyy,tzz,bspx,bspy,bspz,bsdx,bsdy,bsdz,qqq,
     x  bscx,bscy,bscz,ffttable,ww1,ww2,ww3,ahk,zzn,zzd,sss,hon,
     x  dhn,pp,znp,zgc,zgs,crn,lcp,lstout,xxt,yyt,zzt,
     x  lpolar,lthole,athole,dipx,dipy,dipz,emux,emuy,emuz,
     x  efieldkx,efieldky,efieldkz,efdcrecx,efdcrecy,efdcrecz,
     x  efddmurecx,efddmurecy,efddmurecz,potcc,stress,
     x  rcuttb,engtbp,virtbp,listin,lst,lsttbp,ltptbp,prmtbp,
     x  rcutfb,engfbp,virfbp,lstfbp,ltpfbp,prmfbp,rcut4b,
     x  engbnd,virbnd,keybnd,listbnd,prmbnd,xdab,ydab,zdab,
     x  ntangl,engang,virang,keyang,listang,prmang,xdbc,
     x  ydbc,zdbc,ntdihd,engdih,virdih,keydih,listdih,
     x  prmdih,xdcd,ydcd,zdcd,ntshl,engshl,virshl,listshl,
     x  newjob,ntinv,enginv,virinv,keyinv,listinv,prminv,prmshl,
     x  ntteth,keytet,listtet,engtet,virtet,prmtet,xxs,yys,zzs,
     x  keyfld,engfld,virfld,prmfld,weight,elrc,virlrc,
     x  ntpmls,nummols,natms_r,natms_f,lttm,nttm2,listttm2,
     x  nthole,lacs,lads,ascc,ascd,keyres,
     x  uuu,gammattm)

        if (lpolar) then
        write(*,*)'converged dipoles:'
        do ii=1,natms
              write(*,'(I6,3X,3f28.20)')
     x           ii,dipx(ii)*eatd,dipy(ii)*eatd,dipz(ii)*eatd
           enddo
        endif

        call geopt_pt
     x (idnode,imcon,mxnode,natms,ntpatm,engunit,rvdw,volm,
     x  lstvdw,ltpvdw,ltype,numtyp,prmvdw,dens,elrcm,vlrcm,
     x  lmetal,keyens,cell,celprp,newlst,delr,tstep,vxx,vyy,
     x  vzz,xold,yold,zold,engcp3,epsq,rcut,vircp3,lexatm,
     x  nexatm,noxatm,lentry,list,lstfrz,xxx,yyy,zzz,xdf,ydf,
     x  zdf,flx,fly,flz,chge,stresl,lct,link,uxx,uyy,uzz,
     x  buffer,lneut,lms,nneut,neulst,lstneu,lnsq,lgofr,lzeql,
     x  keyfce,multt,nstep,nstbgr,nsteql,numrdf,dlrpot,engcpe,
     x  engsrp,rprim,vircpe,virsrp,ilist,fpx,fpy,fpz,ggg,
     x  rdf,rsqdf,vvv,loglnk,kmax1,kmax2,kmax3,nhko,nlatt,
     x  ntpvdw,nospl,fplan,bplan,alpha,drewd,jlist,nauxfft,
     x  lexatm2,nexatm2,key1,key2,key3,polr,polr2,rho,ewlbuf,toler,
     x  ckc,cks,clm,elc,els,emc,ckr,skr,ercp,ntptbp,rcut3b,
     x  ems,enc,ens,erc,fer,fxx,fyy,fzz,slm,ntbond,ntpfbp,
     x  csp,qqc,txx,tyy,tzz,bspx,bspy,bspz,bsdx,bsdy,bsdz,qqq,
     x  bscx,bscy,bscz,ffttable,ww1,ww2,ww3,ahk,zzn,zzd,sss,hon,
     x  dhn,pp,znp,zgc,zgs,crn,lcp,lstout,xxt,yyt,zzt,
     x  lpolar,lthole,athole,dipx,dipy,dipz,emux,emuy,emuz,
     x  efieldkx,efieldky,efieldkz,efdcrecx,efdcrecy,efdcrecz,
     x  efddmurecx,efddmurecy,efddmurecz,potcc,stress,
     x  rcuttb,engtbp,virtbp,listin,lst,lsttbp,ltptbp,prmtbp,
     x  rcutfb,engfbp,virfbp,lstfbp,ltpfbp,prmfbp,rcut4b,
     x  engbnd,virbnd,keybnd,listbnd,prmbnd,xdab,ydab,zdab,
     x  ntangl,engang,virang,keyang,listang,prmang,xdbc,temp,
     x  ydbc,zdbc,ntdihd,engdih,virdih,keydih,listdih,
     x  prmdih,xdcd,ydcd,zdcd,ntshl,engshl,virshl,listshl,
     x  newjob,ntinv,enginv,virinv,keyinv,listinv,prminv,prmshl,
     x  ntteth,keytet,listtet,engtet,virtet,prmtet,xxs,yys,zzs,
     x  keyfld,engfld,virfld,prmfld,weight,elrc,virlrc,nstrun,
     x  ntpmls,nummols,numsit,wgtsit,numbonds,nums,
     x  listbnd2,listbnd3,listin2,xxx1,yyy1,zzz1,keyres,
     x  xdab2,ydab2,zdab2,xcm,ycm,zcm,vxcm,vycm,vzcm,ascc,ascd,
     x  natms_r,natms_f,lttm,nttm2,listttm2,nthole,lacs,lads,
     x  ooo,p,xi,g,h,pcom,xicom,xt,ftol,npartem,temgap)

c        lcp = .false.

      endif
c
c     test numerical derivatives by finite difference

      if (lfd) then
        icount=0
        ilfd=0
        allocate( fanaly(3,mxatms) )
        fanaly(:,:) = 0.0d0
      endif

! counter for dipole guess
      if( lttm .and. (.not. lcp) ) nstep_for_guess = 0

  100 continue

!
! check system COM velocity and let it stop

      call checkcom(vxx, vyy,vzz,weight,natms,flag_com,nstep)
      
      flag_com = .false.
!
!     M sites are explicitly treated in the case of lfd+lttm

      if( lfd .and. lttm ) then
         natms=natms_f
         n_water = nummols(ntpmls)
         call qdetermine(listttm2,n_water,nttm2,lttm3,imcon,cell,
     $        chge,gammattm,xxx,yyy,zzz)
         goto 900  
      end if
c
c     for ttm2 water model, use fictitious degree of freedom

      if(lttm) then

        ntmp = 1
        call init_ttm_forces(ntmp)
  
        natms=natms_f
c
c     assign coordinates to the M sites jcp 116(2002)5115 Eq. (A5)

        mttm2=listttm2(nttm2)
        do i=1,nummols(ntpmls)
          mttm2=mttm2+1
          ioxy=listttm2(3*i-2)
          ih1=listttm2(3*i-1)
          ih2=listttm2(3*i)
          xdf(1)=xxx(ih1)-xxx(ioxy)
          ydf(1)=yyy(ih1)-yyy(ioxy)
          zdf(1)=zzz(ih1)-zzz(ioxy)
          xdf(2)=xxx(ih2)-xxx(ioxy)
          ydf(2)=yyy(ih2)-yyy(ioxy)
          zdf(2)=zzz(ih2)-zzz(ioxy)
          call images(imcon,0,1,2,cell,xdf,ydf,zdf)
          xxx(mttm2)=xxx(ioxy)+(xdf(1)+xdf(2))*gammattm/2.d0
          yyy(mttm2)=yyy(ioxy)+(ydf(1)+ydf(2))*gammattm/2.d0
          zzz(mttm2)=zzz(ioxy)+(zdf(1)+zdf(2))*gammattm/2.d0
        enddo
!
!     assign partial charges
        n_water = nummols(ntpmls)
        call qdetermine(listttm2,n_water,nttm2,lttm3,imcon,cell,
     $       chge,gammattm,xxx,yyy,zzz)
        
      endif
!
!     lfd + lttm
 900  continue
c     
c     initialize energy and virial accumulators for bonds,angles
c     dihedrals, inversions, short range potential and electrostatics
c     constraint virial,  com virial and field terms
      
      engbnd = 0.d0
      virbnd = 0.d0
      engang = 0.d0
      virang = 0.d0
      engdih = 0.d0
      virdih = 0.d0
      enginv = 0.d0
      virinv = 0.d0
      engtbp = 0.d0
      virtbp = 0.d0
      engfbp = 0.d0
      virfbp = 0.d0
      engsrp = 0.d0
      virsrp = 0.d0
      engcpe = 0.d0
      vircpe = 0.d0
      vircon = 0.d0
      vircom = 0.d0
      engfld = 0.d0
      virfld = 0.d0
      engshl = 0.d0
      virshl = 0.d0
      shlke  = 0.d0
      engtet = 0.d0
      virtet = 0.d0
      virpmf = 0.d0
c     
c     conserved quantity (other than K + U)
      
      consv = 0.d0
      consvd = 0.d0
c     
c     energy accumulators
      
      engke = 0.d0
      engrot = 0.d0
c
c     dipole kinetic energy

      engdke=0.d0
c     
c     zero stress tensor
      
      do i = 1,9
        stress(i) = 0.d0
      enddo
c
c     reset sutton chen long range corrections (constant pressure only)

      if(lmetal) then
        if(keyens.ge.4.and.keyens.le.7) then
          call lrcmetal
     x      (idnode,imcon,mxnode,natms,ntpatm,engunit,rvdw,volm,
     x      lstvdw,ltpvdw,ltype,numtyp,prmvdw,dens,elrcm,vlrcm)

        endif
      endif
c     
c     increase step counter
      
      nstep=nstep+1
c     
c     desired kinetic energy
      
      sigma=temp*boltz*degfre*0.5d0
      if (lpolar.and.lcp.and.ldpts)
     x   sigmad=diptmp*boltz*dble(natms)*1.5d0
c     
c     initialise the force arrays
      
      do i=1,natms
        
        fxx(i)=0.d0
        fyy(i)=0.d0
        fzz(i)=0.d0
        
      enddo
c
c     initialize electric field and dipole force arrays
c     if polarizability is used with car-parrinello

      if (lpolar .and. lcp) then

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

      endif
c     
c     calculate volume of simulation cell
      
      if(imcon.ne.0.and.imcon.ne.6)then
        
        call dcell(cell,celprp)
        volm=celprp(10)
        if(imcon.eq.4)then

          volm=0.5d0*celprp(10)

        elseif(imcon.eq.5)then
        
          volm=0.5d0*celprp(10)

        elseif(imcon.eq.7)then
        
          volm=0.5d0*celprp(10)

        endif

      else
        
        volm=0.d0
        
      endif
c     
c     test for updating of Verlet list
      
      call vertest
     x  (newlst,idnode,mxnode,natms,delr,tstep,vxx,vyy,vzz,
     x  xold,yold,zold)

c     
c     set up nonbonded interaction (verlet) list
      
      newlst=((newjob).or.(newlst))

      if (newlst) then
c     
c     coulombic accumulators for interactions outside rcut
        
        engcp3 = 0.d0
        vircp3 = 0.d0

        if(.not.lneut) then
          
          if(lnsq) then 
c     
c     calculate distant interactions explicitly
            
            call parlst_nsq
     x        (newlst,natms,idnode,mxnode,imcon,engcp3,
     x        epsq,rcut,vircp3,lexatm,nexatm,noxatm,lentry,
     x        list,lstfrz,cell,xxx,yyy,zzz,xdf,ydf,zdf,flx,
     x        fly,flz,chge,stresl)
            
          elseif(loglnk) then
c     
c     ignore real space distant interactions
            
            call parlink
     x        (newlst,natms,idnode,mxnode,imcon,rcut,delr,
     x        lct,link,lexatm,nexatm,lentry,list,lstfrz,
     x        cell,xxx,yyy,zzz,uxx,uyy,uzz,buffer)

            
          else
            
            call parlst
     x        (newlst,natms,idnode,mxnode,imcon,rcut,delr,lexatm,
     x        nexatm,noxatm,lentry,list,lstfrz,cell,xxx,yyy,zzz,
     x        xdf,ydf,zdf)
            
          endif
          
        else
          
          if(.not.loglnk) then 

            call parneulst
     x        (newlst,lneut,lms,nneut,idnode,mxnode,imcon,rcut,delr,
     x        lentry,list,lstfrz,neulst,cell,xxx,yyy,zzz,xdf,ydf,zdf)
          
          else

            call parlinkneu
     x        (newlst,lneut,natms,nneut,idnode,mxnode,imcon,rcut,delr,
     x        lentry,lct,link,lstneu,list,lstfrz,neulst,
     x        cell,xxx,yyy,zzz,uxx,uyy,uzz,buffer)

          endif

        endif

        if (lttm) then
c
c     remove intramolecular m-site interactions for ttm2 water

          ii=0
          do iatm=idnode+1,natms,mxnode

            ii=ii+1

            mttm1=listttm2(1)
            mttm2=listttm2(nttm2)
            ki=1
 1005       continue
            do k=ki,lentry(ii)
              jatm=list(ii,k)
              if (iatm.gt.mttm2 .and. jatm.le.mttm2
     x            .and. jatm.ge.mttm1) then
                idum=iatm-mttm2
                jdum1=listttm2(3*idum-2)
                jdum2=listttm2(3*idum-1)
                jdum3=listttm2(3*idum)
                if (jatm.eq.jdum1 .or. jatm.eq.jdum2
     x              .or. jatm.eq.jdum3) then
                  do l=k,lentry(ii)-1
                    list(ii,l)=list(ii,l+1)
                  enddo
                  lentry(ii)=lentry(ii)-1
                  ki=k
                  goto 1005
                endif
              endif
              if (jatm.gt.mttm2 .and. iatm.le.mttm2
     x            .and. iatm.ge.mttm1) then
                idum=jatm-mttm2
                jdum1=listttm2(3*idum-2)
                jdum2=listttm2(3*idum-1)
                jdum3=listttm2(3*idum)
                if (iatm.eq.jdum1 .or. iatm.eq.jdum2
     x              .or. iatm.eq.jdum3) then
                  do l=k,lentry(ii)-1
                    list(ii,l)=list(ii,l+1)
                  enddo
                  lentry(ii)=lentry(ii)-1
                  ki=k
                  goto 1005
                endif
              endif

            enddo

          enddo

c          ii=0
c          do i=idnode+1,natms,mxnode
c            ii=ii+1
c            write(nrite,*)'i,lentry(ii)',i,lentry(ii)
c            do m=1,lentry(ii)
c              write(nrite,*)i,list(ii,m)
c            enddo
c          enddo

        endif

      endif
c     
c     calculate pair forces, including coulombic forces

      do i = 1,9
        stress(i) = stresl(i)
      enddo

ccc      goto 1003
      
      if(lnsq) then
c     
c     multiple timestep - all-pairs
        
        call multiple_nsq
     x    (lnsq,lgofr,lzeql,newlst,idnode,imcon,keyfce,
     x    multt,mxnode,natms,nstep,nstbgr,nsteql,numrdf,
     x    delr,dlrpot,engcpe,engsrp,engcp3,epsq,rcut,
     x    rprim,rvdw,vircpe,virsrp,vircp3,ilist,lentry,
     x    list,lstvdw,ltpvdw,ltype,buffer,cell,chge,flx,fly,
     x    flz,fxx,fyy,fzz,fpx,fpy,fpz,ggg,rdf,rsqdf,vvv,xdf,
     x    xxx,ydf,yyy,zdf,zzz,stress)
        
      elseif(.not.lneut) then         
c     
c     single timestep
        
        if(multt.eq.1) then
! 
! for good initial guess
!
          if( lttm .and. (.not.lcp) ) then
             
             nstep_for_guess = nstep_for_guess + 1 
             if(laspc)then
                 call aspc_predict(nstep_for_guess,natms,
     x                             dipx,dipy,dipz)
             else
                 call set_guess_dipole(nstep_for_guess,natms,
     x                                 dipx,dipy,dipz)
             endif

          endif

          call forces
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
     x  lpolar,lthole,athole,athole12,athole13,dipx,dipy,dipz,
     x  emux,emuy,emuz,nthole,ascc,ascd,efieldkx,efieldky,efieldkz,
     x  efdcrecx,efdcrecy,efdcrecz,
     x  efddmurecx,efddmurecy,efddmurecz,potcc,keyres,
     x  n_water,gammattm,.false.,laspc,sor_omega,aspc_iter_nstep,
     x  nstep_for_guess)

          if(lpolar .and. laspc .and. (.not.lcp)) then
              call aspc_set(natms, dipx,dipy, dipz)
          endif

        else

          call multiple
     x      (loglnk,lgofr,lzeql,newlst,idnode,imcon,keyfce,nlatt,
     x      kmax1,kmax2,kmax3,nhko,multt,mxnode,natms,nstep,nstbgr,
     x      nsteql,numrdf,nospl,fplan,bplan,alpha,dlrpot,drewd,
     x      engcpe,engsrp,epsq,rcut,rprim,rvdw,vircpe,virsrp,volm,
     x      ilist,jlist,lentry,lexatm,list,lstvdw,ltpvdw,ltype,nexatm,
     x      nauxfft,key1,key2,key3,buffer,cell,chge,ckc,cks,clm,elc,
     x      els,emc,ems,enc,ens,erc,fer,flx,fly,flz,fpx,fpy,fpz,
     x      fxx,fyy,fzz,ggg,rdf,rsqdf,slm,stress,vvv,xdf,xxx,ydf,yyy,
     x      zdf,zzz,ewlbuf,csp,qqc,txx,tyy,tzz,bspx,bspy,bspz,bsdx,
     x      bsdy,bsdz,qqq,bscx,bscy,bscz,ffttable,ww1,ww2,ww3,ahk,
     x      zzn,zzd,sss,hon,dhn,pp,znp,zgc,zgs,crn)

        endif

      elseif(lneut) then
c     
c     neutral groups
        
        if(multt.eq.1) then
          
          call forcesneu
     x      (lgofr,lzeql,idnode,imcon,keyfce,mxnode,natms,
     x      nneut,nstbgr,nstep,nsteql,numrdf,dlrpot,engcpe,
     x      engsrp,epsq,rcut,delr,rvdw,vircpe,virsrp,ilist,jlist,
     x      lentry,lexatm,list,lstfrz,link,lstout,lstvdw,
     x      ltype,nexatm,neulst,buffer,cell,chge,fxx,fyy,
     x      fzz,ggg,rdf,rsqdf,stress,vvv,xxx,yyy,zzz,xdf,
     x      ydf,zdf,txx,tyy,tzz,xxt,yyt,zzt,uxx,uyy,uzz)
          
        else
          
          call multipleneu
     x      (lgofr,lzeql,newlst,idnode,imcon,keyfce,multt,
     x      mxnode,natms,nneut,nstbgr,nstep,nsteql,numrdf,
     x      delr,dlrpot,engcpe,engsrp,epsq,rprim,rcut,rvdw,
     x      vircpe,virsrp,ilist,jlist,lentry,lexatm,list,
     x      lstout,lstvdw,ltype,nexatm,neulst,link,lstfrz,
     x      buffer,cell,chge,fxx,fyy,fzz,flx,fly,flz,fpx,
     x      fpy,fpz,ggg,rdf,rsqdf,vvv,xdf,xxx,ydf,yyy,zdf,
     x      zzz,stress,txx,tyy,tzz,xxt,yyt,zzt,uxx,
     x      uyy,uzz)
          
        endif
        
      endif
c
c     add in long range corrections to energy and pressure

      engsrp = engsrp + elrc + elrcm(1)
      virsrp = virsrp + virlrc + vlrcm(1)
c
c     calculate three body forces

      if (ntptbp.gt.0) call thbfrc
     x  (idnode,mxnode,natms,imcon,rcuttb,engtbp,virtbp,
     x  listin,ltype,lst,lct,link,lsttbp,ltptbp,prmtbp,cell,
     x  xxx,yyy,zzz,fxx,fyy,fzz,rcut3b,stress,buffer)
c     
c     calculate four body forces

      if (ntpfbp.gt.0) call fbpfrc
     x  (idnode,mxnode,natms,imcon,rcutfb,engfbp,virfbp,
     x  listin,ltype,lst,lct,link,lstfbp,ltpfbp,prmfbp,cell,
     x  xxx,yyy,zzz,fxx,fyy,fzz,rcut4b,stress,buffer)
c     
c     calculate bond forces

      if (ntbond.gt.0) call bndfrc
     x  (idnode,imcon,mxnode,ntbond,engbnd,virbnd,keybnd,listbnd,
     x  cell,fxx,fyy,fzz,prmbnd,xxx,yyy,zzz,xdab,ydab,zdab,stress,
     x  buffer,lttm)
c
c     calculate valence angle forces
      
      if (ntangl.gt.0) call angfrc
     x  (idnode,imcon,mxnode,ntangl,engang,virang,keyang,listang,
     x  cell,fxx,fyy,fzz,prmang,xxx,yyy,zzz,xdab,ydab,zdab,xdbc,
     x  ydbc,zdbc,stress,buffer,lttm)
!   
!     intramolecular potential for TTM2 model
!      if( lttm .and. ntangl.gt.0 ) then 
!         open(nintra,file='UINTRA',position='append')
!         engintra = engbnd + engang
!         write(nintra,'(i7,f20.10)') nstep, engintra / engunit
!         close(nintra)
!      endif

c     
c     calculate dihedral forces
      
      if (ntdihd.gt.0) call dihfrc
     x  (idnode,imcon,mxnode,ntdihd,keyfce,dlrpot,epsq,engcpe,
     x  engdih,engsrp,rcut,rvdw,vircpe,virdih,virsrp,keydih,listdih,
     x  ltype,lstvdw,buffer,cell,chge,fxx,fyy,fzz,prmdih,xxx,yyy,
     x  zzz,xdab,ydab,zdab,xdbc,ydbc,zdbc,xdcd,ydcd,zdcd,vvv,ggg,
     x  stress)
c     
c     calculate inversion forces
      
      if (ntinv.gt.0) call invfrc
     x  (idnode,imcon,mxnode,ntinv,enginv,virinv,keyinv,listinv,
     x  buffer,cell,fxx,fyy,fzz,prminv,xxx,yyy,zzz,xdab,ydab,zdab,
     x  xdbc,ydbc,zdbc,xdcd,ydcd,zdcd,stress)

c     
c     calculate tethered atom forces
      
      if(ntteth.gt.0) call tethfrc
     x  (idnode,mxnode,imcon,natms,nstep,ntteth,keytet,listtet,
     x  engtet,virtet,buffer,cell,fxx,fyy,fzz,prmtet,xxx,yyy,
     x  zzz,xxs,yys,zzs,xdab,ydab,zdab,stress)
c
c     calculate shell model forces
      
      if (ntshl.gt.0) call shlfrc
     x   (idnode,imcon,mxnode,ntshl,engshl,virshl,listshl,cell,
     x   fxx,fyy,fzz,prmshl,xxx,yyy,zzz,xdab,ydab,zdab,stress,
     x   buffer)
c     
c     external field
      
      if(keyfld.gt.0) call extnfld
     x  (idnode,imcon,keyfld,mxnode,natms,engfld,virfld,cell,
     x  chge,prmfld,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz,weight)
c     
c     global summation of force arrays (basic replicated data strategy)

      if(mxnode.gt.1) then
#ifdef VAMPIR
      call VTBEGIN(28, ierr)
#endif
        j=0
        do i=1,natms

          buffer(j+1)=fxx(i)
          buffer(j+2)=fyy(i)
          buffer(j+3)=fzz(i)
          j=j+3

        enddo

        call gdsum(buffer(1),3*natms,buffer(3*natms+1))

        j=0
        do i=1,natms
        
          fxx(i)=buffer(j+1)
          fyy(i)=buffer(j+2)
          fzz(i)=buffer(j+3)
          j=j+3

        enddo
#ifdef VAMPIR
      call VTEND(28, ierr)
#endif
      endif

#ifdef STRESS      
      if(mxnode.gt.1) call gdsum(stress,9,buffer)
c
c     add long range correction to diagonal terms of stress tensor

      stress(1) = stress(1) - (virlrc+vlrcm(1))/3.d0
      stress(5) = stress(5) - (virlrc+vlrcm(1))/3.d0
      stress(9) = stress(9) - (virlrc+vlrcm(1))/3.d0
#endif
! 
!     numerical differentiation
      if( lfd ) then

         if (llfirst) then
            
            llfirst = .false.

            do i=1,natms
               
               fanaly(1,i) = fxx(i)
               fanaly(2,i) = fyy(i)
               fanaly(3,i) = fzz(i)
               
            end do

         end if

         engtmp = engsrp+engcpe+engbnd+engang+engdih+engfld+
     $        engtbp+engfbp+engshl+enginv

         if( lfdx ) then 

            call fordebug( lreturn, 1, icount, ilfd, natms,
     $           delta, engtmp, fanaly, xxx )

            if( lreturn ) goto 100 

            lfdx = .false.
            ilfd = 0
         
         end if

         if( lfdy ) then

            call fordebug( lreturn, 2, icount, ilfd, natms,
     $           delta, engtmp, fanaly, yyy )

            if( lreturn ) goto 100 

            lfdy = .false.
            ilfd = 0
            
         end if

         if( lfdz ) then
            
            call fordebug( lreturn, 3, icount, ilfd, natms,
     $           delta, engtmp, fanaly, zzz)

            if( lreturn ) goto 100 

         end if

         write(nrite,*) 'finite difference test is finished'
         stop

      end if  ! lfd

      if (lttm) then
!
! bug check for force decompositions
!         if( mod( nstep, nstbpo ) == 0 ) then
!           do i=1,natms
!              write(6,'(i5,3d30.12)') i, fxx(i), fyy(i), fzz(i)
!              write(6,'(i5,3d30.12)') i,
!     $             fxx(i) * 3.80880d-6 * 0.529177249d0,
!     $             fyy(i) * 3.80880d-6 * 0.529177249d0,
!     $             fzz(i) * 3.80880d-6 * 0.529177249d0
!            enddo
!         endif            
c
c     for ttm2 water use real degree of freedom

        natms=natms_r
c
c     distribute force on the m-site to oxygen and hydrogen

        mttm2=listttm2(nttm2)
        do i=1,nummols(ntpmls)
          mttm2=mttm2+1
          ioxy=listttm2(3*i-2)
          ih1=listttm2(3*i-1)
          ih2=listttm2(3*i)
          fxx(ioxy)=fxx(ioxy)+(1.d0-gammattm)*fxx(mttm2)
          fyy(ioxy)=fyy(ioxy)+(1.d0-gammattm)*fyy(mttm2)
          fzz(ioxy)=fzz(ioxy)+(1.d0-gammattm)*fzz(mttm2)
          fxx(ih1)=fxx(ih1)+fxx(mttm2)*gammattm/2.d0
          fyy(ih1)=fyy(ih1)+fyy(mttm2)*gammattm/2.d0
          fzz(ih1)=fzz(ih1)+fzz(mttm2)*gammattm/2.d0
          fxx(ih2)=fxx(ih2)+fxx(mttm2)*gammattm/2.d0
          fyy(ih2)=fyy(ih2)+fyy(mttm2)*gammattm/2.d0
          fzz(ih2)=fzz(ih2)+fzz(mttm2)*gammattm/2.d0
        enddo

      endif

c     
c     cap forces in equilibration mode
      
      if(nstep.le.nsteql.and.lfcap) 
     x  call fcap(lfcap,idnode,mxnode,natms,fmax,temp,fxx,fyy,fzz)
c     
c     zero Kelvin structure optimisation 
      
      if(loptim) call strucopt
     x  (loptim,idnode,mxnode,imcon,natms,ngrp,ntfree,
     x  lstfre,lstgtp,lstrgd,numgsit,cell,gvxx,gvyy,gvzz,
     x  gcmx,gcmy,gcmz,omx,omy,omz,vxx,vyy,vzz,fxx,fyy,
     x  fzz,xxx,yyy,zzz,xxt,yyt,zzt,q0,q1,q2,q3,rotinx,
     x  rotiny,rotinz)
      
c     
c     frozen atoms option
      
      call freeze(idnode,mxnode,natms,lstfrz,vxx,vyy,vzz,fxx,fyy,fzz)
c     
c     total virial (excluding constraint virial and c.o.m virial)
c      for npt routines     note: virsrp already includes virlrc
      
      virtot = vircpe+virsrp+virbnd+virtbp+virfld+virang+virshl+virtet

 1003 continue
cty
        if (lopt2) then

          write(*,*)'in dlpoly:'
          write(*,*)'engsrp+engcpe+engbnd+engang+engdih+engfld+engtbp+
     x     +engfbp+engshl+enginv'
          write(*,*)engsrp/418.4,engcpe/418.4,engbnd/418.4,
     x     engang/418.4,engdih/418.4

          uuu=engsrp+engcpe+engbnd+engang+engdih+engfld+engtbp+engfbp
     x       +engshl+enginv

          write(*,*)'in dlpoly, uuu =',uuu/engunit
          stop

        endif
c
c     propagate dipole degree of freedom

      if (lpolar .and. lcp .and. ldpts) then

         dipx_old(:) = dipx(:)
         dipy_old(:) = dipy(:)
         dipz_old(:) = dipz(:)

        call dcp_h0
     x  (dipmas,dpms,polr2,dipx,dipy,dipz,vdxx,vdyy,vdzz,
     x  udxx,udyy,udzz,fdxx,fdyy,fdzz,chitd,engdke,emux,emuy,emuz,
     x  efieldkx,efieldky,efieldkz,efdcrecx,efdcrecy,efdcrecz,
     x  efddmurecx,efddmurecy,efddmurecz,sigmad,tautd,consvd,conintd,
     x  idnode,mxnode,natms_f,tstep,buffer)

      elseif (lpolar .and. lcp) then

         dipx_old(:) = dipx(:)
         dipy_old(:) = dipy(:)
         dipz_old(:) = dipz(:)

        call dcp_0
     x  (dipmas,dpms,polr2,dipx,dipy,dipz,vdxx,vdyy,vdzz,
     x  engdke,emux,emuy,emuz,
     x  efieldkx,efieldky,efieldkz,
     x  efdcrecx,efdcrecy,efdcrecz,
     x  efddmurecx,efddmurecy,efddmurecz,
     x  idnode,mxnode,natms,imcon,tstep,buffer)

      endif

ccc      goto 1004
!
!    save old coordinates 
!    (coordinates after integration do not belong to the same time step 
!     as forces)

      xxold(:) = xxx(:)
      yyold(:) = yyy(:)
      zzold(:) = zzz(:)
c
c     integrate equations of motion

      if(ngrp.eq.0) then
        
        if(keyens.eq.0) then
c     
c     verlet leapfrog 

          call nve_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,ntcons,
     x      engke,tolnce,tstep,vircon,lashap,lishap,listcon,
     x      listme,listot,lstfrz,buffer,cell,dxt,dxx,dyt,
     x      dyy,dzt,dzz,fxx,fyy,fzz,prmcon,txx,tyy,tzz,uxx,
     x      uyy,uzz,vxx,vyy,vzz,weight,xdf,xxt,xxx,ydf,yyt,
     x      yyy,zdf,zzt,zzz,stress,lttm)
          
        else if(keyens.eq.1) then
c     
c     Evans Gaussian Temperature constraints
          
          call nvt_e1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,ntcons,
     x      engke,sigma,tolnce,tstep,vircon,lashap,lishap,
     x      listcon,listme,listot,lstfrz,buffer,cell,dxt,
     x      dxx,dyt,dyy,dzt,dzz,fxx,fyy,fzz,prmcon,txx,tyy,
     x      tzz,uxx,uyy,uzz,vxx,vyy,vzz,weight,xdf,xxt,xxx,
     x      ydf,yyt,yyy,zdf,zzt,zzz,stress)
          
        else if(keyens.eq.2) then
c     
c     Berendsen thermostat
          
          call nvt_b1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,ntcons,
     x      engke,taut,sigma,tolnce,tstep,vircon,lashap,
     x      lishap,listcon,listme,lstfrz,listot,buffer,
     x      cell,dxt,dxx,dyt,dyy,dzt,dzz,fxx,fyy,fzz,prmcon,
     x      txx,tyy,tzz,uxx,uyy,uzz,vxx,vyy,vzz,weight,xdf,
     x      xxt,xxx,ydf,yyt,yyy,zdf,zzt,zzz,vx1,vy1,vz1,stress)
          
        else if(keyens.eq.3) then
c     
c     Nose-Hoover thermostat
          
          call nvt_h1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,ntcons,
     x      chit,consv,conint,engke,taut,sigma,tolnce,tstep,vircon,
     x      lashap,lishap,listcon,listme,listot,lstfrz,buffer,
     x      cell,dxt,dxx,dyt,dyy,dzt,dzz,fxx,fyy,fzz,prmcon,txx,
     x      tyy,tzz,uxx,uyy,uzz,vxx,vyy,vzz,weight,xdf,xxt,xxx,
     x      ydf,yyt,yyy,zdf,zzt,zzz,vx1,vy1,vz1,stress,lttm)

        elseif(keyens.eq.4) then
c     
c     Berendsen thermostat and isotropic barostat
          
          call npt_b1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,ntpatm,nscons,
     x      ntcons,elrc,engke,virlrc,press,taup,taut,sigma,tolnce,
     x      tstep,virtot,vircon,volm,lashap,lishap,listcon,listme,
     x      lstfrz,listot,buffer,cell,dens,dxt,dxx,dyt,dyy,dzt,dzz,
     x      fxx,fyy,fzz,prmcon,txx,tyy,tzz,uxx,uyy,uzz,vxx,vyy,vzz,
     x      weight,xdf,xxt,xxx,ydf,yyt,yyy,zdf,zzt,zzz,vx1,vy1,vz1,
     x      stress,eta,dens0)

        else if(keyens.eq.5) then
c     
c     Nose-Hoover thermostat and isotropic barostat 

          call npt_h1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,ntpatm,nscons,
     x      ntcons,chip,chit,conint,consv,elrc,engke,virlrc,press,
     x      taup,taut,sigma,temp,tolnce,tstep,virtot,vircon,volm,
     x      lashap,lishap,listcon,listme,lstfrz,listot,buffer,
     x      cell,dens,dxt,dxx,dyt,dyy,dzt,dzz,fxx,fyy,fzz,prmcon,
     x      txx,tyy,tzz,uxx,uyy,uzz,vxx,vyy,vzz,weight,xdf,xxt,
     x      xxx,ydf,yyt,yyy,zdf,zzt,zzz,vx1,vy1,vz1,stress,eta,dens0)

        else if(keyens.eq.6) then
c     
c     Berendsen thermostat and barostat (cell shape varying)

          call npt_b3
     x      (safe,lshmov,idnode,imcon,mxnode,natms,ntpatm,nscons,
     x      ntcons,elrc,engke,virlrc,press,taup,taut,sigma,tolnce,
     x      tstep,vircon,volm,lashap,lishap,listcon,listme,
     x      lstfrz,listot,buffer,cell,dens,dxt,dxx,dyt,dyy,dzt,dzz,
     x      fxx,fyy,fzz,prmcon,txx,tyy,tzz,uxx,uyy,uzz,vxx,vyy,vzz,
     x      weight,xdf,xxt,xxx,ydf,yyt,yyy,zdf,zzt,zzz,vx1,vy1,vz1,
     x      stress,eta,dens0)

        else if(keyens.eq.7) then
c     
c     Nose-Hoover thermostat and barostat (cell shape varying)
          
          call npt_h3
     x      (safe,lshmov,idnode,imcon,mxnode,natms,ntpatm,nscons,
     x      ntcons,chit,conint,consv,elrc,engke,virlrc,press,
     x      taup,taut,sigma,temp,tolnce,tstep,vircon,volm,
     x      lashap,lishap,listcon,listme,lstfrz,listot,buffer,
     x      cell,dens,dxt,dxx,dyt,dyy,dzt,dzz,fxx,fyy,fzz,prmcon,
     x      txx,tyy,tzz,uxx,uyy,uzz,vxx,vyy,vzz,weight,xdf,xxt,
     x      xxx,ydf,yyt,yyy,zdf,zzt,zzz,vx1,vy1,vz1,stress,eta,dens0)

        elseif(keyens.eq.8) then
c
c     Potential of mean force in NVE

            call pmf_1
     x        (safe,safep,lshmov,idnode,imcon,mxnode,natms,nscons,
     x        ntcons,nspmf,ntpmf,engke,prmpmf,tolnce,tstep,vircon,
     x        virpmf,lstpmf,npmf,lashap,lishap,listcon,listpm,lstpmt,
     x        listme,listot,lstfrz,buffer,cell,dxt,dxx,dyt,
     x        dyy,dzt,dzz,pmfwght,fxx,fyy,fzz,prmcon,txx,tyy,tzz,
     x        uxx,uyy,uzz,vxx,vyy,vzz,weight,xdf,xxt,xxx,ydf,yyt,
     x        yyy,zdf,zzt,zzz,stress,xa,ya,za,dxp,dyp,dzp,dsq)

            vircon = vircon + virpmf

        endif

      elseif(ngrp.gt.0) then
c     
c     apply rigid body equations of motion
        
        if(keyens.eq.0) then
          
          if(.not.lcnb) then

            call nveq_1
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,quattol,tolnce,tstep,vircom,
     x        vircon,lstme,lishap,lashap,lstfrz,listcon,listme,lstfre,
     x        lstgtp,lstrgd,listot,numgsit,cell,gcmx1,gcmy1,gcmz1,dxx,
     x        dyy,dzz,dxt,dyt,dzt,uxx,uyy,uzz,prmcon,omx,omy,omz,opx,
     x        opy,opz,oqx,oqy,oqz,q0,q1,q2,q3,qn0,qn1,qn2,qn3,gxx,gyy,
     x        gzz,xxx,yyy,zzz,vxx,vyy,vzz,gvxx,gvyy,gvzz,fxx,fyy,fzz,
     x        tqx,tqy,tqz,weight,rotinx,rotiny,rotinz,gcmx,gcmy,gcmz,
     x        gmass,buffer,xxt,yyt,zzt,txx,tyy,tzz,vx1,vy1,vz1,stress)

          else

            call nveq_2
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,quattol,tolnce,tstep,vircom,
     x        vircon,lstme,lishap,lashap,lstfrz,listcon,listme,lstfre,
     x        lstcsit,lstbod,lstgtp,lstrgd,listot,numgsit,cell,gcmx1,
     x        gcmy1,gcmz1,dxx,dyy,dzz,dxt,dyt,dzt,uxx,uyy,uzz,prmcon,
     x        omx,omy,omz,opx,opy,opz,oqx,oqy,oqz,q0,q1,q2,q3,qn0,qn1,
     x        qn2,qn3,gxx,gyy,gzz,xxx,yyy,zzz,vxx,vyy,vzz,gvxx,gvyy,
     x        gvzz,fxx,fyy,fzz,tqx,tqy,tqz,weight,rotinx,rotiny,rotinz,
     x        gcmx,gcmy,gcmz,gmass,buffer,xxt,yyt,zzt,txx,tyy,tzz,vx1,
     x        vy1,vz1,stress,gvx1,gvy1,gvz1,omx1,omy1,omz1,redmass,
     x        esig1)

          endif

        elseif(keyens.eq.1) then
c
c     invalid option
          call error(idnode,430)
          
        elseif(keyens.eq.2) then
          
          if(.not.lcnb) then

            call nvtq_b1
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,quattol,sigma,taut,tolnce,
     x        tstep,vircom,vircon,lstme,lishap,lashap,lstfrz,listcon,
     x        listme,lstfre,lstgtp,lstrgd,listot,numgsit,cell,gcmx1,
     x        gcmy1,gcmz1,gvx1,gvy1,gvz1,dxx,dyy,dzz,dxt,dyt,dzt,uxx,
     x        uyy,uzz,prmcon,xdf,ydf,zdf,omx,omy,omz,omx1,omy1,omz1,
     x        opx,opy,opz,oqx,oqy,oqz,q0,q1,q2,q3,qn0,qn1,qn2,qn3,gxx,
     x        gyy,gzz,xxx,yyy,zzz,vxx,vyy,vzz,gvxx,gvyy,gvzz,fxx,fyy,
     x        fzz,tqx,tqy,tqz,weight,rotinx,rotiny,rotinz,gcmx,gcmy,
     x        gcmz,gmass,buffer,xxt,yyt,zzt,txx,tyy,tzz,vx1,vy1,vz1,
     x        stress)

          else

            call nvtq_b2
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,quattol,sigma,taut,tolnce,
     x        tstep,vircom,vircon,lstme,lishap,lashap,lstfrz,listcon,
     x        listme,lstfre,lstcsit,lstbod,lstgtp,lstrgd,listot,
     x        numgsit,cell,gcmx1,gcmy1,gcmz1,dxx,dyy,dzz,dxt,dyt,dzt,
     x        uxx,uyy,uzz,prmcon,omx,omy,omz,opx,opy,opz,oqx,oqy,oqz,
     x        q0,q1,q2,q3,qn0,qn1,qn2,qn3,gxx,gyy,gzz,xxx,yyy,zzz,vxx,
     x        vyy,vzz,gvxx,gvyy,gvzz,fxx,fyy,fzz,tqx,tqy,tqz,weight,
     x        rotinx,rotiny,rotinz,gcmx,gcmy,gcmz,gmass,buffer,xxt,yyt,
     x        zzt,txx,tyy,tzz,vx1,vy1,vz1,stress,gvx1,gvy1,gvz1,omx1,
     x        omy1,omz1,redmass,esig1)
          
          endif

        elseif(keyens.eq.3) then
          
          if(.not.lcnb) then 

            call nvtq_h1
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,chit,consv,conint,engke,engrot,quattol,
     x        sigma,taut,tolnce,tstep,vircom,vircon,lstme,lishap,
     x        lashap,lstfrz,listcon,listme,lstfre,lstgtp,lstrgd,
     x        listot,numgsit,cell,gcmx1,gcmy1,gcmz1,dxx,dyy,dzz,
     x        dxt,dyt,dzt,uxx,uyy,uzz,prmcon,xdf,ydf,zdf,omx,omy,
     x        omz,omx1,omy1,omz1,opx,opy,opz,oqx,oqy,oqz,q0,q1,q2,q3,
     x        qn0,qn1,qn2,qn3,gxx,gyy,gzz,xxx,yyy,zzz,vxx,vyy,vzz,gvxx,
     x        gvyy,gvzz,fxx,fyy,fzz,tqx,tqy,tqz,weight,rotinx,rotiny,
     x        rotinz,gcmx,gcmy,gcmz,gmass,buffer,xxt,yyt,zzt,txx,
     x        tyy,tzz,gvx1,gvy1,gvz1,vx1,vy1,vz1,stress)

          else

            call nvtq_h2
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,conint,consv,chit,engke,engrot,quattol,
     x        sigma,taut,tolnce,tstep,vircom,vircon,lstme,lishap,
     x        lashap,lstfrz,listcon,listme,lstfre,lstcsit,lstbod,
     x        lstgtp,lstrgd,listot,numgsit,cell,gcmx1,gcmy1,gcmz1,dxx,
     x        dyy,dzz,dxt,dyt,dzt,uxx,uyy,uzz,prmcon,omx,omy,omz,opx,
     x        opy,opz,oqx,oqy,oqz,q0,q1,q2,q3,qn0,qn1,qn2,qn3,gxx,gyy,
     x        gzz,xxx,yyy,zzz,vxx,vyy,vzz,gvxx,gvyy,gvzz,fxx,fyy,fzz,
     x        tqx,tqy,tqz,weight,rotinx,rotiny,rotinz,gcmx,gcmy,gcmz,
     x        gmass,buffer,xxt,yyt,zzt,txx,tyy,tzz,vx1,vy1,vz1,stress,
     x        gvx1,gvy1,gvz1,omx1,omy1,omz1,redmass,esig1)

          endif
            
        elseif(keyens.eq.4) then

          if(.not.lcnb) then

            call nptq_b1
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,elrc,engke,engrot,virlrc,press,
     x        quattol,sigma,taup,taut,tolnce,tstep,virtot,vircom,
     x        vircon,volm,lstme,lishap,lashap,lstfrz,listcon,
     x        listme,lstfre,lstgtp,lstrgd,listot,numgsit,cell,
     x        dens,gcmx1,gcmy1,gcmz1,dxx,dyy,dzz,dxt,dyt,dzt,
     x        gvx1,gvy1,gvz1,uxx,uyy,uzz,prmcon,xdf,ydf,zdf,omx,
     x        omy,omz,omx1,omy1,omz1,opx,opy,opz,oqx,oqy,oqz,
     x        q0,q1,q2,q3,qn0,qn1,qn2,qn3,gxx,gyy,gzz,xxx,yyy,
     x        zzz,vxx,vyy,vzz,gvxx,gvyy,gvzz,fxx,fyy,fzz,tqx,tqy,
     x        tqz,weight,rotinx,rotiny,rotinz,gcmx,gcmy,gcmz,gmass,
     x        buffer,xxt,yyt,zzt,txx,tyy,tzz,vx1,vy1,vz1,stress,eta,
     x        dens0)
          
          else

            call nptq_b2
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,elrc,engke,engrot,virlrc,press,
     x        quattol,sigma,taup,taut,tolnce,tstep,vircom,vircon,
     x        virtot,volm,lstme,lishap,lashap,lstfrz,listcon,listme,
     x        lstfre,lstcsit,lstbod,lstgtp,lstrgd,listot,numgsit,cell,
     x        gcmx1,gcmy1,gcmz1,dens,dxx,dyy,dzz,dxt,dyt,dzt,uxx,uyy,
     x        uzz,prmcon,omx,omy,omz,opx,opy,opz,oqx,oqy,oqz,q0,q1,q2,
     x        q3,qn0,qn1,qn2,qn3,gxx,gyy,gzz,xxx,yyy,zzz,vxx,vyy,vzz,
     x        gvxx,gvyy,gvzz,fxx,fyy,fzz,tqx,tqy,tqz,weight,rotinx,
     x        rotiny,rotinz,gcmx,gcmy,gcmz,gmass,buffer,xxt,yyt,zzt,
     x        txx,tyy,tzz,vx1,vy1,vz1,stress,gvx1,gvy1,gvz1,omx1,omy1,
     x        omz1,redmass,eta,dens0,esig1)

          endif

        elseif(keyens.eq.5) then
          
          if(.not.lcnb) then 

            call nptq_h1
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,chip,chit,consv,conint,elrc,engke,
     x        engrot,virlrc,press,quattol,sigma,taup,taut,temp,tolnce,
     x        tstep,virtot,vircom,vircon,volm,lstme,lishap,lashap,
     x        lstfrz,listcon,listme,lstfre,lstgtp,lstrgd,listot,
     x        numgsit,cell,dens,dxx,dyy,dzz,dxt,dyt,dzt,gcmx1,gcmy1,
     x        gcmz1,gvx1,gvy1,gvz1,uxx,uyy,uzz,prmcon,xdf,ydf,zdf,omx,
     x        omy,omz,omx1,omy1,omz1,opx,opy,opz,oqx,oqy,oqz,q0,q1,q2,
     x        q3,qn0,qn1,qn2,qn3,gxx,gyy,gzz,xxx,yyy,zzz,vxx,vyy,vzz,
     x        gvxx,gvyy,gvzz,fxx,fyy,fzz,tqx,tqy,tqz,weight,rotinx,
     x        rotiny,rotinz,gcmx,gcmy,gcmz,gmass,buffer,xxt,yyt,zzt,
     x        txx,tyy,tzz,vx1,vy1,vz1,stress,eta,dens0)

          else

            call nptq_h2
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,chip,chit,consv,conint,elrc,engke,
     x        engrot,virlrc,press,quattol,sigma,taup,taut,temp,tolnce,
     x        tstep,vircom,vircon,virtot,volm,lstme,lishap,lashap,
     x        lstfrz,listcon,listme,lstfre,lstcsit,lstbod,lstgtp,lstrgd,
     x        listot,numgsit,cell,gcmx1,gcmy1,gcmz1,dens,dxx,dyy,dzz,
     x        dxt,dyt,dzt,uxx,uyy,uzz,prmcon,omx,omy,omz,opx,opy,opz,
     x        oqx,oqy,oqz,q0,q1,q2,q3,qn0,qn1,qn2,qn3,gxx,gyy,gzz,xxx,
     x        yyy,zzz,vxx,vyy,vzz,gvxx,gvyy,gvzz,fxx,fyy,fzz,tqx,tqy,
     x        tqz,weight,rotinx,rotiny,rotinz,gcmx,gcmy,gcmz,gmass,
     x        buffer,xxt,yyt,zzt,txx,tyy,tzz,vx1,vy1,vz1,stress,gvx1,
     x        gvy1,gvz1,omx1,omy1,omz1,redmass,eta,dens0,esig1)

          endif
            
        elseif(keyens.eq.6) then

          if(.not.lcnb) then

            call nptq_b3
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,elrc,engke,engrot,virlrc,press,
     x        quattol,sigma,taup,taut,tolnce,tstep,vircom,vircon,volm,
     x        lstme,lishap,lashap,lstfrz,listcon,listme,lstfre,lstgtp,
     x        lstrgd,listot,numgsit,cell,dens,gcmx1,gcmy1,gcmz1,dxx,
     x        dyy,dzz,dxt,dyt,dzt,gvx1,gvy1,gvz1,uxx,uyy,uzz,prmcon,
     x        xdf,ydf,zdf,omx,omy,omz,omx1,omy1,omz1,opx,opy,opz,oqx,
     x        oqy,oqz,q0,q1,q2,q3,qn0,qn1,qn2,qn3,gxx,gyy,gzz,xxx,yyy,
     x        zzz,vxx,vyy,vzz,gvxx,gvyy,gvzz,fxx,fyy,fzz,tqx,tqy,
     x        tqz,weight,rotinx,rotiny,rotinz,gcmx,gcmy,gcmz,gmass,
     x        buffer,xxt,yyt,zzt,txx,tyy,tzz,vx1,vy1,vz1,stress,eta,
     x        dens0)

          else

            call nptq_b4
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,elrc,engke,engrot,virlrc,press,
     x        quattol,sigma,taup,taut,tolnce,tstep,vircom,vircon,volm,
     x        lstme,lishap,lashap,lstfrz,listcon,listme,lstfre,lstgtp,
     x        lstrgd,listot,lstcsit,lstbod,numgsit,cell,dens,gcmx1,
     x        gcmy1,gcmz1,dxx,dyy,dzz,dxt,dyt,dzt,gvx1,gvy1,gvz1,uxx,
     x        uyy,uzz,prmcon,omx,omy,omz,omx1,omy1,omz1,opx,opy,opz,
     x        oqx,oqy,oqz,q0,q1,q2,q3,qn0,qn1,qn2,qn3,gxx,gyy,gzz,xxx,
     x        yyy,zzz,vxx,vyy,vzz,gvxx,gvyy,gvzz,fxx,fyy,fzz,tqx,tqy,
     x        tqz,weight,rotinx,rotiny,rotinz,gcmx,gcmy,gcmz,gmass,
     x        buffer,xxt,yyt,zzt,txx,tyy,tzz,vx1,vy1,vz1,stress,redmass,
     x        eta,dens0,esig1)

          endif

        elseif(keyens.eq.7) then

          if(.not.lcnb) then

            call nptq_h3
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,chit,conint,consv,elrc,engke,engrot,
     x        virlrc,press,quattol,sigma,taup,taut,temp,tolnce,tstep,
     x        vircom,vircon,volm,lstme,lishap,lashap,lstfrz,listcon,
     x        listme,lstfre,lstgtp,lstrgd,listot,numgsit,cell,dens,
     x        gcmx1,gcmy1,gcmz1,dxx,dyy,dzz,dxt,dyt,dzt,gvx1,gvy1,gvz1,
     x        uxx,uyy,uzz,prmcon,xdf,ydf,zdf,omx,omy,omz,omx1,omy1,
     x        omz1,opx,opy,opz,oqx,oqy,oqz,q0,q1,q2,q3,qn0,qn1,qn2,qn3,
     x        gxx,gyy,gzz,xxx,yyy,zzz,vxx,vyy,vzz,gvxx,gvyy,gvzz,fxx,
     x        fyy,fzz,tqx,tqy,tqz,weight,rotinx,rotiny,rotinz,gcmx,
     x        gcmy,gcmz,gmass,buffer,xxt,yyt,zzt,txx,tyy,tzz,vx1,vy1,
     x        vz1,stress,eta,dens0)

          else

            call nptq_h4
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,chit,conint,consv,elrc,engke,engrot,
     x        virlrc,press,quattol,sigma,taup,taut,temp,tolnce,tstep,
     x        vircom,vircon,volm,lstme,lishap,lashap,lstfrz,listcon,
     x        listme,lstfre,lstgtp,lstrgd,listot,lstcsit,lstbod,numgsit,
     x        cell,dens,gcmx1,gcmy1,gcmz1,dxx,dyy,dzz,dxt,dyt,dzt,gvx1,
     x        gvy1,gvz1,uxx,uyy,uzz,prmcon,omx,omy,omz,omx1,
     x        omy1,omz1,opx,opy,opz,oqx,oqy,oqz,q0,q1,q2,q3,qn0,qn1,qn2,
     x        qn3,gxx,gyy,gzz,xxx,yyy,zzz,vxx,vyy,vzz,gvxx,gvyy,gvzz,
     x        fxx,fyy,fzz,tqx,tqy,tqz,weight,rotinx,rotiny,rotinz,gcmx,
     x        gcmy,gcmz,gmass,buffer,xxt,yyt,zzt,txx,tyy,tzz,vx1,vy1,
     x        vz1,stress,eta,redmass,dens0,esig1)

            
          endif

        else
c
c     invalid option
          call error(idnode,430)

        endif

      endif
c     
c    check on convergence of pmf-shake

      if(ntpmf.gt.0) then
        if(mxnode.gt.1) call gstate(safep)
        if(.not.safep) call error(idnode,438)
      endif    
c     
c    check on convergence of shake

      if(ntcons.gt.0) then
        if(mxnode.gt.1) call gstate(safe)
        if(.not.safe) call error(idnode,105)
      endif    
c
c     check on convergence of quaternion algorithm

      if(ngrp.gt.0) then

        if(mxnode.gt.1) call gstate(safeq)
        if(.not.safeq) call error(idnode,321)

      endif
c
c     calculate core-shell kinetic energy

      if(ntshl.gt.0) then

        if(ngrp.gt.0) then

          do i = 1,4
            call shlqnch
     x        (idnode,mxnode,ntshl,temp,listshl,weight,vxx,vyy,vzz,
     x        buffer)

            if(ngrp.gt.0) call quatqnch
     x        (idnode,imcon,mxnode,natms,ngrp,lstgtp,lstrgd,numgsit,
     x        lstme,gxx,gyy,gzz,buffer,cell,xxt,yyt,zzt,gcmx,gcmy,
     x        gcmz,gmass,gvxx,gvyy,gvzz,q0,q1,q2,q3,omx,omy,omz,
     x        rotinx,rotiny,rotinz,vxx,vyy,vzz,weight,xxx,yyy,zzz)

          enddo

        endif

        call corshl
     x    (idnode,mxnode,ntshl,shlke,listshl,weight,vxx,vyy,vzz,buffer)

      endif
c     
c     apply temperature scaling
      
      if((ltscal.and.nstep.le.nsteql).and.
     x  mod(nstep-nsteql,nstbts).eq.0)then
        
        chit = 0.d0
        chip = 0.d0
        do i = 1,9
          eta(i) = 0.d0
        enddo

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
     x      (idnode,mxnode,imcon,natms,ngrp,sigma,lstfrz,lstgtp,lstrgd,
     x      numgsit,lstme,buffer,cell,gcmx,gcmy,gcmz,gmass,
     x      gvxx,gvyy,gvzz,gxx,gyy,gzz,q0,q1,q2,q3,omx,omy,omz,weight,
     x      vxx,vyy,vzz,xxx,yyy,zzz,xxt,yyt,zzt,rotinx,rotiny,rotinz)
        
        endif

      endif
c
c     scale t=0 tether reference positions (constant pressure only)

      if(keyens.ge.4.and.keyens.le.7) then

        call xscale
     x  (idnode,mxnode,natms,keyens,imcon,tstep,
     x    cell,xxs,yys,zzs,weight,eta,buffer)

      endif

 1004 continue
c     
c     calculate physical quantities
      
      call static
     x  (lzeql,cfgname,idnode,intsta,imcon,keyens,natms,nstack,nstep,
     x  nsteql,ntpatm,numacc,mxnode,consv,degfre,degrot,engang,
     x  engbnd,engcpe,engdih,enginv,engke,engrot,engsrp,engunit,stpcfg,
     x  stpeng,stpeth,stpprs,stptmp,stpvir,stpvol,tstep,virbnd,engfbp,
     x  vircom,vircon,vircpe,virsrp,engfld,virfld,engtbp,virtbp,
     x  virpmf,virshl,engshl,engtet,virtet,degshl,shlke,virang,
     x  width,ltype,numtyp,buffer,cell,chge,fxx,fyy,fzz,ravval,
     x  ssqval,stkval,stpval,sumval,vxx,vyy,vzz,xxx,yyy,zzz,zumval,
     x  xx0,yy0,zz0,weight,stress,amsd,engdke,consvd)

c     
c     z density calculation

      if(lzden.and.((.not.lzeql).or.(nstep.gt.nsteql))) then

        zlen=abs(cell(3))+abs(cell(6))+abs(cell(9))
        call zden0(idnode,natms,mxnode,nzden,zlen,ltype,zzz,zdens)

      endif
c
c     center of mass of the molecules for every step

c      call cenmas
c     x  (idnode,imcon,mxnode,ntbond,listbnd,cell,
c     x  natms,ntpmls,nummols,numsit,wgtsit,numbonds,nums,
c     x  listbnd2,listbnd3,listin2,tstep,xmsd,ymsd,zmsd,
c     x  xxx,yyy,zzz,vxx,vyy,vzz,xxx1,yyy1,zzz1,xdab,ydab,zdab,
c     x  xdab2,ydab2,zdab2,xcm,ycm,zcm,vxcm,vycm,vzcm,buffer)
c     
c     terminate program if boundary conditions violated
      
      if(imcon.gt.0.and.rcut.gt.width)then
        
        call revive
     x    (lgofr,lzden,lpolar,cfgname,atmnam,idnode,imcon,mxnode,
     x    natms,nstep,nzden,numacc,numrdf,chip,chit,chitd,conint,
     x    tstep,buffer,cell,fxx,fyy,fzz,ravval,rdf,ssqval,
     x    stkval,stpval,sumval,vxx,vyy,vzz,xxx,yyy,zumval,
     x    zzz,xx0,yy0,zz0,zdens,xxs,yys,zzs,eta,dipx,dipy,dipz,
     x    vdxx,vdyy,vdzz,xmsd,ymsd,zmsd,lcp,ldpts,conintd)
        
        call error(idnode,95)

      endif
c     
c     line-printer output every nstbpo steps
      
#ifdef VAMPIR
      call VTBEGIN(68, ierr)
#endif

      if(lines.eq.0.or.mod(nstep,nstbpo).eq.0)then

        call timchk(0,timelp)
        if(idnode.eq.0)then

          if(mod(lines,npage).eq.0)
     x      write(nrite,"(1x,130('-'),
     x      /,/,1x,'    step',5x,'eng_tot',4x,'temp_tot',5x,'eng_cfg',
     x      5x,'eng_vdw',5x,'eng_cou',5x,'eng_bnd',5x,'eng_ang',5x,
     x      'eng_dih',5x,'eng_tet',/,1x,'time(ps)',5x,' eng_pv',4x,
     x      'temp_rot',5x,'vir_cfg',5x,'vir_vdw',5x,'vir_cou',5x,
     x      'vir_bnd',5x,'vir_ang',5x,'vir_con',5x,'vir_tet',/,
     x      1x,'cpu  (s)',6x,'volume',4x,'temp_shl',5x,'eng_shl',
     x      5x,'vir_shl',7x,'alpha',8x,'beta',7x,'gamma',5x,'vir_pmf',
     x      7x,'press',/,/,
     x      1x,130('-'))")          
          write(nrite,"(1x,i10,1p,9e12.4,/,1x,0p,f10.4,1p,9e12.4,
     x      /,1x,0p,f10.2,1p,9e12.4)")
     x      nstep,(stpval(i),i=1,9),
     x      dble(nstep)*tstep,(stpval(i),i=10,18),
     x      timelp,(stpval(i),i=19,27)
          write(nrite,"(/,1x,' rolling',1p,9e12.4,/,1x,'averages',
     x      1p,9e12.4,/,9x,9e12.4)") (ravval(i),i=1,27)
          write(nrite,"(1x,130('-'))")

!          idum=ntpatm+27
!
!          write(nrite,"(/,/,16x,'Average pressure tensor',
!     x      39x,'  m.s. fluctuations ',/)")
!
!          do i = idum,idum+6,3
!            write(nrite,'(9x,1p,3e12.4,24x,3e12.4))')
!     x        (sumval(i+j),j = 1,3),(ssqval(i+j),j = 1,3)
!          enddo
!
!          write(nrite,'(/,9x,a,24x,a)')
!     x         '    average surface tension (mN / m)',
!     x         'fluctuation of surface tension (mN / m)'
!          surtendum= cell(9)*(2.d0*sumval(idum+9)-sumval(idum+5)
!     x        -sumval(idum+1))/4.d0*10.d0
!          std1=dsqrt(ssqval(idum+9)/dble(nstep))
!          std2=dsqrt(ssqval(idum+5)/dble(nstep))
!          std3=dsqrt(ssqval(idum+1)/dble(nstep))
!          surtenflc= cell(9)*(2.d0*std1+std2+std3)/4.d0*10.d0
!          write(nrite,'(/,21x,e12.4,48x,e12.4,/,/)')
!     x        surtendum,surtenflc
!c     x        cell(9)*(2.d0*ssqval(idum+9)-ssqval(idum+5)
!c     x        -ssqval(idum+1))/4.d0*10.d0
!
!          if (mod(nstep,nstbpo).eq.0) then
!
!            tdd=dble(nstep)*tstep
!
!            open(39,position='append')
!            open(85,file='pressure_tensor',position='append')
!            open(86,file='surface_tension',position='append')
!            if (lpolar) open(87,file='zdipol',position='append')
!            open(88,file='zcoord',position='append')
!            open(89,file='zveloc',position='append')
!ccc            open(413,file='cmcoor',position='append')
!ccc            open(414,file='cmvelc',position='append')
!ccc            open(415,file='cmdisp',position='append')
!
!            write(39,'(f10.4,3x,f14.4)')tdd,stpval(1)
!            pdum=prsunt/stpvol
!            write(85,'(7f12.6,f12.4)')tdd,
!     x        stress(2)*pdum,stress(3)*pdum,stress(6)*pdum,
!     x        stress(1)*pdum,stress(5)*pdum,stress(9)*pdum,
!     x   cell(9)*(2.d0*stress(9)-stress(1)-stress(5))/4.d0*pdum*10.d0
!            write(86,'(2f12.6)')tdd,surtendum
!            if (lpolar) write(87,'(f12.6)')tdd
!            write(88,'(f12.6)')tdd
!            write(89,'(f12.6)')tdd
!            do i=1,natms
!              if (lpolar) write(87,'(3f9.5)')dipx(i),dipy(i),dipz(i)
!ccc     x           dipx(i)*eatd,dipy(i)*eatd,dipz(i)*eatd
!              write(88,'(3f9.4)')xxx(i),yyy(i),zzz(i)
!              write(89,'(3f10.4)')vxx(i),vyy(i),vzz(i)
!            enddo
!            if (lpolar .and. lcp) then
!              dptmp=2.d0*engdke/(boltz*3.d0*natms)
!              write(nrite,'(14x,a,d18.10,/)')'dipole temp: ',dptmp
!              open(36,position='append')
!              write(36,'(f12.4,3x,f10.4)')tdd,dptmp
!              close(36)
!            endif
!ccc            write(413,'(f18.6)')tdd
!ccc            write(414,'(f18.6)')tdd
!ccc            write(415,'(f18.6)')tdd
!ccc            do i=1,nums
!ccc              write(413,'(3f9.4)')xcm(i),ycm(i),zcm(i)
!ccc              write(414,'(3f10.4)')vxcm(i),vycm(i),vzcm(i)
!ccc              write(415,'(3f12.4)')xmsd(i),ymsd(i),zmsd(i)
!ccc            enddo
!
!            close(39)
!            close(85)
!            close(86)
!            if (lpolar) close(87)
!            close(88)
!            close(89)
!ccc            close(413)
!ccc            close(414)
!ccc            close(415)
!          endif
!
!c          write(nrite,'(/,14x,a,/)')'instantaneous stress tensor'
!c          write(nrite,'(9x,1p,3e12.4)')
!c     x        (stress(i)*prsunt/stpvol,i = 1,9)
!c          write(nrite,'(/,12x,a,1p,e12.4,/)') 
!c     x        'instantaneous surface tension (mN / m)',
!c     x        cell(9)*(2.d0*stress(9)-stress(1)-stress(5))
!c     x           /4.d0*prsunt/stpvol*10.d0

        endif
        
        lines=lines+1
        
      endif
#ifdef VAMPIR
      call VTEND(68, ierr)
#endif
c     
c     report end of equilibration period
      
      if(nstep.ge.nsteql)then
        
        if((ltscal.and.idnode.eq.0).and.(nstep.eq.nsteql))
     x    write(nrite,"(/,/,1x,'switching off temperature ',
     x    'scaling at step ',i6,/,/,/,1x,130('-'))") nstep
        ltscal=.false.
        
      endif
c     
c     write trajectory data
!
!     note:
!     xxold, yyold, and zzold are used in order to make 
!     forces and coordinates belong to the same time step.
!     Coordinates and forces written in trajectory 
!     will belong to time step n.
!     Velocities will belong to time step n+1/2
!
! ALSO: output necessary data (TTM2 model)
!
! note1: natms do not include M sites (no dipoles on M sites)
! note2: xxold, yyold, and zzold are used for force output so as to 
!        make both coordinates and forces belong to the same time step
!
      if(ltraj.and.nstep.ge.nstraj) then
        if(idnode.eq.0.and.mod(nstep-nstraj,istraj).eq.0)then
          
          if( lttm ) then 
             ptime  = dble( nstep ) * tstep
             ndummy = nummols(ntpmls)
             if( lcp ) then
                call output_ttm2_dip(listttm2,ndummy,nttm2,
     $               dipx_old,dipy_old,dipz_old,.false.)
             else
                call output_ttm2_dip(listttm2,ndummy,nttm2,
     $               dipx,dipy,dipz,.false.)
             endif
c             call output_ttm2_com(imcon,ndummy,ptime,cell,weight,
c     $            listttm2,xxold,yyold,zzold)
             if(ldecomp) call output_ttm2_force(atmnam,natms,nstep,
     $            nttm2,ndummy,listttm2,gammattm,xxold,yyold,zzold)
!             if(ldecomp) call output_ttm2_force(atmnam,natms_f,nstep,
!     $            ,xxold,yyold,zzold)
!
!            output site dipoles (non TTM water)
!
             open(901,file='SITE_DIPOLE',position='append')
             ndummy=listttm2(1)-1
             do i=1,ndummy
                if(lcp) then
                   write(901,'(i6,1x,a8,3f20.12)') i, atmnam(i),
     $                  dipx_old(i), dipy_old(i), dipz_old(i)
                else
                   write(901,'(i6,1x,a8,3f20.12)') i, atmnam(i),
     $                  dipx(i), dipy(i), dipz(i)
                endif
             enddo
             close(901)

          end if

           if(lunformat) then  ! UNFORMATTED

              call traject_u
     $           (ltraj,cfgname,atmnam,idnode,imcon,istraj,keytrj,natms,
     $            nstraj,nstep,tstep,cell,chge,weight,
     $            xxold,yyold,zzold,vxx,vyy,vzz,fxx,fyy,fzz)

           else   ! FORMATTED
              
              if( lttm ) then

                 call traject
     $                (ltraj,cfgname,atmnam,idnode,imcon,istraj,keytrj,
     $                natms,nstraj,nstep,tstep,cell,chge,weight,
     $                xxold,yyold,zzold,vx_at_t,vy_at_t,vz_at_t,
     $                fxx,fyy,fzz)
c$$$     $                xxold,yyold,zzold,vx0_at_t,vy0_at_t,vz0_at_t,
c$$$     $                fxx,fyy,fzz)
              else 

                 call traject
     $                (ltraj,cfgname,atmnam,idnode,imcon,istraj,keytrj,
     $                natms,nstraj,nstep,tstep,cell,chge,weight,
     $                xxold,yyold,zzold,vxx,vyy,vzz,fxx,fyy,fzz)
                 
              endif

           endif

!          call traject
!     x      (ltraj,cfgname,atmnam,idnode,imcon,istraj,keytrj,natms,
!     x      nstraj,nstep,tstep,cell,chge,weight,xxx,yyy,zzz,vxx,vyy,
!     x      vzz,fxx,fyy,fzz)
          
        endif
        
      endif
!
!     For REVCON: 'call revive' in 'subroutine result'
!     reset natms in case of TTM2 water
!     assign coordinates to the M sites for REVCON
!     forces for M-sites are zero

      if( lttm ) then 
         natms = natms_f        ! including M site
         mttm2=listttm2(nttm2)
         do i=1,nummols(ntpmls)
            mttm2=mttm2+1
            ioxy=listttm2(3*i-2)
            ih1=listttm2(3*i-1)
            ih2=listttm2(3*i)
            xdf(1)=xxx(ih1)-xxx(ioxy)
            ydf(1)=yyy(ih1)-yyy(ioxy)
            zdf(1)=zzz(ih1)-zzz(ioxy)
            xdf(2)=xxx(ih2)-xxx(ioxy)
            ydf(2)=yyy(ih2)-yyy(ioxy)
            zdf(2)=zzz(ih2)-zzz(ioxy)
            call images(imcon,0,1,2,cell,xdf,ydf,zdf)
            xxx(mttm2)=xxx(ioxy)+(xdf(1)+xdf(2))*gammattm/2.d0
            yyy(mttm2)=yyy(ioxy)+(ydf(1)+ydf(2))*gammattm/2.d0
            zzz(mttm2)=zzz(ioxy)+(zdf(1)+zdf(2))*gammattm/2.d0
            fxx(mttm2)=0.0d0
            fyy(mttm2)=0.0d0
            fzz(mttm2)=0.0d0
         enddo
      endif
c     
c     save restart data in event of system crash
      
c      if(mod(nstep,ndump).eq.0.and.nstep.ne.nstrun)then
      if(mod(nstep,nstbpo).eq.0.and.nstep.ne.nstrun)then

         call revive
     x    (lgofr,lzden,lpolar,cfgname,atmnam,idnode,imcon,mxnode,
     x    natms,nstep,nzden,numacc,numrdf,chip,chit,chitd,conint,
     x    tstep,buffer,cell,fxx,fyy,fzz,ravval,rdf,ssqval,
     x    stkval,stpval,sumval,vxx,vyy,vzz,xxx,yyy,zumval,
     x    zzz,xx0,yy0,zz0,zdens,xxs,yys,zzs,eta,dipx,dipy,dipz,
     x    vdxx,vdyy,vdzz,xmsd,ymsd,zmsd,lcp,ldpts,conintd)
        
      endif
c     
c     complete time check
      
      call timchk(0,timelp)
      
      newjob=.false.
     
!
!     output forces: in the case of steepest-descent like opt.

      if( lgeoopt ) then 

         call output_force_optimize( fxx, fyy, fzz,
     $        natms, nstep, atmnam )
         
         call check_conv_optimize( fxx, fyy, fzz, xxx, yyy, zzz,
     $        natms, nstep, nstrun, atmnam )

      end if
 
      if(nstep.lt.nstrun.and.timjob-timelp.gt.timcls) go to 100
      
c***********************************************************************
c     end of molecular dynamics calculations
c***********************************************************************
!
!     deallocation 
      if( lttm ) then
         call end_ttm_forces
         call end_dipole_moments
         if( ldms ) call end_intra_grad_dms
         if( .not.lcp ) then 
            deallocate( dip1x, dip1y, dip1z )
            deallocate( dip2x, dip2y, dip2z )
            deallocate( dip3x, dip3y, dip3z )
         end if
      end if

      if(idnode.eq.0)write(nrite,
     x  "(/,/,1x,'run terminating. elapsed cpu time = ',f13.3,
     x  ', job time = ',f13.3,', close time = ',f13.3,/)")
     x  timelp,timjob,timcls
      
c     
c     produce summary of simulation
      
      call result
     x  (lgofr,lpgr,lzden,cfgname,atmnam,unqatm,idnode,imcon,
     x  keyens,mxnode,natms,nzden,nstep,ntpatm,ntpvdw,numacc,
     x  numrdf,chip,chit,conint,rcut,tstep,volm,lstvdw,buffer,
     x  cell,dens,eta,fxx,fyy,fzz,ravval,rdf,ssqval,stkval,
     x  stpval,sumval,vxx,vyy,vzz,xxx,yyy,zumval,zzz,xx0,yy0,
     x  zz0,zdens,xxs,yys,zzs,dipx,dipy,dipz,lpolar,
     x  vdxx,vdyy,vdzz,xmsd,ymsd,zmsd,lcp,ldpts,chitd,conintd)

c     
c     close output channels
      
      if(idnode.eq.0) then

        close (nrite)
        close (nstats)
        close (nhist)
      
      endif
c
c     terminate job

#ifdef VAMPIR
      call VTEND(99, ierr)
#endif
      call mb_finalize()

      if(laspc) call aspc_fini()

      end
