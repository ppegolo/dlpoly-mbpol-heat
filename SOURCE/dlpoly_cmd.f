ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc, moldipx, moldipy, moldipz
c PIMD/CMD code
c
c March 2015
c

c#define dont_DO_2D_IR yes
c#define dont_DO_RECALC_DIP yes

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
     x     end_dipole_moments, dip1x, dip1y, dip1z, dip2x, dip2y, dip2z,
     x     dip3x, dip3y, dip3z, moldipx, moldipy, moldipz
      use format_history,  only: lunformat
      use param_type_ttm,  only: gammattm
      use ps_type_dms,     only: ldms, init_intra_grad_dms,
     x                           end_intra_grad_dms
      use ttm_forces,      only: lpoloff, ltip4p_geo, ldecomp,
     x                           fx_ttm_per, fy_ttm_per, fz_ttm_per,
     x                           fx_ttm_ind, fy_ttm_ind, fz_ttm_ind,
     x                           init_ttm_forces, end_ttm_forces,
     x                           vx_at_t, vy_at_t, vz_at_t
      use unit_parameters, only: eatd
      use optimize,        only: lgeoopt


      use multibead
      use centroid_dcp
      use centroid_make_cavity
      use mbpol
      use mbnrg
#ifdef HEAT_CURRENT
      use heatcurrent, only:
     x                              total_energy_ang,
     x                              total_stress_ang,
     x                              total_energy_bnd,
     x                              total_stress_bnd,
     x                              total_energy_ew1p,
     x                              total_stress_ew1p,
     x                              total_energy_ew2p,
     x                              total_stress_ew2p,
     x                              total_energy_ew3p,
     x                              total_stress_ew3p,
     x                              total_stress_qd,
     x                              total_energy_mbpol,
     x                              total_stress_mbpol,
     x                              total_energy_srf,
     x                              total_stress_srf,
     x                              total_energy_polar,
     x                              init_heat,
     x                              zero_heat,
     x                              kinetic_energy_per_atom,
     x                              potential_energy_per_atom,
     x                              potential_energy_per_atom2,
     x                              compute_heat_flux,
     x                              compute_heat_flux2,
     x                              write_heat_flux,
     x                              write_heat_flux2,
     x                              stress_per_atom,
     x                              total_energy_cpe,
     x                              total_potential_energy,
     x                              total_kinetic_energy,
     x                              final_sum,
     x                              final_sum_epot,
     x                              final_sum_heat_flux,
     x                              update_kinetic_energy, !PPnote_: remove if unnecessary
     x                              total_stress,
     x                              write_energy,
     x                              write_stress,
     x                              distribute_stress,
     x                              distribute_energy,
     x                              distribute_energy2,
     x                              write_force_matrix,
     x                              distribute_forces,
     x                              deallocate_heat,
     x                              initialize_potential_energy,
     x                              save_last_energy,
     x                              check_forces,
     x                              compute_stress
#endif /* HEAT_CURRENT */
#ifdef HEAT_CHECK
      use heat_check, only:
     x                              write_check
#endif /* HEAT_CHECK */
c
c----------------------------------------------------------------
! PIMD/CMD
      use nose_hoover_module, only:
     x                        thermo_lnv, x_lnv, x_lnv_old,
     x                        v_lnv, f_lnv, c2_lnv, mass_lnv,
     x                        Thermostat_link,
     x                        Thermostat_switch,
     x                        Thermostat_integrate,
     x                        Thermostat_hamiltonian,
     x                        Barostat_thermostat_integrate

      use centroid_module, only:
     x                     conv_energy,
     x                     beta, kT, NkT, tau_vol,
     x                     nbead, natom, nmol,
     x                     hod_in_h2o,hod_in_d2o,
     x                     d2o_in_h2o,h2o_in_d2o,
     x                     nmol_mix,
     x                     activate, nchain, nthermo, thermo,
     x                     cmd_run,
     x                     restart_run, restart_dipole,
     x                     restart_cmd,
     x                     pimd_npt,
     x                     set_vel, iseed, cmd_nve,
     x                     trpmd_run,
     x                     file_out,
     x                     file_nh, file_npt, file_stress,
     x                     file_pimd, file_cmd,
     x                     file_cv, file_force_cv,
     x                     file_pmf, file_ind,
     x                     fac_nm, w_bead,
     x                     adiab_param, fict_mass_nmode,
     x                     Epot_spring, Epot_deriv,
     x                     Equal_part, Equal_part_prim,
     x                     Etot_fict, Epot_fict, Ekin_fict,
     x                     Etot_pimd, Epot_pimd, Ekin_pimd,
     x                     Etot_prim, Ekin_prim,
     x                     Etot_cmd, Epot_cmd, Ekin_cmd,
     x                     temp_fict, temp_cmd,
     x                     Etot_nose,
     x                     P_pimd,
     x                     vir_pimd

      use aspc_module, only:
     x             aspc_init, aspc_predict, aspc_set, aspc_fini

! Langevin
      use Langevin_Thermostat, only:
     x                     langevin,evolve_in_nm,force_scale,
     x                     Langevin_init, Langevin_thermostat_integrate,
     x                     Langevin_evolve_first,Langevin_evolve_second,
     x                     Langevin_Barostat_integrate
! CV
      use CV,                  only:
     x                     CV_Ql,rmin_ql,rmax_ql,nlist_cutoff,ev_pos,
     x                     ev_kappa,nColVar,ev_force,cv_pos,Ql_En,
     x                     CV_FORCE,CV_STRESS,CV_init,cv_eval
! Opt
      use optimization,    only:
     x                     opt_init,lbfgs_init,opt_force_energy,
     x                     opt_finish
      use vibanal,    only:
     x                     vib_init,vibrational_analysis,vib_finish


c----------------------------------------------------------------

#ifdef DO_2D_IR
      use centroid_2D_IR
#endif
#ifdef DO_RECALC_DIP
      use centroid_recalc_dip
#endif

#ifdef IPI
      use ipi
#endif

      parameter (mxmc=4,mxmi=36,mxmr=89)

#include "dl_params.inc"
#include "mpif.h"

      logical lhead,pimd_head,lfirststep
      logical ltscal,lzeql,loptim,lopt2,ltraj,lgofr,lpgr,lfcap,safe,
     x        safeq
      logical newjob,newlst,lneut,loglnk,lnsq,lzden,lshmov,lcnb,lmetal
      logical safep,lpolar,lthole,lacs,lads,lcp,ldpts,lttm,lttm3,lcavity
      logical lfd, llfirst, lfdx, lfdy, lfdz, lreturn  ! for lfd
      logical flag_com          ! for COM velocity check
      logical lmsite,lmbpol,lmbnrg

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

#ifdef HEAT_CURRENT
      logical :: firststep = .true.
      ! PP_:
      real(8), allocatable :: fxx_tmp(:)
#endif /*HEAT_CURRENT*/

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
      real*8, allocatable, target :: fxx(:),fyy(:),fzz(:)
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
      real*8, allocatable, target :: vxx(:),vyy(:),vzz(:)
      real*8, allocatable :: weight(:),chgsit(:),wgtsit(:)
      real*8, allocatable :: polarsit(:),polr(:),potcc(:)
      real*8, allocatable :: polarsit2(:),polr2(:)
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
      real*8, allocatable, target :: xxx(:),yyy(:),zzz(:)
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

c PIMD/CMD: working variables.
      integer :: idim, iatom, ibead
      real*8 :: wfac, dtx, dtx2, time
      real*8 :: centvir, tot_mass
      real*8 :: aa, arg2, poly, e2, e4, e6, e8
      real*8 :: scale_volm, tmp, tmp1
      real*8 :: pos_cm_x, pos_cm_y, pos_cm_z

c PIMD/CMD VB
      real*8, allocatable, target :: npx(:),npy(:),npz(:)
      real*8, allocatable, target :: nvx(:),nvy(:),nvz(:)
      real*8, allocatable, target :: nfx(:),nfy(:),nfz(:)
c PIMD/CMD FP
      real*8 :: stress_pimd(9)

c ASPC GRM
      integer :: aspc_k_param
      integer :: aspc_iter_nstep
      real*8 :: sor_omega
      logical :: laspc

c Opt

      logical          :: lopt,lvib
      integer          :: opt_flag,vib_flag
      real(kind=8),allocatable     :: f(:),r(:)

c mbnrg

      integer, allocatable :: mbnrg_list(:,:),mbnrg_index(:)
      integer :: mbnrg_key



! Langevin
      integer ::  LangThermoType, numS_gle
      real*8 :: omega0,tau_cell,vir_pimd_nospring, alpha_cell
      real*8 :: volume,cell_mass,P_lang
      logical :: fix_com
      real*8, allocatable :: gle_vx(:),gle_vy(:),gle_vz(:)
      character*256  gle_Aunit,gle_Cunit
      character*256  gle_Afile,gle_Cfile
      real*8, pointer, dimension(:) :: ptr_xx,ptr_yy,ptr_zz
      real*8, pointer, dimension(:) :: ptr_vx,ptr_vy,ptr_vz
      real*8, pointer, dimension(:) :: ptr_fx,ptr_fy,ptr_fz
! CV
      real*8 :: cv_energy
! time
      real*8 :: timepb0,timepb1
      real*8 :: timeint,timetot

!SR  for water monomer energies -- mbpol

      real*8,allocatable         ::    wat_monomer_ener(:)


#ifdef PLUMED
      real(8) energyUnits,lengthUnits,timeUnits,dummy_vol,
     x   plumed_stress(9),plumed_KbT
      integer(8) get_comms
      logical, save :: lplumed


#endif /* PLUMED */


#ifdef DO_2D_IR
      logical continue_2D_IR
#endif
#ifdef DO_RECALC_DIP
      logical continue_RECALC_DIP
#endif

#ifdef IPI
      real(kind=8),allocatable :: com3n(:,:)
      character*1024 serveraddr
      real(8),allocatable :: indx(:), indy(:), indz(:),ipi_dipole(:,:)
#endif


#ifdef FFTW
      FFTW_PLAN_TYPE fplan, bplan
#else
      integer fplan, bplan
#endif
      data memc/mxmc*0/,memi/mxmi*0/,memr/mxmr*0/,ifail/0/
      data npage,lines/8,0/,newjob/.true./,safe/.true./safeq/.true./
      data safep/.true./
#ifdef VAMPIR
#include "VT.inc"
#endif
c
c     set up the communications

      call mb_init()
      if (ring_size.lt.1) then
         if (mb_rank.eq.0) print '(/a/)',
     x   ' ** Error ** : PIMD_CMD_NUM_BEADS is too small (or not set)'
         call mb_abort()
      end if

      if (mb_rank.eq.0) print '(/a,i3,a/)',
     x   ' ** CMD/PIMD ** : ',ring_size,' beads requested'

      timepb0=0.
      timepb1=0.
      timeint=0.
      timetot=0.


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

      pimd_head=.false.
      if (mb_rank.eq.0) pimd_head=.true.

c     open main printing file
      if(idnode.eq.0)open(nrite,file='OUTPUT'//bead_suffix)
      if(idnode.eq.0) write (nrite,
     x  "(/,20x,'DL_POLY Version 2.0 (as PIMD/CMD bead ',
     x  i3,' of ',i3,')',
     x  /,/,30x,'Running on ',i3,' nodes',/,/)")
     x  (ring_rank+1),ring_size,mxnode


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

#ifdef HEAT_CURRENT
      call init_heat()
#endif /* HEAT_CURRENT */
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
#ifdef HEAT_CURRENT
      allocate (fxx_tmp(mxatms))
#endif
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
      allocate (xxx(mxatms),yyy(mxatms),zzz(mxatms),
     x          npx(mxatms),npy(mxatms),npz(mxatms),
     x          nvx(mxatms),nvy(mxatms),nvz(mxatms),
     x          nfx(mxatms),nfy(mxatms),nfz(mxatms),stat=memr(55))
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

c  mbnrg

      allocate(mbnrg_list(mxtmls,1),mbnrg_index(mxtmls))


c
c     check integer memory allocation

      nexatm2(:) = 0

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

      call simdef_cmd
     x  (lfcap,lgofr,lnsq,loptim,lpgr,ltraj,ltscal,lzeql,lzden,lpolar,
     x  lthole,sysname,idnode,mxnode,intsta,istraj,keyens,keyfce,lopt2,lopt,lvib,
     x  lfd,keyres,keytrj,keybin,kmax1,kmax2,kmax3,multt,nstack,nstbgr,nstbpo,
     x  nhko,nlatt,nstbts,nsteql,nstraj,nstrun,lplumed,
     x  nospl,fplan,bplan,ftol,
     x  alpha,delr,epsq,fmax,press,quattol,rcut,rprim,rvdw,taup,taut,
     x  athole,athole12,athole13,
     x  athole_ion,ithole,n_ions,athole_ionwat,
     x  temp,timcls,timjob,tolnce,tstep,ffttable,
     x  lcp,ldpts,lads,lacs,lttm,lttm3,nthole,dipmas,diptmp,tautd,toler,
     x  npartem,temgap,delta,ascc,ascd,gammattm,lmsite,lmbpol,lmbnrg,
     x  nbead,restart_run,restart_dipole,restart_cmd,pimd_npt,
     x  set_vel,iseed,cmd_run,cmd_nve,adiab_param,trpmd_run,
     x  nchain,fac_nm,w_bead,
     x  laspc,sor_omega,aspc_k_param,aspc_iter_nstep,
     x  hod_in_h2o,hod_in_d2o,h2o_in_d2o,d2o_in_h2o,
     x  nmol_mix,molmix_nskip,lcavity,
     x  LangThermoType,numS_gle,omega0,
     x  fix_com,gle_Aunit,gle_Cunit,gle_Afile,gle_Cfile,omega_cell,
     x  tau_cell,CV_Ql,rmin_ql,rmax_ql,nlist_cutoff,nColVar,ev_pos,
     x  ev_kappa,serveraddr)
c  SR: add serveraddr variable for ipi
c
c     input the system force field


      call sysdef_cmd
     x  (lneut,lmetal,lnsq,molnam,mbnrg_index,sitnam,
     x  unqatm,idnode,mxnode,keyfce,
     x  keyfld,natms,ngrp,ntpatm,ntpmls,ntpvdw,ntptbp,ntpfbp,nshels,
     x  nhko,nlatt,alpha,dlrpot,drewd,engunit,prmpmf,rcut,rvdw,rcuttb,
     x  rcutfb,indpmf,keyang,keybnd,keydih,keyinv,keytet,lfzsit,listyp,
     x  lstang,lstbnd,lstcon,lstdih,lstinv,lstgst,lstvdw,lsttbp,lstfbp,
     x  lstshl,ltpsit,lsttet,ltpvdw,ltptbp,ltpfbp,npmf,nugrp,numang,
     x  numbonds,numcon,numdih,numinv,numgrp,numgsit,nummols,numsit,
     x  numpmf,numteth,numshl,chgsit,lpolar,polarsit,polarsit2,
     x  erc,fer,ercp,pmfwght,ggg,prmang,prmbnd,
     x  prmcon,prmdih,prminv,prmfld,prmtet,prmvdw,prmtbp,prmfbp,prmshl,
     x  vvv,wgtsit,rcut3b,rcut4b,buffer,ahk,hon,dhn,fon,
     x  keyumb,prmumb,lmbnrg)

      if(lmetal.and.multt.gt.1)call error(idnode,153)
c
c     construct initial configuration of system

      call sysgen_cmd
     x  (lhead,loglnk,lneut,lopt2,atmnam,cfgname,sitnam,idnode,imcon,
     x  mbnrg_list,mbnrg_nmol,mbnrg_index,mbnrg_key,
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


!SR : water monomer mbpol



!     always do iterations for 1D potential computations
#ifdef DO_2D_IR
      lcp = .false.
      continue_2D_IR = .true.
#endif
#ifdef DO_RECALC_DIP
      lcp = .false.
      continue_RECALC_DIP = .true.

        if(laspc) then

          write(nrite,*) 'Change polarization method to iterative'
          write(nrite,*) 'ASPC method should not be used for
     x  recalculating dipoles'
          write(nrite,*) 'terminating the job ............'
          goto 101
        endif


#endif

#ifdef IPI
      if (nbead .gt. 1 ) then
         write(*,*) "You must run with nbead=1"
         call exit(-1)
      endif

cSR:  For convenience, setting keybin=10
      keybin=10   ! to avoid writing output files
#endif


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
      if (lttm.or.lmsite) then

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

      endif ! lttm.or.lmsite

!SR:  for water monomer mbpol
      allocate(wat_monomer_ener(nttm2/3))

!SR:  If /* OPT */ is enabled, it should not be used to run cmd or pimd calculation

#ifdef OPT
      if (nbead .gt. 1 ) then
         write(*,*) "This is executable is only for
     x     optimization or classical MD"
         write(*,*) "You must not use this executable to run PIMD or CMD
     x     simulations"
         call exit(-1)
      endif
#endif /* OPT */


      if(lmbpol) then
        call mbpol_init(nrite,nttm2)
      endif

      if(lmbnrg) then
! MRR & DZ - Add multiple mbnrg ions
!        call mbnrg_init(nrite,nttm2,n_ions,mbnrg_key)\
!! debug
        call mbnrg_init(nrite,nttm2,n_ions,mxtmls,mbnrg_index,nummols)
! END MRR & DZ
      endif
c
c     set initial system temperature

      if(.not.lopt2) call systemp_cmd
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

      call timchk(1,tzero)
!
!    for optimization

      if( loptim ) then

         call input_opt

         if( lhead .and. lgeoopt ) call output_opt

      end if

      ibead = ring_rank + 1

c Morales: To avoid problems, I'm keeping all Nose-Hoover related calls
c          intact even when using langevin
c     NOTE NOTE NOTE: Different normalization in transformation matrix
c                     with langevin thermostat to be consistent with
c                     Michele's matrices.
      if(langevin) then
        call Langevin_init(lhead,pimd_head,ibead,mb_size,mb_rank,
     x          idnode,tstep,nbead,omega0,natms,LangThermoType,numS_gle,
     x          weight,file_out,temp,gle_Aunit,gle_Cunit,
     x          gle_Afile,gle_Cfile,press,tau_cell,omega_cell,cell_mass)

        if(pimd_npt .and. (LangThermoType.ne.4.and.LangThermoType.ne.6))
     x     call error(idnode,712)

        if(.not.pimd_npt.and.(LangThermoType.eq.4
     x                    .or.LangThermoType.eq.6))
     x     call error(idnode,713)

        alpha_cell = 0.00001d0 ! to avoid issues with divisions by 0
        vir_pimd_nospring = 0.d0

        if(LangThermoType.eq.2 .or. LangThermoType.eq.3  .or.
     x     LangThermoType.eq.5 .or. LangThermoType.eq.6) then
           ! Morales: is there a situation where (mxatom != natms) ????.
           ! Assimung no
           allocate(gle_vx((numS_gle+1)*natms),
     x              gle_vy((numS_gle+1)*natms),
     x              gle_vz((numS_gle+1)*natms))
           ! Initializing to zero now, read from file later in restart
           gle_vx(:) = 0.d0
           gle_vy(:) = 0.d0
           gle_vz(:) = 0.d0
        endif

      endif

      if(evolve_in_nm) then
        ptr_xx => npx
        ptr_yy => npy
        ptr_zz => npz
        ptr_vx => nvx
        ptr_vy => nvy
        ptr_vz => nvz
        ptr_fx => nfx
        ptr_fy => nfy
        ptr_fz => nfz
      else
        ptr_xx => xxx
        ptr_yy => yyy
        ptr_zz => zzz
        ptr_vx => vxx
        ptr_vy => vyy
        ptr_vz => vzz
        ptr_fx => fxx
        ptr_fy => fyy
        ptr_fz => fzz
      endif

#ifdef PLUMED

       if((idnode.eq.0).and.(lplumed)) write(*,*) " PLUMED IS ON ",
     x   lplumed

      if(lplumed)then
       plumed_KbT=boltz*temp

c         call plumed_f_installed(plumedavaiable)

c       if (plumedavaiable<=0) then
c        if(idnode.eq.0) write(nrite,*)"PLUMED NOT AVAILABLE"
c        stop
c       else
       call plumed_f_gcreate()
       call plumed_f_gcmd("setMPIFComm"//char(0),MPI_COMM_WORLD)
       call plumed_f_gcmd("setRealPrecision"//char(0),8)
       energyUnits=0.01        ! 10 J
       lengthUnits=0.1         !  nm
       timeUnits=1             !  ps
       call plumed_f_gcmd("setMDEnergyUnits"//char(0),energyUnits)
       call plumed_f_gcmd("setMDLengthUnits"//char(0),lengthUnits)
       call plumed_f_gcmd("setMDTimeUnits"//char(0),timeUnits)
       call plumed_f_gcmd("setPlumedDat"//char(0),"plumed.dat"//char(0))
       call plumed_f_gcmd("setLogFile"//char(0),"PLUMED.OUT"//char(0))
       call plumed_f_gcmd("setNatoms"//char(0),natms)
       call plumed_f_gcmd("setMDEngine"//char(0),"DL_POLY_TTM"//char(0))
       call plumed_f_gcmd("setTimestep"//char(0),tstep)
       call plumed_f_gcmd("setKbT"//char(0),plumed_KbT)
c       call plumed_f_gcmd("setRestart"//char(0),restart_plumed)   ! not needed
       call plumed_f_gcmd("init"//char(0),0)
c       end if
      endif

#endif /* PLUMED */


#ifdef OPTORVIB
! Opt
      if(lopt) then
!         call opt_init(natms_r,nrite,idnode,mxnode,nread)
!         call lbfgs_init()

!        allocate(f(3*natms_r),r(3*natms_r))

        if(laspc) then

          write(nrite,*) 'Change polarization method to iterative'
          write(nrite,*) 'ASPC method should not be used for geometry
     x optimization'
          write(nrite,*) 'terminating the job ............'
          goto 101
        endif

        if(.not.lvib) then
          call opt_init(natms_r,nrite,idnode,mxnode,nread)
          call lbfgs_init()

         else
          call vib_init(natms_r,nrite,idnode,mxnode,nread)
          vib_flag=0

        endif

        allocate(f(3*natms_r),r(3*natms_r))

      endif
#endif

! PIMD/CMD
c     Initialize centroid paramters.
      if (.not.lttm.and..not.lmsite) then
         natms_f = natms
         natms_r = natms
      end if
      call centroid_init_params
     x     (mxatms,mxtmls,ntpmls,nummols,listttm2(1),natms_r,
     x      molmix_nskip,weight, temp)

      if(laspc) then
          !check for errors in CONTROL file
          if(laspc.and.lcp) call error(idnode, 817)
          if(laspc.and.(.not.lpolar)) call error(idnode, 818)

          call aspc_init(aspc_k_param, aspc_iter_nstep)
      endif

c     counter for dipole guess
      if( lpolar .and. (.not. lcp) ) nstep_for_guess = 0

! CV
      if(CV_Ql) then
        call CV_init(lhead,pimd_head,natms)
      endif

! CV
      if(CV_Ql) then
        call CV_init(lhead,pimd_head,natms)
      endif
!SR


#ifdef DEBUG
      open(4444,file='DEBUG',action='write')
#endif


#ifndef IPI
      nchain = max(2, nchain)
      if (pimd_head) then
         open(file_out,file='OUTPUT_PIMD',position='append')
         write(file_out,*) '*****************************************'

#ifdef DO_2D_IR
         write(file_out,*)
     x   'PIMD/CMD (version ' // VERSION // ' : 2D-IR)'
#elif defined(DO_RECALC_DIP)
         write(file_out,*)
     x   'PIMD/CMD (version ' // VERSION // ': recalc_dip)'
#else
         write(file_out,*) 'PAESANI GROUP - DL_TTM'
#endif

         write (file_out,*) '*****************************************'

         if (cmd_run) then
            write (file_out,'(a)') 'this is a CMD run'
            write (file_out,'(a,i5)') 'number of beads =', nbead
            if (.not.restart_run) then
               write (file_out,'(a)') 'CMD needs to start from PIMD'
               stop
            else
               if (restart_cmd) then
                  write (file_out,'(a)') 'continuing a CMD trajectory'
               else
                  write (file_out,'(a)') 'starting a CMD trajectory'
               end if
            end if
            write (file_out,'(a,f8.4)') 'adiabaticity parameter =',
     x                                  adiab_param
            write (file_out,'(a,f8.4)') 'thermostat coupling =', fac_nm
            if (w_bead) then
               write (file_out,'(a)') 'beads coordinates are printed'
            end if
            write (file_out,'(a,i5,a)') 'trajectory saved every',
     x            istraj, ' steps'
         else
            write (file_out,'(a)') 'this is a PIMD run'
            write (file_out,'(a,i5)') 'number of beads =', nbead
            if (restart_run .and. .not.set_vel) then
               write (file_out,'(a)') 'restarting from a previous run'
            else
               write (file_out,'(a)') 'this is a new run'
            end if
            if (pimd_npt) then
                write (file_out,'(a)') 'NPT ensemble'
            else
                write (file_out,'(a)') 'NVT ensemble'
            end if
            if (langevin) then
               write (file_out,'(a)') 'using Langevin thermostat'
               write (file_out,'(a)') '  thermoType = ', LangThermoType
            else
               write (file_out,'(a)') 'using Nose-Hoover chains'
               write (file_out,'(a,i6)') '  number of chains = ', nchain
            end if
            write (file_out,'(a,i5,a)') 'trajectory saved every',
     x            istraj, ' steps'
         end if

         write(file_out,*)

         if (lpolar) then
           if (lcp) then
             write (file_out,'(a)') 'using Car-Parrinello for dipoles'
             if (restart_dipole) then
               write (file_out,'(a)') 'restarting from previous dipoles'
             else
               write (file_out,'(a)') 'determining new dipoles'
             end if
           else if (laspc) then
              write (file_out,'(a)') 'using ASCP method for dipoles'
              write (file_out,'(a,f12.6)') '  omega =', sor_omega
              write (file_out,'(a,i8)') '  extrapolation level =',
     x                               aspc_k_param
           else
              write (file_out,'(a)') 'using iteration for dipoles'
           end if
         end if

         if (hod_in_h2o) then
            write (file_out,*) '          '
            write (file_out,'(i4,a)') nmol_mix, ' HOD in H2O'
            write (file_out,'(a)') 'Atomic masses'
            do iatom = 1, natom
               write (file_out,'(i6,f8.4)') iatom, weight(iatom)
            end do
         else if (hod_in_d2o) then
            write (file_out,*) '          '
            write (file_out,'(i4,a)') nmol_mix, ' HOD in D2O'
            write (file_out,'(a)') 'Atomic masses'
            do iatom = 1, natom
               write (file_out,'(i6,f8.4)') iatom, weight(iatom)
            end do
         else if (h2o_in_d2o) then
            write (file_out,*) '          '
            write (file_out,'(i4,a)') nmol_mix, ' H2O in D2O'
            write (file_out,'(a)') 'Atomic masses'
            do iatom = 1, natom
               write (file_out,'(i6,f8.4)') iatom, weight(iatom)
            end do
         else if (d2o_in_h2o) then
            write (file_out,*) '          '
            write (file_out,'(i4,a)') nmol_mix, ' D2O in H2O'
            write (file_out,'(a)') 'Atomic masses'
            do iatom = 1, natom
               write (file_out,'(i6,f8.4)') iatom, weight(iatom)
            end do
         endif

         if (CV_Ql) then
            write (file_out,*)
            write (file_out,*) 'Using Collective Variables Q4/Q6.'
            write (file_out,*) 'Q4:      ',ev_pos(1)
            write (file_out,*) 'Q6:      ',ev_pos(2)
            write (file_out,*) 'Q4_KAPPA:    ',ev_kappa(1)
            write (file_out,*) 'Q6_KAPPA:    ',ev_kappa(2)
            write (file_out,*) ''
         endif

         close (file_out)

      endif

#ifdef DO_2D_IR
      call centroid_2D_IR_init(natom, continue_2D_IR)
#endif
#ifdef DO_RECALC_DIP
      call centroid_recalc_dip_init(natom, continue_RECALC_DIP)
#endif
      if(idnode.eq.0) flush(nrite)

      if(lcavity) call make_cavity_init(natom)

c     !---------------------------------------

!     timestep for Nose-Hoover chains.
      dtx = tstep
      dtx2 = dtx / 2.d0
      dtx4 = dtx / 4.d0

c     !---------------------------------------
c     using pointers to avoid multiple if/else statements below
c      if(langevin ) then
c        ptr_xx => xxx
c        ptr_yy => yyy
c        ptr_zz => zzz
c        ptr_vx => vxx
c        ptr_vy => vyy
c        ptr_vz => vzz
c        ptr_fx => fxx
c        ptr_fy => fyy
c        ptr_fz => fzz
c      else
c        ptr_xx => npx
c        ptr_yy => npy
c        ptr_zz => npz
c        ptr_vx => nvx
c        ptr_vy => nvy
c        ptr_vz => nvz
c        ptr_fx => nfx
c        ptr_fy => nfy
c        ptr_fz => nfz
c      endif

c     !---------------------------------------
      lfirststep=.true.

c     Initial PIMD/CMD configuration.
      if (.not.restart_run) then

         time = 0.d0

c        If a PIMD/CMD run starting from scratch.
c        ignore forces/velocities

         fxx(1:natom) = 0.d0
         fyy(1:natom) = 0.d0
         fzz(1:natom) = 0.d0

         dipx(1:natms_f) = 0.d0
         dipy(1:natms_f) = 0.d0
         dipz(1:natms_f) = 0.d0

         if (lpolar.and.lcp) then
            vdxx(1:natms_f) = 0.d0
            vdyy(1:natms_f) = 0.d0
            vdzz(1:natms_f) = 0.d0
            fdxx(1:natms_f) = 0.d0
            fdyy(1:natms_f) = 0.d0
            fdzz(1:natms_f) = 0.d0
         endif

      else

c        If a PIMD/CMD run starting from a previous PIMD/CMD run.

         call centroid_read_config
     x           (cfgname,nstep,time,levcfg,imcon,cell,
     x            mxatms,natom,atmnam,
     x            xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz,buffer)

         if (levcfg.ne.2) then
            if (pimd_head) then
               open(file_out,file='OUTPUT_PIMD',position='append')
               write(file_out,'(/,/,1x,a)')
     x               '** Error ** : levcfg.ne.2 for restart_run'
               close(file_out)
            end if ! pimd_head
            call mb_abort()
         end if ! levcfg.ne.2

         if (imcon.gt.0) then
            call dcell(cell,celprp)

            width = min(celprp(7),celprp(8),celprp(9))/2.d0
            if (imcon.eq.4) width = sqrt(3.d0)*cell(1)/4.d0
            if (imcon.eq.5) width = cell(1)/2.d0
            if (imcon.eq.6) width = min(celprp(7),celprp(8))/2.d0
c
c           halt program if potential cutoff exceeds cell width

            if (rcut.gt.width) then
               if (pimd_head) then
                  open(file_out,file='OUTPUT_PIMD',position='append')
                  write(file_out,'(/,/,1x,a,f6.2,a,f6.2)')
     x                  '** Error ** : potential cutoff ',rcut,
     x                  ' exceeds half-cell width ',width
                  close(file_out)
               endif ! pimd_head
               call mb_abort()
            endif ! rcut.gt.width
         endif ! imcon.gt.0

         if (lpolar.and.lcp.and.restart_dipole) then
            lfirststep = .false.
            call centroid_read_dipole
     x                (natms_f, atmnam, dipx, dipy, dipz,
     x                 vdxx, vdyy, vdzz, fdxx, fdyy, fdzz, buffer)
         endif

#ifdef HEAT_CURRENT
         nstrun=nstrun+nstep
#endif /* HEAT_CURRENT */
      endif ! .not.restart_run

      if (nstep>=nstrun) then
         if (pimd_head) then
            open(file_out,file='OUTPUT_PIMD',position='append')
            write(file_out,*) 'Nstep > nstrun: no more steps required'
            write(file_out,*) 'Stop.'
            close(file_out)
         endif
         call mb_finalize()
      endif

c     !---------------------------------------

      if (evolve_in_nm) then
c        Transform cartesian positions into normal mode positions.
         call transform_position_from_cart_to_nmode
     x        (xxx,yyy,zzz,npx,npy,npz)

c       If a restart run:
c       Transform cartesian positions into normal mode positions.
        if (restart_run)
     x     call transform_velocity_from_cart_to_nmode
     x          (vxx,vyy,vzz,nvx,nvy,nvz)

c       Transform cartesian positions into normal mode positions.
        call transform_force_from_cart_to_nmode(fxx,fyy,fzz,nfx,nfy,nfz)

c       !---------------------------------------
c       If a PIMD/CMD run starting from scratch:
c       Initialize normal mode velocities.
        if ( .not.restart_run .or. set_vel ) then
           call setup_vel_nmode(iseed,nvx,nvy,nvz)
           call fix_vel_nmode(nvx,nvy,nvz,temp)
        end if

      else

c       !---------------------------------------
c       If a PIMD/CMD run starting from scratch:
c       Initialize normal mode velocities.
        if (.not.restart_run .or. set_vel) then
           call setup_vel_nmode(iseed,vxx,vyy,vzz)
           call fix_vel_nmode(vxx,vyy,vzz,temp)
        end if
      endif

c     !---------------------------------------

      !if(langevin) then
      !  if(cmd_run) call error(idnode,711)
      !else
c       Initialize Nose'-Hoover chains.
        call centroid_init_nose_hoover(iseed,nvx,nvy,nvz)

c       For NPT PIMD.
        if ( pimd_npt ) then

           call setup_npt_pimd(iseed,tau_vol)
           e2 = 1.d0 / ( 2.d0 * 3.d0 )
           e4 = e2 / ( 4.d0 * 5.d0 )
           e6 = e4 / ( 6.d0 * 7.d0 )
           e8 = e6 / ( 8.d0 * 9.d0 )
           x_lnv = log( cell(1) * cell(5) * cell(9) ) / 3.d0
        end if
      !endif

      P_pimd = press

c     For CMD.
      if ( cmd_run .and. cmd_nve ) then

c        De-activate thermostat for path-centroid.
         do ithermo=1,size(thermo,2)
            if ( ring_rank==0 )  then
               activate = .false.
            else
               activate = .true.
            end if
            do idim = 1, 3
               call Thermostat_switch(thermo(idim,ithermo),activate)
            enddo
         enddo

c        Scale path-centroid velocity and set total momentum equal to zero.
         if (.not.restart_cmd.and.cmd_nve) then
            call scale_vel_centroid(nvx,nvy,nvz,temp)
            if (pimd_head) then
               open(file_out,file='OUTPUT_PIMD',position='append')
               write(file_out,*) 'Initial velocities have been scaled'
               write(file_out,*) '                  '
               close(file_out)
            endif
         endif

      endif

      Ekin_fict = 0.d0

      do iatom = (idnode*natom)/mxnode+1, ((idnode+1)*natom)/mxnode
         Ekin_fict = Ekin_fict
     x              + fict_mass_nmode(iatom,ibead)
     x      * ( ptr_vx(iatom)**2 + ptr_vy(iatom)**2 + ptr_vz(iatom)**2 )
      enddo

      call MPI_ALLREDUCE(Ekin_fict, buffer, 1,
     x        MPI_DOUBLE_PRECISION, MPI_SUM, comm_mb, iatom)
      Ekin_fict = buffer(1)

      tot_mass = 0.d0
      do iatom = 1, natom
         tot_mass = tot_mass + fict_mass_nmode(iatom,1)
      end do

      if ( cmd_run .and. cmd_nve ) then
         temp_fict = Ekin_fict / dble( 3*natom*nbead-3 ) / boltz
      else if ( langevin .and. fix_com ) then
         temp_fict = Ekin_fict / dble( 3*natom*nbead-3 ) / boltz
      else
         temp_fict = Ekin_fict / dble( 3*natom*nbead ) / boltz
      end if

      if (pimd_head) then
         open(file_out,file='OUTPUT_PIMD',position='append')
         write(file_out,*) 'Initial temperature = ', temp_fict
         write(file_out,*) '                      '
         close(file_out)
      end if

      call centroid_vertest
     x  (newlst,idnode,mxnode,natms_r,delr,imcon,cell,
     x   xxx,yyy,zzz,xold,yold,zold)


      if (lpolar.and.lcp) then
         if (ldpts) then
            call centroid_dcp_init(natms_f,dipmas,diptmp,tautd,polr2)
         else
            diptmp=1.0d-1
            call centroid_dcp_init
     x              (natms_f,dipmas,diptmp,diptau_max+1.0d0,polr2)
         end if ! ldpts
      end if ! lpolar.and.lcp

      if (lpolar.and..not.lcp) then
         ndipole = 0
         do i=1,natms_f
            if (polr2(i).gt.1.d-6) ndipole=ndipole+1
         end do
      end if ! lpolar.and..not.lcp
      if(idnode.eq.0) flush(nrite)
#endif  /* IPI */
c***********************************************************************
c     start of molecular dynamics calculations
c***********************************************************************
!
!     for check COM velocity
!      flag_com = .true.
!      if( nstep > 0 ) flag_com = .false.

#ifdef DO_2D_IR

      ! read first frame
      call centroid_2D_IR_read_frame(natom, nttm2, listttm2,
     x                               continue_2D_IR)

      call MPI_BCAST(continue_2D_IR, 1,
     x               MPI_LOGICAL, 0, comm_bead, ierr)

      do while(continue_2D_IR)
      do ihod = listttm2(1)+3*(nihod-1), listttm2(1)+3*(nfhod-1)!pass OW index
      !do ihod = listttm2(1)+nihod-1, listttm2(1)+nfhod-1!pass OW index
      do ivib = 1, n_vib
#endif /* DO_2D_IR */
#ifdef DO_RECALC_DIP
      ! read first frame
      call centroid_recalc_dip_read_frame(natom, nttm2, listttm2,
     x                               continue_RECALC_DIP)

      do while(continue_RECALC_DIP)
#endif /* DO_RECALC_DIP */

#ifdef IPI
        nstep=0      !for safety

! for dipole
      ndummy = nummols(ntpmls)

      allocate(indx(ndummy),indy(ndummy),indz(ndummy))
      allocate(ipi_dipole(ndummy,3))
#endif /* IPI */


#ifdef OPTORVIB
c Opt
      if(lopt) then
         nstep=0
!        ndim=3*natms_r+2000
!        MSAVE=7
!        nwork=ndim*(2*msave +1)+2*msave
!
!        allocate(diag(ndim),w(nwork))
      endif
#endif

  100 continue


!PIMD/CMD

c     increase step counter
      nstep = nstep + 1
      time = time + tstep

      iatm1 = (idnode*natms)/mxnode+1
      iatm2 = ((idnode+1)*natms)/mxnode

#ifdef HEAT_CURRENT
      call zero_heat()
#endif

#ifdef OPTORVIB
c vibrational analysis

      if(nstep>1) then


        do i=1,natms_r
          r(3*i-2)=xxx(i)
          r(3*i-1)=yyy(i)
          r(3*i)=zzz(i)

          f(3*i-2)=-fxx(i)
          f(3*i-1)=-fyy(i)
          f(3*i)=-fzz(i)
        enddo

!        write(*,*) "EPOT",Epot_pimd/418.d0
        pot=Epot_pimd/engunit    ! sending value in kcal/mol

c  vibrational analysis

        if(lopt.and.lvib) then

          call vibrational_analysis(atmnam,weight,r,f,pot,vib_flag)

           if(vib_flag<0) goto 101

c   Opt
         elseif(lopt.and..not.lvib) then

          call opt_force_energy(atmnam,r,f,pot,opt_flag)

c writing out the stress tensor useful to cell parameter optimization

          if(opt_flag<=0) then
             P_pimd = NkT/volm - vir_pimd/volm/3.d0
             open(file_stress,file='STATIS_STRESS',position='append')
             if(idnode==0) then
             write(file_stress,*) 'writing stress_tensor
     x  at last optimization step'
             write(file_stress,'(f12.5,10f10.5)')
     x               time, P_pimd*prsunt, stress_pimd(1:9)
             endif
            close(file_stress)


          endif

          if(opt_flag<0)  goto 101     ! termination error in optimization
          if(opt_flag==0) goto 101     ! geometry optimization successful

        endif

        do i=1,natms_r
          xxx(i)=r(3*i-2)
          yyy(i)=r(3*i-1)
          zzz(i)=r(3*i)
        enddo

          call transform_position_from_cart_to_nmode
     x        (xxx,yyy,zzz,npx,npy,npz)
      endif

#endif /* OPTORVIB */


#ifdef IPI

      ! OPENS SOCKET
       CALL IPI_INIT(serveraddr, idnode)
      ! DOES ONE "STEP"
       pot=Epot_pimd
       stress=stress_pimd
!       pot=(engsrp+engcpe+engbnd+engang+engdih+engfld+engtbp+
!     x         engfbp+engshl+enginv+engter+engmet)

       ! you must get the forces and put them in a 3,n array

       if (hasdata) then
!      write(2599,*) natms_r
!      write(2599,*) "FRAME",nstep
!      do iloop=1,natms_r
!        write(2599,*)   fxx(iloop),fyy(iloop),fzz(iloop)
!      enddo

       do iloop=1,natms_r
           com3n(1,iloop) = fxx(iloop)
           com3n(2,iloop) = fyy(iloop)
           com3n(3,iloop) = fzz(iloop)
       enddo
       endif
       CALL IPI_STEP(volm, pot, stress, cell, com3n, ipi_dipole)

!       if(idnode==0) write(*,*) 'dlpoly','output step',nstep
       if (.not. imglock) then
!          if (ionode) write(*,*) "Resetting ASPC history"
          nstep_for_guess = 0  ! if the image has changed, we must recompute stuff as this was the first step
       endif

       do iloop=1,natms_r
           xxx(iloop)=com3n(1,iloop)
           yyy(iloop)=com3n(2,iloop)
           zzz(iloop)=com3n(3,iloop)
       enddo
        call transform_position_from_cart_to_nmode
     x        (xxx,yyy,zzz,npx,npy,npz)

!      write(2699,*) natms_r
!      write(2699,*) "FRAME",nstep
!      do iloop=1,natms_r
!        write(2699,*)   xxx(iloop),yyy(iloop),zzz(iloop)
!      enddo


#endif /* IPI */


#ifdef OPTORVIB

! Opt
      if(lopt) goto 2501

#endif /* OPTORVIB */

#ifndef IPI
c
c Velocity-Verlet propagation + Nose'-Hoover chains.
c
c     Thermostat integration.
c     N.B. If nstep > nequilib_cmd the thermostat for path-centroid
c     ( ibead = 1 ) is de-activated since it must strictly follow
c     Newton equation of motions.


      if(langevin) then

        call Langevin_thermostat_integrate(idnode,ibead,iatm1,iatm2
     x       ,nvx,nvy,nvz,vxx,vyy,vzz,
     x       gle_vx,gle_vy,gle_vz,fix_com,fict_mass_nmode,alpha_cell)

      else
        if (pimd_npt.and.mb_rank.eq.0) then
           call Barostat_thermostat_integrate(nchain,dtx2) ! T-baro
        endif
        call Thermostat_integrate(nchain,thermo,nthermo,dtx2)
      endif

c     Velocity-Verlet, 1st step:
c     update normal mode full timestep positions and
c     half timestep velocities.

      if(langevin) then

        volume = volm

        call Langevin_evolve_first(idnode,ibead,iatm1,iatm2,
     x       xxx,yyy,zzz,npx,npy,npz,nvx,nvy,nvz,vxx,vyy,vzz,
     x       fxx,fyy,fzz,nfx,nfy,nfz,alpha_cell,volume,P_lang
     x             ,vir_pimd_nospring,Epot_deriv )

        if(pimd_npt) then
          cell(1:9) = cell(1:9) * (volume / volm)**(1.0/3.0)

          if (nbead.gt.1) then
             call MPI_BCAST(cell, 9, MPI_DOUBLE_PRECISION,
     x                     0, comm_mb, ierr)
          endif ! nbead.gt.1

          volm = cell(1) * cell(5) * cell(9)
        endif

      elseif ( pimd_npt ) then

         Ekin2=0.d0
         do iatom=iatm1,iatm2
            Ekin2=EKin2+fict_mass_nmode(iatom,1)
     x   *(nvx(iatom)**2+nvy(iatom)**2+nvz(iatom)**2)
         end do
         call gdsum(Ekin2,1,buffer)

         if (mb_rank.eq.0) then
            ! wrong at the 1st step as P_pimd=press set above
            f_lnv = (c2_lnv-1.d0)*Ekin2+3.d0*volm*(P_pimd-press)
            v_lnv = v_lnv + dtx2*f_lnv/mass_lnv
         endif

         if (mxnode.gt.1.and.ibead.eq.1)
     x      call MPI_BCAST(v_lnv,1,MPI_DOUBLE_PRECISION,
     x                     0,comm_bead,ierr)

         ! for position
         aa = exp( dtx2 * v_lnv )
         arg2 = v_lnv * dtx2 * v_lnv * dtx2
         poly = 1.d0 + arg2 *
     x      ( e2 + arg2 * ( e4 + arg2 * ( e6 + arg2 * e8 ) ) )

         ! for velocity
         aav = exp(-c2_lnv*dtx4*v_lnv)
         arg2v = (c2_lnv*dtx4*v_lnv)*(c2_lnv*dtx4*v_lnv)
         polyv = 1.d0 + arg2v *
     x      ( e2 + arg2v * ( e4 + arg2v * ( e6 + arg2v * e8 ) ) )

         if (ibead.eq.1) then
            do iatom = iatm1, iatm2
               wfac = polyv*dtx2 / fict_mass_nmode(iatom,1)
               nvx(iatom) = aav*(nvx(iatom)*aav + nfx(iatom)*wfac)
               npx(iatom) = aa*(npx(iatom)*aa+nvx(iatom)*poly*dtx)
               nvy(iatom) = aav*(nvy(iatom)*aav + nfy(iatom)*wfac)
               npy(iatom) = aa*(npy(iatom)*aa+nvy(iatom)*poly*dtx)
               nvz(iatom) = aav*(nvz(iatom)*aav + nfz(iatom)*wfac)
               npz(iatom) = aa*(npz(iatom)*aa+nvz(iatom)*poly*dtx)
            enddo

         else ! ibead.ne.1

            do iatom = iatm1, iatm2
               wfac = dtx2 / fict_mass_nmode(iatom,ibead)
               nvx(iatom) = nvx(iatom) + nfx(iatom)*wfac
               npx(iatom) = npx(iatom) + nvx(iatom)*dtx
               nvy(iatom) = nvy(iatom) + nfy(iatom)*wfac
               npy(iatom) = npy(iatom) + nvy(iatom)*dtx
               nvz(iatom) = nvz(iatom) + nfz(iatom)*wfac
               npz(iatom) = npz(iatom) + nvz(iatom)*dtx
            enddo

         endif ! ibead.eq.1

         x_lnv_old = x_lnv
         x_lnv = x_lnv_old + v_lnv * dtx
         scale_volm = exp( x_lnv - x_lnv_old )
         cell(1:9) = cell(1:9) * scale_volm

         if (nbead.gt.1) then
            call MPI_BCAST(cell, 9, MPI_DOUBLE_PRECISION,
     x                     0, comm_mb, ierr)
         endif ! nbead.gt.1

         volm = cell(1) * cell(5) * cell(9)

      else ! not pimd_npt

         if (ibead.eq.1) then
            do iatom = iatm1, iatm2
               wfac = dtx2 / fict_mass_nmode(iatom,1)
               nvx(iatom) = nvx(iatom) + nfx(iatom)*wfac
               npx(iatom) = npx(iatom) + nvx(iatom)*dtx
               nvy(iatom) = nvy(iatom) + nfy(iatom)*wfac
               npy(iatom) = npy(iatom) + nvy(iatom)*dtx
               nvz(iatom) = nvz(iatom) + nfz(iatom)*wfac
               npz(iatom) = npz(iatom) + nvz(iatom)*dtx
            enddo

         else ! ibead.ne.1

            do iatom = iatm1, iatm2
               wfac = dtx2 / fict_mass_nmode(iatom,ibead)
               nvx(iatom) = nvx(iatom) + nfx(iatom)*wfac
               npx(iatom) = npx(iatom) + nvx(iatom)*dtx
               nvy(iatom) = nvy(iatom) + nfy(iatom)*wfac
               npy(iatom) = npy(iatom) + nvy(iatom)*dtx
               nvz(iatom) = nvz(iatom) + nfz(iatom)*wfac
               npz(iatom) = npz(iatom) + nvz(iatom)*dtx
            enddo
         endif ! ibead.eq.1
      endif ! pimd_npt

      if ( .not.cmd_run ) then

!         if(langevin.and..not.evolve_in_nm) then
!           call transform_position_from_cart_to_nmode
!     x          (xxx,yyy,zzz,npx,npy,npz)
!         endif
         if (ibead.eq.1) then
           if(evolve_in_nm.and..not.langevin) then
            call images(imcon,idnode,mxnode,natom,cell,npx,npy,npz)
!           else
! M Morales: Keeping this until I understand why it is here!
!            npx(:) = npx(:)/sqrt(dble(nbead))
!            npy(:) = npy(:)/sqrt(dble(nbead))
!            npz(:) = npz(:)/sqrt(dble(nbead))
!            call images(imcon,idnode,mxnode,natom,cell,npx,npy,npz)
!            npx(:) = npx(:)*sqrt(dble(nbead))
!            npy(:) = npy(:)*sqrt(dble(nbead))
!            npz(:) = npz(:)*sqrt(dble(nbead))
           endif
         endif ! ibead.eq.1
!         if(langevin.and..not.evolve_in_nm) then
!           call transform_position_from_nmode_to_cart
!     x          (npx,npy,npz,xxx,yyy,zzz)
!         endif
      endif


#endif  /* IPI */

! Opt
2501     continue

!FP: check
!      write(*,*) 'after COM'

c     !---------------------------------------

      if(evolve_in_nm) then
c       Transform normal mode positions into cartesian positions.
        call transform_position_from_nmode_to_cart
     x                   (npx,npy,npz,xxx,yyy,zzz)
      else
        ! not sure this is needed
        call transform_position_from_cart_to_nmode
     x          (xxx,yyy,zzz,npx,npy,npz)
      endif

      call merge(idnode,mxnode,natom,mxbuff,xxx,yyy,zzz,buffer)

!FP: check
!      write(*,*) 'after POS NMODE_TO_CART'

c     Update Verlet list
      call centroid_vertest
     x     (newlst,idnode,mxnode,natms_r,delr,imcon,cell,
     x      xxx,yyy,zzz,xold,yold,zold)

#ifdef DO_2D_IR
      newlst = .true.
      call centroid_2D_IR_set(cell,xxx,yyy,zzz,ihod,ivib)
#endif /* DO_2D_IR */
#ifdef DO_RECALC_DIP
      newlst = .true.
      call centroid_recalc_dip_set(cell,xxx,yyy,zzz)
#endif /* DO_RECALC_DIP */

!FP: check
!      write(*,*) 'after VERTEST'


!FP: check
!         write(*,*) 'before dcp_h0_1st'

c        1st step of velocity-Verlet + Nose'-Hoover for dipoles.
      if (lpolar.and.lcp) then
         call centroid_dcp_1st(polr2,dipx,dipy,dipz,
     x         vdxx,vdyy,vdzz,fdxx,fdyy,fdzz,tstep)
         call merge(idnode,mxnode,natms_f,
     x              mxbuff,dipx,dipy,dipz,buffer)
      endif
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

      if(lttm.or.lmsite) then

#ifdef TTM_FORCE_DECOMPOSITION
        ntmp = 1
        call init_ttm_forces(ntmp)
#endif /* TTM_FORCE_DECOMPOSITION */

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

!SR : for water monomer mbpol
      wat_monomer_ener=0.d0

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

c     zero long range component of stress
      do i = 1,9
        stresl(i) = 0.d0
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

         do i=1,natms_f

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

      endif ! lpolar.and.lcp
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

      endif ! lnewlst

c
c     calculate pair forces, including coulombic forces

#ifdef DEBUG
      call timchk(0,timeint)
      timetot=timeint
#endif
      do i = 1,9
        stress(i) = stresl(i)
      enddo

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
          if( lpolar .and. (.not.lcp) ) then
             nstep_for_guess = nstep_for_guess + 1

             if(laspc)then
                 call aspc_predict(nstep_for_guess,natms,
     x                             dipx,dipy,dipz)
             else
                 ! for some reason, the dipole guess is disabled
                 call set_guess_dipole(1,natms,
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
     x  lpolar,lthole,athole,athole12,athole13,
     x  athole_ion,ithole,n_ions,athole_ionwat,
     x  dipx,dipy,dipz,
     x  emux,emuy,emuz,nthole,ascc,ascd,efieldkx,efieldky,efieldkz,
     x  efdcrecx,efdcrecy,efdcrecz,
     x  efddmurecx,efddmurecy,efddmurecz,potcc,keyres,
     x  n_water,gammattm,restart_dipole,laspc,sor_omega,aspc_iter_nstep,
     x  nstep_for_guess)
c    Uncomment to print dipoles for each bead
c       if(lhead) then
c         write(4000+ibead,'(2i4,f15.6)') ibead,
c    x        idnode, dipx(natms)
c         write(2000+ibead,'(2i4,f15.6)') ibead,
c    x        idnode, dipy(natms)
c         write(3000+ibead,'(2i4,f15.6)') ibead,
c    x        idnode, dipz(natms)
c       endif

          if(lpolar .and. laspc .and. (.not.lcp)) then
              call aspc_set(natms, dipx, dipy, dipz)
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
#ifdef DEBUG
      call timchk(0,timeint)
      write(4444,"(a,F10.3,a)") 'Time taken for electrostatics    = ',timeint-timetot,' sec'
      write(4444,*)
#endif

      ! for dipoles output in CP case
      if (lfirststep.and.lttm.and.lcp) then
         lfirststep = .false.
         do i = 1, natms_f
            dipx_old(i) = dipx(i)
            dipy_old(i) = dipy(i)
            dipz_old(i) = dipz(i)
         end do
      endif

c
c     add in long range corrections to energy and pressure
      call lrcorrect
     x  (idnode,imcon,keyfce,mxnode,natms,ntpatm,elrc,engunit,
     x   virlrc,rvdw,volm,lstvdw,ltpvdw,ltype,numtyp,numfrz,
     x   lstfrz,prmvdw,dens)
      engsrp = engsrp + elrc + elrcm(1)
      virsrp = virsrp + virlrc + vlrcm(1)

!VB: at that point forces differ from pimdpol because here
!    qdforces()'s was already invoked; comment it in forces.f for
!    check

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

#ifdef DEBUG
      call timchk(0,timeint)
      timetot=timeint
#endif


      if (ntbond.gt.0) call bndfrc
     x  (idnode,imcon,mxnode,ntbond,engbnd,virbnd,keybnd,listbnd,
     x  cell,fxx,fyy,fzz,prmbnd,xxx,yyy,zzz,xdab,ydab,zdab,stress,
     x  buffer,lttm,wat_monomer_ener)

c
c     calculate valence angle forces

      if (ntangl.gt.0) call angfrc
     x  (idnode,imcon,mxnode,ntangl,engang,virang,keyang,listang,
     x  cell,fxx,fyy,fzz,prmang,xxx,yyy,zzz,xdab,ydab,zdab,xdbc,
     x  ydbc,zdbc,stress,buffer,lttm,wat_monomer_ener)


!
!     intramolecular potential for TTM2 model
!     if( lttm .and. ntangl.gt.0 ) then
!        open(nintra,file='UINTRA',position='append')
!        engintra = engbnd + engang
!        if(idnode==0) write(nintra,'(i7,f20.10)') nstep, engintra / engunit
!        close(nintra)
!     endif

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

#ifdef DEBUG

      call timchk(0,timeint)

      write(4444,"(a,F10.3,a)") 'Time taken for intramolecular part (bond, angle,
     x     dihedral, improper, inversion)    =   ',timeint-timetot,' sec'
      write(4444,*)

#endif
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

      if(lcavity) call make_cavity_frc
     x  (natom, nttm2, listttm2, cell,
     x   xxx,yyy,zzz,fxx,fyy,fzz,engsrp,virsrp)

c      if(idnode==0) write(*,*) 'Nframe',nstep

      ! NOTE NOTE NOTE FP: added stress
      if(lmbpol) then
#ifdef DEBUG
      call timchk(0,timeint)
      timetot=timeint
#endif
         call mbpol_forces(imcon,cell,nttm2,listttm2,wat_monomer_ener,
     x                     xxx,yyy,zzz,fxx,fyy,fzz,engsrp,virsrp,stress)
#ifdef DEBUG
      call timchk(0,timeint)
      write(4444,"(a,F10.3,a)") 'Time taken for MB-pol 2B+3B    = ',timeint-timetot,' sec'
      write(4444,*)
#endif
      endif



c mbnrg temporary

      if(lmbnrg) then
#ifdef DEBUG
      call timchk(0,timeint)
      timetot=timeint
#endif
         call mbnrg_forces(imcon,cell,nttm2,listttm2,n_ions,
     x                     mbnrg_list(1,:),
     x                     xxx,yyy,zzz,fxx,fyy,fzz,engsrp,virsrp,stress)
#ifdef DEBUG
      call timchk(0,timeint)
      write(4444,"(a,F10.3,a)") 'Time taken for MB-nrg,
     x  both 2B+3B    = ',timeint-timetot,' sec'
      write(4444,*)
#endif
      endif


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

      if (lttm.or.lmsite) then
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
#ifdef HEAT_CURRENT
          call distribute_stress(ioxy,ih1,ih2,mttm2,gammattm)
          call distribute_energy(ioxy,ih1,ih2,mttm2,gammattm)
          call distribute_forces(ioxy,ih1,ih2,mttm2,gammattm)
#endif /* HEAT_CURRENT */

        enddo

      endif

!VB: at that point forces match with pimdpol

c
c     cap forces in equilibration mode

      if(nstep.le.nsteql.and.lfcap)
     x  call fcap(lfcap,idnode,mxnode,natms,fmax,temp,fxx,fyy,fzz)
c
c     total virial (excluding constraint virial and c.o.m virial)
c      for npt routines     note: virsrp already includes virlrc

      virtot = vircpe+virsrp+virbnd+virtbp+virfld+virang+virshl+virtet

c
c        2nd step of velocity-Verlet + Nose'-Hoover for dipoles.
      if (lpolar.and.lcp) then
         call centroid_dcp_2nd
     x      (polr2,dipx,dipy,dipz,vdxx,vdyy,vdzz,
     x             fdxx,fdyy,fdzz,emux,emuy,emuz,
     x             efieldkx,efieldky,efieldkz,
     x             efdcrecx,efdcrecy,efdcrecz,
     x             efddmurecx,efddmurecy,efddmurecz,engdke,tstep)
      endif
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
     x  xx0,yy0,zz0,weight,stress,amsd,engdke,consvd,stress_pimd)

#if 0
      write(100+mb_rank,*) 'engbnd+engang=', (engbnd+engang)/engunit
      write(100+mb_rank,*) 'engsrp=', engsrp/engunit
      write(100+mb_rank,*) 'engcpe=', engcpe/engunit
      write(100+mb_rank,*) 'engcfg=', stpcfg/engunit
      write(100+mb_rank,*) 'virtot=', virtot

      if (bead_rank.eq.0) then
      omyvir=0.d0
!     derivatives:
      do iatom=1,natms_r
         write(100+mb_rank,'(i3,1x,3(1x,ES16.9))'),iatom,
     x -fxx(iatom)/engunit,-fyy(iatom)/engunit,-fzz(iatom)/engunit
      omyvir=omyvir-fxx(iatom)*xxx(iatom)-fyy(iatom)*yyy(iatom)
     x -fzz(iatom)*zzz(iatom)
      end do
      write(100+mb_rank,*)'omyvir=',omyvir
      end if
#endif

c CV
      if(CV_Ql) then
        call CV_EVAL(cell,xxx,yyy,zzz,cv_energy)
        do iatom = 1, natms_r
           fxx(iatom) = fxx(iatom) +  CV_FORCE(1,iatom)
           fyy(iatom) = fyy(iatom) +  CV_FORCE(2,iatom)
           fzz(iatom) = fzz(iatom) +  CV_FORCE(3,iatom)
        enddo
        ! must modify virial! and check stress convention  (divide by V
        ! or not)
        ! also need to add cv_energy to total energy
        stress(1:9) = stress(1:9) + CV_STRESS(1:9)
      endif

! PIMD/CMD

      do iatom = 1, natms_r
         fxx(iatom) = fxx(iatom) / dble(nbead) * force_scale
         fyy(iatom) = fyy(iatom) / dble(nbead) * force_scale
         fzz(iatom) = fzz(iatom) / dble(nbead) * force_scale
      enddo
c
c     Quantum potential energy.
      Epot_pimd = stpcfg / dble(nbead)

c     Virial and pressure.
      vir_pimd = stpvir / dble(nbead)

      buffer(1) = Epot_pimd
      buffer(2) = vir_pimd

      call MPI_ALLREDUCE(buffer(1),buffer(3),2,MPI_DOUBLE_PRECISION,
     x                   MPI_SUM, comm_ring, ierr)

      Epot_pimd = buffer(3)
      vir_pimd = buffer(4)
      vir_pimd_nospring = vir_pimd


!FP: added
c     Stress tensor
      stress_pimd(1:9) = stress_pimd(1:9) / dble(nbead)

      buffer(1:9) = stress_pimd(1:9)

      call MPI_ALLREDUCE(buffer(1),buffer(10),9,MPI_DOUBLE_PRECISION,
     x                   MPI_SUM, comm_ring, ierr)

      stress_pimd(1:9) = buffer(10:18)


#ifdef IPI
!SR
! passing the dipole moments of water molecules
!      --- initialize induced dipole moments ---
!

      indx(:) = 0.0d0;  indy(:) = 0.0d0;  indz(:) = 0.0d0


!
! --- output water dipole moments ---
!
      call merge(idnode,mxnode,natms_f,mxbuff,dipx,dipy,dipz,buffer)

      mttm2 = listttm2(nttm2)

      do i=1,ndummy

         mttm2 = mttm2 + 1

         iox = listttm2(3*i-2)
         ih1 = listttm2(3*i-1)
         ih2 = listttm2(3*i)
!        imm = n_water*3 + i
         imm = listttm2(3*ndummy) + i

       indx(i) = indx(i) + dipx(iox) + dipx(ih1) + dipx(ih2) + dipx(imm)
       indy(i) = indy(i) + dipy(iox) + dipy(ih1) + dipy(ih2) + dipy(imm)
       indz(i) = indz(i) + dipz(iox) + dipz(ih1) + dipz(ih2) + dipz(imm)

         ipi_dipole(i,1)=moldipx(i)+indx(i)
         ipi_dipole(i,2)=moldipy(i)+indy(i)
         ipi_dipole(i,3)=moldipz(i)+indz(i)

         ipi_dipole(i,:)=ipi_dipole(i,:)*eatd
      end do

!       write(327222,*) 'dip',ipi_dipole(i,1),ipi_dipole(i,2),ipi_dipole(i,3)


      if(stpvol<1.d-3) stpvol=1.d0          ! takes care of cluster simulations
      stress_pimd=stress_pimd*stpvol/prsunt  ! converting back to dlpoly internal units

      hasdata=.true.
#endif /* IPI */


c  Opt
       if(lopt) goto 100

!SR: Don't need dynamic part for IPI
#ifndef IPI

!FP: added
#if 0
      open(100+mb_rank, position='append')
      write(100+mb_rank,'(3f15.8)') stress_pimd(1:3)
      write(100+mb_rank,'(3f15.8)') stress_pimd(4:6)
      write(100+mb_rank,'(3f15.8)') stress_pimd(7:9)
      close(100+mb_rank)
!FP: end
#endif


c     Add harmonic forces between beads.
c     Compute quantum energies.
c     Compute fictitious potential energy.
      if ( nbead>1 ) then
         call centroid_spring_force
     x        (xxx,yyy,zzz,npx,npy,npz,fxx,fyy,fzz,
     x         Epot_spring,Epot_deriv)

         Ekin_pimd = Equal_part + Epot_deriv
         Ekin_prim = Equal_part_prim - Epot_spring
         Etot_pimd = Ekin_pimd + Epot_pimd
         Etot_prim = Ekin_prim + Epot_pimd
         Epot_fict = Epot_pimd + Epot_spring
      else
         Epot_fict = Epot_pimd
      endif


#ifdef PLUMED


      if(lplumed) then
c     sandeep
c       if (pimd_head) then
c        write(501,*) "FORCES -- Before PLUMED"
c        do i=1,6
c          write(501,"(3F15.8)") fxx(i),fyy(i),fzz(i)
c        enddo
c
c        write(502,*) "Energy -- Before PLUMED"
c          write(502,*) Epot_pimd
c
c        write(503,*) "Stress -- Before PLUMED"
c          write(503,"(3F15.8)") stress_pimd
c       endif


c  sandeep: Modifies forces and stress tensor in general
c  if 'ENERGY' is used as CV, then it modifies PE value as well
c sandeep: Change stress tensor for PLUMED

         dummy_vol=stpvol
         if(stpvol<1.d-3) dummy_vol=1.d0            ! sandeep: actual value should be zero, but changed it to account for numerical errors
c      stress_pimd=stress_pimd/prsunt
         stress_pimd=stress_pimd*dummy_vol/prsunt      ! converting back to internal units
         plumed_stress=stress_pimd
         stress_pimd=-stress_pimd                    ! PLUMED expects it with -ve sign

         call plumed_f_gcmd("setStep"//char(0),nstep)
         call plumed_f_gcmd("setMasses"//char(0),weight)
         call plumed_f_gcmd("setCharges"//char(0),chge)
         call plumed_f_gcmd("setPositionsX"//char(0),xxx)
         call plumed_f_gcmd("setPositionsY"//char(0),yyy)
         call plumed_f_gcmd("setPositionsZ"//char(0),zzz)
         call plumed_f_gcmd("setBox"//char(0),cell)
         call plumed_f_gcmd("setEnergy"//char(0),Epot_pimd)   ! stpcfg or Epot_pimd  -sandeep
         call plumed_f_gcmd("setForcesX"//char(0),fxx)
         call plumed_f_gcmd("setForcesY"//char(0),fyy)
         call plumed_f_gcmd("setForcesZ"//char(0),fzz)
         call plumed_f_gcmd("setVirial"//char(0),stress_pimd)   ! changed -sandeep

         call plumed_f_gcmd("calc"//char(0) )
c    correct for virial -- to get correct pressure -- virtot and stpvir
c    ???   virtot  or vir_pimd
      plumed_stress=-plumed_stress-stress_pimd  ! new - old

      vir_pimd=vir_pimd-(plumed_stress(1)+plumed_stress(5)+
     x   plumed_stress(9))       !  not diving by 3V   ! vir contribution added
c     here, negative sign (though we are adding) used because vir_pimd=-ve of actual vir_pimd
c      virtot=virtot-((plumed_stress(1)+plumed_stress(5)+
c     x   plumed_stress(9))/3.d0)       !  assuming 3 D of the system

c   rescaling the stress back to original    : also to katm
      stress_pimd=-stress_pimd*prsunt/(dummy_vol)



c     if (pimd_head) then
c      write(501,*) "Forces -- After PLUMED"
c      do i=1,6
c        write(501,"(3F15.8)") fxx(i),fyy(i),fzz(i)
c      enddo
c
c      write(502,*) "Energy -- After PLUMED"
c        write(502,*) Epot_pimd
c
c      write(503,*) "Stress -- After PLUMED"
c        write(503,"(3F15.8)") stress_pimd
c     endif

      endif
c     PLUMED

#endif /* PLUMED */




      if(evolve_in_nm) then
c       Transform cartesian forces into normal mode forces.
        call transform_force_from_cart_to_nmode
     x   (fxx,fyy,fzz,nfx,nfy,nfz)
      endif

c     !---------------------------------------

c
c     Velocity-Verlet, 2nd step:
c     Update normal mode full timestep velocities.

      P_pimd = NkT/volm - vir_pimd/volm/3.d0
!     if (bead_rank.eq.0.and.mod(nstep,50).eq.0) then
!        write(ibead+99,*) nstep, P_pimd*prsunt
!     end if

      if(langevin) then

        ! to get the pressure in NVT
        Ekin2 = 2.0*Equal_part/mxnode

        ! testing
         Ekin2 = 0.d0
         if (ibead.eq.1) then
           do iatom = iatm1, iatm2
             Ekin2 = Ekin2 + fict_mass_nmode(iatom,1)
     x   *(nvx(iatom)**2+nvy(iatom)**2+nvz(iatom)**2)
           enddo
         endif
         Ekin2 = Ekin2 / force_scale

        volume = volm

        call Langevin_evolve_second(idnode,ibead,iatm1,iatm2,
     x       xxx,yyy,zzz,npx,npy,npz,nvx,nvy,nvz,vxx,vyy,vzz,
     x       fxx,fyy,fzz,nfx,nfy,nfz,alpha_cell,volume,P_lang
     x             ,vir_pimd_nospring,Epot_deriv )

         if(pimd_npt) then
           cell(1:9) = cell(1:9) * (volume / volm)**(1.0/3.0)

           if (nbead.gt.1) then
              call MPI_BCAST(cell, 9, MPI_DOUBLE_PRECISION,
     x                     0, comm_mb, ierr)
           endif ! nbead.gt.1

           volm = cell(1) * cell(5) * cell(9)
        endif

      elseif (pimd_npt) then

         Ekin2 = 0.d0

         ! for velocity
         aav = exp(-c2_lnv*dtx4*v_lnv)
         arg2v = (c2_lnv*dtx4*v_lnv)*(c2_lnv*dtx4*v_lnv)
         polyv = 1.d0 + arg2v *
     x      ( e2 + arg2v * ( e4 + arg2v * ( e6 + arg2v * e8 ) ) )

         if (ibead.eq.1) then
            do iatom = iatm1, iatm2
               wfac = polyv*dtx2 / fict_mass_nmode(iatom,1)
               nvx(iatom) = aav*(nvx(iatom)*aav + nfx(iatom)*wfac)
               nvy(iatom) = aav*(nvy(iatom)*aav + nfy(iatom)*wfac)
               nvz(iatom) = aav*(nvz(iatom)*aav + nfz(iatom)*wfac)
               Ekin2 = Ekin2 + fict_mass_nmode(iatom,1)
     x               *(nvx(iatom)**2+nvy(iatom)**2+nvz(iatom)**2)
            enddo
         else ! ibead.ne.1
            do iatom = iatm1, iatm2
               wfac = dtx2 / fict_mass_nmode(iatom,ibead)
               nvx(iatom) = nvx(iatom) + nfx(iatom)*wfac
               nvy(iatom) = nvy(iatom) + nfy(iatom)*wfac
               nvz(iatom) = nvz(iatom) + nfz(iatom)*wfac
            enddo
         end if ! ibead.eq.1

      else ! not NPT

         ! to get the pressure in NVT
         Ekin2 = 2.0*Equal_part/mxnode
         do iatom = iatm1, iatm2
            wfac = dtx2 / fict_mass_nmode(iatom,ibead)
            nvx(iatom) = nvx(iatom) + nfx(iatom)*wfac
            nvy(iatom) = nvy(iatom) + nfy(iatom)*wfac
            nvz(iatom) = nvz(iatom) + nfz(iatom)*wfac
         enddo

      endif ! pimd_npt/langevin


      P_pimd = 0.d0

c     PIMD: pressure.
      if (nbead.gt.1) then
         ! VB: this is the centroid virial; see vir_pimd above
         ! VB: and Epot_deriv in centroid_spring_force.f
         vir_pimd = vir_pimd - 2.d0*Epot_deriv
      end if ! nbead.gt.1

c Morales: I still want to know the pressure even if NVT
      if (ibead.eq.1) then
         call gdsum(Ekin2,1,buffer)

         ! VB: this is the correct estimator for this integrator
         P_pimd = (Ekin2 - vir_pimd)/volm/3.d0
      endif

      if (pimd_npt.and.ibead.eq.1) then

         if (pimd_npt) then
            if (mb_rank.eq.0) then
               f_lnv = (c2_lnv-1.d0)*Ekin2+3.d0*volm*(P_pimd-press)
               v_lnv = v_lnv + dtx2*f_lnv/mass_lnv
            end if
         endif


c         if (mxnode.gt.1)
c     x      call MPI_BCAST(v_lnv,1,MPI_DOUBLE_PRECISION,
c     x                     0,comm_bead,ierr)

      endif ! pimd_npt


c     Thermostat integration.
c     N.B. If nstep > nequilib_cmd the thermostat path-centroid
c     ( ibead = 1 ) is de-activated since it must strictly follow
c     Newton equation of motions.


      if(langevin) then

        call Langevin_thermostat_integrate(idnode,ibead,iatm1,iatm2
     x       ,nvx,nvy,nvz,vxx,vyy,vzz,
     x       gle_vx,gle_vy,gle_vz,fix_com,fict_mass_nmode,alpha_cell)

      else
        call Thermostat_integrate(nchain,thermo,nthermo,dtx2)
        if (pimd_npt.and.mb_rank.eq.0) then
           call Barostat_thermostat_integrate(nchain,dtx2)
c        no need to broadcast v_lnv here as it is needed only on master
        endif
      endif


c     For PIMD.
      if ( .not.cmd_run ) then

         Ekin_fict = 0.d0

         do iatom = iatm1, iatm2
            Ekin_fict = Ekin_fict
     x        + fict_mass_nmode(iatom,ibead)
     x        * (ptr_vx(iatom)**2 + ptr_vy(iatom)**2 + ptr_vz(iatom)**2)

#ifdef HEAT_CURRENT
        call update_kinetic_energy(iatom,0.5d0*fict_mass_nmode
     x       (iatom,ibead)*(ptr_vx(iatom)**2+ptr_vy(iatom)**2
     x          +ptr_vz(iatom)**2))
#endif /*HEAT_CURRENT*/

         enddo

         if(langevin) then
           Etot_nose = 0.d0
         ! NPT thermostat lives at mb_rank==0 bead only
         elseif (mb_rank.eq.0) then
            Etot_nose =
     x         Thermostat_hamiltonian(nchain,thermo,nthermo,pimd_npt)
            if (pimd_npt) Etot_nose = Etot_nose + press*volm
         else
            Etot_nose =
     x         Thermostat_hamiltonian(nchain,thermo,nthermo,.false.)
         endif

         buffer(1) = Ekin_fict
         buffer(2) = Etot_nose

         call MPI_ALLREDUCE(buffer(1),buffer(3),2,
     x        MPI_DOUBLE_PRECISION,MPI_SUM,comm_mb,i)

         Ekin_fict = buffer(3)
         Etot_nose = buffer(4)

c        Temperature for the fictitious classical system.
         if ( langevin .and. fix_com ) then
            temp_fict = Ekin_fict / dble( 3*natms_r*nbead-3 ) / boltz
     x                            / force_scale
         else
            temp_fict = Ekin_fict / dble( 3*natms_r*nbead ) / boltz
     x                            / force_scale
         endif

c        Energetics for the fictitious classical system.
         Ekin_fict = 0.5d0 * Ekin_fict
         Etot_fict = Ekin_fict + Epot_fict

         if(langevin) then
            ! add pseudo-conserved quantity
         else
c          Conserved energy for the extended Lagrangian.
           Etot_nose = Etot_nose + Etot_fict
         endif

c        For nbead = 1 (classical simulation).
         if ( nbead==1 ) then
            Ekin_pimd = Ekin_fict
            Etot_pimd = Ekin_pimd + Epot_pimd
            Etot_prim = Etot_pimd
         endif

c        Output energies and positions for PIMD.
         if ( mod( nstep, istraj )==0 ) then

            if (lpolar .and. lcp) then
               call induced_error(natms_f, polr2, dipx, dipy, dipz,
     x              efieldkx, efieldky, efieldkz, emux, emuy, emuz,
     x              efddmurecx, efddmurecy, efddmurecz,
     x              efdcrecx, efdcrecy, efdcrecz, extx, exty, extz,
     x              err_L2, err_Linf)

               err_L2 = err_L2/sqrt(dble(ndipole))

               if (lhead) then
                  open(file_ind,file='STATIS_IND'//bead_suffix,
     x                  position='append')
                  if (lcp) then
                     write(file_ind,'(i12,4e16.8)')
     x                  nstep, time, err_L2, err_Linf,
     x                  engdke/boltz/dble(ndipole)/1.5d0
                  else
                     write(file_ind,'(i12,3e16.8)')
     x                  nstep, time, err_L2, err_Linf
                  end if ! lcp
                  close(file_ind)
               end if ! lhead
            end if ! lpolar

            if (pimd_head) then

               call timchk(0,timepb1)

               open(file_pimd,file='STATIS_PIMD',position='append')
               open(file_nh,file='STATIS_NH',position='append')
               open(file_stress,file='STATIS_STRESS',position='append')

               write(file_pimd,'(f12.5,5f15.5)')
     x               time, Etot_pimd / conv_energy,
     x                     Epot_pimd / conv_energy,
     x                     Ekin_pimd / conv_energy,
     x                     Etot_prim / conv_energy,
     x                     timepb1-timepb0

               write(file_stress,'(f12.5,10f10.5)')
     x               time, P_pimd*prsunt, stress_pimd(1:9)

               write(file_nh,'(i12,3f18.5)')
     x              nstep, time, temp_fict, Etot_nose / conv_energy

               close(file_pimd)
               close(file_nh)
               close(file_stress)

               timepb0 = timepb1

               if (pimd_npt) then
                  open(file_npt,file='STATIS_NPT',position='append')
                  if(langevin) then
                    write(file_npt,'(i12,8f18.5)')
     x                 nstep, P_pimd * prsunt,
     x                 volm,
     x                 alpha_cell/cell_mass,
     x                 P_lang * prsunt,
     x                 cell(1), cell(5), cell(9)
                  else
                    write(file_npt,'(i12,12f14.5)')
     x                 nstep, time, P_pimd * prsunt,
     x                 volm, cell(1:9)
                  endif
                  close(file_npt)
               end if

               if(CV_Ql) then
                  open(file_cv,file='STATIS_CV',position='append')
                  write(file_cv,'(i12,15f18.5)')
     x                 nstep,Ql_En(1),Ql_En(2),cv_pos(1),cv_pos(2),
     x                 ev_force(1),ev_force(2),
     x                 (CV_STRESS(i),i=1,9)
                  close(file_cv)
                  ! write CV_FORCE here if needed
               endif

            end if ! pimd_head

         end if ! mod( nstep, istraj ) == 0

         if (ltraj) then
         if (mod(nstep-nstraj,istraj).eq.0) then
            if(evolve_in_nm) then
c             For restart.
              call transform_velocity_from_nmode_to_cart
     x              (nvx,nvy,nvz,vxx,vyy,vzz)
           else
c             Transform cartesian positions into normal mode positions.
              call transform_position_from_cart_to_nmode
     x          (xxx,yyy,zzz,npx,npy,npz)

c             Transform cartesian positions into normal mode positions.
              call transform_velocity_from_cart_to_nmode
     x          (vxx,vyy,vzz,nvx,nvy,nvz)

c             Transform cartesian positions into normal mode positions.
              call transform_force_from_cart_to_nmode
     x          (fxx,fyy,fzz,nfx,nfy,nfz)
           endif

           call merge(idnode,mxnode,natms_r,mxbuff,npx,npy,npz,buffer)
           call merge(idnode,mxnode,natms_r,mxbuff,nvx,nvy,nvz,buffer)
           call merge(idnode,mxnode,natms_r,mxbuff,nfx,nfy,nfz,buffer)

           call merge(idnode,mxnode,natms_r,mxbuff,xxx,yyy,zzz,buffer)
           call merge(idnode,mxnode,natms_r,mxbuff,vxx,vyy,vzz,buffer)
           call merge(idnode,mxnode,natms_r,mxbuff,fxx,fyy,fzz,buffer)

           if ( lpolar .and. lcp ) then
        call merge(idnode,mxnode,natms_f,mxbuff,dipx,dipy,dipz,buffer)
        call merge(idnode,mxnode,natms_f,mxbuff,vdxx,vdyy,vdzz,buffer)
        call merge(idnode,mxnode,natms_f,mxbuff,fdxx,fdyy,fdzz,buffer)
           end if

#ifdef DEBUG
      write(4444,*) 'writing forces in kcal/mol '

      do iatom=1,natms_r
         write(4444,'(i3,1x,3(1x,ES16.9))'),iatom,
     x fxx(iatom)/engunit,fyy(iatom)/engunit,fzz(iatom)/engunit
      end do

#endif

#ifdef HEAT_CURRENT

      call potential_energy_per_atom()
#ifdef HEAT_STRESS
      call stress_per_atom()
      call compute_heat_flux(vxx,vyy,vzz)
#endif /* HEAT_STRESS */
      /* if (firststep) then
        call initialize_potential_energy(Epot_pimd,vxx,vyy,vzz)
        firststep = .false.
      else
        call potential_energy_per_atom2(tstep,vxx,vyy,vzz)
      end if
       */

      call compute_heat_flux2(xxx,yyy,zzz,vxx,vyy,vzz,gammattm)

      call final_sum_heat_flux()
      call final_sum()
#endif /* HEAT_CURRENT */

            if (lhead) then
               call centroid_traject_pimd
     x           (keytrj,nstep,keybin,time,mxatms,natms_r,atmnam,
     x            imcon,cell,xxx,yyy,zzz,vxx,vyy,vzz,
     x            fxx,fyy,fzz,lcavity)

#ifdef HEAT_CURRENT
            /* if (idnode==0) call write_heat_flux(time,
     x           bead_suffix) */
            /* if (idnode==0) call write_force_matrix
     x                     (time,bead_suffix) */
            if ((.not. firststep) .and. idnode==0) then
              call write_heat_flux2(time,bead_suffix)
            end if
            /* if (idnode==0) then
              call save_last_energy(nstep,time,
     x          bead_suffix)
            end if    */

#ifdef HEAT_CHECK
            if (idnode==0) call write_check(time,
     x            conv_energy,prsunt/stpvol,bead_suffix)
#endif /* HEAT_CHECK */

#endif /* HEAT_CURRENT */

               call centroid_write_config
     x                 (cfgname,nstep,time,levcfg,imcon,cell,
     x                  mxatms,natms_r,atmnam,nfict,
     x                  xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)


               if ( lpolar .and. lcp ) then
                  call centroid_write_dipole
     x                 (natms_f,atmnam,dipx,dipy,dipz,
     x                  vdxx,vdyy,vdzz,fdxx,fdyy,fdzz)
               endif

! VB: instead of the following pimdpol line
!               call centroid_traj_dip_pimd(nstep,time)
! VB: proceed as non-PIMD/CMD version does
! (there is no need to merge moldip? -- computed on every PE)
              if (lttm) then
                  ndummy = nummols(ntpmls)
c                  call output_ttm2_com(imcon,ndummy,time,cell,
c     x                                 weight,listttm2,xxx,yyy,zzz)
                  if (lcp) then
                      call output_ttm2_dip(listttm2,keybin,ndummy,nttm2,
     x                   dipx_old,dipy_old,dipz_old,.false.)
                      do i = 1, natms_f
                         dipx_old(i) = dipx(i)
                         dipy_old(i) = dipy(i)
                         dipz_old(i) = dipz(i)
                      end do
                  else
                      call output_ttm2_dip(listttm2,keybin,ndummy,nttm2,
     x                                     dipx,dipy,dipz,.false.)
                  endif ! lcp
               endif ! ttm

            end if ! lhead

         end if ! mod(nstep, istraj).eq.0
         end if ! ltraj

         ! PP_NOTE

      else ! cmd_run

c        Output energies and positions for CMD.
         Ekin_fict = 0.d0
         Ekin_cmd  = 0.d0

         do iatom = iatm1, iatm2
            if (ibead.gt.1) then
               Ekin_fict = Ekin_fict
     x              + fict_mass_nmode(iatom,ibead)
     x              * (nvx(iatom)**2 + nvy(iatom)**2 + nvz(iatom)**2)
#ifdef HEAT_CURRENT
               call update_kinetic_energy(iatom,0.5d0*
     x              fict_mass_nmode(iatom,ibead)
     x              * (nvx(iatom)**2 + nvy(iatom)**2 + nvz(iatom)**2))
#endif /* HEAT_CURRENT */
            else
               Ekin_cmd = Ekin_cmd
     x              + fict_mass_nmode(iatom,ibead)
     x              * (nvx(iatom)**2 + nvy(iatom)**2 + nvz(iatom)**2)
#ifdef HEAT_CURRENT
              call update_kinetic_energy(iatom,0.5d0*
     x          fict_mass_nmode(iatom,ibead)
     x              * (nvx(iatom)**2 + nvy(iatom)**2 + nvz(iatom)**2))
#endif /* HEAT_CURRENT */
            endif
         enddo

         buffer(1) = Ekin_fict
         buffer(2) = Ekin_cmd

         call MPI_ALLREDUCE(buffer(1), buffer(3), 2,
     x        MPI_DOUBLE_PRECISION, MPI_SUM, comm_mb, i)

         Ekin_fict = buffer(3)
         Ekin_cmd = buffer(4)

c        Energies and temperature for path-centroid.
         if ( cmd_nve ) then
            temp_cmd = Ekin_cmd / dble( 3*natms_r-3 ) / boltz
         else
            temp_cmd = Ekin_cmd / dble( 3*natms_r ) / boltz
         end if
         Ekin_cmd = 0.5d0 * Ekin_cmd
         Etot_cmd = Ekin_cmd + Epot_pimd

c        Kinetic energy and temperature for the fictitious classical system.
c        The path-centroid is running in NVE.
c        Conserved energy for the extended Lagrangian.
         if ( nbead>1 ) then
            temp_fict = Ekin_fict / dble( 3*natms_r*(nbead-1) ) / boltz
            Ekin_fict = 0.5d0 * Ekin_fict
            Etot_nose =
     x          Thermostat_hamiltonian(nchain,thermo,nthermo,.false.)
            call MPI_ALLREDUCE(Etot_nose,buffer(1),1,
     x           MPI_DOUBLE_PRECISION, MPI_SUM, comm_mb, i)
            Etot_nose = Ekin_fict + Etot_cmd + buffer(1)
         else
            temp_fict = temp_cmd
            Etot_nose = Etot_cmd
         endif

c        Output energies and positions & velocities for path-centroid.
         if ( mod( nstep, istraj )==0 ) then

            if ( pimd_head ) then

               open(file_cmd,file='STATIS_CMD', position='append')
               open(file_nh,file='STATIS_NH', position='append')
               open(file_stress,file='STATIS_STRESS',position='append')

               write(file_stress,'(f12.5,10f10.5)')
     x               time, P_pimd*prsunt, stress_pimd(1:9)

               write(file_cmd,'(i12,5f18.5)')
     x              nstep, time,
     x              Etot_cmd  / conv_energy,
     x              Epot_pimd / conv_energy,
     x              Ekin_cmd  / conv_energy,
     x              temp_cmd

               if ( nbead == 1 ) then
                  write( file_nh, '(i12,2f18.5)' )
     x                 nstep, temp_cmd, Etot_cmd  / conv_energy
               else
                  write( file_nh, '(i12,2f18.5)' )
     x                 nstep, temp_fict, Etot_nose / conv_energy
               end if

               close( file_cmd )
               close( file_nh )

            endif ! pimd_lhead

           if(evolve_in_nm) then
c             For restart.
              call transform_velocity_from_nmode_to_cart
     x              (nvx,nvy,nvz,vxx,vyy,vzz)
           else
c             Transform cartesian positions into normal mode positions.
              call transform_position_from_cart_to_nmode
     x          (xxx,yyy,zzz,npx,npy,npz)

c             Transform cartesian positions into normal mode positions.
              call transform_velocity_from_cart_to_nmode
     x          (vxx,vyy,vzz,nvx,nvy,nvz)

c             Transform cartesian positions into normal mode positions.
              call transform_force_from_cart_to_nmode
     x          (fxx,fyy,fzz,nfx,nfy,nfz)
           endif

           call merge(idnode,mxnode,natms_r,mxbuff,npx,npy,npz,buffer)
           call merge(idnode,mxnode,natms_r,mxbuff,nvx,nvy,nvz,buffer)
           call merge(idnode,mxnode,natms_r,mxbuff,nfx,nfy,nfz,buffer)

           call merge(idnode,mxnode,natms_r,mxbuff,xxx,yyy,zzz,buffer)
           call merge(idnode,mxnode,natms_r,mxbuff,vxx,vyy,vzz,buffer)
           call merge(idnode,mxnode,natms_r,mxbuff,fxx,fyy,fzz,buffer)

           if ( lpolar .and. lcp ) then
        call merge(idnode,mxnode,natms_f,mxbuff,dipx,dipy,dipz,buffer)
        call merge(idnode,mxnode,natms_f,mxbuff,vdxx,vdyy,vdzz,buffer)
        call merge(idnode,mxnode,natms_f,mxbuff,fdxx,fdyy,fdzz,buffer)
           end if

           if (pimd_head) then
              call centroid_traject_cmd
     x           (keytrj,keybin,nstep,time,mxatms,natms_r,atmnam,cell,
     x            npx,npy,npz,nvx,nvy,nvz,nfx,nfy,nfz,xxx,yyy,zzz,
     x            vxx,vyy,vzz,fxx,fyy,fzz)
           end if

#ifdef HEAT_CURRENT

          call potential_energy_per_atom()
#ifdef HEAT_STRESS
          call stress_per_atom()
          call compute_heat_flux(vxx,vyy,vzz)
#endif /* HEAT_STRESS */
          /* if (firststep) then
            call initialize_potential_energy(Epot_pimd,vxx,vyy,vzz)
          else
            call potential_energy_per_atom2(tstep,vxx,vyy,vzz)
          end if
           */
          call compute_heat_flux2(xxx,yyy,zzz,vxx,vyy,vzz,gammattm)

          /* if (firststep) then
            call final_sum_epot()
            firststep = .false.
          end if
           */

          call final_sum_heat_flux()
          call final_sum()
#endif /* HEAT_CURRENT */
#ifdef DEBUG
      write(4444,*) 'writing forces in kcal/mol '

      do iatom=1,natms_r
         write(4444,'(i3,1x,3(1x,ES16.9))'),iatom,
     x fxx(iatom)/engunit,fyy(iatom)/engunit,fzz(iatom)/engunit
      end do

#endif

           if (lhead) then

#ifdef HEAT_CURRENT
            /* if (idnode==0) call write_heat_flux(time,
     x           bead_suffix)
            if (idnode==0) call write_force_matrix
     x                     (time,bead_suffix) */
            if (idnode==0) call compute_stress(time,
     x        xxx, yyy, zzz, prsunt/stpvol)
            if (idnode==0) then
              call write_heat_flux2(time,bead_suffix)
            end if
            /* if (idnode==0) then
              call save_last_energy(nstep,time,
     x          bead_suffix)
            end if */
#ifdef HEAT_CHECK
            if (idnode==0) call write_check(time,
     x            conv_energy,prsunt/stpvol,bead_suffix)
#endif /* HEAT_CHECK */
#endif /* HEAT_CURRENT */

              if (w_bead) then
                 call centroid_traject_bead
     x              (nstep,keybin,time,mxatms,natms_r,atmnam,imcon,cell,
     x               xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)
              end if

               call centroid_write_config
     x                 (cfgname,nstep,time,levcfg,imcon,cell,
     x                  mxatms,natms_r,atmnam,nfict,
     x                  xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)

               if ( lpolar .and. lcp ) then
                   call centroid_write_dipole
     x                 (natms_f,atmnam,dipx,dipy,dipz,
     x                  vdxx,vdyy,vdzz,fdxx,fdyy,fdzz)
               endif
               if (lttm) then
                   ndummy = nummols(ntpmls)
c                  call output_ttm2_com(imcon,ndummy,time,cell,
c    x                                  weight,listttm2,xxx,yyy,zzz)
                   if (lcp) then
                       call output_ttm2_dip(listttm2,keybin,ndummy,nttm2,
     x                   dipx_old,dipy_old,dipz_old,.true.)
                       do i = 1, natms_f
                       dipx_old(i) = dipx(i)
                       dipy_old(i) = dipy(i)
                       dipz_old(i) = dipz(i)
                       end do
                   else
                       call output_ttm2_dip(listttm2,keybin,ndummy,nttm2,
     x                   dipx,dipy,dipz,.true.)
                   endif ! lcp
               endif ! ttm
           end if ! lhead

         endif ! mod(nstep, istraj) == 0

      endif ! cmd_run

c
c     complete time check

      call timchk(0,timelp)

#endif /* IPI */

      newjob=.false.

#ifdef DO_2D_IR
      call merge(idnode,mxnode,natms_f,mxbuff,dipx,dipy,dipz,buffer)
      if (lhead)
     x call centroid_2D_IR_report(ihod,ivib,
     x          stpcfg/engunit,xxx,yyy,zzz,
     x          listttm2,nttm2,weight,
     x          dipx,dipy,dipz)

      enddo
      enddo
      call centroid_2D_IR_read_frame(natom, nttm2, listttm2,
     x                               continue_2D_IR)
      call MPI_BCAST(continue_2D_IR, 1,
     x               MPI_LOGICAL, 0, comm_bead, ierr)
      enddo
#elif defined(DO_RECALC_DIP)
      call merge(idnode,mxnode,natms_f,mxbuff,dipx,dipy,dipz,buffer)

      call centroid_recalc_dip_read_frame(natom, nttm2, listttm2,
     x                                    continue_RECALC_DIP)
      write(22,*) "Finished reading frame, cont: ", continue_RECALC_DIP
      enddo
#elif IPI
      go to 100
#else
      if(nstep.lt.nstrun.and.timjob-timelp.gt.timcls) go to 100
#endif

101   continue

#ifdef OPTORVIB
c Opt
      if(lopt.and..not.lvib) call opt_finish()
      if(lopt.and.lvib) call vib_finish()

#endif /* OPTORVIB */

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

#ifdef HEAT_CURRENT
      call deallocate_heat()
#endif /* HEAT_CURRENT */

#ifdef DO_2D_IR
      call centroid_2D_IR_fini()
#endif /* DO_2D_IR */
#ifdef DO_RECALC_DIP
      call centroid_recalc_dip_fini()
#endif /* DO_RECALC_DIP */

      if(idnode.eq.0)write(nrite,
     x  "(/,/,1x,'run terminating. elapsed cpu time = ',f13.3,
     x  ', job time = ',f13.3,', close time = ',f13.3,/)")
     x  timelp,timjob,timcls
c
c     close output channels

      if(idnode.eq.0) then

        close (nrite)
        close (nstats)
        close (nhist)

      endif
c
c     terminate job

#ifdef PLUMED
      if(lplumed) then
         call plumed_f_gfinalize()
      endif
#endif /* PLUMED */


#ifdef VAMPIR
      call VTEND(99, ierr)
#endif

#ifdef DO_2D_IR
      call centroid_2D_IR_fini()
#endif
#ifdef DO_RECALC_DIP
      call centroid_recalc_dip_fini()
#endif /* DO_RECALC_DIP */
      if(laspc) call aspc_fini()

      ! NOTE NOTE NOTE
      if(lmbpol) call mbpol_fini()
      if(lmbnrg) call mbnrg_fini()

      call mb_finalize()

      end
