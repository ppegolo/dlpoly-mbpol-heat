! *****************************************************************************************
!   Module to accumulate the stress per atom, the energy per atom
!   and to calculate the macroscopic heat current.
!
!   23-03-2018: PP_: first version of update_stress and update_energy.
!   10-06-2018: PP_: qtip4p: Completed force distribution (no ewald1)
! *****************************************************************************************
module heatcurrent
  !
  use global_variables, only: mxatms    ! number of atoms
  use multibead, only: comm_bead, bead_rank

  implicit none
  !
  real(8), dimension(1:3)             :: j1,jkin   ! atomic energy contribution to the heat current
  real(8), dimension(1:3)             :: j2,jpot,jpot_kspace   ! stress contribution to the heat current
  real(8), dimension(1:3)             :: jtot, jtot2   ! stress contribution to the heat current
  real(8), dimension(:,:), allocatable :: f_dot_v_old
  real(8), parameter :: junit=1.0d3  ! conversion form internal units

  real(8), dimension(:,:,:), allocatable :: atomic_stress_ang,atomic_stress_bnd,&
  atomic_stress_ew1p, atomic_stress_ew2p, atomic_stress_ew3p, &
  atomic_stress_qd, atomic_stress_mbpol, atomic_stress_srf, atomic_stress_lrcorrect,atomic_stress  ! Array with stress on each atom

  real(8), dimension(:), allocatable     ::  atomic_energy_ang,atomic_energy_bnd,&
  atomic_energy_ew1p, atomic_energy_ew2p, atomic_energy_ew3p, &
  atomic_energy_mbpol, atomic_energy_srf, atomic_energy_lrcorrect, &
  atomic_energy_polar, atomic_kinetic_energy, atomic_potential_energy, atomic_potential_energy2,&  ! Array with the energy assigned to each atom
  delta_potential_energy(:)
  real(8), dimension(:,:,:), allocatable     :: force_matrix
  real(8), dimension(:,:), allocatable :: force_ew1p

  integer :: ierror
  ! For the Finsl_sum
  real(8), allocatable :: tmpf(:), tmpf3(:), tmpFM3(:), tmpFEW3(:)
  real(8), pointer :: tmpf2(:), tmpFM2(:), tmpFEW2(:)
  real(8), contiguous, pointer :: tmpf4(:,:,:)
  real(8), contiguous, pointer :: tmpf5(:,:,:)
  real(8), contiguous, pointer :: tmpFM4(:,:,:),tmpFM1(:,:,:), tmpFEW1(:,:), tmpFEW4(:,:)
  !
contains
  !
  ! ====================================================================================================================
  ! ============================================ Initialization ========================================================
  ! ====================================================================================================================
  !
  subroutine init_heat()
    allocate(atomic_stress(1:mxatms,1:3,1:3), &
    atomic_stress_ang(1:mxatms,1:3,1:3), &
    atomic_stress_bnd(1:mxatms,1:3,1:3), &
    atomic_stress_ew1p(1:mxatms,1:3,1:3), &
    atomic_stress_ew2p(1:mxatms,1:3,1:3), &
    atomic_stress_ew3p(1:mxatms,1:3,1:3), &
    atomic_stress_qd(1:mxatms,1:3,1:3), &
    atomic_stress_mbpol(1:mxatms,1:3,1:3), &
    atomic_stress_srf(1:mxatms,1:3,1:3), &
    atomic_stress_lrcorrect(1:mxatms,1:3,1:3), &
    atomic_energy_ang(1:mxatms), &
    atomic_energy_bnd(1:mxatms), &
    atomic_energy_ew1p(1:mxatms), &
    atomic_energy_ew2p(1:mxatms), &
    atomic_energy_ew3p(1:mxatms), &
    atomic_energy_mbpol(1:mxatms), &
    atomic_energy_srf(1:mxatms), &
    atomic_energy_lrcorrect(1:mxatms), &
    atomic_energy_polar(1:mxatms), &
    atomic_kinetic_energy(1:mxatms), &
    atomic_potential_energy(1:mxatms),&
    atomic_potential_energy2(1:mxatms),&
    delta_potential_energy(1:mxatms),&
    force_matrix(1:mxatms,1:mxatms,1:3),&
    force_ew1p(1:mxatms,1:3),&
    f_dot_v_old(1:mxatms,1:mxatms))
    atomic_stress=0.d0
    atomic_stress_ang=0.d0
    atomic_stress_bnd=0.d0
    atomic_stress_ew1p=0.d0
    atomic_stress_ew2p=0.d0
    atomic_stress_ew3p=0.d0
    atomic_stress_qd=0.d0
    atomic_stress_mbpol=0.d0
    atomic_stress_srf=0.d0
    atomic_stress_lrcorrect=0.d0
    atomic_energy_bnd=0.d0
    atomic_energy_ang=0.d0
    atomic_energy_ew1p=0.d0
    atomic_energy_ew2p=0.d0
    atomic_energy_ew3p=0.d0
    atomic_energy_mbpol=0.d0
    atomic_energy_srf=0.d0
    atomic_energy_lrcorrect=0.d0
    atomic_energy_polar=0.d0
    atomic_potential_energy2 = 0.d0
    delta_potential_energy = 0.d0
    j1=0.d0
    j2=0.d0
    jkin=0.d0
    jpot=0.d0
    force_matrix=0.d0
    force_ew1p = 0.d0
    !
    allocate(tmpf2(9*mxatms),tmpf4(mxatms,3,3))
    allocate(tmpf(2*mxatms), tmpf3(9*mxatms))
    allocate(tmpFM1(mxatms,mxatms,3))
    allocate(tmpFM2(mxatms*mxatms*3))
    allocate(tmpFM3(mxatms*mxatms*3))
    allocate(tmpFEW1(mxatms,3))
    allocate(tmpFEW2(mxatms*3))
    allocate(tmpFEW3(mxatms*3))
  end subroutine init_heat
  !
  subroutine zero_heat()
    atomic_stress=0.d0
    atomic_stress_ang=0.d0
    atomic_stress_bnd=0.d0
    atomic_stress_ew1p=0.d0
    atomic_stress_ew2p=0.d0
    atomic_stress_ew3p=0.d0
    atomic_stress_qd=0.d0
    atomic_stress_mbpol=0.d0
    atomic_stress_srf=0.d0
    atomic_stress_lrcorrect=0.d0
    atomic_energy_ang=0.d0
    atomic_energy_bnd=0.d0
    atomic_energy_ew1p=0.d0
    atomic_energy_ew2p=0.d0
    atomic_energy_ew3p=0.d0
    atomic_energy_mbpol=0.d0
    atomic_energy_srf=0.d0
    atomic_energy_lrcorrect=0.d0
    atomic_energy_polar=0.d0
    atomic_kinetic_energy=0.d0
    atomic_potential_energy=0.d0
    !atomic_potential_energy2 = 0.d0
    delta_potential_energy = 0.d0
    j1=0.d0
    j2=0.d0
    jkin=0.d0
    jpot=0.d0
    force_matrix=0.d0
    force_ew1p = 0.d0
  end subroutine zero_heat
  !
  ! ====================================================================================================================
  ! ====================================================================================================================
  !
  ! ====================================================================================================================
  ! ========================= LAMPPS-like method: "Naive" energy and stress distribution ===============================
  ! ====================================================================================================================
  !
  ! ____________________________________________________________________________________________________________________
  ! _______________-_____________________ Update energy & stress _______________________________________________________
  ! ____________________________________________________________________________________________________________________
  !
  subroutine update_stress_ang(iatom,icart1,icart2,stress_per_atom)
    integer, intent(in) :: iatom, icart1, icart2  ! atomic index, cartesian indeces
    real(8), intent(in) :: stress_per_atom!, velocity_i_cart
    atomic_stress_ang(iatom,icart1,icart2) = atomic_stress_ang(iatom,icart1,icart2) + stress_per_atom
  end subroutine update_stress_ang
  !
  subroutine update_energy_ang(iatom,energy_per_atom)
    integer, intent(in) :: iatom! atomic index, cartesian index
    real(8), intent(in) :: energy_per_atom
    atomic_energy_ang(iatom) = atomic_energy_ang(iatom) + energy_per_atom
  end subroutine update_energy_ang
  !
  subroutine update_stress_bnd(iatom,icart1,icart2,stress_per_atom)
    integer, intent(in) :: iatom, icart1, icart2  ! atomic index, cartesian indeces
    real(8), intent(in) :: stress_per_atom!, velocity_i_cart
    atomic_stress_bnd(iatom,icart1,icart2) = atomic_stress_bnd(iatom,icart1,icart2) + stress_per_atom
  end subroutine update_stress_bnd
  !
  subroutine update_energy_bnd(iatom,energy_per_atom)
    integer, intent(in) :: iatom! atomic index, cartesian index
    real(8), intent(in) :: energy_per_atom
    atomic_energy_bnd(iatom) = atomic_energy_bnd(iatom) + energy_per_atom
  end subroutine update_energy_bnd
  !
  subroutine update_stress_ew1p(iatom,icart1,icart2,stress_per_atom)
    integer, intent(in) :: iatom, icart1, icart2  ! atomic index, cartesian indeces
    real(8), intent(in) :: stress_per_atom!, velocity_i_cart
    atomic_stress_ew1p(iatom,icart1,icart2) = atomic_stress_ew1p(iatom,icart1,icart2) + stress_per_atom
  end subroutine update_stress_ew1p
  !
  subroutine update_energy_ew1p(iatom,energy_per_atom)
    integer, intent(in) :: iatom! atomic index, cartesian index
    real(8), intent(in) :: energy_per_atom
    atomic_energy_ew1p(iatom) = atomic_energy_ew1p(iatom) + energy_per_atom
  end subroutine update_energy_ew1p
  !
  subroutine update_stress_ew2p(iatom,icart1,icart2,stress_per_atom)
    integer, intent(in) :: iatom, icart1, icart2  ! atomic index, cartesian indeces
    real(8), intent(in) :: stress_per_atom!, velocity_i_cart
    atomic_stress_ew2p(iatom,icart1,icart2) = atomic_stress_ew2p(iatom,icart1,icart2) + stress_per_atom
  end subroutine update_stress_ew2p
  !
  subroutine update_energy_ew2p(iatom,energy_per_atom)
    integer, intent(in) :: iatom! atomic index, cartesian index
    real(8), intent(in) :: energy_per_atom
    atomic_energy_ew2p(iatom) = atomic_energy_ew2p(iatom) + energy_per_atom
  end subroutine update_energy_ew2p
  !
  subroutine update_stress_ew3p(iatom,icart1,icart2,stress_per_atom)
    integer, intent(in) :: iatom, icart1, icart2  ! atomic index, cartesian indeces
    real(8), intent(in) :: stress_per_atom!, velocity_i_cart
    atomic_stress_ew3p(iatom,icart1,icart2) = atomic_stress_ew3p(iatom,icart1,icart2) + stress_per_atom
  end subroutine update_stress_ew3p
  !
  subroutine update_energy_ew3p(iatom,energy_per_atom)
    integer, intent(in) :: iatom! atomic index, cartesian index
    real(8), intent(in) :: energy_per_atom
    atomic_energy_ew3p(iatom) = atomic_energy_ew3p(iatom) + energy_per_atom
  end subroutine update_energy_ew3p
  !
  subroutine update_stress_qd(iatom,icart1,icart2,stress_per_atom)
    integer, intent(in) :: iatom, icart1, icart2  ! atomic index, cartesian indeces
    real(8), intent(in) :: stress_per_atom!, velocity_i_cart
    atomic_stress_qd(iatom,icart1,icart2) = atomic_stress_qd(iatom,icart1,icart2) + stress_per_atom
  end subroutine update_stress_qd
  !
  subroutine update_stress_mbpol(iatom,icart1,icart2,stress_per_atom)
    integer, intent(in) :: iatom, icart1, icart2  ! atomic index, cartesian indeces
    real(8), intent(in) :: stress_per_atom!, velocity_i_cart
    atomic_stress_mbpol(iatom,icart1,icart2) = atomic_stress_mbpol(iatom,icart1,icart2) + stress_per_atom
  end subroutine update_stress_mbpol
  !
  subroutine update_energy_mbpol(iatom,energy_per_atom)
    integer, intent(in) :: iatom! atomic index, cartesian index
    real(8), intent(in) :: energy_per_atom
    atomic_energy_mbpol(iatom) = atomic_energy_mbpol(iatom) + energy_per_atom
  end subroutine update_energy_mbpol
  !
  subroutine update_stress_srf(iatom,icart1,icart2,stress_per_atom)
    integer, intent(in) :: iatom, icart1, icart2  ! atomic index, cartesian indeces
    real(8), intent(in) :: stress_per_atom!, velocity_i_cart
    atomic_stress_srf(iatom,icart1,icart2) = atomic_stress_srf(iatom,icart1,icart2) + stress_per_atom
  end subroutine update_stress_srf
  !
  subroutine update_energy_srf(iatom,energy_per_atom)
    integer, intent(in) :: iatom! atomic index, cartesian index
    real(8), intent(in) :: energy_per_atom
    atomic_energy_srf(iatom) = atomic_energy_srf(iatom) + energy_per_atom
  end subroutine update_energy_srf
  !
  subroutine update_energy_lrcorrect(iatom,energy_per_atom)
    integer, intent(in) :: iatom! atomic index, cartesian index
    real(8), intent(in) :: energy_per_atom
    atomic_energy_lrcorrect(iatom) = atomic_energy_lrcorrect(iatom) + energy_per_atom
  end subroutine update_energy_lrcorrect
  !
  subroutine update_stress_lrcorrect(iatom,icart1,icart2,stress_per_atom)
    integer, intent(in) :: iatom, icart1, icart2  ! atomic index, cartesian indeces
    real(8), intent(in) :: stress_per_atom!, velocity_i_cart
    atomic_stress_lrcorrect(iatom,icart1,icart2) = atomic_stress_lrcorrect(iatom,icart1,icart2) + stress_per_atom
  end subroutine update_stress_lrcorrect
  !
  subroutine update_energy_polar(iatom,energy_per_atom)
    integer, intent(in) :: iatom! atomic index, cartesian index
    real(8), intent(in) :: energy_per_atom
    atomic_energy_polar(iatom) = atomic_energy_polar(iatom) + energy_per_atom
  end subroutine update_energy_polar
  !
  subroutine update_kinetic_energy(iatm,kinetic_energy)
    real(8) :: kinetic_energy
    integer :: iatm
    atomic_kinetic_energy(iatm) = atomic_kinetic_energy(iatm) + kinetic_energy
  end subroutine update_kinetic_energy
  !
  ! ____________________________________________________________________________________________________________________
  ! _________________________ Atomic kinetic and potential energy, atomic stress _______________________________________
  ! ____________________________________________________________________________________________________________________
  subroutine kinetic_energy_per_atom(vxx,vyy,vzz,weight)
    real(8), dimension(1:mxatms) :: vxx,vyy,vzz,weight
    real(8) :: vtmp(3)
    integer :: iatm
    do iatm=1,mxatms
      vtmp = (/ vxx(iatm), vyy(iatm), vzz(iatm) /)
      atomic_kinetic_energy(iatm) = 0.5d0*weight(iatm)*dot_product(vtmp,vtmp)
    end do
  end subroutine kinetic_energy_per_atom
  !
  subroutine potential_energy_per_atom()
    integer :: iatm
    do iatm=1,mxatms
      atomic_potential_energy(iatm) = atomic_energy_ang(iatm)+atomic_energy_bnd(iatm) + &
      atomic_energy_ew1p(iatm)+atomic_energy_ew2p(iatm)+atomic_energy_ew3p(iatm)+ &
      atomic_energy_mbpol(iatm)+atomic_energy_polar(iatm)+atomic_energy_srf(iatm)+&
      atomic_energy_lrcorrect(iatm)
    end do
  end subroutine potential_energy_per_atom
  !
  subroutine stress_per_atom()
    atomic_stress = atomic_stress_ang+atomic_stress_bnd + &
    atomic_stress_ew1p+atomic_stress_ew2p+atomic_stress_ew3p+ &
    atomic_stress_qd+atomic_stress_mbpol+atomic_stress_srf+atomic_stress_lrcorrect
  end subroutine stress_per_atom
  !
  ! ____________________________________________________________________________________________________________________
  ! ______________________________ Distribution of energy and stress in the M sites ____________________________________
  ! ____________________________________________________________________________________________________________________
  !
  subroutine distribute_stress(ioxy,ih1,ih2,mttm2,gammattm)
    integer :: ioxy, ih1, ih2, mttm2
    real(8) :: gammattm
    !
    atomic_stress_ang(ioxy,:,1)=atomic_stress_ang(ioxy,:,1)+(1.d0-gammattm)*atomic_stress_ang(mttm2,:,1)
    atomic_stress_ang(ioxy,:,2)=atomic_stress_ang(ioxy,:,2)+(1.d0-gammattm)*atomic_stress_ang(mttm2,:,2)
    atomic_stress_ang(ioxy,:,3)=atomic_stress_ang(ioxy,:,3)+(1.d0-gammattm)*atomic_stress_ang(mttm2,:,3)
    atomic_stress_ang(ih1,:,1)=atomic_stress_ang(ih1,:,1)+gammattm/2.d0*atomic_stress_ang(mttm2,:,1)
    atomic_stress_ang(ih1,:,2)=atomic_stress_ang(ih1,:,2)+gammattm/2.d0*atomic_stress_ang(mttm2,:,2)
    atomic_stress_ang(ih1,:,3)=atomic_stress_ang(ih1,:,3)+gammattm/2.d0*atomic_stress_ang(mttm2,:,3)
    atomic_stress_ang(ih2,:,1)=atomic_stress_ang(ih2,:,1)+gammattm/2.d0*atomic_stress_ang(mttm2,:,1)
    atomic_stress_ang(ih2,:,2)=atomic_stress_ang(ih2,:,2)+gammattm/2.d0*atomic_stress_ang(mttm2,:,2)
    atomic_stress_ang(ih2,:,3)=atomic_stress_ang(ih2,:,3)+gammattm/2.d0*atomic_stress_ang(mttm2,:,3)
    atomic_stress_ang(mttm2,:,:)=0
    !
    atomic_stress_bnd(ioxy,:,1)=atomic_stress_bnd(ioxy,:,1)+(1.d0-gammattm)*atomic_stress_bnd(mttm2,:,1)
    atomic_stress_bnd(ioxy,:,2)=atomic_stress_bnd(ioxy,:,2)+(1.d0-gammattm)*atomic_stress_bnd(mttm2,:,2)
    atomic_stress_bnd(ioxy,:,3)=atomic_stress_bnd(ioxy,:,3)+(1.d0-gammattm)*atomic_stress_bnd(mttm2,:,3)
    atomic_stress_bnd(ih1,:,1)=atomic_stress_bnd(ih1,:,1)+gammattm/2.d0*atomic_stress_bnd(mttm2,:,1)
    atomic_stress_bnd(ih1,:,2)=atomic_stress_bnd(ih1,:,2)+gammattm/2.d0*atomic_stress_bnd(mttm2,:,2)
    atomic_stress_bnd(ih1,:,3)=atomic_stress_bnd(ih1,:,3)+gammattm/2.d0*atomic_stress_bnd(mttm2,:,3)
    atomic_stress_bnd(ih2,:,1)=atomic_stress_bnd(ih2,:,1)+gammattm/2.d0*atomic_stress_bnd(mttm2,:,1)
    atomic_stress_bnd(ih2,:,2)=atomic_stress_bnd(ih2,:,2)+gammattm/2.d0*atomic_stress_bnd(mttm2,:,2)
    atomic_stress_bnd(ih2,:,3)=atomic_stress_bnd(ih2,:,3)+gammattm/2.d0*atomic_stress_bnd(mttm2,:,3)
    atomic_stress_bnd(mttm2,:,:)=0
    !
    atomic_stress_ew1p(ioxy,:,1)=atomic_stress_ew1p(ioxy,:,1)+(1.d0-gammattm)*atomic_stress_ew1p(mttm2,:,1)
    atomic_stress_ew1p(ioxy,:,2)=atomic_stress_ew1p(ioxy,:,2)+(1.d0-gammattm)*atomic_stress_ew1p(mttm2,:,2)
    atomic_stress_ew1p(ioxy,:,3)=atomic_stress_ew1p(ioxy,:,3)+(1.d0-gammattm)*atomic_stress_ew1p(mttm2,:,3)
    atomic_stress_ew1p(ih1,:,1)=atomic_stress_ew1p(ih1,:,1)+gammattm/2.d0*atomic_stress_ew1p(mttm2,:,1)
    atomic_stress_ew1p(ih1,:,2)=atomic_stress_ew1p(ih1,:,2)+gammattm/2.d0*atomic_stress_ew1p(mttm2,:,2)
    atomic_stress_ew1p(ih1,:,3)=atomic_stress_ew1p(ih1,:,3)+gammattm/2.d0*atomic_stress_ew1p(mttm2,:,3)
    atomic_stress_ew1p(ih2,:,1)=atomic_stress_ew1p(ih2,:,1)+gammattm/2.d0*atomic_stress_ew1p(mttm2,:,1)
    atomic_stress_ew1p(ih2,:,2)=atomic_stress_ew1p(ih2,:,2)+gammattm/2.d0*atomic_stress_ew1p(mttm2,:,2)
    atomic_stress_ew1p(ih2,:,3)=atomic_stress_ew1p(ih2,:,3)+gammattm/2.d0*atomic_stress_ew1p(mttm2,:,3)
    atomic_stress_ew1p(mttm2,:,:)=0
    !
    atomic_stress_ew2p(ioxy,:,1)=atomic_stress_ew2p(ioxy,:,1)+(1.d0-gammattm)*atomic_stress_ew2p(mttm2,:,1)
    atomic_stress_ew2p(ioxy,:,2)=atomic_stress_ew2p(ioxy,:,2)+(1.d0-gammattm)*atomic_stress_ew2p(mttm2,:,2)
    atomic_stress_ew2p(ioxy,:,3)=atomic_stress_ew2p(ioxy,:,3)+(1.d0-gammattm)*atomic_stress_ew2p(mttm2,:,3)
    atomic_stress_ew2p(ih1,:,1)=atomic_stress_ew2p(ih1,:,1)+gammattm/2.d0*atomic_stress_ew2p(mttm2,:,1)
    atomic_stress_ew2p(ih1,:,2)=atomic_stress_ew2p(ih1,:,2)+gammattm/2.d0*atomic_stress_ew2p(mttm2,:,2)
    atomic_stress_ew2p(ih1,:,3)=atomic_stress_ew2p(ih1,:,3)+gammattm/2.d0*atomic_stress_ew2p(mttm2,:,3)
    atomic_stress_ew2p(ih2,:,1)=atomic_stress_ew2p(ih2,:,1)+gammattm/2.d0*atomic_stress_ew2p(mttm2,:,1)
    atomic_stress_ew2p(ih2,:,2)=atomic_stress_ew2p(ih2,:,2)+gammattm/2.d0*atomic_stress_ew2p(mttm2,:,2)
    atomic_stress_ew2p(ih2,:,3)=atomic_stress_ew2p(ih2,:,3)+gammattm/2.d0*atomic_stress_ew2p(mttm2,:,3)
    atomic_stress_ew2p(mttm2,:,:)=0
    !
    atomic_stress_ew3p(ioxy,:,1)=atomic_stress_ew3p(ioxy,:,1)+(1.d0-gammattm)*atomic_stress_ew3p(mttm2,:,1)
    atomic_stress_ew3p(ioxy,:,2)=atomic_stress_ew3p(ioxy,:,2)+(1.d0-gammattm)*atomic_stress_ew3p(mttm2,:,2)
    atomic_stress_ew3p(ioxy,:,3)=atomic_stress_ew3p(ioxy,:,3)+(1.d0-gammattm)*atomic_stress_ew3p(mttm2,:,3)
    atomic_stress_ew3p(ih1,:,1)=atomic_stress_ew3p(ih1,:,1)+gammattm/2.d0*atomic_stress_ew3p(mttm2,:,1)
    atomic_stress_ew3p(ih1,:,2)=atomic_stress_ew3p(ih1,:,2)+gammattm/2.d0*atomic_stress_ew3p(mttm2,:,2)
    atomic_stress_ew3p(ih1,:,3)=atomic_stress_ew3p(ih1,:,3)+gammattm/2.d0*atomic_stress_ew3p(mttm2,:,3)
    atomic_stress_ew3p(ih2,:,1)=atomic_stress_ew3p(ih2,:,1)+gammattm/2.d0*atomic_stress_ew3p(mttm2,:,1)
    atomic_stress_ew3p(ih2,:,2)=atomic_stress_ew3p(ih2,:,2)+gammattm/2.d0*atomic_stress_ew3p(mttm2,:,2)
    atomic_stress_ew3p(ih2,:,3)=atomic_stress_ew3p(ih2,:,3)+gammattm/2.d0*atomic_stress_ew3p(mttm2,:,3)
    atomic_stress_ew3p(mttm2,:,:)=0
    !
    atomic_stress_qd(ioxy,:,1)=atomic_stress_qd(ioxy,:,1)+(1.d0-gammattm)*atomic_stress_qd(mttm2,:,1)
    atomic_stress_qd(ioxy,:,2)=atomic_stress_qd(ioxy,:,2)+(1.d0-gammattm)*atomic_stress_qd(mttm2,:,2)
    atomic_stress_qd(ioxy,:,3)=atomic_stress_qd(ioxy,:,3)+(1.d0-gammattm)*atomic_stress_qd(mttm2,:,3)
    atomic_stress_qd(ih1,:,1)=atomic_stress_qd(ih1,:,1)+gammattm/2.d0*atomic_stress_qd(mttm2,:,1)
    atomic_stress_qd(ih1,:,2)=atomic_stress_qd(ih1,:,2)+gammattm/2.d0*atomic_stress_qd(mttm2,:,2)
    atomic_stress_qd(ih1,:,3)=atomic_stress_qd(ih1,:,3)+gammattm/2.d0*atomic_stress_qd(mttm2,:,3)
    atomic_stress_qd(ih2,:,1)=atomic_stress_qd(ih2,:,1)+gammattm/2.d0*atomic_stress_qd(mttm2,:,1)
    atomic_stress_qd(ih2,:,2)=atomic_stress_qd(ih2,:,2)+gammattm/2.d0*atomic_stress_qd(mttm2,:,2)
    atomic_stress_qd(ih2,:,3)=atomic_stress_qd(ih2,:,3)+gammattm/2.d0*atomic_stress_qd(mttm2,:,3)
    atomic_stress_qd(mttm2,:,:)=0
    !
    atomic_stress_mbpol(ioxy,:,1)=atomic_stress_mbpol(ioxy,:,1)+(1.d0-gammattm)*atomic_stress_mbpol(mttm2,:,1)
    atomic_stress_mbpol(ioxy,:,2)=atomic_stress_mbpol(ioxy,:,2)+(1.d0-gammattm)*atomic_stress_mbpol(mttm2,:,2)
    atomic_stress_mbpol(ioxy,:,3)=atomic_stress_mbpol(ioxy,:,3)+(1.d0-gammattm)*atomic_stress_mbpol(mttm2,:,3)
    atomic_stress_mbpol(ih1,:,1)=atomic_stress_mbpol(ih1,:,1)+gammattm/2.d0*atomic_stress_mbpol(mttm2,:,1)
    atomic_stress_mbpol(ih1,:,2)=atomic_stress_mbpol(ih1,:,2)+gammattm/2.d0*atomic_stress_mbpol(mttm2,:,2)
    atomic_stress_mbpol(ih1,:,3)=atomic_stress_mbpol(ih1,:,3)+gammattm/2.d0*atomic_stress_mbpol(mttm2,:,3)
    atomic_stress_mbpol(ih2,:,1)=atomic_stress_mbpol(ih2,:,1)+gammattm/2.d0*atomic_stress_mbpol(mttm2,:,1)
    atomic_stress_mbpol(ih2,:,2)=atomic_stress_mbpol(ih2,:,2)+gammattm/2.d0*atomic_stress_mbpol(mttm2,:,2)
    atomic_stress_mbpol(ih2,:,3)=atomic_stress_mbpol(ih2,:,3)+gammattm/2.d0*atomic_stress_mbpol(mttm2,:,3)
    atomic_stress_mbpol(mttm2,:,:)=0
    !
    atomic_stress_srf(ioxy,:,1)=atomic_stress_srf(ioxy,:,1)+(1.d0-gammattm)*atomic_stress_srf(mttm2,:,1)
    atomic_stress_srf(ioxy,:,2)=atomic_stress_srf(ioxy,:,2)+(1.d0-gammattm)*atomic_stress_srf(mttm2,:,2)
    atomic_stress_srf(ioxy,:,3)=atomic_stress_srf(ioxy,:,3)+(1.d0-gammattm)*atomic_stress_srf(mttm2,:,3)
    atomic_stress_srf(ih1,:,1)=atomic_stress_srf(ih1,:,1)+gammattm/2.d0*atomic_stress_srf(mttm2,:,1)
    atomic_stress_srf(ih1,:,2)=atomic_stress_srf(ih1,:,2)+gammattm/2.d0*atomic_stress_srf(mttm2,:,2)
    atomic_stress_srf(ih1,:,3)=atomic_stress_srf(ih1,:,3)+gammattm/2.d0*atomic_stress_srf(mttm2,:,3)
    atomic_stress_srf(ih2,:,1)=atomic_stress_srf(ih2,:,1)+gammattm/2.d0*atomic_stress_srf(mttm2,:,1)
    atomic_stress_srf(ih2,:,2)=atomic_stress_srf(ih2,:,2)+gammattm/2.d0*atomic_stress_srf(mttm2,:,2)
    atomic_stress_srf(ih2,:,3)=atomic_stress_srf(ih2,:,3)+gammattm/2.d0*atomic_stress_srf(mttm2,:,3)
    atomic_stress_srf(mttm2,:,:)=0
    !
    atomic_stress_lrcorrect(ioxy,:,1)=atomic_stress_lrcorrect(ioxy,:,1)+&
    (1.d0-gammattm)*atomic_stress_lrcorrect(mttm2,:,1)
    atomic_stress_lrcorrect(ioxy,:,2)=atomic_stress_lrcorrect(ioxy,:,2)+&
    (1.d0-gammattm)*atomic_stress_lrcorrect(mttm2,:,2)
    atomic_stress_lrcorrect(ioxy,:,3)=atomic_stress_lrcorrect(ioxy,:,3)+&
    (1.d0-gammattm)*atomic_stress_lrcorrect(mttm2,:,3)
    atomic_stress_lrcorrect(ih1,:,1)=atomic_stress_lrcorrect(ih1,:,1)+gammattm/2.d0*atomic_stress_lrcorrect(mttm2,:,1)
    atomic_stress_lrcorrect(ih1,:,2)=atomic_stress_lrcorrect(ih1,:,2)+gammattm/2.d0*atomic_stress_lrcorrect(mttm2,:,2)
    atomic_stress_lrcorrect(ih1,:,3)=atomic_stress_lrcorrect(ih1,:,3)+gammattm/2.d0*atomic_stress_lrcorrect(mttm2,:,3)
    atomic_stress_lrcorrect(ih2,:,1)=atomic_stress_lrcorrect(ih2,:,1)+gammattm/2.d0*atomic_stress_lrcorrect(mttm2,:,1)
    atomic_stress_lrcorrect(ih2,:,2)=atomic_stress_lrcorrect(ih2,:,2)+gammattm/2.d0*atomic_stress_lrcorrect(mttm2,:,2)
    atomic_stress_lrcorrect(ih2,:,3)=atomic_stress_lrcorrect(ih2,:,3)+gammattm/2.d0*atomic_stress_lrcorrect(mttm2,:,3)
    atomic_stress_lrcorrect(mttm2,:,:)=0

  end subroutine distribute_stress
  !
  subroutine distribute_energy(ioxy,ih1,ih2,mttm2,gammattm)
    integer :: ioxy, ih1, ih2, mttm2
    real(8) :: gammattm
    !
    atomic_energy_ang(ioxy)=atomic_energy_ang(ioxy)+(1.d0-gammattm)*atomic_energy_ang(mttm2)
    atomic_energy_ang(ih1)=atomic_energy_ang(ih1)+gammattm/2.d0*atomic_energy_ang(mttm2)
    atomic_energy_ang(ih2)=atomic_energy_ang(ih2)+gammattm/2.d0*atomic_energy_ang(mttm2)
    atomic_energy_ang(mttm2)=0
    !
    atomic_energy_bnd(ioxy)=atomic_energy_bnd(ioxy)+(1.d0-gammattm)*atomic_energy_bnd(mttm2)
    atomic_energy_bnd(ih1)=atomic_energy_bnd(ih1)+gammattm/2.d0*atomic_energy_bnd(mttm2)
    atomic_energy_bnd(ih2)=atomic_energy_bnd(ih2)+gammattm/2.d0*atomic_energy_bnd(mttm2)
    atomic_energy_bnd(mttm2)=0
    !
    atomic_energy_ew1p(ioxy)=atomic_energy_ew1p(ioxy)+(1.d0-gammattm)*atomic_energy_ew1p(mttm2)
    atomic_energy_ew1p(ih1)=atomic_energy_ew1p(ih1)+gammattm/2.d0*atomic_energy_ew1p(mttm2)
    atomic_energy_ew1p(ih2)=atomic_energy_ew1p(ih2)+gammattm/2.d0*atomic_energy_ew1p(mttm2)
    atomic_energy_ew1p(mttm2)=0
    !
    atomic_energy_ew2p(ioxy)=atomic_energy_ew2p(ioxy)+(1.d0-gammattm)*atomic_energy_ew2p(mttm2)
    atomic_energy_ew2p(ih1)=atomic_energy_ew2p(ih1)+gammattm/2.d0*atomic_energy_ew2p(mttm2)
    atomic_energy_ew2p(ih2)=atomic_energy_ew2p(ih2)+gammattm/2.d0*atomic_energy_ew2p(mttm2)
    atomic_energy_ew2p(mttm2)=0
    !
    atomic_energy_ew3p(ioxy)=atomic_energy_ew3p(ioxy)+(1.d0-gammattm)*atomic_energy_ew3p(mttm2)
    atomic_energy_ew3p(ih1)=atomic_energy_ew3p(ih1)+gammattm/2.d0*atomic_energy_ew3p(mttm2)
    atomic_energy_ew3p(ih2)=atomic_energy_ew3p(ih2)+gammattm/2.d0*atomic_energy_ew3p(mttm2)
    atomic_energy_ew3p(mttm2)=0
    !
    atomic_energy_polar(ioxy)=atomic_energy_polar(ioxy)+(1.d0-gammattm)*atomic_energy_polar(mttm2)
    atomic_energy_polar(ih1)=atomic_energy_polar(ih1)+gammattm/2.d0*atomic_energy_polar(mttm2)
    atomic_energy_polar(ih2)=atomic_energy_polar(ih2)+gammattm/2.d0*atomic_energy_polar(mttm2)
    atomic_energy_polar(mttm2)=0
    !
    atomic_energy_mbpol(ioxy)=atomic_energy_mbpol(ioxy)+(1.d0-gammattm)*atomic_energy_mbpol(mttm2)
    atomic_energy_mbpol(ih1)=atomic_energy_mbpol(ih1)+gammattm/2.d0*atomic_energy_mbpol(mttm2)
    atomic_energy_mbpol(ih2)=atomic_energy_mbpol(ih2)+gammattm/2.d0*atomic_energy_mbpol(mttm2)
    atomic_energy_mbpol(mttm2)=0
    !
    atomic_energy_srf(ioxy)=atomic_energy_srf(ioxy)+(1.d0-gammattm)*atomic_energy_srf(mttm2)
    atomic_energy_srf(ih1)=atomic_energy_srf(ih1)+gammattm/2.d0*atomic_energy_srf(mttm2)
    atomic_energy_srf(ih2)=atomic_energy_srf(ih2)+gammattm/2.d0*atomic_energy_srf(mttm2)
    atomic_energy_srf(mttm2)=0
    !
    atomic_energy_lrcorrect(ioxy)=atomic_energy_lrcorrect(ioxy)+(1.d0-gammattm)*atomic_energy_lrcorrect(mttm2)
    atomic_energy_lrcorrect(ih1)=atomic_energy_lrcorrect(ih1)+gammattm/2.d0*atomic_energy_lrcorrect(mttm2)
    atomic_energy_lrcorrect(ih2)=atomic_energy_lrcorrect(ih2)+gammattm/2.d0*atomic_energy_lrcorrect(mttm2)
    atomic_energy_lrcorrect(mttm2)=0

  end subroutine distribute_energy
  ! ____________________________________________________________________________________________________________________
  ! _______________________________________ Compute and write the heat current _________________________________________
  ! ____________________________________________________________________________________________________________________
  !
  subroutine compute_heat_flux(vxx,vyy,vzz)
    real(8), dimension(1:mxatms), intent(in) :: vxx, vyy, vzz
    integer :: iatm, icart, jcart
    real(8) :: vtmp(3), jmpi(6), eng
    j1 = 0.d0
    j2 = 0.d0
    jtot = 0.d0
    do iatm=1,mxatms
      vtmp = (/ vxx(iatm), vyy(iatm), vzz(iatm) /)
      eng = (atomic_kinetic_energy(iatm)+atomic_potential_energy(iatm))
      j1 = j1 + eng*vtmp
      j2(1) = j2(1) + atomic_stress(iatm,1,1)*vtmp(1) + atomic_stress(iatm,1,2)*vtmp(2)&
      + atomic_stress(iatm,1,3)*vtmp(3)
      j2(2) = j2(2) + atomic_stress(iatm,1,2)*vtmp(1) + atomic_stress(iatm,2,2)*vtmp(2)&
      + atomic_stress(iatm,2,3)*vtmp(3)
      j2(3) = j2(3) + atomic_stress(iatm,1,3)*vtmp(1) + atomic_stress(iatm,2,3)*vtmp(2)&
      + atomic_stress(iatm,3,3)*vtmp(3)
    enddo
    jtot = j1 - j2
  end subroutine compute_heat_flux
  !
  subroutine write_heat_flux(time,bead_suffix)
    real(8), intent(in) :: time
    character(8), intent(in) :: bead_suffix
    integer :: flux_file = 7000
    open(flux_file, file='heatflux'//bead_suffix,position='append')
    write(flux_file,'(es12.5)') time
    write(flux_file,'(3es15.7)') j1*junit
    write(flux_file,'(3es15.7)') j2*junit
    write(flux_file,'(3es15.7)') jtot*junit
    close(flux_file)
  end subroutine write_heat_flux
  !
  ! ====================================================================================================================
  ! ======================================== END OF THE LAMMPS-LIKE METHOD =============================================
  ! ====================================================================================================================
  !
  ! ====================================================================================================================
  ! ====================== SECOND METHOD: FORCE DECOMPOSITION & TIME INTEGRATION OF POTENTIAL ENERGY ===================
  ! ====================================================================================================================
  !
  ! ____________________________________________________________________________________________________________________
  ! _______________________________________ Initialize potential energy ________________________________________________
  ! ____________________________________________________________________________________________________________________
  subroutine initialize_potential_energy(epot_ini,vxx,vyy,vzz)
    real(8), intent(in) :: epot_ini
    real(8), dimension(1:mxatms), intent(in) :: vxx, vyy, vzz
    real(8) :: eng_ew1, force_ij(3), velocity_ij(3)
    integer :: iatm,jatm
    call total_energy_ew1p(eng_ew1)
    atomic_potential_energy2(1:3*mxatms/4) = atomic_potential_energy(1:3*mxatms/4)-atomic_energy_ew1p(1:3*mxatms/4)
    atomic_potential_energy2(3*mxatms/4+1:mxatms) = 0.d0
    do iatm = 1, 3*mxatms/4
      do jatm = 1, 3*mxatms/4
        if (jatm == iatm) cycle
        force_ij = force_matrix(iatm,jatm,1:3)
        velocity_ij = (/vxx(iatm)-vxx(jatm), vyy(iatm)-vyy(jatm), vzz(iatm)-vzz(jatm)/)
        f_dot_v_old(iatm,jatm) = dot_product(force_ij,velocity_ij)
      end do
    end do
  end subroutine initialize_potential_energy
  ! ____________________________________________________________________________________________________________________
  ! _______________________________________ Build force matrix _________________________________________________________
  ! ____________________________________________________________________________________________________________________
  !
  subroutine update_forces(iatom,jatom,forceij)
    integer, intent(in) :: iatom,jatom
    real(8), intent(in) :: forceij(3)
    force_matrix(iatom,jatom,:) = force_matrix(iatom,jatom,:) + forceij
  end subroutine update_forces
  !
  subroutine update_force_ew1p(iatm,force)
    integer, intent(in) :: iatm
    real(8), intent(in) :: force(3)
    force_ew1p(iatm,:) = force_ew1p(iatm,:) + force
  end subroutine update_force_ew1p
  !
  subroutine distribute_forces(ioxy,ih1,ih2,mttm2,gammattm)
    integer, intent(in) :: ioxy,ih1,ih2,mttm2
    real(8), intent(in) :: gammattm
    integer :: j
    do j = 1, mxatms
      force_matrix(ioxy,j,:) = force_matrix(ioxy,j,:) + (1-gammattm)*force_matrix(mttm2,j,:)
      force_matrix(ih1,j,:) = force_matrix(ih1,j,:) + 0.5d0*gammattm*force_matrix(mttm2,j,:)
      force_matrix(ih2,j,:) = force_matrix(ih2,j,:) + 0.5d0*gammattm*force_matrix(mttm2,j,:)
      force_matrix(mttm2,j,:) = 0.d0
    end do
    force_ew1p(ioxy,:) = force_ew1p(ioxy,:) + (1-gammattm)*force_ew1p(mttm2,:)
    force_ew1p(ih1,:) = force_ew1p(ih1,:) + 0.5d0*gammattm*force_ew1p(mttm2,:)
    force_ew1p(ih2,:) = force_ew1p(ih2,:) + 0.5d0*gammattm*force_ew1p(mttm2,:)
    force_ew1p(mttm2,:) = 0.d0
  end subroutine distribute_forces
  !
  ! ____________________________________________________________________________________________________________________
  ! ________________________________________ Distribute the potential energy ___________________________________________
  ! ____________________________________________________________________________________________________________________
  !
  subroutine integrate_force_velocity(iatm,jatm,f_dot_v,int_force_velocity)
    integer, intent(in) :: iatm,jatm
    real(8), intent(in) :: f_dot_v
    real(8), intent(inout) :: int_force_velocity
    int_force_velocity = 0.5d0*(f_dot_v + f_dot_v_old(iatm,jatm))
  end subroutine integrate_force_velocity
  !
  subroutine potential_energy_per_atom2(dt,vxx,vyy,vzz)
    real(8), intent(in) :: dt
    real(8), dimension(1:mxatms), intent(in) :: vxx, vyy, vzz
    integer :: iatm, jatm
    real(8) :: force_ij(3), velocity_ij(3), f_dot_v, int_force_velocity
    delta_potential_energy = 0.d0
    do iatm = 1, 3*mxatms/4
      do jatm = 1, 3*mxatms/4
        if (jatm == iatm) cycle
        force_ij = force_matrix(iatm,jatm,:)
        velocity_ij = (/vxx(iatm)-vxx(jatm), vyy(iatm)-vyy(jatm), vzz(iatm)-vzz(jatm)/)
        f_dot_v = dot_product(force_ij,velocity_ij)
        call integrate_force_velocity(iatm,jatm,f_dot_v,int_force_velocity)
        delta_potential_energy(iatm) = delta_potential_energy(iatm) + int_force_velocity
        f_dot_v_old(iatm,jatm) = f_dot_v
      end do
      delta_potential_energy(iatm) = -0.5d0*dt*delta_potential_energy(iatm)
      !atomic_potential_energy2(iatm) = atomic_potential_energy2(iatm) + delta_potential_energy(iatm) !+ &
      !atomic_energy_ew1p(iatm)
    end do
  end subroutine potential_energy_per_atom2
  !
  ! PP_: SHOULD BE DELETED
  subroutine distribute_energy2(ioxy,ih1,ih2,mttm2,gammattm)
      integer :: ioxy, ih1, ih2, mttm2
      real(8) :: gammattm
      !
      delta_potential_energy(ioxy)=delta_potential_energy(ioxy)+(1.d0-gammattm)*delta_potential_energy(mttm2)
      delta_potential_energy(ih1)=delta_potential_energy(ih1)+gammattm/2.d0*delta_potential_energy(mttm2)
      delta_potential_energy(ih2)=delta_potential_energy(ih2)+gammattm/2.d0*delta_potential_energy(mttm2)
      delta_potential_energy(mttm2)=0
  end subroutine distribute_energy2
  !
  ! ____________________________________________________________________________________________________________________
  ! _____________________________________ Compute and write the heat current ___________________________________________
  ! ____________________________________________________________________________________________________________________

  ! subroutine compute_heat_flux2(xxx,yyy,zzz,vxx,vyy,vzz,gammattm)
  !   real(8), dimension(1:mxatms), intent(in) :: xxx,yyy,zzz,vxx,vyy,vzz
  !   real(8), intent(in) :: gammattm
  !   integer :: imol,jmol,ioxy,ih1,ih2,mttm2,joxy,jh1,jh2
  !   real(8) :: vtmp(3), eng, fij(3), rij(3)
  !   jkin = 0.d0
  !   jpot = 0.d0
  !   do imol = 1, mxatms/4
  !     ioxy = 3*imol-2
  !     ih1 = 3*imol-1
  !     ih2 = 3*imol
  !     mttm2 = 3*mxatms/4+imol
  !     !call distribute_energy2(ioxy,ih1,ih2,mttm2,gammattm)
  !
  !     vtmp = (/ vxx(ioxy), vyy(ioxy), vzz(ioxy) /)
  !     !eng = (atomic_kinetic_energy(ioxy)+atomic_potential_energy2(ioxy)+delta_potential_energy(ioxy)+&
  !     !atomic_energy_ew1p(ioxy))
  !     eng = atomic_kinetic_energy(ioxy)+atomic_potential_energy(ioxy)
  !     jkin = jkin + eng*vtmp
  !     do jmol = 1, mxatms/4
  !       joxy = 3*jmol-2
  !       jh1 = 3*jmol-1
  !       jh2 = 3*jmol
  !       if (joxy == ioxy) cycle
  !       fij = force_matrix(ioxy,joxy,:)
  !       rij = (/xxx(ioxy)-xxx(joxy),yyy(ioxy)-yyy(joxy),zzz(ioxy)-zzz(joxy)/)
  !       jpot = jpot + 0.5d0*rij*dot_product(fij,vtmp)
  !       fij = force_matrix(ioxy,jh1,:)
  !       rij = (/xxx(ioxy)-xxx(jh1),yyy(ioxy)-yyy(jh1),zzz(ioxy)-zzz(jh1)/)
  !       jpot = jpot + 0.5d0*rij*dot_product(fij,vtmp)
  !       fij = force_matrix(ioxy,jh2,:)
  !       rij = (/xxx(ioxy)-xxx(jh2),yyy(ioxy)-yyy(jh2),zzz(ioxy)-zzz(jh2)/)
  !       jpot = jpot + 0.5d0*rij*dot_product(fij,vtmp)
  !     end do
  !     jpot(1) = jpot(1) + atomic_stress_ew1p(ioxy,1,1)*vtmp(1) + atomic_stress_ew1p(ioxy,1,2)*vtmp(2)&
  !     + atomic_stress_ew1p(ioxy,1,3)*vtmp(3)
  !     jpot(2) = jpot(2) + atomic_stress_ew1p(ioxy,1,2)*vtmp(1) + atomic_stress_ew1p(ioxy,2,2)*vtmp(2)&
  !     + atomic_stress_ew1p(ioxy,2,3)*vtmp(3)
  !     jpot(3) = jpot(3) + atomic_stress_ew1p(ioxy,1,3)*vtmp(1) + atomic_stress_ew1p(ioxy,2,3)*vtmp(2)&
  !     + atomic_stress_ew1p(ioxy,3,3)*vtmp(3)
  !
  !     vtmp = (/ vxx(ih1), vyy(ih1), vzz(ih1) /)
  !     !eng = (atomic_kinetic_energy(ih1)+atomic_potential_energy2(ih1)+delta_potential_energy(ih1)&
  !     !+atomic_energy_ew1p(ih1))
  !     eng = atomic_kinetic_energy(ih1)+atomic_potential_energy(ih1)
  !     jkin = jkin + eng*vtmp
  !     do jmol = 1, mxatms/4
  !       joxy = 3*jmol-2
  !       jh1 = 3*jmol-1
  !       jh2 = 3*jmol
  !       if (jh1 == ih1) cycle
  !       fij = force_matrix(ih1,joxy,:)
  !       rij = (/xxx(ih1)-xxx(joxy),yyy(ih1)-yyy(joxy),zzz(ih1)-zzz(joxy)/)
  !       jpot = jpot + 0.5d0*rij*dot_product(fij,vtmp)
  !       fij = force_matrix(ih1,jh1,:)
  !       rij = (/xxx(ih1)-xxx(jh1),yyy(ih1)-yyy(jh1),zzz(ih1)-zzz(jh1)/)
  !       jpot = jpot + 0.5d0*rij*dot_product(fij,vtmp)
  !       fij = force_matrix(ih1,jh2,:)
  !       rij = (/xxx(ih1)-xxx(jh2),yyy(ih1)-yyy(jh2),zzz(ih1)-zzz(jh2)/)
  !       jpot = jpot + 0.5d0*rij*dot_product(fij,vtmp)
  !     end do
  !     jpot(1) = jpot(1) + atomic_stress_ew1p(ih1,1,1)*vtmp(1) + atomic_stress_ew1p(ih1,1,2)*vtmp(2)&
  !     + atomic_stress_ew1p(ih1,1,3)*vtmp(3)
  !     jpot(2) = jpot(2) + atomic_stress_ew1p(ih1,1,2)*vtmp(1) + atomic_stress_ew1p(ih1,2,2)*vtmp(2)&
  !     + atomic_stress_ew1p(ih1,2,3)*vtmp(3)
  !     jpot(3) = jpot(3) + atomic_stress_ew1p(ih1,1,3)*vtmp(1) + atomic_stress_ew1p(ih1,2,3)*vtmp(2)&
  !     + atomic_stress_ew1p(ih1,3,3)*vtmp(3)
  !
  !     vtmp = (/ vxx(ih2), vyy(ih2), vzz(ih2) /)
  !     !eng = atomic_kinetic_energy(ih2)+atomic_potential_energy2(ih2)+delta_potential_energy(ih2)&
  !     !+atomic_energy_ew1p(ih2)
  !     eng = atomic_kinetic_energy(ih2)+atomic_potential_energy(ih2)
  !     jkin = jkin + eng*vtmp
  !     do jmol = 1, mxatms/4
  !       joxy = 3*jmol-2
  !       jh1 = 3*jmol-1
  !       jh2 = 3*jmol
  !       if (jh2 == ih2) cycle
  !       fij = force_matrix(ih2,joxy,:)
  !       rij = (/xxx(ih2)-xxx(joxy),yyy(ih2)-yyy(joxy),zzz(ih2)-zzz(joxy)/)
  !       jpot = jpot + 0.5d0*rij*dot_product(fij,vtmp)
  !       fij = force_matrix(ih2,jh1,:)
  !       rij = (/xxx(ih2)-xxx(jh1),yyy(ih2)-yyy(jh1),zzz(ih2)-zzz(jh1)/)
  !       jpot = jpot + 0.5d0*rij*dot_product(fij,vtmp)
  !       fij = force_matrix(ih2,jh2,:)
  !       rij = (/xxx(ih2)-xxx(jh2),yyy(ih2)-yyy(jh2),zzz(ih2)-zzz(jh2)/)
  !       jpot = jpot + 0.5d0*rij*dot_product(fij,vtmp)
  !     end do
  !     jpot(1) = jpot(1) + atomic_stress_ew1p(ih2,1,1)*vtmp(1) + atomic_stress_ew1p(ih2,1,2)*vtmp(2) &
  !     + atomic_stress_ew1p(ih2,1,3)*vtmp(3)
  !     jpot(2) = jpot(2) + atomic_stress_ew1p(ih2,1,2)*vtmp(1) + atomic_stress_ew1p(ih2,2,2)*vtmp(2) &
  !     + atomic_stress_ew1p(ih2,2,3)*vtmp(3)
  !     jpot(3) = jpot(3) + atomic_stress_ew1p(ih2,1,3)*vtmp(1) + atomic_stress_ew1p(ih2,2,3)*vtmp(2) &
  !     + atomic_stress_ew1p(ih2,3,3)*vtmp(3)
  !   end do
  !   jtot2 = jkin + jpot
  ! end subroutine compute_heat_flux2

  subroutine compute_heat_flux2(xxx,yyy,zzz,vxx,vyy,vzz,gammattm)
    real(8), dimension(1:mxatms), intent(in) :: xxx,yyy,zzz,vxx,vyy,vzz
    real(8), intent(in) :: gammattm
    integer :: iatm, jatm
    real(8) :: vtmp(3), eng, fij(3), rij(3)
    jkin = 0.d0
    jpot = 0.d0
    jpot_kspace = 0.d0
    do iatm = 1, 3*mxatms/4
      vtmp = (/ vxx(iatm), vyy(iatm), vzz(iatm) /)
      eng = atomic_kinetic_energy(iatm)+atomic_potential_energy(iatm)
      jkin = jkin + eng*vtmp
      do jatm = 1, 3*mxatms/4
        if (iatm == jatm) cycle
        fij = force_matrix(iatm,jatm,:)
        rij = (/xxx(iatm)-xxx(jatm),yyy(iatm)-yyy(jatm),zzz(iatm)-zzz(jatm)/)
        jpot = jpot + rij*dot_product(fij,vtmp)
      end do
      jpot_kspace(1) = jpot_kspace(1) - (atomic_stress_ew1p(iatm,1,1)*vtmp(1) + atomic_stress_ew1p(iatm,1,2)*vtmp(2)&
      + atomic_stress_ew1p(iatm,1,3)*vtmp(3))
      jpot_kspace(2) = jpot_kspace(2) - (atomic_stress_ew1p(iatm,1,2)*vtmp(1) + atomic_stress_ew1p(iatm,2,2)*vtmp(2)&
      + atomic_stress_ew1p(iatm,2,3)*vtmp(3))
      jpot_kspace(3) = jpot_kspace(3) - (atomic_stress_ew1p(iatm,1,3)*vtmp(1) + atomic_stress_ew1p(iatm,2,3)*vtmp(2)&
      + atomic_stress_ew1p(iatm,3,3)*vtmp(3))
    end do
    jtot2 = jkin + jpot + jpot_kspace
  end subroutine compute_heat_flux2
  !
  subroutine write_heat_flux2(time,bead_suffix)
    real(8), intent(in) :: time
    character(8), intent(in) :: bead_suffix
    integer :: flux_file = 7001
    open(flux_file, file='heatflux2'//bead_suffix,position='append')
    write(flux_file,'(es12.5)') time
    write(flux_file,'(3es15.7)') jkin*junit
    write(flux_file,'(3es15.7)') jpot*junit
    write(flux_file,'(3es15.7)') jpot_kspace*junit
    write(flux_file,'(3es15.7)') jtot2*junit
    close(flux_file)
  end subroutine write_heat_flux2
  !
  subroutine compute_stress(time, xxx, yyy, zzz, unit)
    real(8), intent(in) :: xxx(1:mxatms), yyy(1:mxatms), zzz(1:mxatms), time, unit
    real(8) :: stress_hardy(1:9), rij(3), st_ew1p(1:3,1:3), st_ew1p19(1:9)
    integer :: i, j, file_stress = 45
    stress_hardy = 0.d0
    do i = 1, 3*mxatms/4
      do j = 1, 3*mxatms/4
        if (i == j) cycle
        rij = (/xxx(i) - xxx(j), yyy(i) - yyy(j), zzz(i) - zzz(j)/)
        ! stress_hardy(1) = stress_hardy(1) + force_matrix(i, j, 1)*rij(1)
        ! stress_hardy(2) = stress_hardy(2) + force_matrix(i, j, 1)*rij(2)
        ! stress_hardy(3) = stress_hardy(3) + force_matrix(i, j, 1)*rij(3)
        ! stress_hardy(4) = stress_hardy(4) + force_matrix(i, j, 2)*rij(1)
        ! stress_hardy(5) = stress_hardy(5) + force_matrix(i, j, 2)*rij(2)
        ! stress_hardy(6) = stress_hardy(6) + force_matrix(i, j, 2)*rij(3)
        ! stress_hardy(7) = stress_hardy(7) + force_matrix(i, j, 3)*rij(1)
        ! stress_hardy(8) = stress_hardy(8) + force_matrix(i, j, 3)*rij(2)
        ! stress_hardy(9) = stress_hardy(9) + force_matrix(i, j, 3)*rij(3)
        stress_hardy(1) = stress_hardy(1) + force_matrix(i, j, 1)*rij(1)
        stress_hardy(2) = stress_hardy(2) + force_matrix(i, j, 2)*rij(1)
        stress_hardy(3) = stress_hardy(3) + force_matrix(i, j, 3)*rij(1)
        stress_hardy(4) = stress_hardy(4) + force_matrix(i, j, 1)*rij(2)
        stress_hardy(5) = stress_hardy(5) + force_matrix(i, j, 2)*rij(2)
        stress_hardy(6) = stress_hardy(6) + force_matrix(i, j, 3)*rij(2)
        stress_hardy(7) = stress_hardy(7) + force_matrix(i, j, 1)*rij(3)
        stress_hardy(8) = stress_hardy(8) + force_matrix(i, j, 2)*rij(3)
        stress_hardy(9) = stress_hardy(9) + force_matrix(i, j, 3)*rij(3)
      end do
    end do
    stress_hardy = -0.5d0*stress_hardy
    call total_stress_ew1p(st_ew1p)
    st_ew1p19(1:3) = st_ew1p(1:3,1)
    st_ew1p19(4:6) = st_ew1p(1:3,2)
    st_ew1p19(7:9) = st_ew1p(1:3,3)
    open(file_stress,file='STRESS_HARDY', position='append')
    write(file_stress,'(f12.5,9f10.5)') time, stress_hardy*unit
    close(file_stress)
  end subroutine compute_stress
  !
  ! ____________________________________________________________________________________________________________________
  ! _______________________________ Save the last energy configuration for restart _____________________________________
  ! ____________________________________________________________________________________________________________________
  ! PP_: Is this really useful?
  subroutine save_last_energy(step,time,bead_suffix)
    character(8), intent(in) :: bead_suffix
    integer, intent(in) :: step
    real(8), intent(in) :: time
    integer :: energy_file = 1020, i, iox, ih1, ih2
    open(energy_file,file='ENERGY_CONFIG'//bead_suffix)
    write(energy_file,*) ''
    write(energy_file,*) step, time
    do i = 1, mxatms/4
      iox = 3*i-2
      ih1 = 3*i-1
      ih2 = 3*i
      write(energy_file,*) atomic_potential_energy2(iox)
      write(energy_file,*) atomic_potential_energy2(ih1)
      write(energy_file,*) atomic_potential_energy2(ih2)
    end do
    close(energy_file)
  end subroutine save_last_energy
  !
  ! ====================================================================================================================
  ! ========================================= END OF THE SECOND METHOD =================================================
  ! ====================================================================================================================

  ! ====================================================================================================================
  ! ============================================== Check and debug =====================================================
  ! ====================================================================================================================
  !
  subroutine write_force_matrix(time,bead_suffix)
    real(8), intent(in) :: time
    character(8), intent(in) :: bead_suffix
    integer :: FMfile = 1001, i, natoms
    open(FMfile,file='FORCE_MATRIX'//bead_suffix,position='append')
    write(FMfile,*) time
    do i = 1, 3*mxatms/4
      write(FMfile,*) sum(force_matrix(i,:,1))+force_ew1p(i,1),sum(force_matrix(i,:,2))+force_ew1p(i,2),&
      sum(force_matrix(i,:,3))+force_ew1p(i,3)

      !write(FMfile,*) sum(force_matrix(i,:,1)),sum(force_matrix(i,:,2)),sum(force_matrix(i,:,3))
    end do
    ! open(FMfile  ,file='FORCE_MATRIX_X'//bead_suffix,position='append')
    ! open(FMfile+1,file='FORCE_MATRIX_Y'//bead_suffix,position='append')
    ! open(FMfile+2,file='FORCE_MATRIX_Z'//bead_suffix,position='append')
    !
    ! write(FMfile,*) time
    ! write(FMfile+1,*) time
    ! write(FMfile+2,*) time
    !
    ! natoms = 3*mxatms/4
    ! do i = 1, natoms
    !   write(FMfile,*) force_matrix(i,1:natoms,1)
    !   write(FMfile+1,'(768f15.5)') force_matrix(i,1:natoms,2)
    !   write(FMfile+2,'(768f15.5)') force_matrix(i,1:natoms,3)
    ! end do

  end subroutine write_force_matrix
  !
  subroutine check_forces(time,fxx)
    real(8), intent(in) :: fxx(mxatms), time
    integer :: natom
    open(1010,file='CHECK_FORCE.dat',position='append')
    write(1010,*) time
    do natom = 1, mxatms
      write(1010,*) natom, fxx(natom)
      write(1010,*) natom, sum(force_matrix(natom,:,1))! + force_ew1p(natom,1)
      write(1010,*) natom, force_ew1p(natom,1)
    end do
    close(1010)
  end subroutine check_forces
  !
  subroutine write_stress(time)
    real(8) :: time
    open(5000,file='atomic_stress.dat',position='append')
    !
    write(5000,'(A)') 'ang'
    write(5000,*) atomic_stress_ang
    !
    write(5000,*) 'bnd'
    write(5000,'(4608f15.5)') atomic_stress_bnd
    !
    write(5000,*) 'ew1p'
    write(5000,'(4608f15.5)') atomic_stress_ew1p
    !
    write(5000,*) 'ew2p'
    write(5000,'(4608f15.5)') atomic_stress_ew2p
    !
    write(5000,*) 'ew3p'
    write(5000,'(4608f15.5)') atomic_stress_ew3p
    !
    write(5000,*) 'qd'
    write(5000,'(4608f15.5)') atomic_stress_qd
    !
    write(5000,*) 'mbpol'
    write(5000,'(4608f15.5)') atomic_stress_mbpol
    !
    write(5000,*) 'srf'
    write(5000,'(4608f15.5)') atomic_stress_srf
    !
    write(5000,*) 'lrcorrect'
    write(5000,'(4608f15.5)') atomic_stress_lrcorrect
    close(5000)
  end subroutine write_stress
  !
  subroutine write_energy(time)
    real(8) :: time
    open(5000,file='atomic_energy.dat',position='append')
    !
    write(5000,*) 'ang'
    write(5000,'(512f15.5)') atomic_energy_ang
    !
    write(5000,*) 'bnd'
    write(5000,'(512f15.5)') atomic_energy_bnd
    !
    write(5000,*) 'ew1p'
    write(5000,'(512f15.5)') atomic_energy_ew1p
    !
    write(5000,*) 'ew2p'
    write(5000,'(512f15.5)') atomic_energy_ew2p
    !
    write(5000,*) 'ew3p'
    write(5000,'(512f15.5)') atomic_energy_ew3p
    !
    write(5000,*) 'polar'
    write(5000,'(512f15.5)') atomic_energy_polar
    !
    write(5000,*) 'mbpol'
    write(5000,'(512f15.5)') atomic_energy_mbpol
    !
    write(5000,*) 'srf'
    write(5000,'(512f15.5)') atomic_energy_srf
    !
    write(5000,*) 'lrcorrect'
    write(5000,'(512f15.5)') atomic_energy_lrcorrect
    close(5000)
  end subroutine write_energy
  !
  subroutine total_energy_ang(entot)
    real(8) :: entot
    entot = sum(atomic_energy_ang)
  end subroutine total_energy_ang
  !
  subroutine total_stress_ang(stresstot)
    real(8) :: stresstot(1:3,1:3)
    stresstot(1,1) = sum(atomic_stress_ang(:,1,1))
    stresstot(1,2) = sum(atomic_stress_ang(:,1,2))
    stresstot(1,3) = sum(atomic_stress_ang(:,1,3))
    stresstot(2,1) = stresstot(1,2)
    stresstot(2,2) = sum(atomic_stress_ang(:,2,2))
    stresstot(2,3) = sum(atomic_stress_ang(:,2,3))
    stresstot(3,1) = stresstot(1,3)
    stresstot(3,2) = stresstot(2,3)
    stresstot(3,3) = sum(atomic_stress_ang(:,3,3))
  end subroutine total_stress_ang
  !
  subroutine total_energy_bnd(entot)
    real(8) :: entot
    entot = sum(atomic_energy_bnd)
  end subroutine total_energy_bnd
  !
  subroutine total_stress_bnd(stresstot)
    real(8) :: stresstot(1:3,1:3)
    stresstot(1,1) = sum(atomic_stress_bnd(:,1,1))
    stresstot(1,2) = sum(atomic_stress_bnd(:,1,2))
    stresstot(1,3) = sum(atomic_stress_bnd(:,1,3))
    stresstot(2,1) = stresstot(1,2)
    stresstot(2,2) = sum(atomic_stress_bnd(:,2,2))
    stresstot(2,3) = sum(atomic_stress_bnd(:,2,3))
    stresstot(3,1) = stresstot(1,3)
    stresstot(3,2) = stresstot(2,3)
    stresstot(3,3) = sum(atomic_stress_bnd(:,3,3))
  end subroutine total_stress_bnd
  !
  subroutine total_energy_ew1p(entot)
    real(8) :: entot
    entot = sum(atomic_energy_ew1p)
  end subroutine total_energy_ew1p
  !
  subroutine total_stress_ew1p(stresstot)
    real(8) :: stresstot(1:3,1:3)
    stresstot(1,1) = sum(atomic_stress_ew1p(:,1,1))
    stresstot(1,2) = sum(atomic_stress_ew1p(:,1,2))
    stresstot(1,3) = sum(atomic_stress_ew1p(:,1,3))
    stresstot(2,1) = stresstot(1,2)
    stresstot(2,2) = sum(atomic_stress_ew1p(:,2,2))
    stresstot(2,3) = sum(atomic_stress_ew1p(:,2,3))
    stresstot(3,1) = stresstot(1,3)
    stresstot(3,2) = stresstot(2,3)
    stresstot(3,3) = sum(atomic_stress_ew1p(:,3,3))
  end subroutine total_stress_ew1p
  !
  subroutine total_energy_ew2p(entot)
    real(8) :: entot
    entot = sum(atomic_energy_ew2p)
  end subroutine total_energy_ew2p
  !
  subroutine total_stress_ew2p(stresstot)
    real(8) :: stresstot(1:3,1:3)
    stresstot(1,1) = sum(atomic_stress_ew2p(:,1,1))
    stresstot(1,2) = sum(atomic_stress_ew2p(:,1,2))
    stresstot(1,3) = sum(atomic_stress_ew2p(:,1,3))
    stresstot(2,1) = stresstot(1,2)
    stresstot(2,2) = sum(atomic_stress_ew2p(:,2,2))
    stresstot(2,3) = sum(atomic_stress_ew2p(:,2,3))
    stresstot(3,1) = stresstot(1,3)
    stresstot(3,2) = stresstot(2,3)
    stresstot(3,3) = sum(atomic_stress_ew2p(:,3,3))
  end subroutine total_stress_ew2p
  !
  subroutine total_energy_ew3p(entot)
    real(8) :: entot
    entot = sum(atomic_energy_ew3p)
  end subroutine total_energy_ew3p
  !
  subroutine total_stress_ew3p(stresstot)
    real(8) :: stresstot(1:3,1:3)
    stresstot(1,1) = sum(atomic_stress_ew3p(:,1,1))
    stresstot(1,2) = sum(atomic_stress_ew3p(:,1,2))
    stresstot(1,3) = sum(atomic_stress_ew3p(:,1,3))
    stresstot(2,1) = stresstot(1,2)
    stresstot(2,2) = sum(atomic_stress_ew3p(:,2,2))
    stresstot(2,3) = sum(atomic_stress_ew3p(:,2,3))
    stresstot(3,1) = stresstot(1,3)
    stresstot(3,2) = stresstot(2,3)
    stresstot(3,3) = sum(atomic_stress_ew3p(:,3,3))
  end subroutine total_stress_ew3p
  !
  subroutine total_energy_cpe(entot)
    real(8) :: entot
    real(8) :: entmp
    entot=0
    call total_energy_ew1p(entmp)
    entot=entot+entmp
    call total_energy_ew2p(entmp)
    entot=entot+entmp
    call total_energy_ew3p(entmp)
    entot=entot+entmp
    call total_energy_polar(entmp)
    entot=entot+entmp
  end subroutine total_energy_cpe
  !
  subroutine total_stress_qd(stresstot)
    real(8) :: stresstot(1:3,1:3)
    stresstot(1,1) = sum(atomic_stress_qd(:,1,1))
    stresstot(1,2) = sum(atomic_stress_qd(:,1,2))
    stresstot(1,3) = sum(atomic_stress_qd(:,1,3))
    stresstot(2,1) = stresstot(1,2)
    stresstot(2,2) = sum(atomic_stress_qd(:,2,2))
    stresstot(2,3) = sum(atomic_stress_qd(:,2,3))
    stresstot(3,1) = stresstot(1,3)
    stresstot(3,2) = stresstot(2,3)
    stresstot(3,3) = sum(atomic_stress_qd(:,3,3))
  end subroutine total_stress_qd
  !
  subroutine total_energy_mbpol(entot)
    real(8) :: entot
    entot = sum(atomic_energy_mbpol)
  end subroutine total_energy_mbpol
  !
  subroutine total_stress_mbpol(stresstot)
    real(8) :: stresstot(1:3,1:3)
    stresstot(1,1) = sum(atomic_stress_mbpol(:,1,1))
    stresstot(1,2) = sum(atomic_stress_mbpol(:,1,2))
    stresstot(1,3) = sum(atomic_stress_mbpol(:,1,3))
    stresstot(2,1) = stresstot(1,2)
    stresstot(2,2) = sum(atomic_stress_mbpol(:,2,2))
    stresstot(2,3) = sum(atomic_stress_mbpol(:,2,3))
    stresstot(3,1) = stresstot(1,3)
    stresstot(3,2) = stresstot(2,3)
    stresstot(3,3) = sum(atomic_stress_mbpol(:,3,3))
  end subroutine total_stress_mbpol
  !
  subroutine total_energy_srf(entot)
    real(8) :: entot
    entot = sum(atomic_energy_srf)
  end subroutine total_energy_srf
  !
  subroutine total_energy_lrcorrect(entot)
    real(8) :: entot
    entot = sum(atomic_energy_lrcorrect)
  end subroutine total_energy_lrcorrect
  !
  subroutine total_stress_srf(stresstot)
    real(8) :: stresstot(1:3,1:3)
    stresstot(1,1) = sum(atomic_stress_srf(:,1,1))
    stresstot(1,2) = sum(atomic_stress_srf(:,1,2))
    stresstot(1,3) = sum(atomic_stress_srf(:,1,3))
    stresstot(2,1) = stresstot(1,2)
    stresstot(2,2) = sum(atomic_stress_srf(:,2,2))
    stresstot(2,3) = sum(atomic_stress_srf(:,2,3))
    stresstot(3,1) = stresstot(1,3)
    stresstot(3,2) = stresstot(2,3)
    stresstot(3,3) = sum(atomic_stress_srf(:,3,3))
  end subroutine total_stress_srf

  subroutine total_stress_lrcorrect(stresstot)
    real(8) :: stresstot(1:3,1:3)
    stresstot(1,1) = sum(atomic_stress_lrcorrect(:,1,1))
    stresstot(1,2) = sum(atomic_stress_lrcorrect(:,1,2))
    stresstot(1,3) = sum(atomic_stress_lrcorrect(:,1,3))
    stresstot(2,1) = stresstot(1,2)
    stresstot(2,2) = sum(atomic_stress_lrcorrect(:,2,2))
    stresstot(2,3) = sum(atomic_stress_lrcorrect(:,2,3))
    stresstot(3,1) = stresstot(1,3)
    stresstot(3,2) = stresstot(2,3)
    stresstot(3,3) = sum(atomic_stress_lrcorrect(:,3,3))
  end subroutine total_stress_lrcorrect
  !
  subroutine total_energy_polar(entot)
    real(8) :: entot
    entot = sum(atomic_energy_polar)
  end subroutine total_energy_polar
  !
  subroutine total_kinetic_energy(entot)
    real(8) :: entot
    entot = sum(atomic_kinetic_energy)
  end subroutine total_kinetic_energy
  !
  subroutine total_potential_energy(entot)
    real(8) :: entot, atomic_energy
    call potential_energy_per_atom()
    entot = sum(atomic_potential_energy)
  end subroutine total_potential_energy
  !
  subroutine total_stress(stresstot)
    real(8) :: stresstot(1:9)
    integer :: i,j
    call stress_per_atom()
    stresstot(1) = sum(atomic_stress(:,1,1))
    stresstot(2) = sum(atomic_stress(:,1,2))
    stresstot(3) = sum(atomic_stress(:,1,3))
    stresstot(4) = stresstot(2)!sum(atomic_stress(:,2,1))
    stresstot(5) = sum(atomic_stress(:,2,2))
    stresstot(6) = sum(atomic_stress(:,2,3))
    stresstot(7) = stresstot(3)!sum(atomic_stress(:,3,1))
    stresstot(8) = stresstot(6)!sum(atomic_stress(:,3,2))
    stresstot(9) = sum(atomic_stress(:,3,3))
  end subroutine total_stress
  !
  ! ====================================================================================================================
  ! ====================================================================================================================
  !
  ! ====================================================================================================================
  ! ============================================== Parallelization =====================================================
  ! ====================================================================================================================
  !
  subroutine final_sum()
    tmpf(1:mxatms) = atomic_energy_ang
    call gdsum(tmpf(1),mxatms,tmpf(mxatms+1))
    atomic_energy_ang = tmpf(1:mxatms)
    !
    tmpf(1:mxatms) = atomic_energy_bnd
    call gdsum(tmpf(1),mxatms,tmpf(mxatms+1))
    atomic_energy_bnd = tmpf(1:mxatms)
    !
    tmpf(1:mxatms) = atomic_energy_ew1p
    call gdsum(tmpf(1),mxatms,tmpf(mxatms+1))
    atomic_energy_ew1p = tmpf(1:mxatms)
    !
    tmpf(1:mxatms) = atomic_energy_ew2p
    call gdsum(tmpf(1),mxatms,tmpf(mxatms+1))
    atomic_energy_ew2p = tmpf(1:mxatms)
    !
    tmpf(1:mxatms) = atomic_energy_ew3p
    call gdsum(tmpf(1),mxatms,tmpf(mxatms+1))
    atomic_energy_ew3p = tmpf(1:mxatms)
    !
    tmpf(1:mxatms) = atomic_energy_mbpol
    call gdsum(tmpf(1),mxatms,tmpf(mxatms+1))
    atomic_energy_mbpol = tmpf(1:mxatms)
    !
    tmpf(1:mxatms) = atomic_energy_srf
    call gdsum(tmpf(1),mxatms,tmpf(mxatms+1))
    atomic_energy_srf = tmpf(1:mxatms)
    !
    tmpf(1:mxatms) = atomic_energy_polar
    call gdsum(tmpf(1),mxatms,tmpf(mxatms+1))
    atomic_energy_polar = tmpf(1:mxatms)
    !
    tmpf4 = atomic_stress_ang
    tmpf2(1:9*mxatms) => tmpf4
    call gdsum(tmpf2(1),9*mxatms,tmpf3(1))
    tmpf5(1:mxatms,1:3,1:3) => tmpf2
    atomic_stress_ang = tmpf5
    !
    tmpf4 = atomic_stress_bnd
    tmpf2(1:9*mxatms) => tmpf4
    call gdsum(tmpf2(1),9*mxatms,tmpf3(1))
    tmpf5(1:mxatms,1:3,1:3) => tmpf2
    atomic_stress_bnd = tmpf5
    !
    tmpf4 = atomic_stress_ew1p
    tmpf2(1:9*mxatms) => tmpf4
    call gdsum(tmpf2(1),9*mxatms,tmpf3(1))
    tmpf5(1:mxatms,1:3,1:3) => tmpf2
    atomic_stress_ew1p = tmpf5
    !
    tmpf4 = atomic_stress_ew2p
    tmpf2(1:9*mxatms) => tmpf4
    call gdsum(tmpf2(1),9*mxatms,tmpf3(1))
    tmpf5(1:mxatms,1:3,1:3) => tmpf2
    atomic_stress_ew2p = tmpf5
    !
    tmpf4 = atomic_stress_ew3p
    tmpf2(1:9*mxatms) => tmpf4
    call gdsum(tmpf2(1),9*mxatms,tmpf3(1))
    tmpf5(1:mxatms,1:3,1:3) => tmpf2
    atomic_stress_ew3p = tmpf5
    !
    tmpf4 = atomic_stress_qd
    tmpf2(1:9*mxatms) => tmpf4
    call gdsum(tmpf2(1),9*mxatms,tmpf3(1))
    tmpf5(1:mxatms,1:3,1:3) => tmpf2
    atomic_stress_qd = tmpf5
    !
    tmpf4 = atomic_stress_mbpol
    tmpf2(1:9*mxatms) => tmpf4
    call gdsum(tmpf2(1),9*mxatms,tmpf3(1))
    tmpf5(1:mxatms,1:3,1:3) => tmpf2
    atomic_stress_mbpol = tmpf5
    !
    tmpf4 = atomic_stress_srf
    tmpf2(1:9*mxatms) => tmpf4
    call gdsum(tmpf2(1),9*mxatms,tmpf3(1))
    tmpf5(1:mxatms,1:3,1:3) => tmpf2
    atomic_stress_srf = tmpf5
    !
    tmpf4 = atomic_stress
    tmpf2(1:9*mxatms) => tmpf4
    call gdsum(tmpf2(1),9*mxatms,tmpf3(1))
    tmpf5(1:mxatms,1:3,1:3) => tmpf2
    atomic_stress = tmpf5
    !
    tmpFM1 = force_matrix
    tmpFM2(1:3*mxatms**2) => tmpFM1
    call gdsum(tmpFM2(1),3*mxatms**2,tmpFM3(1))
    tmpFM4(1:mxatms,1:mxatms,1:3) => tmpFM2
    force_matrix = tmpFM4
    !
    tmpFEW1 = force_ew1p
    tmpFEW2(1:3*mxatms) => tmpFEW1
    call gdsum(tmpFEW2(1),3*mxatms,tmpFEW3(1))
    tmpFEW4(1:mxatms,1:3) => tmpFEW2
    force_ew1p = tmpFEW4
    !
    tmpf(1:mxatms) = atomic_potential_energy
    call gdsum(tmpf(1),mxatms,tmpf(mxatms+1))
    atomic_potential_energy = tmpf(1:mxatms)

    ! tmpf(1:mxatms) = atomic_potential_energy2
    ! call gdsum(tmpf(1),mxatms,tmpf(mxatms+1))
    ! atomic_potential_energy2 = tmpf(1:mxatms)

    tmpf(1:mxatms) = delta_potential_energy
    call gdsum(tmpf(1),mxatms,tmpf(mxatms+1))
    delta_potential_energy = tmpf(1:mxatms)

    ! tmpf6(1:3) = j1
    ! call gdsum(tmpf6(1),3,tmpf6(4))
    ! j1 = tmpf6(1:3)
    ! !
    ! tmpf6(1:3) = j2
    ! call gdsum(tmpf6(1),3,tmpf6(4))
    ! j2 = tmpf6(1:3)
    ! !
    !
    tmpf(1:mxatms) = atomic_kinetic_energy
    call gdsum(tmpf(1),mxatms,tmpf(mxatms+1))
    atomic_kinetic_energy = tmpf(1:mxatms)
    !
    ! tmpf6(1:3) = jtot
    ! call gdsum(tmpf6(1),3,tmpf6(4))
    ! jtot = tmpf6(1:3)
    !deallocate(tmpf2,tmpf4)
    !deallocate(tmpf, tmpf3, tmpf6)

  end subroutine final_sum
  !
  subroutine final_sum_epot()
    tmpf(1:mxatms) = atomic_potential_energy2
    call gdsum(tmpf(1),mxatms,tmpf(mxatms+1))
    atomic_potential_energy2 = tmpf(1:mxatms)
  end subroutine final_sum_epot
  !
  subroutine final_sum_heat_flux()
    real(8) :: tmp(6)
    !
    tmp(1:3) = j1
    call gdsum(tmp(1),3,tmp(4))
    j1 = tmp(1:3)

    tmp(1:3) = j2
    call gdsum(tmp(1),3,tmp(4))
    j2 = tmp(1:3)

    tmp(1:3) = jtot
    call gdsum(tmp(1),3,tmp(4))
    jtot = tmp(1:3)

    tmp(1:3) = jkin
    call gdsum(tmp(1),3,tmp(4))
    jkin = tmp(1:3)

    tmp(1:3) = jpot
    call gdsum(tmp(1),3,tmp(4))
    jpot = tmp(1:3)

    tmp(1:3) = jtot2
    call gdsum(tmp(1),3,tmp(4))
    jtot2 = tmp(1:3)

  end subroutine final_sum_heat_flux
  !
  ! ====================================================================================================================
  ! ====================================================================================================================
  !
  ! ====================================================================================================================
  ! =============================================== Deallocation =======================================================
  ! ====================================================================================================================
  !
  subroutine deallocate_heat()
    deallocate(&
    atomic_stress, &
    atomic_stress_ang, &
    atomic_stress_bnd, &
    atomic_stress_ew1p, &
    atomic_stress_ew2p, &
    atomic_stress_ew3p, &
    atomic_stress_qd, &
    atomic_stress_mbpol, &
    atomic_stress_srf, &
    atomic_stress_lrcorrect, &
    atomic_energy_ang, &
    atomic_energy_bnd, &
    atomic_energy_ew1p, &
    atomic_energy_ew2p, &
    atomic_energy_ew3p, &
    atomic_energy_mbpol, &
    atomic_energy_srf, &
    atomic_energy_lrcorrect, &
    atomic_energy_polar, &
    atomic_kinetic_energy, &
    atomic_potential_energy,&
    force_matrix)
  end subroutine deallocate_heat
  !
  ! ====================================================================================================================
  ! ====================================================================================================================
  ! ====================================================================================================================
  !
end module heatcurrent
