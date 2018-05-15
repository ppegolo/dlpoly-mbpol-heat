! *****************************************************************************************
!   Module to accumulate the stress per atom, the energy per atom
!   and to calculate the macroscopic heat current.
!
!   23-03-2018: PP_: first version of update_stress and update_energy.
! *****************************************************************************************
module heatcurrent
  !
  use global_variables, only: mxatms    ! number of atoms
  use multibead, only: comm_bead, bead_rank

  implicit none
  !
  real(8), dimension(1:3)             :: j1   ! atomic energy contribution to the heat current
  real(8), dimension(1:3)             :: j2   ! stress contribution to the heat current
  real(8), dimension(1:3)             :: jtot, jtot0 = 0.d0   ! stress contribution to the heat current
  real(8), parameter :: junit=1.0d3  ! conversion form internal units
  real(8), dimension(:,:,:), allocatable :: atomic_stress_ang,atomic_stress_bnd,&
  atomic_stress_ew1p, atomic_stress_ew2p, atomic_stress_ew3p, &
  atomic_stress_qd, atomic_stress_mbpol, atomic_stress_srf, atomic_stress_lrcorrect,atomic_stress  ! Array with stress on each atom
  real(8), dimension(:), allocatable         ::  atomic_energy_ang,atomic_energy_bnd,&
  atomic_energy_ew1p, atomic_energy_ew2p, atomic_energy_ew3p, &
  atomic_energy_mbpol, atomic_energy_srf, atomic_energy_lrcorrect, &
  atomic_energy_polar, atomic_kinetic_energy, atomic_potential_energy  ! Array with the energy assigned to each atom
  integer :: ierror
  ! For the Finsl_sum
  real(8), allocatable :: tmpf(:), tmpf3(:), tmpf6(:)
  real(8), pointer :: tmpf2(:)
  real(8), contiguous, pointer :: tmpf4(:,:,:)
  real(8), contiguous, pointer :: tmpf5(:,:,:)
  !
contains
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
    atomic_potential_energy(1:mxatms))
    atomic_stress=0
    atomic_stress_ang=0
    atomic_stress_bnd=0
    atomic_stress_ew1p=0
    atomic_stress_ew2p=0
    atomic_stress_ew3p=0
    atomic_stress_qd=0
    atomic_stress_mbpol=0
    atomic_stress_srf=0
    atomic_stress_lrcorrect=0
    atomic_energy_ang=0
    atomic_energy_bnd=0
    atomic_energy_ew1p=0
    atomic_energy_ew2p=0
    atomic_energy_ew3p=0
    atomic_energy_mbpol=0
    atomic_energy_srf=0
    atomic_energy_lrcorrect=0
    atomic_energy_polar=0
    j1=0
    j2=0
    !
    allocate(tmpf2(9*mxatms),tmpf4(mxatms,3,3))
    allocate(tmpf(2*mxatms), tmpf3(9*mxatms), tmpf6(6))
  end subroutine init_heat
  !
  subroutine zero_heat()
    atomic_stress=0
    atomic_stress_ang=0
    atomic_stress_bnd=0
    atomic_stress_ew1p=0
    atomic_stress_ew2p=0
    atomic_stress_ew3p=0
    atomic_stress_qd=0
    atomic_stress_mbpol=0
    atomic_stress_srf=0
    atomic_stress_lrcorrect=0
    atomic_energy_ang=0
    atomic_energy_bnd=0
    atomic_energy_ew1p=0
    atomic_energy_ew2p=0
    atomic_energy_ew3p=0
    atomic_energy_mbpol=0
    atomic_energy_srf=0
    atomic_energy_lrcorrect=0
    atomic_energy_polar=0
    atomic_kinetic_energy=0
    atomic_potential_energy=0
    j1=0
    j2=0
  end subroutine zero_heat
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

  ! PP_: calculates the atomic kinetic energy: TO BE USED IN THE HEAT FLUX CALCULATION
  subroutine kinetic_energy_per_atom(vxx,vyy,vzz,weight)
    real(8), dimension(1:mxatms) :: vxx,vyy,vzz,weight
    real(8) :: vtmp(3)
    integer :: iatm
    do iatm=1,mxatms
      vtmp = (/ vxx(iatm), vyy(iatm), vzz(iatm) /)
      atomic_kinetic_energy(iatm) = 0.5d0*weight(iatm)*dot_product(vtmp,vtmp)
    end do
  end subroutine kinetic_energy_per_atom

  subroutine update_kinetic_energy(iatm,kinetic_energy)
    real(8) :: kinetic_energy
    integer :: iatm
    atomic_kinetic_energy(iatm) = atomic_kinetic_energy(iatm) + kinetic_energy
  end subroutine update_kinetic_energy
  !
  subroutine total_kinetic_energy(entot)
    real(8) :: entot
    entot = sum(atomic_kinetic_energy)
  end subroutine total_kinetic_energy
  ! PP_: calculate the total atomic potential energy: TO BE USED IN THE HEAT FLUX CALCULATION
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
    jtot = j1 + j2
  end subroutine compute_heat_flux
  !
  subroutine write_heat_flux(time)
    real(8), intent(in) :: time
    integer :: flux_file = 7000
    open(flux_file, file='heatflux.dat',position='append')
    write(flux_file,'(es12.5)') time
    write(flux_file,'(3es15.7)') j1*junit
    write(flux_file,'(3es15.7)') j2*junit
    write(flux_file,'(3es15.7)') jtot*junit
    close(flux_file)
  end subroutine write_heat_flux
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
    tmpf(1:mxatms) = atomic_potential_energy
    call gdsum(tmpf(1),mxatms,tmpf(mxatms+1))
    atomic_potential_energy = tmpf(1:mxatms)

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

  end subroutine final_sum_heat_flux
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
    atomic_potential_energy, &
    tmpf2, &
    tmpf4, &
    tmpf, &
    tmpf3, &
    tmpf6)
  end subroutine deallocate_heat
  !
  ! PPnote_: deallocate arrays
end module heatcurrent
