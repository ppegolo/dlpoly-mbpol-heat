#include "assert.h"

module centroid_dcp

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

implicit none

public :: centroid_dcp_init
public :: centroid_dcp_fini

public :: centroid_dcp_1st
public :: centroid_dcp_2nd

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

real(8), save, private :: dipmass
real(8), save, private :: dipsigma
real(8), save, private :: diptau

! diptau > diptau_max means "no thermostat"
real(8), public, parameter :: diptau_max = 1000.d0

real(8), private, parameter :: boltz = 8.31451115d-1
real(8), private, parameter :: r4pie0 = 138935.4835d0

integer, save, public :: ndipole
integer, save, private :: natom, iatm1, iatm2

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine centroid_dcp_init(n_atoms, dip_mass, dip_temp, dip_tau, polr2)

   use multibead

   implicit none

   integer, intent(in) :: n_atoms

   real(8), intent(in) :: dip_mass, dip_temp, dip_tau, polr2(*)

   integer :: i

   assert(n_atoms.gt.0)
   assert(dip_mass.gt.0.d0)
   assert(dip_temp.gt.0.d0)
   assert(dip_tau.gt.0.d0)

   natom = n_atoms

   iatm1 = (bead_rank*natom)/bead_size+1
   iatm2 = ((bead_rank+1)*natom)/bead_size

   ndipole = 0
   do i = 1, natom
      if (polr2(i).gt.1.d-6) &
         ndipole = ndipole + 1
   end do

   assert(ndipole.gt.0)

   dipmass = dip_mass
   dipsigma = dip_temp*boltz*dble(ndipole)*1.5d0

   diptau  = dip_tau

end subroutine centroid_dcp_init

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine centroid_dcp_fini()

   implicit none

   natom = 0
   iatm1 = 0
   iatm2 = 0

end subroutine centroid_dcp_fini

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine centroid_dcp_1st &
   (polr2, dipx, dipy, dipz, vdxx, vdyy, vdzz, &
    fdxx, fdyy, fdzz, tstep)

   implicit none

   real(8), intent(in) :: polr2(*), tstep

   real(8), intent(inout) :: dipx(*), dipy(*), dipz(*)
   real(8), intent(inout) :: vdxx(*), vdyy(*), vdzz(*)

   real(8), intent(in) :: fdxx(*), fdyy(*), fdzz(*)

   real(8) :: dpms1
   integer :: i

   do i = iatm1, iatm2

      if (.not.polr2(i).gt.1.d-6) &
         cycle

      dpms1 = 0.5d0*tstep*polr2(i)/r4pie0/dipmass

!     update dipole velocities
      vdxx(i) = vdxx(i) + fdxx(i)*dpms1
      vdyy(i) = vdyy(i) + fdyy(i)*dpms1
      vdzz(i) = vdzz(i) + fdzz(i)*dpms1

!     advance dipoles using velocity-verlet
      dipx(i) = dipx(i) + tstep*vdxx(i)
      dipy(i) = dipy(i) + tstep*vdyy(i)
      dipz(i) = dipz(i) + tstep*vdzz(i)

   end do

end subroutine centroid_dcp_1st

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine centroid_dcp_2nd &
   (polr2, dipx, dipy, dipz, vdxx, vdyy, vdzz, &
    fdxx, fdyy, fdzz, emux, emuy, emuz, &
    efieldkx, efieldky, efieldkz, efdcrecx, efdcrecy, efdcrecz, &
    efddmurecx, efddmurecy, efddmurecz, engdke, tstep)

   use multibead, only : bead_size

   implicit none

   real(8), intent(in) :: polr2(*)

   real(8), intent(inout) :: dipx(*), dipy(*), dipz(*)
   real(8), intent(inout) :: vdxx(*), vdyy(*), vdzz(*)
   real(8), intent(inout) :: fdxx(*), fdyy(*), fdzz(*)

   real(8), intent(in) :: emux(*), emuy(*), emuz(*)
   real(8), intent(in) :: efieldkx(*), efieldky(*), efieldkz(*)
   real(8), intent(in) :: efdcrecx(*), efdcrecy(*), efdcrecz(*)
   real(8), intent(in) :: efddmurecx(*), efddmurecy(*), efddmurecz(*)

   real(8), intent(out) :: engdke
   real(8), intent(in) :: tstep

   real(8) :: dpms1, chi
   integer :: i

   engdke = 0.d0

!  calculate forces applied on dipoles
   do i = iatm1, iatm2

      if (.not.polr2(i).gt.1.d-6) &
         cycle

      fdxx(i) = - dipx(i)/polr2(i)*r4pie0 + efieldkx(i) + efdcrecx(i) &
                 + emux(i) + efddmurecx(i)
      fdyy(i) = - dipy(i)/polr2(i)*r4pie0 + efieldky(i) + efdcrecy(i) &
                 + emuy(i) + efddmurecy(i)
      fdzz(i) = - dipz(i)/polr2(i)*r4pie0 + efieldkz(i) + efdcrecz(i) &
                 + emuz(i) + efddmurecz(i)

      dpms1 = 0.5*tstep*polr2(i)/r4pie0/dipmass

!     update dipole velocities
      vdxx(i) = vdxx(i) + fdxx(i)*dpms1
      vdyy(i) = vdyy(i) + fdyy(i)*dpms1
      vdzz(i) = vdzz(i) + fdzz(i)*dpms1

!     kinetic energy at half timestep
      engdke = engdke + (vdxx(i)**2 + vdyy(i)**2 + vdzz(i)**2)/polr2(i)

   end do

   engdke = 0.5d0*engdke*dipmass*r4pie0

   if (bead_size.gt.1) call gdsum(engdke, 1, chi)

   if (diptau.le.diptau_max) then
      chi = sqrt(1.d0 + tstep*(dipsigma/engdke - 1.d0)/diptau)

!     scale the velocities
      if (chi.lt.1.d0) then
         do i = iatm1, iatm2

            vdxx(i) = chi*vdxx(i)
            vdyy(i) = chi*vdyy(i)
            vdzz(i) = chi*vdzz(i)

         end do
      end if ! chi.lt.1.d0

   end if ! diptau.le.diptau_max

end subroutine centroid_dcp_2nd

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module centroid_dcp
