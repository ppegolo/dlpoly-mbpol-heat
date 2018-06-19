#include "assert.h"

module centroid_make_cavity

!===============================================================================

implicit none

!===============================================================================

!
! parameters
!

real(8), public, parameter :: a_cavity = 627.509d0
real(8), public            :: b_cavity = 2.d0
real(8), public, parameter :: c_cavity = 0.05d0

real(8), public, parameter :: k_sphere = 224.237692118d0 ! kcal/mol/(\AA^2)
! k_sphere = 0.1 h/bohr^2

real(8), public, parameter :: hole_x = 0.0d0
real(8), public, parameter :: hole_y = 0.0d0
real(8), public, parameter :: hole_z = 0.0d0

real(8), private, parameter :: engunit = 418.4d0

real(8), dimension(9), private :: my_cell
integer, private :: my_imcon, my_natom

!
! subroutines
!

public :: make_cavity_init, make_cavity_frc
private :: vsphere, gsphere, vcavity, gcavity

! private

integer, parameter, private :: output_unit = 11232
integer, parameter, private :: trajectory_unit = 33113

!===============================================================================

contains

!===============================================================================

subroutine make_cavity_init(natom)

   use multibead

   implicit none

#include "mpif.h"

   integer, intent(in)    :: natom

   integer, parameter :: param_file = 33114
   integer :: ierr

   real(8) :: rtmp

   ! load INP_CRD

   if (is_bead_head()) then
      open (param_file, file = 'INP_CAVITY')
         read (unit = param_file, fmt = *, iostat = ierr)  rtmp ! start
         if (ierr.eq.0) b_cavity = rtmp
      close (param_file)

      ! allocate & load
      if (natom.lt.1) then
         print '(/a/)', ' ** Error ** : zero atoms in INP_CRD for 2D-IR run'
         call mb_abort()
      end if

      open (output_unit, file = 'OUTPUT.CAVITY')
      write (output_unit, *) ' a_cavity = ', a_cavity
      write (output_unit, *) ' b_cavity = ', b_cavity
      write (output_unit, *) ' c_cavity = ', c_cavity
      write (output_unit, *) ' hole_x   = ', hole_x
      write (output_unit, *) ' hole_y   = ', hole_y
      write (output_unit, *) ' hole_z   = ', hole_z
      close (output_unit)

   end if

   call MPI_BCAST(b_cavity, 1, &
                  MPI_DOUBLE_PRECISION, 0, comm_bead, ierr)

end subroutine make_cavity_init

!------------------------------------------------------------------------------!

subroutine make_cavity_frc(natom, nttm2, listttm2, cell, &
                           xxx, yyy, zzz, fxx, fyy, fzz, engacc, viracc)

   use multibead

   implicit none

#include "mpif.h"

   integer, intent(in) :: natom
   integer, intent(in) :: nttm2, listttm2(:)

   real(8), intent(in) :: cell(*), xxx(*), yyy(*), zzz(*)
   real(8), intent(inout) :: fxx(*), fyy(*), fzz(*)
   real(8), intent(inout) :: engacc, viracc

   integer :: itmp, ierr, iatom, nmol, io, t, j
   real(8) :: rtmp, dx(3), rcell(9), det
   real(8) :: ssx, ssy, ssz, xss, yss, zss
   character(8) :: ctmp

   real(8) :: rOOsq, rOO, crd(6), hole_cut2, grad(3)
   real(8) :: hole_cut, e2, crd_com(3), sphere_force

   ! image the coordinates

   nmol = nttm2/3

   hole_cut = b_cavity + 1.0d0

   hole_cut2 = hole_cut*hole_cut
   e2 = 0.d0
   crd_com(:) = 0.0d0

   call invert(cell,rcell,det)
   do t = 1, nmol
      io = listttm2(3*(t-1)+1) !oxygen index

      ! image: water relative to the cavity, get the O-hole distance

      dx(1) = xxx(io) - hole_x
      dx(2) = yyy(io) - hole_y
      dx(3) = zzz(io) - hole_z

      ssx=(rcell(1)*dx(1)+rcell(4)*dx(2)+rcell(7)*dx(3))
      ssy=(rcell(2)*dx(1)+rcell(5)*dx(2)+rcell(8)*dx(3))
      ssz=(rcell(3)*dx(1)+rcell(6)*dx(2)+rcell(9)*dx(3))

      xss=ssx-nint(ssx)
      yss=ssy-nint(ssy)
      zss=ssz-nint(ssz)

      dx(1)=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
      dx(2)=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
      dx(3)=(cell(3)*xss+cell(6)*yss+cell(9)*zss)

      !imaging done, now calc distance 

      rOOsq = dx(1)**2 + dx(2)**2 + dx(3)**2
      rOO = sqrt(rOOsq)

      crd(1) = hole_x
      crd(2) = hole_y
      crd(3) = hole_z

      !O
      crd(4) = crd(1) + dx(1)
      crd(5) = crd(2) + dx(2)
      crd(6) = crd(3) + dx(3)

      do j=1,3
          crd_com(j) = crd_com(j) + crd(3+j)
      end do

      ! Now, add force to oxygen atoms that have strayed within the cavity
      ! The 'hard sphere' potential is a*[1 - tanh((r-b)/c)]

      if (rOOsq.lt.(hole_cut2)) then
          e2 = e2 + vcavity(rOO)
          call gcavity(dx,rOO,grad)
          fxx(io) = fxx(io) + grad(1)
          fyy(io) = fyy(io) + grad(2)
          fzz(io) = fzz(io) + grad(3)
      
!           write(87,'(i4,a,4f12.8)') io, ' xyz:', xxx(io),yyy(io),zzz(io),rOO
!           write(87,'(i4,a,3f12.8)') io, ' crd:', crd(4),crd(5),crd(6)
!           write(87,'(i4,a,3f12.8)') io, ' frc:', grad(1),grad(2),grad(3)
!           write(87,*)
      
          do j=1,3
              viracc = viracc + crd(3+j)*grad(j)
          enddo
      end if

   end do

   do j=1,3
       crd_com(j) = crd_com(j)/nmol
   enddo
   sphere_force = k_sphere

   e2 = e2 + vsphere(crd_com, k_sphere)
   call gsphere(crd_com, grad, k_sphere, nmol)

   ! Now distribute the center-of-mass constraint force

   do t = 1, nmol
      io = listttm2(3*(t-1)+1) !oxygen index

      ! image: water relative to the cavity, get the O-hole distance

      dx(1) = xxx(io) - hole_x
      dx(2) = yyy(io) - hole_y
      dx(3) = zzz(io) - hole_z

      ssx=(rcell(1)*dx(1)+rcell(4)*dx(2)+rcell(7)*dx(3))
      ssy=(rcell(2)*dx(1)+rcell(5)*dx(2)+rcell(8)*dx(3))
      ssz=(rcell(3)*dx(1)+rcell(6)*dx(2)+rcell(9)*dx(3))

      xss=ssx-nint(ssx)
      yss=ssy-nint(ssy)
      zss=ssz-nint(ssz)

      dx(1)=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
      dx(2)=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
      dx(3)=(cell(3)*xss+cell(6)*yss+cell(9)*zss)

      !imaging done, now calc distance 

      rOOsq = dx(1)**2 + dx(2)**2 + dx(3)**2
      rOO = sqrt(rOOsq)

      crd(1) = hole_x
      crd(2) = hole_y
      crd(3) = hole_z

      !O
      crd(4) = crd(1) + dx(1)
      crd(5) = crd(2) + dx(2)
      crd(6) = crd(3) + dx(3)

      fzz(io) = fzz(io) + grad(3)
      do j=1,3
          viracc = viracc + crd(3+j)*grad(j)
      end do

   end do

   engacc = engacc + e2*engunit

end subroutine make_cavity_frc

   real(8) function vsphere(r, k)
       implicit none
       real(8) :: r(3),k

       !vsphere = k*(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
       vsphere = k*(r(3)*r(3))

   end function vsphere

   subroutine gsphere(rij,g, k,nw)
       implicit none
       real(8), intent(in) :: rij(3), k
       integer, intent(in) :: nw
       real(8), intent(out) :: g(3)

!       g(:) = -2.0d0*k/real(nw)*rij(:)*engunit
       g(1) = 0.0d0
       g(2) = 0.0d0
       g(3) =-2.0d0*k/real(nw)*rij(3)*engunit

   end subroutine gsphere

   real(8) function vcavity(r)
       implicit none
       real(8) :: r

       vcavity = a_cavity*(1.0d0 - tanh((r - b_cavity)/c_cavity))

   end function vcavity

   subroutine gcavity(rij,r,g)
       implicit none
       real(8), intent(in) :: rij(3),r
       real(8), intent(out) :: g(3)
       real(8) :: tmp

       tmp = 1.0d0/cosh((r-b_cavity)/c_cavity) !sech(x)
       tmp = tmp*tmp !sech^2(x)
       tmp = a_cavity*tmp/c_cavity/r
       g(:) = tmp*rij(:)*engunit

   end subroutine gcavity

end module centroid_make_cavity

!===============================================================================
