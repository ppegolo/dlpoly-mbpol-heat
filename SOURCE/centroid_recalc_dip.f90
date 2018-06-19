#include "assert.h"

module centroid_recalc_dip

!===============================================================================

implicit none

!===============================================================================

!
! trajectory
!

integer, public, save :: nframe = 0
real(8), public, save :: centroid_time = 0.d0
real(8), public, save :: start_time = 0.d0

real(8), public, save, allocatable :: centroid_pos(:) ! 3*natom
real(8), public, save, allocatable :: centroid_box(:) ! 9

!
! subroutines
!

public :: centroid_recalc_dip_init
public :: centroid_recalc_dip_fini

public :: centroid_recalc_dip_read_frame
public :: centroid_recalc_dip_set

! private

integer, parameter, private :: trajectory_unit = 33113

!===============================================================================

contains

!===============================================================================

subroutine centroid_recalc_dip_init(natom, continue_recalc_dip)

   use multibead

   implicit none

   integer, intent(in)    :: natom
   logical, intent(inout) :: continue_recalc_dip

   integer, parameter :: start_file = 33114
   integer :: ierr

   real(8) :: rtmp

   if (ring_size.ne.1) then
      if (is_bead_head()) &
         print '(/a/)', ' ** Error ** : more than one bead for recalc-dip run'
      call mb_abort()
   end if

   ! load INP_CRD

   if (is_bead_head()) then
      open (start_file, file = 'INP_START')
         read (unit = start_file, fmt = *, iostat = ierr)  rtmp ! start
         if (ierr.eq.0) start_time = rtmp

write(*,*) "start_time = ", start_time

      open (trajectory_unit, file = 'INP_CRD')

      ! allocate & load
      if (natom.lt.1) then
         print '(/a/)', ' ** Error ** : zero atoms in INP_CRD for recalc-dip run'
         call mb_abort()
      end if

      allocate(centroid_pos(3*natom), centroid_box(9))

   end if

   continue_recalc_dip = .true.

   if (bead_rank.ne.0) then
      allocate(centroid_pos(3*natom), centroid_box(9))
   end if

end subroutine centroid_recalc_dip_init

subroutine centroid_recalc_dip_read_frame(natom, nttm2, listttm2, continue_recalc_dip)

   use multibead

   implicit none

#  include "mpif.h"

   integer, intent(in) :: natom
   integer, intent(in) :: nttm2, listttm2(:)
   logical, intent(inout) :: continue_recalc_dip

   integer :: itmp, ierr, iatom, nmol, imol, io3, ih3, t
   real(8) :: rtmp, dx(3), rcell(9), det
   real(8) :: ssx, ssy, ssz, xss, yss, zss
   character(8) :: ctmp

   if (is_bead_head()) then
      if(continue_recalc_dip) then
91687 continue
         read (unit = trajectory_unit, fmt = *, iostat = ierr) &
            itmp, centroid_time, rtmp ! step time natm
         if (ierr.ne.0) then
             continue_recalc_dip = .false.
                call MPI_BCAST(continue_recalc_dip, 1, &
                    MPI_LOGICAL, 0, comm_bead, ierr)
             return
         end if

         read (unit = trajectory_unit, fmt = *, iostat = ierr) &
            centroid_box(1), centroid_box(2), centroid_box(3) 
         if (ierr.ne.0) then
             continue_recalc_dip = .false.
                call MPI_BCAST(continue_recalc_dip, 1, &
                    MPI_LOGICAL, 0, comm_bead, ierr)
             return
         end if

         read (unit = trajectory_unit, fmt = *, iostat = ierr) &
            centroid_box(4), centroid_box(5), centroid_box(6) 
         if (ierr.ne.0) then
             continue_recalc_dip = .false.
                call MPI_BCAST(continue_recalc_dip, 1, &
                    MPI_LOGICAL, 0, comm_bead, ierr)
             return
         end if

         read (unit = trajectory_unit, fmt = *, iostat = ierr) &
            centroid_box(7), centroid_box(8), centroid_box(9) 
         if (ierr.ne.0) then
             continue_recalc_dip = .false.
                call MPI_BCAST(continue_recalc_dip, 1, &
                    MPI_LOGICAL, 0, comm_bead, ierr)
             return
         end if

         do iatom = 1, natom
            read (unit = trajectory_unit, fmt = *, iostat = ierr) ctmp, &
               centroid_pos(3*iatom - 2), &
               centroid_pos(3*iatom - 1), &
               centroid_pos(3*iatom - 0)
            if (ierr.ne.0) then
                continue_recalc_dip = .false.
                call MPI_BCAST(continue_recalc_dip, 1, &
                    MPI_LOGICAL, 0, comm_bead, ierr)
                return
            end if
         end do
        if(centroid_time .lt. start_time) goto 91687
      end if
   end if

   ! broadcast to all PEs

   call MPI_BCAST(centroid_pos, 3*natom, &
      MPI_DOUBLE_PRECISION, 0, comm_bead, ierr)

   call MPI_BCAST(centroid_box, 9, &
      MPI_DOUBLE_PRECISION, 0, comm_bead, ierr)

   ! image the coordinates

   nmol = nttm2/3

   call invert(centroid_box,rcell,det)
   do t = 1, nmol
      imol = listttm2(3*(t-1)+1) !oxygen index
      io3 = 3*(imol-1)+1
      do iatom = 2, 3
         ih3 = io3 + 3*(iatom - 1)! + 1
         dx(:) = centroid_pos(ih3:ih3+2) - centroid_pos(io3:io3+2)

         ssx=(rcell(1)*dx(1)+rcell(4)*dx(2)+rcell(7)*dx(3))
         ssy=(rcell(2)*dx(1)+rcell(5)*dx(2)+rcell(8)*dx(3))
         ssz=(rcell(3)*dx(1)+rcell(6)*dx(2)+rcell(9)*dx(3))

         xss=ssx-nint(ssx)
         yss=ssy-nint(ssy)
         zss=ssz-nint(ssz)

         dx(1)=(centroid_box(1)*xss+centroid_box(4)*yss+centroid_box(7)*zss)
         dx(2)=(centroid_box(2)*xss+centroid_box(5)*yss+centroid_box(8)*zss)
         dx(3)=(centroid_box(3)*xss+centroid_box(6)*yss+centroid_box(9)*zss)

         centroid_pos(ih3:ih3+2) = centroid_pos(io3:io3+2) + dx(:)
      end do
   end do

end subroutine centroid_recalc_dip_read_frame

!===============================================================================

subroutine centroid_recalc_dip_fini()

   implicit none

   deallocate(centroid_pos, centroid_box)

end subroutine centroid_recalc_dip_fini

!===============================================================================

subroutine centroid_recalc_dip_set(cell,xxx, yyy, zzz)

   implicit none

   real(8), intent(inout) :: cell(*), xxx(*), yyy(*), zzz(*)

   integer :: iatm, natom, io3, ih3, ih, i
   real(8) :: R0, th, csth, snth, snph, csph, x, y, z, Rvib

   assert(allocated(centroid_pos))

   do i = 1,9
      cell(i) = centroid_box(i)
   end do

   natom = size(centroid_pos, 1)/3

   do iatm = 1, natom
      xxx(iatm) = centroid_pos(3*iatm - 2)
      yyy(iatm) = centroid_pos(3*iatm - 1)
      zzz(iatm) = centroid_pos(3*iatm - 0)
   end do

end subroutine centroid_recalc_dip_set

!===============================================================================

end module centroid_recalc_dip

!===============================================================================
