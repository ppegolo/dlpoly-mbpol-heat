#include "assert.h"

module centroid_2D_IR

!===============================================================================

implicit none

!===============================================================================

!
! parameters
!

! number of points, r_min and r_max for the 1D potential
integer, parameter, public :: n_vib = 11
real(8), parameter, public :: r_vmin = 0.5d0
real(8), parameter, public :: r_vmax = 1.5d0

! initial/final indices of HOD
integer, public :: nihod = 1
integer, public :: nfhod = 1

!
! trajectory
!

real(8), public, save :: centroid_time = 0.d0
real(8), public, save :: start_time = 0.d0
real(8), public, save :: final_time = 40000.d0

real(8), public, save, allocatable :: centroid_pos(:) ! 3*natom
real(8), public, save, allocatable :: centroid_box(:) ! 9

!
! subroutines
!

public :: centroid_2D_IR_init
public :: centroid_2D_IR_fini

public :: centroid_2D_IR_read_frame
public :: centroid_2D_IR_set
public :: centroid_2D_IR_report

! private

integer, parameter, private :: output_unit = 11232
integer, parameter, private :: trajectory_unit = 33113

!===============================================================================

contains

!===============================================================================

subroutine centroid_2D_IR_init(natom, continue_2D_IR)

   use multibead

   implicit none

#  include "mpif.h"

   integer, intent(in)    :: natom
   logical, intent(inout) :: continue_2D_IR

   integer, parameter :: start_file = 33114
   integer :: ierr

   real(8) :: rtmp

   if (ring_size.ne.1) then
      if (is_bead_head()) &
         print '(/a/)', ' ** Error ** : more than one bead for 2D-IR run'
      call mb_abort()
   end if

   ! load INP_CRD

   if (is_bead_head()) then
      open (start_file, file = 'INP_START')
         read (unit = start_file, fmt = *, iostat = ierr)  rtmp ! start
         if (ierr.eq.0) start_time = rtmp
         read (unit = start_file, fmt = *, iostat = ierr)  rtmp ! finish
         if (ierr.eq.0) final_time = rtmp
      close (start_file)

write(*,*) "start_time = ", start_time
write(*,*) "final_time = ", final_time

      open (start_file, file = 'INP_IHOD')
         read (unit = start_file, fmt = *, iostat = ierr)  rtmp ! start
         if (ierr.eq.0) then 
            nihod = rtmp
            nfhod = rtmp
         endif   
      close (start_file)

      open (trajectory_unit, file = 'INP_CRD')
!      ! find out the number of frames
!      do while (.true.)
!         read (unit = trajectory_unit, fmt = *, iostat = ierr) &
!            itmp, rtmp(1), rtmp(2) ! step time box
!         if (ierr.ne.0) exit
!
!         do iatom = 1, 3 !dummy index
!            read (unit = trajectory_unit, fmt = *, iostat = ierr) &
!               rtmp(1), rtmp(2), rtmp(3) ! box
!            if (ierr.ne.0) exit
!         end do
!
!         do iatom = 1, natom
!            read (unit = trajectory_unit, fmt = *, iostat = ierr) &
!               ctmp, rtmp(1), rtmp(2), rtmp(3) ! atm x y z
!            if (ierr.ne.0) exit
!         end do
!
!         if (ierr.ne.0) exit
!         nframe = nframe + 1
!      end do ! while (.true.)

      ! allocate & load
      if (natom.lt.1) then
         print '(/a/)', ' ** Error ** : zero atoms in INP_CRD for 2D-IR run'
         call mb_abort()
      end if

      allocate(centroid_pos(3*natom), centroid_box(9))

      open (output_unit, file = 'OUTPUT.2D-IR')
      write (output_unit, *) '  n_vib = ', n_vib
      write (output_unit, *) ' r_vmin = ', r_vmin
      write (output_unit, *) ' r_vmax = ', r_vmax
      write (output_unit, *) '  nihod = ', nihod
      write (output_unit, *) '  nfhod = ', nfhod
      close (output_unit)

   end if

   call MPI_BCAST(nihod, 1, &
       MPI_INTEGER, 0, comm_bead, ierr)
   call MPI_BCAST(nfhod, 1, &
       MPI_INTEGER, 0, comm_bead, ierr)

   continue_2D_IR = .true.

   if (bead_rank.ne.0) then
      allocate(centroid_pos(3*natom), centroid_box(9))
   end if

end subroutine centroid_2D_IR_init

subroutine centroid_2D_IR_read_frame(natom, nttm2, listttm2, continue_2D_IR)

   use multibead

   implicit none

#  include "mpif.h"

   integer, intent(in) :: natom
   integer, intent(in) :: nttm2, listttm2(:)
   logical, intent(inout) :: continue_2D_IR

   integer :: itmp, ierr, iatom, nmol, imol, io3, ih3, t
   real(8) :: rtmp, dx(3), rcell(9), det
   real(8) :: ssx, ssy, ssz, xss, yss, zss
   character(8) :: ctmp

   if (is_bead_head()) then
      if(continue_2D_IR .eqv. .true.) then
91687 continue
         read (unit = trajectory_unit, fmt = *, iostat = ierr) &
            itmp, centroid_time, rtmp ! step time natm
         if (ierr.ne.0) then
             continue_2D_IR = .false.
             return
         end if
        if(centroid_time .gt. final_time) then
            call centroid_2D_IR_fini()
            stop
        endif

         read (unit = trajectory_unit, fmt = *, iostat = ierr) &
            centroid_box(1), centroid_box(2), centroid_box(3) 
         if (ierr.ne.0) then
             continue_2D_IR = .false.
             return
         end if

         read (unit = trajectory_unit, fmt = *, iostat = ierr) &
            centroid_box(4), centroid_box(5), centroid_box(6) 
         if (ierr.ne.0) then
             continue_2D_IR = .false.
             return
         end if

         read (unit = trajectory_unit, fmt = *, iostat = ierr) &
            centroid_box(7), centroid_box(8), centroid_box(9) 
         if (ierr.ne.0) then
             continue_2D_IR = .false.
             return
         end if

         do iatom = 1, natom
            read (unit = trajectory_unit, fmt = *, iostat = ierr) ctmp, &
               centroid_pos(3*iatom - 2), &
               centroid_pos(3*iatom - 1), &
               centroid_pos(3*iatom - 0)
            if (ierr.ne.0) then
                continue_2D_IR = .false.
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

end subroutine centroid_2D_IR_read_frame

!===============================================================================

subroutine centroid_2D_IR_fini()

   implicit none

   deallocate(centroid_pos, centroid_box)

end subroutine centroid_2D_IR_fini

!===============================================================================

subroutine centroid_2D_IR_set(cell,xxx, yyy, zzz, ihod, ivib)

   implicit none

   real(8), intent(inout) :: cell(*), xxx(*), yyy(*), zzz(*)
   integer, intent(in) :: ihod, ivib

   integer :: iatm, natom, io3, ih3, ih, i
   real(8) :: R0, th, csth, snth, snph, csph, x, y, z, Rvib

   assert(allocated(centroid_pos))

   !assert(ihod.ge.nihod.and.ihod.le.nfhod)
   assert(ivib.ge.1.and.ivib.le.n_vib)

   do i = 1,9
      cell(i) = centroid_box(i)
   end do

   natom = size(centroid_pos, 1)/3

   do iatm = 1, natom
      xxx(iatm) = centroid_pos(3*iatm - 2)
      yyy(iatm) = centroid_pos(3*iatm - 1)
      zzz(iatm) = centroid_pos(3*iatm - 0)
   end do

   
!   ! "deform" the molecule #ihod

   io3 = 3*(ihod - 1) + 1 ! Oxygen
   ih3 = io3 + 3

   ! O-H distance
   R0 = sum((centroid_pos(io3:io3+2) &
           - centroid_pos(ih3:ih3+2))**2)
   R0 = sqrt(R0)

   csth = (centroid_pos(ih3+2) - centroid_pos(io3+2))/R0
   th = acos(csth)
   snth = sin(th)
   snph = (centroid_pos(ih3+1) - centroid_pos(io3+1))/R0/snth
   csph = (centroid_pos(ih3) - centroid_pos(io3))/R0/snth

   Rvib = r_vmin + dble(ivib - 1)*((r_vmax - r_vmin)/(n_vib - 1))

   x = Rvib*snth*csph
   y = Rvib*snth*snph
   z = Rvib*csth

   ih = ihod + 1

   xxx(ih) = centroid_pos(io3    ) + x
   yyy(ih) = centroid_pos(io3 + 1) + y
   zzz(ih) = centroid_pos(io3 + 2) + z

!   write(81,*) 3
!   write(81,*) io3, ih3, ihod, ih
!   write(81,'(a8,3f15.8)') "OW", xxx(ihod), yyy(ihod), zzz(ihod)
!   write(81,'(a8,3f15.8)') "HW", xxx(ihod+1), yyy(ihod+1), zzz(ihod+1)
!   write(81,'(a8,3f15.8)') "HW", xxx(ihod+2), yyy(ihod+2), zzz(ihod+2)

end subroutine centroid_2D_IR_set

!===============================================================================

subroutine centroid_2D_IR_report &
   (ihod,ivib,Epot,xxx,yyy,zzz, &
    listttm2,nttm2,weight,dipx,dipy,dipz)

   use unit_parameters, only : eatd
   use dipole_moments, only : moldipx, moldipy, moldipz

   implicit none

   integer, intent(in) :: ihod, ivib, nttm2, listttm2(*)
   real(8), intent(in) :: Epot
   real(8), intent(in) :: xxx(*), yyy(*), zzz(*), weight(*)
   real(8), intent(in) :: dipx(*), dipy(*), dipz(*)

   integer :: iox, ih1, ih2, imm, iwat
   real(8) :: Rvib, com_pos(3), weight_sum, moldip
   real(8) :: indx, indy, indz, ind

   Rvib = r_vmin + dble(ivib - 1)*((r_vmax - r_vmin)/(n_vib - 1))

   open(10000+ihod,position='append')
   write(10000+ihod,'(3f18.5)') centroid_time, Rvib, Epot
   if (ivib.eq.n_vib) write(10000+ihod,'(2x)')
   close(10000+ihod)

   iox = ihod  
   ih1 = ihod+1
   ih2 = ihod+2

   iwat = (ihod - listttm2(1))/3 + 1

   open(20000+ihod,position='append')
   write(20000+ihod,'(11f10.5)') centroid_time, Rvib, &
      xxx(iox), yyy(iox), zzz(iox), &
      xxx(ih1), yyy(ih1), zzz(ih1), &
      xxx(ih2), yyy(ih2), zzz(ih2)
   if (ivib.eq.n_vib) write(20000+ihod,'(2x)')
   close(20000+ihod)

   weight_sum = weight(iox) + weight(ih1) + weight(ih2)

   com_pos(1) = (weight(iox)*xxx(iox) &
      + weight(ih1)*xxx(ih1) + weight(ih2)*xxx(ih2))/weight_sum
   com_pos(2) = (weight(iox)*yyy(iox) &
      + weight(ih1)*yyy(ih1) + weight(ih2)*yyy(ih2))/weight_sum
   com_pos(3) = (weight(iox)*zzz(iox) &
      + weight(ih1)*zzz(ih1) + weight(ih2)*zzz(ih2))/weight_sum

   !moldipx is nwat long
   moldip = sqrt(moldipx(iwat)*moldipx(iwat) &
               + moldipy(iwat)*moldipy(iwat) &
               + moldipz(iwat)*moldipz(iwat))*eatd

   open(30000+ihod,position='append')
   write(30000+ihod,'(9f10.5)') centroid_time, Rvib, com_pos, &
      moldipx(iwat)*eatd, moldipy(iwat)*eatd, moldipz(iwat)*eatd, moldip
   if (ivib.eq.n_vib) write(30000+ihod,'(2x)')
   close(30000+ihod)

   imm = listttm2(nttm2) + iwat

   !dipx is natm + nfict_atm long
   indx = (dipx(iox) + dipx(ih1) + dipx(ih2) + dipx(imm))*eatd
   indy = (dipy(iox) + dipy(ih1) + dipy(ih2) + dipy(imm))*eatd
   indz = (dipz(iox) + dipz(ih1) + dipz(ih2) + dipz(imm))*eatd

   ind = sqrt(indx*indx + indy*indy + indz*indz)

   open(40000+ihod,position='append')
   write(40000+ihod,'(9f10.5)') centroid_time, Rvib, com_pos, &
      indx, indy, indz, ind
   if (ivib.eq.n_vib) write(40000+ihod,'(2x)')
   close(40000+ihod)

end subroutine centroid_2D_IR_report

!===============================================================================

end module centroid_2D_IR

!===============================================================================
