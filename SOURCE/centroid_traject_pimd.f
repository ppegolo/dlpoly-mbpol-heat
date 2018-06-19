c==============================================================================
c
      subroutine centroid_traject_pimd
     x  (keytrj,nstep,keybin,t,mxatms,natms,atmnam,imcon,cell,
     x   xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz,lcavity)
c
c==============================================================================

      use multibead
      use centroid_make_cavity
      use centroid_module, only: file_history_pimd

c     !..................................................

      implicit none

      integer, intent(in) :: keytrj

      integer :: nstep, mxatms, natms, keybin
      logical :: lcavity

      character(8),target :: atmnam(mxatms)
      character(*), parameter :: cfgname = ' >>> PIMD >>>'

      real(8) :: t
      real(8) :: cell(9)
      integer :: imcon, natms_cavity
      real(8),target :: xxx(*), yyy(*), zzz(*)
      real(8) :: vxx(*), vyy(*), vzz(*)
      real(8) :: fxx(*), fyy(*), fzz(*)
      integer :: file_pos_pimd=33, file_vel_pimd=34
      integer :: file_force_pimd=35

      integer,allocatable :: tmp_atms(:)

      real(8),pointer :: tmp_xxx(:), tmp_yyy(:), tmp_zzz(:)
      character(8),pointer :: tmp_atmnam(:)

c     !..................................................

      logical, save :: newjob = .true.
      integer :: i

c     !..................................................

      allocate(tmp_atms(natms))

      if(keybin==0) then   ! binary or ascii file format
        ! PP_:begin
        open (file_pos_pimd,file='POSITION_PIMD',position='APPEND')

        if (keytrj >= 1) then
           open (file_vel_pimd,file='VELOCITY_PIMD',position='APPEND')
           write(file_vel_pimd,'(i10,1f15.7,i10)') nstep, t, natms
        endif

       if ( keytrj >= 2 ) then
          open(file_force_pimd,file='FORCE_PIMD',position='APPEND')
          write(file_force_pimd,'(i10,1f15.7,i10)') nstep, t, natms
       endif

        write(file_pos_pimd,'(i10,1f15.7,i10)') nstep, t, natms
        write(file_pos_pimd,'(3f20.10)') cell(1),cell(2),cell(3)
        write(file_pos_pimd,'(3f20.10)') cell(4),cell(5),cell(6)
        write(file_pos_pimd,'(3f20.10)') cell(7),cell(8),cell(9)     ! PP_:end


        open(file_history_pimd,file='HISTORY_PIMD'//bead_suffix,
     x     position='APPEND')

        if(lcavity) then
            natms_cavity = natms + 1
        else
            natms_cavity = natms
        end if


!FP_fix
! Not necessary and makes it difficult to combine the trajectories.
!      if (newjob) then
!         write(file_history_pimd,'(80a1)') cfgname
!         write(file_history_pimd,'(3i10)') keytrj,imcon,natms_cavity
!         newjob = .false.
!      end if ! newjob

        write(file_history_pimd,'(a8,4i10,f12.6)') 'timestep',
     x       nstep,natms_cavity,keytrj,imcon,t


        if(imcon.gt.0) write(file_history_pimd,'(3g12.4)') cell

        do i = 1,natms
          write(file_history_pimd,'(a8,i10,2f12.6)') atmnam(i),i,
     x       0.d0,0.d0
          write(file_history_pimd,'(1p,3e12.4)') xxx(i),yyy(i),zzz(i)

          ! PP_:begin
          write (file_pos_pimd,'(a8,3f15.7)')
     x           atmnam(i),xxx(i),yyy(i),zzz(i)

          if (keytrj >= 1) then
             write (file_vel_pimd,'(a8,3f15.7)')
     x              atmnam(i),vxx(i),vyy(i),vzz(i)
          endif

          if (keytrj >= 2) then
             write (file_force_pimd,'(a8,3f15.7)')
     x              atmnam(i),fxx(i),fyy(i),fzz(i)
          endif ! PP_:end
        enddo

        if(lcavity) then
          write(file_history_pimd,'(a8,i10,2f12.6)') "CAVITY",
     x       natms_cavity, 0.d0,0.d0
          write(file_history_pimd,'(1p,3e12.4)') hole_x, hole_y, hole_z
        end if

      elseif(keybin==1) then

        open(file_history_pimd,file='HISTORY_PIMD'//bead_suffix,
     x     position='APPEND',form='unformatted' )

        if(lcavity) then
            natms_cavity = natms + 1
        else
            natms_cavity = natms
        end if

!FP_fix
! Not necessary and makes it difficult to combine the trajectories.
!      if (newjob) then
!         write(file_history_pimd,'(80a1)') cfgname
!         write(file_history_pimd,'(3i10)') keytrj,imcon,natms_cavity
!         newjob = .false.
!      end if ! newjob

        tmp_xxx => xxx(1:natms)
        tmp_yyy => yyy(1:natms)
        tmp_zzz => zzz(1:natms)
        tmp_atmnam => atmnam(1:natms)
        tmp_atms=(/(i,i=1,natms)/)

        write(file_history_pimd) 'timestep',
     x       nstep,natms_cavity,keytrj,imcon,t


        if(imcon.gt.0) then
           write(file_history_pimd) cell(1:3)
           write(file_history_pimd) cell(4:6)
           write(file_history_pimd) cell(7:9)
        endif

          write(file_history_pimd) tmp_atmnam,tmp_atms,
     x       0.d0,0.d0
          write(file_history_pimd) tmp_xxx,tmp_yyy,tmp_zzz

!        do i = 1,natms
!          write(file_history_pimd) atmnam(i),i,
!     x       0.d0,0.d0
!          write(file_history_pimd) xxx(i),yyy(i),zzz(i)
!        enddo

        if(lcavity) then
          write(file_history_pimd) "CAVITY",
     x       natms_cavity, 0.d0,0.d0
          write(file_history_pimd) hole_x, hole_y, hole_z
        end if

      endif

      close (file_history_pimd)
      deallocate(tmp_atms)

      end subroutine
