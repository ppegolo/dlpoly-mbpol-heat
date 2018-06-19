c==============================================================================
c
      subroutine centroid_traject_bead 
     x        (nstep,keybin,t,mxatms,natms,atmnam,imcon,cell,
     x         xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)
c
c==============================================================================

      use multibead
      use centroid_module, only: file_history_bead

c     !..................................................

      implicit none

      integer, parameter :: keytrj = 1

      integer :: nstep, mxatms, natms,keybin

      character(8) :: atmnam(mxatms)
      character(*), parameter :: cfgname = ' >>> BEAD >>>'

      real(8) :: t
      real(8) :: cell(9)
      integer :: imcon
      real(8) :: xxx(*), yyy(*), zzz(*)
      real(8) :: vxx(*), vyy(*), vzz(*)
      real(8) :: fxx(*), fyy(*), fzz(*)

c     !..................................................

      logical, save :: newjob = .true.
      integer :: i

c     !..................................................

      if(keybin==0) then 

         open(file_history_bead,file='HISTORY_BEAD'//bead_suffix,
     x     position='APPEND')

c     if (newjob) then
c        write(file_history_bead,'(80a1)') cfgname
c        write(file_history_bead,'(3i10)') keytrj,imcon,natms
c        newjob = .false.
c     end if ! newjob

         write(file_history_bead,'(a8,4i10,f12.6)') 'timestep',
     x       nstep,natms,keytrj,imcon,t


         if(imcon.gt.0) write(file_history_bead,'(3g12.4)') cell

         do i = 1,natms
            write(file_history_bead,'(a8,i10,2f12.6)') atmnam(i),i,
     x         0.d0,0.d0
            write(file_history_bead,'(1p,3e12.4)') xxx(i),yyy(i),zzz(i)
            write(file_history_bead,'(1p,3e12.4)') vxx(i),vyy(i),vzz(i)
         enddo

       elseif (keybin==1) then 

         open(file_history_bead,file='HISTORY_BEAD'//bead_suffix,
     x     position='APPEND',form='unformatted' )

c     if (newjob) then
c        write(file_history_bead,'(80a1)') cfgname
c        write(file_history_bead,'(3i10)') keytrj,imcon,natms
c        newjob = .false.
c     end if ! newjob

         write(file_history_bead) 'timestep',
     x       nstep,natms,keytrj,imcon,t


         if(imcon.gt.0) write(file_history_bead) cell

         do i = 1,natms
            write(file_history_bead) atmnam(i),i,
     x         0.d0,0.d0
            write(file_history_bead) xxx(i),yyy(i),zzz(i)
            write(file_history_bead) vxx(i),vyy(i),vzz(i)
         enddo

      endif 

      close (file_history_bead)

      end subroutine
