c==============================================================================
c
      subroutine centroid_traject_cmd
     x           (keytrj,keybin,nstep,t,mxatms,natms,atmnam,cell,
     x            npx,npy,npz,nvx,nvy,nvz,nfx,nfy,nfz,xxx,yyy,zzz,
     x            vxx,vyy,vzz,fxx,fyy,fzz)
c
c==============================================================================

      use centroid_module, only:
     x                     file_pos_cmd, file_vel_cmd, file_force_cmd

c     !..................................................

      implicit none

      integer :: keytrj,keybin, nstep, mxatms, natms

      character(8) :: atmnam(mxatms)

      real(8) :: t
      real(8) :: cell(9)
      real(8) :: npx(*), npy(*), npz(*)
      real(8) :: nvx(*), nvy(*), nvz(*)
      real(8) :: nfx(*), nfy(*), nfz(*)
      real(8) :: xxx(*), yyy(*), zzz(*)
      real(8) :: vxx(*), vyy(*), vzz(*)
      real(8) :: fxx(*), fyy(*), fzz(*)

c     !..................................................

      integer :: iatom

c     !..................................................

      if(keybin==0) then
         open (file_pos_cmd,file='POSITION_CMD',position='APPEND')

         if (keytrj >= 1) then
            open (file_vel_cmd,file='VELOCITY_CMD',position='APPEND')
            write(file_vel_cmd,'(i10,1f15.7,i10)') nstep, t, natms
         endif

        if ( keytrj >= 2 ) then
           open(file_force_cmd,file='FORCE_CMD',position='APPEND')
           write(file_force_cmd,'(i10,1f15.7,i10)') nstep, t, natms
        endif

         write(file_pos_cmd,'(i10,1f15.7,i10)') nstep, t, natms
         write(file_pos_cmd,'(3f20.10)') cell(1),cell(2),cell(3)
         write(file_pos_cmd,'(3f20.10)') cell(4),cell(5),cell(6)
         write(file_pos_cmd,'(3f20.10)') cell(7),cell(8),cell(9)

         do iatom = 1, natms
           write (file_pos_cmd,'(a8,3f15.7)')
!     x           atmnam(iatom),npx(iatom),npy(iatom),npz(iatom)
     x           atmnam(iatom),xxx(iatom),yyy(iatom),zzz(iatom)

           if (keytrj >= 1) then
              write (file_vel_cmd,'(a8,3f15.7)')
!     x              atmnam(iatom),nvx(iatom), nvy(iatom),nvz(iatom)
     x              atmnam(iatom),vxx(iatom),vyy(iatom),vzz(iatom)

           endif

           if (keytrj >= 2) then
              write (file_force_cmd,'(a8,3f15.7)')
     x              atmnam(iatom),fxx(iatom),fyy(iatom),fzz(iatom)
           endif


c        if ( ltraj_cmd>=2 ) then
c           write(file_force_cmd,'(a8,3f22.7)')
c     x            atmnam(iatom),nfx(iatom),nfy(iatom),nfz(iatom)
c        endif

         enddo

       elseif (keybin==1) then

         open (file_pos_cmd,file='POSITION_CMD',position='APPEND',
     x     form='unformatted')

         if (keytrj .eq. 1) then
            open (file_vel_cmd,file='VELOCITY_CMD',position='APPEND',
     x     form='unformatted')
            write(file_vel_cmd) nstep, t, natms
         endif

c        if ( ltraj_cmd>=2 ) then
c           open(file_force_cmd,file='FORCE_CMD',position='APPEND')
c        endif

         write(file_pos_cmd) nstep, t, natms
         write(file_pos_cmd) cell(1),cell(2),cell(3)
         write(file_pos_cmd) cell(4),cell(5),cell(6)
         write(file_pos_cmd) cell(7),cell(8),cell(9)

!         do iatom = 1, natms
           write (file_pos_cmd)
     x           atmnam(1:natms),npx(1:natms),npy(1:natms),npz(1:natms)

           if (keytrj .eq. 1) then
              write (file_vel_cmd)
     x          atmnam(1:natms),nvx(1:natms),nvy(1:natms),nvz(1:natms)
          endif

c        if ( ltraj_cmd>=2 ) then
c           write(file_force_cmd,'(a8,3f15.7)')
c    x            atmnam(natms),nfx(natms),nfy(natms),nfz(natms)
c        endif

c         enddo
      endif   ! binary or ascill file

      if(keybin==0.or.keybin==1) then
        close (file_pos_cmd)

        if (keytrj >= 1) then
           close(file_vel_cmd)
        endif
        if (keytrj >= 2) then
           close(file_force_cmd)
        endif
      endif

c     if ( ltraj_cmd>=2 ) then
c        close( file_force_cmd )
c     endif

      end subroutine
