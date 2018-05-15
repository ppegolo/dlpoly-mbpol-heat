c==============================================================================
c
      subroutine centroid_read_config
     x           (cfgname,nstep,t,levcfg,imcon,cell,
     x            mxatms,natms,atmnam,
     x            xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz,buffer)
c
c==============================================================================
      use multibead
      use centroid_module, only: file_initial_config_pimd

c     !..................................................

      implicit none

      integer :: levcfg, nstep, imcon, mxatms, natms
      character(80) :: cfgname
      character(8) :: atmnam(mxatms)
      real(8) :: t
      real(8) :: cell(9)
      real(8) :: xxx(*), yyy(*), zzz(*)
      real(8) :: vxx(*), vyy(*), vzz(*)
      real(8) :: fxx(*), fyy(*), fzz(*)
      real(8) :: buffer(*)

#include "mpif.h"

      integer :: iatom, j, ierr, iatm

c     !..................................................

      if (bead_rank.ne.0) goto 100

      open(file_initial_config_pimd, file='FINAL_CONFIG'//bead_suffix)

c     !..................................................
      read(file_initial_config_pimd,"(a80)") cfgname
      read(file_initial_config_pimd,*)
     x     levcfg, imcon, nstep, t

      if ( imcon>0 ) then
         read(file_initial_config_pimd,*) cell(1:3)
         read(file_initial_config_pimd,*) cell(4:6)
         read(file_initial_config_pimd,*) cell(7:9)
      end if

      j = 0
      do iatom = 1, natms
         read(file_initial_config_pimd,*) atmnam(iatom), iatm
         read(file_initial_config_pimd,*)
     x             xxx(iatom), yyy(iatom), zzz(iatom)
         buffer(j + 1) = xxx(iatom)
         buffer(j + 2) = yyy(iatom)
         buffer(j + 3) = zzz(iatom)
         if(levcfg > 0) then
             read(file_initial_config_pimd,*)
     x             vxx(iatom), vyy(iatom), vzz(iatom)
             buffer(j + 4) = vxx(iatom)
             buffer(j + 5) = vyy(iatom)
             buffer(j + 6) = vzz(iatom)
         endif
         if(levcfg > 1) then
             read(file_initial_config_pimd,*)
     x             fxx(iatom), fyy(iatom), fzz(iatom)
             buffer(j + 7) = fxx(iatom)
             buffer(j + 8) = fyy(iatom)
             buffer(j + 9) = fzz(iatom)
         endif
         j = j + 9
      enddo

 100  continue

      call MPI_BCAST(t,1,MPI_DOUBLE_PRECISION,
     x               0, comm_mb, ierr)
      call MPI_BCAST(nstep,1,MPI_INTEGER,
     x               0, comm_mb, ierr)
      call MPI_BCAST(buffer,9*natms,MPI_DOUBLE_PRECISION,
     x               0, comm_bead, ierr)

      j = 0
      do iatom = 1, natms
         xxx(iatom) = buffer(j+1)
         yyy(iatom) = buffer(j+2)
         zzz(iatom) = buffer(j+3)

         vxx(iatom) = buffer(j+4)
         vyy(iatom) = buffer(j+5)
         vzz(iatom) = buffer(j+6)

         fxx(iatom) = buffer(j+7)
         fyy(iatom) = buffer(j+8)
         fzz(iatom) = buffer(j+9)
         j = j + 9
      enddo

      end subroutine
