c==============================================================================
c
      subroutine centroid_read_dipole 
     x           (natms, atmnam, dipx, dipy, dipz,
     x            vdxx, vdyy, vdzz, fdxx, fdyy, fdzz,
     x            buffer)
c
c==============================================================================

      use multibead
      use centroid_module, only: file_initial_dipole_pimd

c     !..................................................

      implicit none

      integer ::  natms
      character(8) :: atmnam(*)

      real(8) :: dipx(*), dipy(*), dipz(*)
      real(8) :: vdxx(*), vdyy(*), vdzz(*)
      real(8) :: fdxx(*), fdyy(*), fdzz(*)

      real(8) :: buffer(*)

#include "mpif.h"

c     !..................................................

      integer :: iatom, j, iatm, ierr

c     !..................................................

      if (bead_rank.gt.0) goto 100

      open(file_initial_dipole_pimd, file='INITIAL_DIPOLE'//bead_suffix)

c     !..................................................

      j = 0
      do iatom = 1, natms
         read(file_initial_dipole_pimd,'(a8,i10)') atmnam(iatom), iatm
         read(file_initial_dipole_pimd,'(3g20.10)')
     x        dipx(iatom),dipy(iatom),dipz(iatom)
         buffer(j+1) = dipx(iatom)
         buffer(j+2) = dipy(iatom)
         buffer(j+3) = dipz(iatom)
         read(file_initial_dipole_pimd,'(3g20.10)')
     x        vdxx(iatom),vdyy(iatom),vdzz(iatom)
         buffer(j+4) = vdxx(iatom)
         buffer(j+5) = vdyy(iatom)
         buffer(j+6) = vdzz(iatom)
         read(file_initial_dipole_pimd,'(3g20.10)')
     x        fdxx(iatom),fdyy(iatom),fdzz(iatom)
         buffer(j+7) = fdxx(iatom)
         buffer(j+8) = fdyy(iatom)
         buffer(j+9) = fdzz(iatom)
         j = j + 9
      enddo

      close(file_initial_dipole_pimd)

 100  continue

      call MPI_BCAST(buffer,9*natms,MPI_DOUBLE_PRECISION,
     x               0, comm_bead, ierr)

      j = 0
      do iatom = 1, natms
         dipx(iatom) = buffer(j+1)
         dipy(iatom) = buffer(j+2)
         dipz(iatom) = buffer(j+3)

         vdxx(iatom) = buffer(j+4)
         vdyy(iatom) = buffer(j+5)
         vdzz(iatom) = buffer(j+6)

         fdxx(iatom) = buffer(j+7)
         fdyy(iatom) = buffer(j+8)
         fdzz(iatom) = buffer(j+9)
         j = j + 9
      enddo

      end subroutine
