c==============================================================================
c
      subroutine centroid_write_dipole 
     x           (natms, atmnam, dipx, dipy, dipz,
     x            vdxx, vdyy, vdzz, fdxx, fdyy, fdzz)
c
c==============================================================================

      use multibead
      use centroid_module, only: file_final_dipole_pimd

c     !..................................................

      implicit none

      integer ::  natms
      character(8) :: atmnam(*)

      real(8) :: dipx(*), dipy(*), dipz(*)
      real(8) :: vdxx(*), vdyy(*), vdzz(*)
      real(8) :: fdxx(*), fdyy(*), fdzz(*)

c     !..................................................

      integer :: iatom

c     !..................................................

      open(file_final_dipole_pimd,file='FINAL_DIPOLE'//bead_suffix)

c     !..................................................

      do iatom = 1, natms
         write(file_final_dipole_pimd,'(a8,i10)') 
     x      atmnam(iatom), iatom
         write(file_final_dipole_pimd,'(3g20.10)') 
     x      dipx(iatom), dipy(iatom), dipz(iatom)
         write(file_final_dipole_pimd,'(3g20.10)') 
     x      vdxx(iatom), vdyy(iatom), vdzz(iatom)
         write(file_final_dipole_pimd,'(3g20.10)') 
     x      fdxx(iatom), fdyy(iatom), fdzz(iatom)
      enddo

      close(file_final_dipole_pimd)
 

      end subroutine
