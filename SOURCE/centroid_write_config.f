c==============================================================================
c
      subroutine centroid_write_config 
     x           (cfgname,nstep,t,levcfg,imcon,cell, 
     x            mxatms,natms,atmnam,nfict,
     x            xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)
      use multibead
      use centroid_module, only: file_final_config_pimd
c
c==============================================================================

      implicit none

      real(8), parameter :: zero = 0.d0

      integer :: levcfg, imcon, nstep, mxatms, natms, nfict

      character(80) :: cfgname
      character(80) :: command

      character(8) :: atmnam(mxatms)

      character(8) :: Msite = 'MW'

      real(8) :: t
      real(8) :: cell(9)

      real(8) :: xxx(*),yyy(*),zzz(*)
      real(8) :: vxx(*),vyy(*),vzz(*)
      real(8) :: fxx(*),fyy(*),fzz(*)

c     !..................................................

      integer :: iatom, ifict

c     !..................................................

      open(file_final_config_pimd,file='new_FINAL_CONFIG'//bead_suffix,
     x     form='formatted' )

c     !..................................................
 
      write(file_final_config_pimd,'(a80)') cfgname

      write(file_final_config_pimd,'(3i10,es25.16)') 
     x   2, imcon, nstep, t

      if ( imcon>0 ) 
     x   write(file_final_config_pimd ,'(3es25.16)') cell(:)

      do iatom = 1, natms

         write(file_final_config_pimd,'(a8,i10)') 
     x      atmnam(iatom), iatom

         write(file_final_config_pimd,'(3es25.16)') 
     x      xxx(iatom), yyy(iatom), zzz(iatom)

         write(file_final_config_pimd,'(3es25.16)') 
     x      vxx(iatom), vyy(iatom), vzz(iatom)

         write(file_final_config_pimd,'(3es25.16)') 
     x      fxx(iatom), fyy(iatom), fzz(iatom)

      enddo

      iatom = natms

      do ifict = 1, Nfict

         iatom = iatom + 1

         write(file_final_config_pimd,'(a8,i10)')
     x      Msite, iatom

         write(file_final_config_pimd,'(3es25.16)')
     x      zero, zero, zero

         write(file_final_config_pimd,'(3es25.16)')
     x      zero, zero, zero

         write(file_final_config_pimd,'(3es25.16)')
     x      zero, zero, zero

      enddo

      close(file_final_config_pimd)
      command = 'mv new_FINAL_CONFIG'//bead_suffix//' FINAL_CONFIG'//
     x   bead_suffix
      call system(command)

      end subroutine
