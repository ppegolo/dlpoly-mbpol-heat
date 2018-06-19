!
! === output docomposed forces (TTM2 model) === 
!
! Last updated: Mar 17, 2007 by S. Iuchi
!

subroutine output_ttm2_force( atmnam, natm, nstep, nttm2, n_water, listttm2, gamma, xxx, yyy, zzz ) 

  use global_variables, only: nforce1, nforce2, nforce3, nforce4, nforce5, mxatms
  use ttm_forces,       only: fx_ttm_per,   fy_ttm_per,   fz_ttm_per,   &
                            & fx_ttm_ind,   fy_ttm_ind,   fz_ttm_ind,   &
                            & fx_ttm_intra, fy_ttm_intra, fz_ttm_intra, &
                            & fx_ttm_lj,    fy_ttm_lj,    fz_ttm_lj
  use unit_parameters,  only: bohr

  implicit none

! arguments
  
  character(len=8), intent(in) :: atmnam(mxatms)

  integer, intent(in) :: natm, nstep, n_water, nttm2
  integer, intent(in) :: listttm2(mxatms)
  
  real(8), intent(in) :: gamma
  real(8), intent(in) :: xxx(mxatms), yyy(mxatms), zzz(mxatms)

! local variables

  integer :: i, iox, ih1, ih2
  integer :: mttm2
  
  real(8), parameter :: engconv = 3.80880d-6  ! a.u. / (10J/mol)
  
  real(8) :: fx(natm), fy(natm), fz(natm)

  real(8) :: fx_ind(natm+n_water),   fy_ind(natm+n_water),   fz_ind(natm+n_water)
  real(8) :: fx_per(natm+n_water),   fy_per(natm+n_water),   fz_per(natm+n_water)
  real(8) :: fx_lj(natm+n_water),    fy_lj(natm+n_water),    fz_lj(natm+n_water)
  real(8) :: fx_intra(natm+n_water), fy_intra(natm+n_water), fz_intra(natm+n_water)
  real(8) :: fx_tot(natm+n_water),   fy_tot(natm+n_water),   fz_tot(natm+n_water)

!
! --- save 4-sites data ---
!

  fx_ind(:) = fx_ttm_ind(:)
  fy_ind(:) = fy_ttm_ind(:)
  fz_ind(:) = fz_ttm_ind(:)

  fx_per(:) = fx_ttm_per(:)
  fy_per(:) = fy_ttm_per(:)
  fz_per(:) = fz_ttm_per(:)

  fx_lj(:) = fx_ttm_lj(:)
  fy_lj(:) = fy_ttm_lj(:)
  fz_lj(:) = fz_ttm_lj(:)

  fx_intra(:) = fx_ttm_intra(:)
  fy_intra(:) = fy_ttm_intra(:)
  fz_intra(:) = fz_ttm_intra(:)

!
! --- distribute forces on the m-site to oxygen and hydrogen sites ---
!

  mttm2 = listttm2(nttm2)

  do i=1,n_water
     
     mttm2 = mttm2 + 1

     iox = listttm2(3*i-2)
     ih1 = listttm2(3*i-1)
     ih2 = listttm2(3*i)

! CC part

     fx_ttm_per(iox) = fx_ttm_per(iox) + ( 1.0d0 - gamma ) * fx_ttm_per(mttm2)
     fy_ttm_per(iox) = fy_ttm_per(iox) + ( 1.0d0 - gamma ) * fy_ttm_per(mttm2)
     fz_ttm_per(iox) = fz_ttm_per(iox) + ( 1.0d0 - gamma ) * fz_ttm_per(mttm2)

     fx_ttm_per(ih1) = fx_ttm_per(ih1) + fx_ttm_per(mttm2) * gamma / 2.0d0
     fy_ttm_per(ih1) = fy_ttm_per(ih1) + fy_ttm_per(mttm2) * gamma / 2.0d0
     fz_ttm_per(ih1) = fz_ttm_per(ih1) + fz_ttm_per(mttm2) * gamma / 2.0d0

     fx_ttm_per(ih2) = fx_ttm_per(ih2) + fx_ttm_per(mttm2) * gamma / 2.0d0
     fy_ttm_per(ih2) = fy_ttm_per(ih2) + fy_ttm_per(mttm2) * gamma / 2.0d0
     fz_ttm_per(ih2) = fz_ttm_per(ih2) + fz_ttm_per(mttm2) * gamma / 2.0d0

! DC-DD part 

     fx_ttm_ind(iox) = fx_ttm_ind(iox) + ( 1.0d0 - gamma ) * fx_ttm_ind(mttm2)
     fy_ttm_ind(iox) = fy_ttm_ind(iox) + ( 1.0d0 - gamma ) * fy_ttm_ind(mttm2)
     fz_ttm_ind(iox) = fz_ttm_ind(iox) + ( 1.0d0 - gamma ) * fz_ttm_ind(mttm2)

     fx_ttm_ind(ih1) = fx_ttm_ind(ih1) + fx_ttm_ind(mttm2) * gamma / 2.0d0
     fy_ttm_ind(ih1) = fy_ttm_ind(ih1) + fy_ttm_ind(mttm2) * gamma / 2.0d0
     fz_ttm_ind(ih1) = fz_ttm_ind(ih1) + fz_ttm_ind(mttm2) * gamma / 2.0d0

     fx_ttm_ind(ih2) = fx_ttm_ind(ih2) + fx_ttm_ind(mttm2) * gamma / 2.0d0
     fy_ttm_ind(ih2) = fy_ttm_ind(ih2) + fy_ttm_ind(mttm2) * gamma / 2.0d0
     fz_ttm_ind(ih2) = fz_ttm_ind(ih2) + fz_ttm_ind(mttm2) * gamma / 2.0d0

! LJ part

     fx_ttm_lj(iox) = fx_ttm_lj(iox) + ( 1.0d0 - gamma ) * fx_ttm_lj(mttm2)
     fy_ttm_lj(iox) = fy_ttm_lj(iox) + ( 1.0d0 - gamma ) * fy_ttm_lj(mttm2)
     fz_ttm_lj(iox) = fz_ttm_lj(iox) + ( 1.0d0 - gamma ) * fz_ttm_lj(mttm2)

     fx_ttm_lj(ih1) = fx_ttm_lj(ih1) + fx_ttm_lj(mttm2) * gamma / 2.0d0
     fy_ttm_lj(ih1) = fy_ttm_lj(ih1) + fy_ttm_lj(mttm2) * gamma / 2.0d0
     fz_ttm_lj(ih1) = fz_ttm_lj(ih1) + fz_ttm_lj(mttm2) * gamma / 2.0d0

     fx_ttm_lj(ih2) = fx_ttm_lj(ih2) + fx_ttm_lj(mttm2) * gamma / 2.0d0
     fy_ttm_lj(ih2) = fy_ttm_lj(ih2) + fy_ttm_lj(mttm2) * gamma / 2.0d0
     fz_ttm_lj(ih2) = fz_ttm_lj(ih2) + fz_ttm_lj(mttm2) * gamma / 2.0d0
     
! intra. part

     fx_ttm_intra(iox) = fx_ttm_intra(iox) + ( 1.0d0 - gamma ) * fx_ttm_intra(mttm2)
     fy_ttm_intra(iox) = fy_ttm_intra(iox) + ( 1.0d0 - gamma ) * fy_ttm_intra(mttm2)
     fz_ttm_intra(iox) = fz_ttm_intra(iox) + ( 1.0d0 - gamma ) * fz_ttm_intra(mttm2)

     fx_ttm_intra(ih1) = fx_ttm_intra(ih1) + fx_ttm_intra(mttm2) * gamma / 2.0d0
     fy_ttm_intra(ih1) = fy_ttm_intra(ih1) + fy_ttm_intra(mttm2) * gamma / 2.0d0
     fz_ttm_intra(ih1) = fz_ttm_intra(ih1) + fz_ttm_intra(mttm2) * gamma / 2.0d0

     fx_ttm_intra(ih2) = fx_ttm_intra(ih2) + fx_ttm_intra(mttm2) * gamma / 2.0d0
     fy_ttm_intra(ih2) = fy_ttm_intra(ih2) + fy_ttm_intra(mttm2) * gamma / 2.0d0
     fz_ttm_intra(ih2) = fz_ttm_intra(ih2) + fz_ttm_intra(mttm2) * gamma / 2.0d0
     
  end do

!
! --- set total forces ---
!

  do i=1,natm
     
     fx(i) = fx_ttm_per(i) + fx_ttm_ind(i) + fx_ttm_intra(i) + fx_ttm_lj(i)
     fy(i) = fy_ttm_per(i) + fy_ttm_ind(i) + fy_ttm_intra(i) + fy_ttm_lj(i)
     fz(i) = fz_ttm_per(i) + fz_ttm_ind(i) + fz_ttm_intra(i) + fz_ttm_lj(i)
     
  end do

  do i=1,natm+n_water
     
     fx_tot(i) = fx_per(i) + fx_ind(i) + fx_intra(i) + fx_lj(i)
     fy_tot(i) = fy_per(i) + fy_ind(i) + fy_intra(i) + fy_lj(i)
     fz_tot(i) = fz_per(i) + fz_ind(i) + fz_intra(i) + fz_lj(i)
     
  end do

!
! --- convert from 10J/mol/ang to a.u. ---
!

  do i=1,natm

     fx(i) = fx(i) * engconv * bohr
     fy(i) = fy(i) * engconv * bohr
     fz(i) = fz(i) * engconv * bohr

     fx_ttm_per(i) = fx_ttm_per(i) * engconv * bohr
     fy_ttm_per(i) = fy_ttm_per(i) * engconv * bohr
     fz_ttm_per(i) = fz_ttm_per(i) * engconv * bohr

     fx_ttm_ind(i) = fx_ttm_ind(i) * engconv * bohr
     fy_ttm_ind(i) = fy_ttm_ind(i) * engconv * bohr
     fz_ttm_ind(i) = fz_ttm_ind(i) * engconv * bohr
     
     fx_ttm_intra(i) = fx_ttm_intra(i) * engconv * bohr
     fy_ttm_intra(i) = fy_ttm_intra(i) * engconv * bohr
     fz_ttm_intra(i) = fz_ttm_intra(i) * engconv * bohr

     fx_ttm_lj(i) = fx_ttm_lj(i) * engconv * bohr
     fy_ttm_lj(i) = fy_ttm_lj(i) * engconv * bohr
     fz_ttm_lj(i) = fz_ttm_lj(i) * engconv * bohr

  end do

  do i=1,natm+n_water

     fx_ind(i) = fx_ind(i) * engconv * bohr
     fy_ind(i) = fy_ind(i) * engconv * bohr
     fz_ind(i) = fz_ind(i) * engconv * bohr

     fx_per(i) = fx_per(i) * engconv * bohr
     fy_per(i) = fy_per(i) * engconv * bohr
     fz_per(i) = fz_per(i) * engconv * bohr

     fx_lj(i) = fx_lj(i) * engconv * bohr
     fy_lj(i) = fy_lj(i) * engconv * bohr
     fz_lj(i) = fz_lj(i) * engconv * bohr

     fx_intra(i) = fx_intra(i) * engconv * bohr
     fy_intra(i) = fy_intra(i) * engconv * bohr
     fz_intra(i) = fz_intra(i) * engconv * bohr

     fx_tot(i) = fx_tot(i) * engconv * bohr
     fy_tot(i) = fy_tot(i) * engconv * bohr
     fz_tot(i) = fz_tot(i) * engconv * bohr

  end do
  
!
! --- open files ---
!

!  open(nforce1,file='FORCE_CC'   ,position='append')
!  open(nforce2,file='FORCE_DC_DD',position='append')
!  open(nforce3,file='FORCE_TOTAL',position='append')
!  open(nforce4,file='FORCE_INTRA',position='append')
!  open(nforce5,file='FORCE_LJ'   ,position='append')
  open(896,file='FORCE_LJ_RAW',   position='append')
  open(897,file='FORCE_CC_RAW',   position='append')
  open(898,file='FORCE_TOTAL_RAW',position='append')
  open(899,file='FORCE_DC_DD_RAW',position='append')
  open(951,file='FORCE_LJ_ALL',   position='append')
  open(952,file='FORCE_CC_ALL',   position='append')
  open(953,file='FORCE_TOTAL_ALL',position='append')
  open(954,file='FORCE_DC_DD_ALL',position='append')

!
! --- output components of forces ---
!
!
!  write(nforce1,'(a10,i10)') 'timestep= ', nstep
!  write(nforce2,'(a10,i10)') 'timestep= ', nstep
!  write(nforce3,'(a10,i10)') 'timestep= ', nstep
!  write(nforce4,'(a2,f20.10)') '# ', time
!  write(nforce5,'(a2,f20.10)') '# ', time
  write(896,'(a10,i10)') 'timestep= ', nstep
  write(897,'(a10,i10)') 'timestep= ', nstep
  write(898,'(a10,i10)') 'timestep= ', nstep
  write(899,'(a10,i10)') 'timestep= ', nstep
  write(951,'(a10,i10)') 'timestep= ', nstep
  write(952,'(a10,i10)') 'timestep= ', nstep
  write(953,'(a10,i10)') 'timestep= ', nstep
  write(954,'(a10,i10)') 'timestep= ', nstep

!  do i=1,natm
!     
!     write(nforce1,'(a2,6d25.12)') trim(atmnam(i)), &
!          & xxx(i) / bohr, yyy(i) / bohr, zzz(i) / bohr, &
!          & fx_ttm_per(i), fy_ttm_per(i), fz_ttm_per(i)
!     write(nforce2,'(a2,6d25.12)') trim(atmnam(i)), &
!          & xxx(i) / bohr, yyy(i) / bohr, zzz(i) / bohr, &
!          & fx_ttm_ind(i), fy_ttm_ind(i), fz_ttm_ind(i)
!     write(nforce3,'(a2,6d25.12)') trim(atmnam(i)), &
!          & xxx(i) / bohr, yyy(i) / bohr, zzz(i) / bohr, &
!          & fx(i), fy(i), fz(i)
!     write(nforce4,'(a2,6d25.12)') trim(atmnam(i)), &
!          & xxx(i) / bohr, yyy(i) / bohr, zzz(i) / bohr, &
!          & fx_ttm_intra(i), fy_ttm_intra(i), fz_ttm_intra(i)
!     write(nforce5,'(a2,6d25.12)') trim(atmnam(i)), &
!          & xxx(i) / bohr, yyy(i) / bohr, zzz(i) / bohr, &
!          & fx_ttm_lj(i), fy_ttm_lj(i), fz_ttm_lj(i)
!     
!  end do

  mttm2 = listttm2(nttm2)

  do i=1,n_water
     
     mttm2 = mttm2 + 1
     
     iox = listttm2(3*i-2)
     ih1 = listttm2(3*i-1)
     ih2 = listttm2(3*i)

! polarization part
     
     write(899,'(a2,6d25.12)') trim(atmnam(iox)), &
          & xxx(iox) / bohr, yyy(iox) / bohr, zzz(iox) / bohr, &
          & fx_ind(iox), fy_ind(iox), fz_ind(iox)
     write(899,'(a2,6d25.12)') trim(atmnam(ih1)), &
          & xxx(ih1) / bohr, yyy(ih1) / bohr, zzz(ih1) / bohr, &
          & fx_ind(ih1), fy_ind(ih1), fz_ind(ih1)
     write(899,'(a2,6d25.12)') trim(atmnam(ih2)), &
          & xxx(ih2) / bohr, yyy(ih2) / bohr, zzz(ih2) / bohr, &
          & fx_ind(ih2), fy_ind(ih2), fz_ind(ih2)
     write(899,'(a2,6d25.12)') 'MW', &
          & xxx(mttm2) / bohr, yyy(mttm2) / bohr, zzz(mttm2) / bohr, &
          & fx_ind(mttm2), fy_ind(mttm2), fz_ind(mttm2)

! total

     write(898,'(a2,6d25.12)') trim(atmnam(iox)), &
          & xxx(iox) / bohr, yyy(iox) / bohr, zzz(iox) / bohr, &
          & fx_tot(iox), fy_tot(iox), fz_tot(iox)
     write(898,'(a2,6d25.12)') trim(atmnam(ih1)), &
          & xxx(ih1) / bohr, yyy(ih1) / bohr, zzz(ih1) / bohr, &
          & fx_tot(ih1), fy_tot(ih1), fz_tot(ih1)
     write(898,'(a2,6d25.12)') trim(atmnam(ih2)), &
          & xxx(ih2) / bohr, yyy(ih2) / bohr, zzz(ih2) / bohr, &
          & fx_tot(ih2), fy_tot(ih2), fz_tot(ih2)
     write(898,'(a2,6d25.12)') 'MW', &
          & xxx(mttm2) / bohr, yyy(mttm2) / bohr, zzz(mttm2) / bohr, &
          & fx_tot(mttm2), fy_tot(mttm2), fz_tot(mttm2)

!     
! permanent part
     
     write(897,'(a2,6d25.12)') trim(atmnam(iox)), &
          & xxx(iox) / bohr, yyy(iox) / bohr, zzz(iox) / bohr, &
          & fx_per(iox), fy_per(iox), fz_per(iox)
     write(897,'(a2,6d25.12)') trim(atmnam(ih1)), &
          & xxx(ih1) / bohr, yyy(ih1) / bohr, zzz(ih1) / bohr, &
          & fx_per(ih1), fy_per(ih1), fz_per(ih1)
     write(897,'(a2,6d25.12)') trim(atmnam(ih2)), &
          & xxx(ih2) / bohr, yyy(ih2) / bohr, zzz(ih2) / bohr, &
          & fx_per(ih2), fy_per(ih2), fz_per(ih2)
     write(897,'(a2,6d25.12)') 'MW', &
          & xxx(mttm2) / bohr, yyy(mttm2) / bohr, zzz(mttm2) / bohr, &
          & fx_per(mttm2), fy_per(mttm2), fz_per(mttm2)

! LJ part
     
     write(896,'(a2,6d25.12)') trim(atmnam(iox)), &
          & xxx(iox) / bohr, yyy(iox) / bohr, zzz(iox) / bohr, &
          & fx_lj(iox), fy_lj(iox), fz_lj(iox)
     write(896,'(a2,6d25.12)') trim(atmnam(ih1)), &
          & xxx(ih1) / bohr, yyy(ih1) / bohr, zzz(ih1) / bohr, &
          & fx_lj(ih1), fy_lj(ih1), fz_lj(ih1)
     write(896,'(a2,6d25.12)') trim(atmnam(ih2)), &
          & xxx(ih2) / bohr, yyy(ih2) / bohr, zzz(ih2) / bohr, &
          & fx_lj(ih2), fy_lj(ih2), fz_lj(ih2)
     write(896,'(a2,6d25.12)') 'MW', &
          & xxx(mttm2) / bohr, yyy(mttm2) / bohr, zzz(mttm2) / bohr, &
          & fx_lj(mttm2), fy_lj(mttm2), fz_lj(mttm2)

  end do

  do i=1,natm+n_water
     
     write(951,'(a2,6d25.12)') trim(atmnam(i)), &
          & xxx(i) / bohr, yyy(i) / bohr, zzz(i) / bohr, &
          & fx_lj(i), fy_lj(i), fz_lj(i)
     write(952,'(a2,6d25.12)') trim(atmnam(i)), &
          & xxx(i) / bohr, yyy(i) / bohr, zzz(i) / bohr, &
          & fx_per(i), fy_per(i), fz_per(i)
     write(953,'(a2,6d25.12)') trim(atmnam(i)), &
          & xxx(i) / bohr, yyy(i) / bohr, zzz(i) / bohr, &
          & fx_tot(i), fy_tot(i), fz_tot(i)
     write(954,'(a2,6d25.12)') trim(atmnam(i)), &
          & xxx(i) / bohr, yyy(i) / bohr, zzz(i) / bohr, &
          & fx_ind(i), fy_ind(i), fz_ind(i)

  end do
     
!
! --- close files ---
!
!
!  close(nforce1)
!  close(nforce2)
!  close(nforce3)
!  close(nforce4)
!  close(nforce5)
  close(896)
  close(897)
  close(898)
  close(899)
  close(951)
  close(952)
  close(953)
  close(954)

  return

end subroutine output_ttm2_force
