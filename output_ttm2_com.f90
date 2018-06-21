!
! === output COM data (TTM2 model) ===
!
! Last updated: Oct 14, 2005 by S. Iuchi
!
! Note:
! N_atom do not include M sites; no dipoles on M sites
!

subroutine output_ttm2_com( imcon, N_wat, time, cell, weight, listttm2, &
                          & xxx, yyy, zzz )
 
  use multibead, only: bead_suffix
  use global_variables, only: ncomwat, mxatms

  implicit none

! arguments
  
  integer, intent(in) :: imcon
  integer, intent(in) :: listttm2(mxatms)
  integer, intent(in) :: N_wat
  
  real(8), intent(in) :: cell(9)
  real(8), intent(in) :: time
  real(8), intent(in) :: xxx(mxatms), yyy(mxatms), zzz(mxatms)
  real(8), intent(in) :: weight(mxatms)

! local variables

  integer :: i, iox, ih1, ih2
  
  real(8) :: comx(N_wat), comy(N_wat), comz(N_wat)
  real(8) :: dumx, dumy, dumz
  real(8) :: summass
  real(8) :: xh1, yh1, zh1, xh2, yh2, zh2

!
! --- COM of TTM2 water ---
!

  comx(:) = 0.0d0;  comy(:) = 0.0d0;  comz(:) = 0.0d0
  
  do i=1,N_wat
     
     iox = listttm2(3*i-2)
     ih1 = listttm2(3*i-1)
     ih2 = listttm2(3*i)

     dumx = xxx(iox) - xxx(ih1)
     dumy = yyy(iox) - yyy(ih1)
     dumz = zzz(iox) - zzz(ih1)

     call images( imcon, 0, 1, 1, cell, dumx, dumy, dumz )

     xh1 = xxx(iox) - dumx
     yh1 = yyy(iox) - dumy
     zh1 = zzz(iox) - dumz

     dumx = xxx(iox) - xxx(ih2)
     dumy = yyy(iox) - yyy(ih2)
     dumz = zzz(iox) - zzz(ih2)

     call images( imcon, 0, 1, 1, cell, dumx, dumy, dumz )

     xh2 = xxx(iox) - dumx
     yh2 = yyy(iox) - dumy
     zh2 = zzz(iox) - dumz
     
     comx(i) = weight(iox) * xxx(iox) + weight(ih1) * xh1 + weight(ih2) * xh2 
     comy(i) = weight(iox) * yyy(iox) + weight(ih1) * yh1 + weight(ih2) * yh2
     comz(i) = weight(iox) * zzz(iox) + weight(ih1) * zh1 + weight(ih2) * zh2
     
     summass = weight(iox) + weight(ih1) + weight(ih2)

     comx(i) = comx(i) / summass
     comy(i) = comy(i) / summass
     comz(i) = comz(i) / summass

 end do
  
!
! --- open files ---
!

  open(ncomwat,file='COM'//bead_suffix,position='append')

!
! --- output COM of waters ---
!

  write(ncomwat,'(f20.10)') time

  do i=1,N_wat
     
     write(ncomwat,'(i5,3f20.10)') i, comx(i), comy(i), comz(i)

  end do

!
! --- close files ---
!

  close(ncomwat)

  return

end subroutine output_ttm2_com
