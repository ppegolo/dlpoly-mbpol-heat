!============================================================================

module mbpol_utils

!============================================================================

implicit none

!----------------------------------------------------------------------------

public :: mbpol_image
public :: mbpol_image_dimer
public :: mbpol_utils_set_cell

!----------------------------------------------------------------------------

private

integer, save :: my_imcon = 0
real(8), save :: my_cell(9)

!============================================================================

contains

!============================================================================

subroutine mbpol_utils_set_cell(imcon, cell)

   implicit none

   integer, intent(in) :: imcon
   real(8), intent(in) :: cell(*)

   my_imcon = imcon
   my_cell(1:9) = cell(1:9)

end subroutine mbpol_utils_set_cell

!----------------------------------------------------------------------------

subroutine mbpol_image(n, xyz)

   implicit none

   integer, intent(in) :: n
   real(8), intent(inout) :: xyz(*)

   integer :: k, k3
   real(8) :: a1, b1, c1

   if (my_imcon.eq.0) then
      return
   else if (my_imcon.eq.1) then
      a1 = 1.d0/my_cell(1)
      do k = 1, 3*n
         xyz(k) = xyz(k) - my_cell(1)*nint(xyz(k)*a1)
      end do
   else if (my_imcon.eq.2) then
      a1 = 1.d0/my_cell(1)
      b1 = 1.d0/my_cell(5)
      c1 = 1.d0/my_cell(9)

      do k = 1, n
         k3 = 3*k
         xyz(k3-2) = xyz(k3-2) - my_cell(1)*nint(xyz(k3-2)/my_cell(1))
         xyz(k3-1) = xyz(k3-1) - my_cell(5)*nint(xyz(k3-1)/my_cell(5))
         xyz(k3-0) = xyz(k3-0) - my_cell(9)*nint(xyz(k3-0)/my_cell(9))
      end do
   else
      print *, ' ** Error ** : '//__FILE__//': not implemented'
   end if

end subroutine mbpol_image

!----------------------------------------------------------------------------

subroutine mbpol_image_dimer(crd)

   implicit none

   real(8), intent(inout) :: crd(*) ! O H H O H H

   real(8) :: tmp(15) ! o1h11, o1h12, o2h21, o2h22, o1o2
   integer :: k

   do k = 1, 3
      tmp(0 + k) = crd(3 + k) - crd(k)
      tmp(3 + k) = crd(6 + k) - crd(k)
      tmp(6 + k) = crd(12 + k) - crd(k + 9)
      tmp(9 + k) = crd(15 + k) - crd(k + 9)
      tmp(12 + k) = crd(k + 9) - crd(k)
   end do

   call mbpol_image(5, crd)

   do k = 1, 3
      crd(k + 3) = crd(k) + tmp(k + 0) ! H11
      crd(k + 6) = crd(k) + tmp(k + 3) ! H12
      crd(k + 9) = crd(k) + tmp(k + 12) ! O2
      crd(k + 12) = crd(k + 9) + tmp(k + 6) ! H21
      crd(k + 15) = crd(k + 9) + tmp(k + 9) ! H22
   end do

end subroutine mbpol_image_dimer

!============================================================================

end module mbpol_utils

!============================================================================
