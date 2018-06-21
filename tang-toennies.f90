!
!  T(n, x) = 1 - exp(-x)*(1 + x/1! + x^2/2! + x^3/3! + ... + x^n/n!)
!  diff(T(n, x), x) = T(n - 1, x) - T(n, x) = exp(-x)*x^n/n!
!

real(8) function tang_toennies(n, x)

   implicit none

   integer, intent(in) :: n
   real(8), intent(in) :: x

   integer :: nn
   real(8) :: sum, term

   sum = 1.d0
   do nn = n, 1, -1
      sum = 1.d0 + sum*x/nn
   end do

   tang_toennies = 1.d0 - sum*dexp(-x)

   if (abs(tang_toennies).lt.1.0d-8) then

      term = 1.d0
      do nn = n, 1, -1
         term = term*x/nn
      end do

      sum = 0.d0
      do nn = n + 1, 1000
         term = term*x/nn
         sum = sum + term
         if (abs(term/sum) < 1.0d-8) exit
      end do

      tang_toennies = sum*dexp(-x)
   end if

end function tang_toennies

real(8) function tang_toennies_lr(n, x)

   implicit none

   integer, intent(in) :: n
   real(8), intent(in) :: x

   real(8) :: tang_toennies

   tang_toennies_lr = tang_toennies(n, x)/x**(n - 3) &
   & + exp(-x)*(x*x*(3.d0 + x) + 6*(1.d0 + x))/factorial(n)

   tang_toennies_lr = tang_toennies_lr/(n - 3)

contains

real(8) function factorial(z)

   implicit none

   integer, intent(in) :: z

   integer :: i

   factorial = 1.d0

   do i = z, 1, -1
      factorial = factorial*z
   end do

end function factorial

end function tang_toennies_lr
