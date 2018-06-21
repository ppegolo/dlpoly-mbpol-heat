!
! === guess dipole using second(or third)-order extrapolation ===
!
! Last updated: 22 Nov 2006 by S. Iuchi
! 
! (Comments updated: 13 Feb 2007)
!

subroutine set_guess_dipole( nstep, natms, dipx, dipy, dipz )

  use dipole_moments,   only: dip1x, dip1y, dip1z, dip2x, dip2y, dip2z, &
                            & dip3x, dip3y, dip3z
  use global_variables, only: mxatms

  implicit none

! arguments 

  integer, intent(in) :: nstep, natms

  real(8) :: dipx(mxatms), dipy(mxatms), dipz(mxatms)
  
! local variables
 
  integer :: i

  real(8) :: dip0x(mxatms), dip0y(mxatms), dip0z(mxatms)
  
!
! --- initialization ---
!

  if( nstep == 1 ) then

     dip1x(:) = 0.0d0;  dip1y(:) = 0.0d0;  dip1z(:) = 0.0d0
     dip2x(:) = 0.0d0;  dip2y(:) = 0.0d0;  dip2z(:) = 0.0d0
     dip3x(:) = 0.0d0;  dip3y(:) = 0.0d0;  dip3z(:) = 0.0d0

  end if

!
! --- set dipole at (n-1)dt ---
!

  do i=1,natms

     dip0x(i) = dipx(i)
     dip0y(i) = dipy(i)
     dip0z(i) = dipz(i)

  end do

!!$!
!!$! --- set guess dipole using second-order extrapolation ---
!!$!  
!!$
!!$  do i=1,natms
!!$
!!$     dipx(i) = 3.0d0 * dip0x(i) - 3.0d0 * dip1x(i) + dip2x(i)
!!$     dipy(i) = 3.0d0 * dip0y(i) - 3.0d0 * dip1y(i) + dip2y(i)
!!$     dipz(i) = 3.0d0 * dip0z(i) - 3.0d0 * dip1z(i) + dip2z(i)
!!$
!!$  end do

!
! --- set guess dipole using third-order extrapolation ---
!  

  if( nstep >= 5 ) then

     do i=1,natms
        
        dipx(i) = 4.0d0 * dip0x(i) - 6.0d0 * dip1x(i) + 4.0d0 * dip2x(i) - dip3x(i)
        dipy(i) = 4.0d0 * dip0y(i) - 6.0d0 * dip1y(i) + 4.0d0 * dip2y(i) - dip3y(i)
        dipz(i) = 4.0d0 * dip0z(i) - 6.0d0 * dip1z(i) + 4.0d0 * dip2z(i) - dip3z(i)
        
     end do

  end if

!
! --- set dipole at (n-2)dt, (n-3)dt, and (n-4)dt ---
!

  do i=1,natms

     dip3x(i) = dip2x(i)  ! (n-4)dt
     dip3y(i) = dip2y(i)  ! (n-4)dt
     dip3z(i) = dip2z(i)  ! (n-4)dt

     dip2x(i) = dip1x(i)  ! (n-3)dt
     dip2y(i) = dip1y(i)  ! (n-3)dt
     dip2z(i) = dip1z(i)  ! (n-3)dt

     dip1x(i) = dip0x(i)  ! (n-2)dt
     dip1y(i) = dip0y(i)  ! (n-2)dt
     dip1z(i) = dip0z(i)  ! (n-2)dt

  end do

end subroutine set_guess_dipole
