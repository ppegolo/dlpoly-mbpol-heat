!
! === check system COM velocity: let it stop ===
!
! Written by Satoru Iuchi
!
! Last updated: 28 Nov 2005
!
! flag: .true.  --- let com stop
!       .false. --- just check com velocity
! 

subroutine checkcom( vxx, vyy, vzz, weight, natms, flag, nstep )

  use global_variables, only: mxatms
  
  implicit none

! arguments

  logical, intent(inout) :: flag

  integer, intent(in) :: natms, nstep
    
  real(8), intent(in)    :: weight(mxatms)
  real(8), intent(inout) :: vxx(mxatms), vyy(mxatms), vzz(mxatms)

! local variables

  integer :: i

  real(8), parameter :: tol = 1.0d-10

  real(8) :: comx, comy, comz
  real(8) :: mass

!!$  real(8) :: dumx, dumy, dumz

!
! --- com velocity of the system ---
!

  comx = 0.0d0;   comy = 0.0d0;   comz = 0.0d0
  
  mass = 0.0d0

!  write(6,*) 'CHECK', natms, mxatms

  do i=1,natms
     
     comx = comx + weight(i) * vxx(i)
     comy = comy + weight(i) * vyy(i)
     comz = comz + weight(i) * vzz(i)
     
     mass = mass + weight(i)

!!$     if( mod(i,3) == 0 ) then
!!$
!!$        dumx = weight(i-2) * vxx(i-2) + weight(i-1) * vxx(i-1) + weight(i) * vxx(i) 
!!$        dumy = weight(i-2) * vyy(i-2) + weight(i-1) * vyy(i-1) + weight(i) * vyy(i) 
!!$        dumz = weight(i-2) * vzz(i-2) + weight(i-1) * vzz(i-1) + weight(i) * vzz(i) 
!!$ 
!!$        dumx = dumx / ( weight(i-2) + weight(i-1) + weight(i) )
!!$        dumy = dumy / ( weight(i-2) + weight(i-1) + weight(i) )
!!$        dumz = dumz / ( weight(i-2) + weight(i-1) + weight(i) )
!!$
!!$        write(6,'(i5,3f20.10)') i-2,dumx/vxx(i-2), dumy/vyy(i-2),dumz/vzz(i-2)
!!$
!!$     end if

  end do

  comx = comx / mass
  comy = comy / mass
  comz = comz / mass

  if( ( abs( comx ) > tol ) .or. ( abs( comy ) > tol ) .or. &
                              &  ( abs( comz ) > tol ) ) then
 
     write(6,*) 'Warning: system COM velocity is not 0'
     write(6,'(a5,i7,3f20.15)') 'COM V', nstep + 1, comx, comy, comz

  end if
  
  if( .not. flag ) return
  
! write(6,'(a28,3f20.15)') 'system COM velocity (before)', comx, comy, comz

!
! --- let system COM velocity stop ---
!     

  do i=1,natms
     
     vxx(i) = vxx(i) - comx
     vyy(i) = vyy(i) - comy
     vzz(i) = vzz(i) - comz
     
  end do

!
! --- com velocity of the system: again ---
!

  comx = 0.0d0;   comy = 0.0d0;   comz = 0.0d0
  
  mass = 0.0d0
  
  do i=1,natms
     
     comx = comx + weight(i) * vxx(i)
     comy = comy + weight(i) * vyy(i)
     comz = comz + weight(i) * vzz(i)
     
     mass = mass + weight(i)
     
  end do
  
  comx = comx / mass
  comy = comy / mass
  comz = comz / mass
  
! write(6,'(a28,3f20.15)') 'system COM velocity (after) ', comx, comy, comz

  return

end subroutine checkcom
