!
! === Always Stable Predictor-Corrector (ASPC) method ===
! 
!    following Kolafa, J. Computat. Chem., 25 (2004)
!              http://dx.doi.org/10.1002/jcc.10385
!

!===============================================================================

module aspc_module

!===============================================================================

implicit none

! Parameters

real(8), private, save, allocatable :: dip4x(:)
real(8), private, save, allocatable :: dip4y(:)
real(8), private, save, allocatable :: dip4z(:)
real(8), private, save, allocatable :: dip5x(:)
real(8), private, save, allocatable :: dip5y(:)
real(8), private, save, allocatable :: dip5z(:)
real(8), private, save, allocatable :: dip6x(:)
real(8), private, save, allocatable :: dip6y(:)
real(8), private, save, allocatable :: dip6z(:)

integer, private, save              :: aspc_k
integer, private, save              :: aspc_niter
real(8), private, save, allocatable :: aspc_B(:) ! k

! Subroutines

public :: aspc_init
public :: aspc_predict
public :: aspc_set
public :: aspc_fini

!===============================================================================

contains

!===============================================================================

subroutine aspc_init( k , niter )

  use dipole_moments,   only: dip1x, dip1y, dip1z, &
                              dip2x, dip2y, dip2z, &
                              dip3x, dip3y, dip3z
  use global_variables, only: mxatms

  implicit none

! arguments 

  integer, intent(in) :: k
  integer, intent(in) :: niter

! local variables
 
  integer :: i

  real(8) :: rk, tmp, w1

  aspc_k = k
  aspc_niter = niter
  if(aspc_k .ne. 2 .and. aspc_k .ne. 4) then
      write(*,*) 'ASPC currently only supports k=2 and k=4'
      stop
  endif
  if(aspc_niter .ne. aspc_k+2) then
      write(*,*) 'For some reason, aspc_niter != k+2'
      write(*,*) "aspc_niter = ", aspc_niter
      write(*,*) "aspc_k     = ", aspc_k
  endif

  allocate(aspc_B(k+2))
  allocate(dip4x(mxatms))
  allocate(dip4y(mxatms))
  allocate(dip4z(mxatms))

  if(aspc_k .gt. 2) then
      allocate(dip5x(mxatms))
      allocate(dip5y(mxatms))
      allocate(dip5z(mxatms))
      allocate(dip6x(mxatms))
      allocate(dip6y(mxatms))
      allocate(dip6z(mxatms))
  endif


  ! initialize dipole history
  dip1x(:) = 0.0d0;  dip1y(:) = 0.0d0;  dip1z(:) = 0.0d0
  dip2x(:) = 0.0d0;  dip2y(:) = 0.0d0;  dip2z(:) = 0.0d0
  dip3x(:) = 0.0d0;  dip3y(:) = 0.0d0;  dip3z(:) = 0.0d0
  dip4x(:) = 0.0d0;  dip4y(:) = 0.0d0;  dip4z(:) = 0.0d0

  if(aspc_k .gt. 2) then
      dip5x(:) = 0.0d0;  dip5y(:) = 0.0d0;  dip5z(:) = 0.0d0
      dip6x(:) = 0.0d0;  dip6y(:) = 0.0d0;  dip6z(:) = 0.0d0
  endif

  rk = real(k)
  tmp = 1.0d0/(rk+3.0d0)
  w1 = (4.0d0*rk+6.0d0)

  tmp = 1.0d0 / (rk + 3.0d0)
  do i=1,k+2
    aspc_B(i) = ((-1.0d0)**(i+1)) * real(i) * w1 * tmp
    tmp = tmp * (rk + 2.0d0 - real(i)) / (rk + 3.0d0 + real(i))
  enddo

end subroutine aspc_init

!
! predicted dipole (eq. 6 in Kofala):
! 
! \mu^{predicted}(t+h) = \sum_{j=0}^{k+1} B_{j+1} * \mu(t - j*h)
!

subroutine aspc_predict(istep, natms, dipx, dipy, dipz)

  use dipole_moments,   only: dip1x, dip1y, dip1z, &
                              dip2x, dip2y, dip2z, &
                              dip3x, dip3y, dip3z
  use global_variables, only: mxatms


  implicit none                              

  integer, intent(in)    :: istep, natms
  real(8), intent(inout) :: dipx(mxatms), dipy(mxatms), dipz(mxatms)


  integer                :: i

  ! collect \mu^{predicted}(t+h) in dipx, etc.
  ! here, input dipx  = \mu(t)
  !             dip1x = \mu(t-h)
  !             dip2x = \mu(t-2h)
  !             dip3x = \mu(t-3h), etc.

  if( istep > aspc_niter ) then

      if(aspc_k .eq. 4)then
          do i=1,natms
             
             dipx(i) = aspc_B(1)*dip1x(i) + aspc_B(2)*dip2x(i) &
                     + aspc_B(3)*dip3x(i) + aspc_B(4)*dip4x(i) &
                     + aspc_B(5)*dip5x(i) + aspc_B(6)*dip6x(i)
             dipy(i) = aspc_B(1)*dip1y(i) + aspc_B(2)*dip2y(i) &
                     + aspc_B(3)*dip3y(i) + aspc_B(4)*dip4y(i) &
                     + aspc_B(5)*dip5y(i) + aspc_B(6)*dip6y(i)
             dipz(i) = aspc_B(1)*dip1z(i) + aspc_B(2)*dip2z(i) &
                     + aspc_B(3)*dip3z(i) + aspc_B(4)*dip4z(i) &
                     + aspc_B(5)*dip5z(i) + aspc_B(6)*dip6z(i)
             
          end do
      else
          do i=1,natms
             
             dipx(i) = aspc_B(1)*dip1x(i) + aspc_B(2)*dip2x(i) &
                     + aspc_B(3)*dip3x(i) + aspc_B(4)*dip4x(i)
             dipy(i) = aspc_B(1)*dip1y(i) + aspc_B(2)*dip2y(i) &
                     + aspc_B(3)*dip3y(i) + aspc_B(4)*dip4y(i)
             dipz(i) = aspc_B(1)*dip1z(i) + aspc_B(2)*dip2z(i) &
                     + aspc_B(3)*dip3z(i) + aspc_B(4)*dip4z(i)
             
          end do
      endif

  end if

  ! now, dipx = \mu^{predicted}(t+h)

end subroutine aspc_predict

!
! Advance dipole histories to new timestep
! 

subroutine aspc_set(natms, dipx, dipy, dipz)

  use dipole_moments,   only: dip1x, dip1y, dip1z, &
                              dip2x, dip2y, dip2z, &
                              dip3x, dip3y, dip3z
  use global_variables, only: mxatms


  implicit none                              

  integer, intent(in)    :: natms
  real(8), intent(inout) :: dipx(mxatms), dipy(mxatms), dipz(mxatms)

  integer                :: i

  if(aspc_k .eq. 4)then
      do i=1,natms

          dip6x(i) = dip5x(i)  ! (n-6)dt
          dip6y(i) = dip5y(i)  ! (n-6)dt
          dip6z(i) = dip5z(i)  ! (n-6)dt

          dip5x(i) = dip4x(i)  ! (n-5)dt
          dip5y(i) = dip4y(i)  ! (n-5)dt
          dip5z(i) = dip4z(i)  ! (n-5)dt

          dip4x(i) = dip3x(i)  ! (n-4)dt
          dip4y(i) = dip3y(i)  ! (n-4)dt
          dip4z(i) = dip3z(i)  ! (n-4)dt

          dip3x(i) = dip2x(i)  ! (n-3)dt
          dip3y(i) = dip2y(i)  ! (n-3)dt
          dip3z(i) = dip2z(i)  ! (n-3)dt

          dip2x(i) = dip1x(i)  ! (n-2)dt
          dip2y(i) = dip1y(i)  ! (n-2)dt
          dip2z(i) = dip1z(i)  ! (n-2)dt

          dip1x(i) = dipx(i)   ! (n-1)dt
          dip1y(i) = dipy(i)   ! (n-1)dt
          dip1z(i) = dipz(i)   ! (n-1)dt

      end do
  else ! aspc_k .eq. 2
      do i=1,natms

          dip4x(i) = dip3x(i)  ! (n-4)dt
          dip4y(i) = dip3y(i)  ! (n-4)dt
          dip4z(i) = dip3z(i)  ! (n-4)dt

          dip3x(i) = dip2x(i)  ! (n-3)dt
          dip3y(i) = dip2y(i)  ! (n-3)dt
          dip3z(i) = dip2z(i)  ! (n-3)dt

          dip2x(i) = dip1x(i)  ! (n-2)dt
          dip2y(i) = dip1y(i)  ! (n-2)dt
          dip2z(i) = dip1z(i)  ! (n-2)dt

          dip1x(i) = dipx(i)   ! (n-1)dt
          dip1y(i) = dipy(i)   ! (n-1)dt
          dip1z(i) = dipz(i)   ! (n-1)dt

      end do
  endif

end subroutine aspc_set

! ---------------------------------------------------------------------------- !

subroutine aspc_fini()

  implicit none

  deallocate(aspc_B)
  deallocate(dip4z)
  deallocate(dip4y)
  deallocate(dip4x)

  if(aspc_k .gt. 2) then
      deallocate(dip5z)
      deallocate(dip5y)
      deallocate(dip5x)
      deallocate(dip6z)
      deallocate(dip6y)
      deallocate(dip6x)
  endif


end subroutine aspc_fini

!===============================================================================

end module aspc_module

!===============================================================================

