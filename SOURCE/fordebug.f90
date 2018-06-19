!
! === three point numerical differentiation for debug ===
!
! Written by S. Iuchi
!
! Last updated: Oct 25, 2005
!

subroutine fordebug( lreturn, ixyz, icount, ilfd, natms, delta, eng, fanaly, xx )

  use global_variables, only: mxatms

  implicit none

! arguments
  
  logical :: lreturn

  integer, intent(in) :: ixyz, natms

  integer :: icount, ilfd
  
  real(8), intent(in) :: delta, eng
  real(8), intent(in) :: fanaly(3,mxatms)

  real(8) :: xx(mxatms)

! local variables

  real(8), save :: eng_plus, eng_minus
  real(8) :: fdum

!
! --- diff. ---
!

  lreturn = .false.
  
  icount = icount + 1
  
  if( mod( icount, 3 ) == 1 ) then  ! move +delta
            
     ilfd = ilfd + 1
     xx(ilfd) = xx(ilfd) + delta
     
     lreturn=.true.
     
  else if( mod( icount, 3 ) == 2 ) then ! store plus data
     
     eng_plus = eng
        
     xx(ilfd) = xx(ilfd) - 2.0d0 * delta ! move -delta

     lreturn=.true.
            
  else if( mod( icount, 3 ) == 0 ) then  ! store minus data

     eng_minus = eng
            
     xx(ilfd) = xx(ilfd) + delta  ! initialize

  end if

  if( lreturn ) return

  fdum = -( eng_plus - eng_minus ) / delta / 2.0d0
         
  write(6,'(i6,i3,3f20.10)') ilfd, ixyz, fdum, fanaly(ixyz,ilfd),&
       & abs( fdum - fanaly(ixyz,ilfd) )
  
  if( ilfd /= natms ) lreturn=.true.

  return

end subroutine fordebug
  
