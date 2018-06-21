!
! === input for optimization ===
!
! Last update: 28 Dec 2005 by S. Iuchi
!

subroutine input_opt

  use global_variables, only: nread
  use optimize,         only: lgeoopt, tol_opt

  implicit none

! local variables

  logical :: optim
  
  real(8) :: tol
  
  namelist /geoopt/ optim, tol
  
!
! --- read namelist ---
!

  optim = .false.
  
  tol = 1.0d+00 ! DLPOLY internal unit in force
  
  open( nread, file='CONTROL', status='old' )

  rewind nread

!  read( unit=nread, nml=geoopt )
  read( nread, geoopt, end=200 )
  
200 continue
  
  lgeoopt = optim;   tol_opt = tol
  
  close(nread)
  
  return
  
end subroutine input_opt

!
! === output optimization conditions ===
!

subroutine output_opt

  use global_variables, only: nrite
  use optimize,         only: tol_opt
  use unit_parameters,  only: bohr
  
  implicit none

!
! --- output conditions ---
!

  write(nrite,'(/)') 
  write(nrite,'(a21)') 'Geometry Optimization'
  write(nrite,'(a44)') 'Using steepest descent like method in DLPOLY'
  write(nrite,'(a18,f20.10)') 'tolerance in force', tol_opt
  write(nrite,'(a18,f20.10)')  '(in a.u.)        ', &
       & tol_opt * 3.80880d-4 * 1.0d-2 * bohr  
  write(nrite,'(/)')
  
  return
  
end subroutine output_opt

!
! === output forces in the case of optimization method ===
!

subroutine output_force_optimize( fxx, fyy, fzz, natms, nstep, atmnam )
  
  use global_variables, only: mxatms, noutopt

  implicit none
  
! arguments

  character(len=8), intent(in) :: atmnam(mxatms)

  integer, intent(in) :: natms, nstep
  
  real(8), intent(in) :: fxx(mxatms), fyy(mxatms), fzz(mxatms)
  
! local variables

  integer :: i
  
!
! --- output forces ---
! 
 
  open(noutopt,file='FORCES',position='append')
  
  write(noutopt,'(a6,i8)') '# step', nstep
  
  do i=1,natms
     
     write(noutopt,'(a8,3f20.10)') atmnam(i), fxx(i), fyy(i), fzz(i)
     
  end do
  
  close(noutopt)
  
  return

end subroutine output_force_optimize
  
!
! === convergence check ===
!

subroutine check_conv_optimize( fxx, fyy, fzz, xxx, yyy, zzz,&
                              & natms, nstep, nstrun, atmnam )

  use global_variables, only: mxatms, noutopt
  use optimize,         only: tol_opt
  
  implicit none

! arguments

  character(len=8), intent(in) :: atmnam(mxatms)

  integer, intent(in)    :: natms, nstrun
  integer, intent(inout) :: nstep

  real(8), intent(in) :: fxx(mxatms), fyy(mxatms), fzz(mxatms)
  real(8), intent(in) :: xxx(mxatms), yyy(mxatms), zzz(mxatms)

! local variables

  integer :: i, j
  
  real(8) :: diff, tmp
  real(8) :: rxij, ryij, rzij, rij

!
! --- set maximum forces ---
!
  
  diff = 0.0d0
         
  do i=1,natms
            
     tmp = max( abs( fxx(i) ), abs( fyy(i) ), abs( fzz(i) ) )
            
     if( tmp > diff ) diff = tmp

  end do
  
!  write(6,*) 'CHECK', diff
  
  if( diff < tol_opt ) then   ! convergence 
     
     write(6,*) 'geometry is converged'
     
     open(noutopt,file='FORCES',position='append')
     
     write(noutopt,'(/)')
     write(noutopt,*) '--- final geometry ---'
     write(noutopt,*) natms
     write(noutopt,*) 'timestep:     1'

     do i=1,natms
        
        write(noutopt,'(a8,i4,3f20.10)') atmnam(i), i, xxx(i), yyy(i), zzz(i)
        
     end do
     
     write(noutopt,'(/)') 
     
! note: periodic boundary is not applied. 
! In the case with Ewald, box should be large enough to ignore the effect

     do i=1,natms-1
        do j=i+1,natms
           
           rxij = xxx(i) - xxx(j)
           ryij = yyy(i) - yyy(j)
           rzij = zzz(i) - zzz(j)
           
           rij = sqrt( rxij * rxij + ryij * ryij + rzij * rzij )
           
           write(noutopt,'(2a8,2i4,f15.8)') atmnam(i), atmnam(j), i, j, rij
           
        end do
     end do

     close(noutopt)
     
     nstep = nstrun + 1
     
  end if
  
  if( nstep == nstrun ) then
     
     write(6,*) 'geometry is not converged'

  end if
  
  return

end subroutine check_conv_optimize
