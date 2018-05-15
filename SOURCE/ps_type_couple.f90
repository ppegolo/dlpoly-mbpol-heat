!
! === Coupling terms of PES for TTM2-F water model ===
!
! Last updated: June 13 2006 by S. Iuchi
!
! Ref. 1: C. J. Burnham and S. S. Xantheas, JCP 116, 5115 (2002)
! Ref. 2: G. S. Fanourgakis and S. S. Xantheas, JPCA 110, 4100 (2006)
! Ref. 3: H. Partridge and D. W. Schwenke, JCP 106, 4618 (1997)
!
! Some part of this routine were written based on the fortran codes mentioned 
! in the above references. See README file for the details.  
! Important comments are included into that file. 
!

subroutine ps_type_couple( roh1, xoh1, yoh1, zoh1, roh2, xoh2, yoh2, zoh2, &
                         & theta, pterm, fdum, vterm )
  
  use ps_type_pes, only: b1, c5z, cmass, cose, idx, idxm, reoh, xm1, xm2
  
  implicit none

! arguments

  real(8), intent(in)  :: roh1, roh2, theta
  real(8), intent(in)  :: xoh1, yoh1, zoh1, xoh2, yoh2, zoh2 

  real(8) :: fdum(3,2), pterm, vterm

! note: xoh1, xoh2 etc. are already divided by roh1 and roh2 in angfrc

! local variables

  integer :: j

  real(8) :: cost
  real(8) :: dum, dumx, dumy, dumz
  real(8) :: dumcorrx, dumcorry, dumcorrz
  real(8) :: dum1x, dum1y, dum1z
  real(8) :: dum2x, dum2y, dum2z
  real(8) :: dum3x, dum3y, dum3z
  real(8) :: dummy1, dummy2
  real(8) :: efac, fmat(0:15,3)
  real(8) :: fxa, fya, fza, fxc, fyc, fzc
  real(8) :: term
  real(8) :: v1, v2
  real(8) :: x1, x2, x3
  
!
! --- set (r1-re/re)^j, (r2-re/re)^j, and (cos(theta)-cos(thetae))^j ---
! 
! Based on the reference code. 
!
  
  cost = cos( theta )

  x1 = ( roh1 - reoh ) / reoh
  x2 = ( roh2 - reoh ) / reoh
  x3 = cost - cose
  
  fmat(:,:) = 0.0d0

  fmat(1,1) = 1.0d0
  fmat(1,2) = 1.0d0
  fmat(1,3) = 1.0d0
  
  do j=2,15
     
     fmat(j,1) = fmat(j-1,1) * x1 
     fmat(j,2) = fmat(j-1,2) * x2
     fmat(j,3)=  fmat(j-1,3) * x3
     
  end do

! note: fmat dimension is set to be (0:15) for force evaluations.
!       fmat(j,xk) denote xk^j-1, so idx(j,*)-1 is correct for 
!       dum*xyz below   

!
! --- potential (pterm), forces, and virial (vterm) terms ---
!
! Potential and force terms were written based on the reference code. 
! Modifications were made to be consistent to the current code.
!

  pterm = 0.0d0  ! for potential
  vterm = 0.0d0  ! for virial

  fxa = 0.0d0;  fya = 0.0d0;  fza = 0.0d0
  fxc = 0.0d0;  fyc = 0.0d0;  fzc = 0.0d0
  
  coeff: do j=2,245
     
     dum =  fmat(idx(j,1),1) * fmat(idx(j,2),2) + &
          & fmat(idx(j,2),1) * fmat(idx(j,1),2)
 
     term  = c5z(j) * dum * fmat(idx(j,3),3)
     
     pterm = pterm + term
     
     dum1x = ( idx(j,1) - 1 ) * fmat(idx(j,1)-1,1) * xoh1 / reoh
     dum1y = ( idx(j,1) - 1 ) * fmat(idx(j,1)-1,1) * yoh1 / reoh
     dum1z = ( idx(j,1) - 1 ) * fmat(idx(j,1)-1,1) * zoh1 / reoh
     
     dum2x = ( idx(j,2) - 1 ) * fmat(idx(j,2)-1,1) * xoh1 / reoh
     dum2y = ( idx(j,2) - 1 ) * fmat(idx(j,2)-1,1) * yoh1 / reoh
     dum2z = ( idx(j,2) - 1 ) * fmat(idx(j,2)-1,1) * zoh1 / reoh
 
     dum3x = ( idx(j,3) - 1 ) * fmat(idx(j,3)-1,3) * ( xoh2 - cost * xoh1 ) / roh1
     dum3y = ( idx(j,3) - 1 ) * fmat(idx(j,3)-1,3) * ( yoh2 - cost * yoh1 ) / roh1
     dum3z = ( idx(j,3) - 1 ) * fmat(idx(j,3)-1,3) * ( zoh2 - cost * zoh1 ) / roh1

     dumx = dum1x * fmat(idx(j,2),2) * fmat(idx(j,3),3) + &
          & dum2x * fmat(idx(j,1),2) * fmat(idx(j,3),3) + dum * dum3x
     dumy = dum1y * fmat(idx(j,2),2) * fmat(idx(j,3),3) + &
          & dum2y * fmat(idx(j,1),2) * fmat(idx(j,3),3) + dum * dum3y
     dumz = dum1z * fmat(idx(j,2),2) * fmat(idx(j,3),3) + &
          & dum2z * fmat(idx(j,1),2) * fmat(idx(j,3),3) + dum * dum3z

     fxa = fxa - c5z(j) * dumx  
     fya = fya - c5z(j) * dumy  
     fza = fza - c5z(j) * dumz  

     dum1x = ( idx(j,2) - 1 ) * fmat(idx(j,2)-1,2) * xoh2 / reoh
     dum1y = ( idx(j,2) - 1 ) * fmat(idx(j,2)-1,2) * yoh2 / reoh
     dum1z = ( idx(j,2) - 1 ) * fmat(idx(j,2)-1,2) * zoh2 / reoh
     
     dum2x = ( idx(j,1) - 1 ) * fmat(idx(j,1)-1,2) * xoh2 / reoh
     dum2y = ( idx(j,1) - 1 ) * fmat(idx(j,1)-1,2) * yoh2 / reoh
     dum2z = ( idx(j,1) - 1 ) * fmat(idx(j,1)-1,2) * zoh2 / reoh

     dum3x = ( idx(j,3) - 1 ) * fmat(idx(j,3)-1,3) * ( xoh1 - cost * xoh2 ) / roh2
     dum3y = ( idx(j,3) - 1 ) * fmat(idx(j,3)-1,3) * ( yoh1 - cost * yoh2 ) / roh2
     dum3z = ( idx(j,3) - 1 ) * fmat(idx(j,3)-1,3) * ( zoh1 - cost * zoh2 ) / roh2

     dumx = fmat(idx(j,1),1) * dum1x * fmat(idx(j,3),3) + &
          & fmat(idx(j,2),1) * dum2x * fmat(idx(j,3),3) + dum * dum3x
     dumy = fmat(idx(j,1),1) * dum1y * fmat(idx(j,3),3) + &
          & fmat(idx(j,2),1) * dum2y * fmat(idx(j,3),3) + dum * dum3y
     dumz = fmat(idx(j,1),1) * dum1z * fmat(idx(j,3),3) + &
          & fmat(idx(j,2),1) * dum2z * fmat(idx(j,3),3) + dum * dum3z

     fxc = fxc - c5z(j) * dumx
     fyc = fyc - c5z(j) * dumy
     fzc = fzc - c5z(j) * dumz

     vterm = vterm + c5z(j) * ( &
          & roh1 * fmat(idx(j,2),2) * ( idx(j,1) - 1 ) * fmat(idx(j,1)-1,1) / reoh + &
          & roh2 * fmat(idx(j,1),1) * ( idx(j,2) - 1 ) * fmat(idx(j,2)-1,2) / reoh + &
          & roh1 * fmat(idx(j,1),2) * ( idx(j,2) - 1 ) * fmat(idx(j,2)-1,1) / reoh + &
          & roh2 * fmat(idx(j,2),1) * ( idx(j,1) - 1 ) * fmat(idx(j,1)-1,2) / reoh ) * fmat(idx(j,3),3)

  end do coeff

!
! --- remove mass correction: potential (pterm), forces, and virial (vterm) terms ---
!
! This part was written based on the reference code. 
! Modifications were made to be consistent to the current implementation.
!

  v1 = 0.0d0;  v2 = 0.0d0

  coeff_mass: do j=1,9

     v1 = v1 + cmass(j) * fmat(idxm(j,1),1) * fmat(idxm(j,2),2) * fmat(idxm(j,3),3)
     v2 = v2 + cmass(j) * fmat(idxm(j,2),1) * fmat(idxm(j,1),2) * fmat(idxm(j,3),3)
     
     dum1x = ( idxm(j,1) - 1 ) * fmat(idxm(j,1)-1,1) * xoh1 / reoh
     dum1y = ( idxm(j,1) - 1 ) * fmat(idxm(j,1)-1,1) * yoh1 / reoh
     dum1z = ( idxm(j,1) - 1 ) * fmat(idxm(j,1)-1,1) * zoh1 / reoh
     
     dum2x = ( idxm(j,2) - 1 ) * fmat(idxm(j,2)-1,1) * xoh1 / reoh
     dum2y = ( idxm(j,2) - 1 ) * fmat(idxm(j,2)-1,1) * yoh1 / reoh
     dum2z = ( idxm(j,2) - 1 ) * fmat(idxm(j,2)-1,1) * zoh1 / reoh
 
     dum3x = ( idxm(j,3) - 1 ) * fmat(idxm(j,3)-1,3) * ( xoh2 - cost * xoh1 ) / roh1
     dum3y = ( idxm(j,3) - 1 ) * fmat(idxm(j,3)-1,3) * ( yoh2 - cost * yoh1 ) / roh1
     dum3z = ( idxm(j,3) - 1 ) * fmat(idxm(j,3)-1,3) * ( zoh2 - cost * zoh1 ) / roh1

     dumcorrx = xm1 * ( dum1x * fmat(idxm(j,2),2) * fmat(idxm(j,3),3)   &
          &           + fmat(idxm(j,1),1) * fmat(idxm(j,2),2) * dum3x ) &
          &   + xm2 * ( dum2x * fmat(idxm(j,1),2) * fmat(idxm(j,3),3)   &
          &           + fmat(idxm(j,2),1) * fmat(idxm(j,1),2) * dum3x ) 
     dumcorry = xm1 * ( dum1y * fmat(idxm(j,2),2) * fmat(idxm(j,3),3)   &
          &           + fmat(idxm(j,1),1) * fmat(idxm(j,2),2) * dum3y ) &
          &   + xm2 * ( dum2y * fmat(idxm(j,1),2) * fmat(idxm(j,3),3)   &
          &           + fmat(idxm(j,2),1) * fmat(idxm(j,1),2) * dum3y ) 
     dumcorrz = xm1 * ( dum1z * fmat(idxm(j,2),2) * fmat(idxm(j,3),3)   &
          &           + fmat(idxm(j,1),1) * fmat(idxm(j,2),2) * dum3z ) &
          &   + xm2 * ( dum2z * fmat(idxm(j,1),2) * fmat(idxm(j,3),3)   &
          &           + fmat(idxm(j,2),1) * fmat(idxm(j,1),2) * dum3z ) 

     fxa = fxa - cmass(j) * dumcorrx
     fya = fya - cmass(j) * dumcorry
     fza = fza - cmass(j) * dumcorrz

     dum1x = ( idxm(j,2) - 1 ) * fmat(idxm(j,2)-1,2) * xoh2 / reoh
     dum1y = ( idxm(j,2) - 1 ) * fmat(idxm(j,2)-1,2) * yoh2 / reoh
     dum1z = ( idxm(j,2) - 1 ) * fmat(idxm(j,2)-1,2) * zoh2 / reoh
     
     dum2x = ( idxm(j,1) - 1 ) * fmat(idxm(j,1)-1,2) * xoh2 / reoh
     dum2y = ( idxm(j,1) - 1 ) * fmat(idxm(j,1)-1,2) * yoh2 / reoh
     dum2z = ( idxm(j,1) - 1 ) * fmat(idxm(j,1)-1,2) * zoh2 / reoh

     dum3x = ( idxm(j,3) - 1 ) * fmat(idxm(j,3)-1,3) * ( xoh1 - cost * xoh2 ) / roh2
     dum3y = ( idxm(j,3) - 1 ) * fmat(idxm(j,3)-1,3) * ( yoh1 - cost * yoh2 ) / roh2
     dum3z = ( idxm(j,3) - 1 ) * fmat(idxm(j,3)-1,3) * ( zoh1 - cost * zoh2 ) / roh2

     dumcorrx = xm1 * ( fmat(idxm(j,1),1) * dum1x * fmat(idxm(j,3),3)   &
          &           + fmat(idxm(j,1),1) * fmat(idxm(j,2),2) * dum3x ) &
          &   + xm2 * ( fmat(idxm(j,2),1) * dum2x * fmat(idxm(j,3),3)   &
          &           + fmat(idxm(j,2),1) * fmat(idxm(j,1),2) * dum3x ) 
     dumcorry = xm1 * ( fmat(idxm(j,1),1) * dum1y * fmat(idxm(j,3),3)   &
          &           + fmat(idxm(j,1),1) * fmat(idxm(j,2),2) * dum3y ) &
          &   + xm2 * ( fmat(idxm(j,2),1) * dum2y * fmat(idxm(j,3),3)   &
          &           + fmat(idxm(j,2),1) * fmat(idxm(j,1),2) * dum3y ) 
     dumcorrz = xm1 * ( fmat(idxm(j,1),1) * dum1z * fmat(idxm(j,3),3)   &
          &           + fmat(idxm(j,1),1) * fmat(idxm(j,2),2) * dum3z ) &
          &   + xm2 * ( fmat(idxm(j,2),1) * dum2z * fmat(idxm(j,3),3)   &
          &           + fmat(idxm(j,2),1) * fmat(idxm(j,1),2) * dum3z ) 

     fxc = fxc - cmass(j) * dumcorrx
     fyc = fyc - cmass(j) * dumcorry
     fzc = fzc - cmass(j) * dumcorrz

     vterm = vterm + cmass(j) * ( &
          & xm1 * roh1 * fmat(idxm(j,2),2) * ( idxm(j,1) - 1 ) * fmat(idxm(j,1)-1,1) / reoh + &
          & xm1 * roh2 * fmat(idxm(j,1),1) * ( idxm(j,2) - 1 ) * fmat(idxm(j,2)-1,2) / reoh + &
          & xm2 * roh1 * fmat(idxm(j,1),2) * ( idxm(j,2) - 1 ) * fmat(idxm(j,2)-1,1) / reoh + &
          & xm2 * roh2 * fmat(idxm(j,2),1) * ( idxm(j,1) - 1 ) * fmat(idxm(j,1)-1,2) / reoh ) * fmat(idxm(j,3),3)

  end do coeff_mass

  pterm = pterm + xm1 * v1 + xm2 * v2

!
! --- multiply prefactor ---
!
! note: xoh1, xoh2 etc. are already divided by roh1 and roh2 in angfrc
!
! Based on the reference code. 
! Modifications were made to be consistent to the parent subroutine. 
!

  efac = exp( -b1 * ( ( roh1 - reoh )**2 + ( roh2 - reoh )**2 ) )

  term  = pterm
  pterm = term * efac + c5z(1) 

! correction term with respect to isolated monomer(cjb)
  pterm = pterm + 2.0384848702093006d-6 * 2625.50d0 / 100.0d0
  
  fxa = fxa * efac + 2.0d0 * b1 * ( roh1 - reoh ) * xoh1 * efac * term
  fya = fya * efac + 2.0d0 * b1 * ( roh1 - reoh ) * yoh1 * efac * term
  fza = fza * efac + 2.0d0 * b1 * ( roh1 - reoh ) * zoh1 * efac * term

  fxc = fxc * efac + 2.0d0 * b1 * ( roh2 - reoh ) * xoh2 * efac * term
  fyc = fyc * efac + 2.0d0 * b1 * ( roh2 - reoh ) * yoh2 * efac * term
  fzc = fzc * efac + 2.0d0 * b1 * ( roh2 - reoh ) * zoh2 * efac * term

  dummy1 = -2.0d0 * b1 * ( roh1 - reoh ) * roh1 * efac
  dummy2 = -2.0d0 * b1 * ( roh2 - reoh ) * roh2 * efac
  vterm = dummy1 * term + dummy2 * term + efac * vterm
  
!
! --- final forces ---
!

  fdum(1,1) = fxa;  fdum(2,1) = fya;  fdum(3,1) = fza
  fdum(1,2) = fxc;  fdum(2,2) = fyc;  fdum(3,2) = fzc

  return

end subroutine ps_type_couple
