!
! === molecular dipole and partial charges ===
!
! Last updated: 10 July 2006 by S. Iuchi
!
! Ref. 1: C. J. Burnham and S. S. Xantheas, JCP 116, 5115 (2002)
! Ref. 2: G. S. Fanourgakis and S. S. Xantheas, JPCA 110, 4100 (2006)
! Ref. 3: H. Partridge and D. W. Schwenke, JCP 106, 4618 (1997)
!
! Some part of this routine were written based on the fortran codes mentioned 
! in the above references. See README file for the details.  
! Important comments are included into that file. 
!

subroutine qdetermine(listttm2,n_water,nttm2,lttm3,imcon,cell,chge,gamma,xxx,yyy,zzz)

  use dipole_moments,   only: moldipx, moldipy, moldipz, chkdipx, chkdipy, chkdipz
  use global_variables, only: mxatms
  use ps_type_dms,      only: ldms, lttmfold, adms, bdms, b1dms, coef, &
                            & c0, c1, c2, isump, idx, &
                            & dp1dr1, dp1dr2, dp2dr1, dp2dr2, dp1dcos, dp2dcos, &
                            & vesp, vesp_k, vesp_r, vesp_s, vesp_c, &
                            & vesp_dc_k, vesp_dc_r, vesp_dc_c, &
                            & reoh_dip, thetae_dip, &
                            & fac_r, fac_theta
  use ps_type_pes,      only: cose, reoh, thetae

  implicit none

! arguments 

  logical, intent(in) :: lttm3

  integer, intent(in) :: imcon
  integer, intent(in) :: listttm2(mxatms)
  integer, intent(in) :: nttm2, n_water

  real(8), intent(in) :: cell(9)
  real(8)             :: chge(mxatms)
  real(8), intent(in) :: gamma
  real(8), intent(in) :: xxx(mxatms), yyy(mxatms), zzz(mxatms)

! local variables 

  integer :: i, j, iox, ih1, ih2
  integer :: mttm2

  real(8) :: cost, sint, damp
  real(8) :: dumx, dumy, dumz
  real(8) :: ddampdr1, ddampdr2
  real(8) :: dpc0dr1, dpc0dr2, dpc0dcos
  real(8) :: factor
  real(8) :: fmat(0:10,3)
  real(8) :: p1, p2, pl1, pl2, pc0
  real(8) :: r1x, r1y, r1z, r2x, r2y, r2z, r1, r2
  real(8) :: sumchge
  real(8) :: x1, x2, x3
  real(8) :: xh1, yh1, zh1, xh2, yh2, zh2
  real(8) :: theta, tmp_x, tmp_y, tmp_z, tmp

!
! --- molecular dipole ---
!
  
  mttm2 = listttm2(nttm2) 

  do i=1,n_water

     mttm2 = mttm2 + 1

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
     
     moldipx(i) = chge(mttm2) * xxx(mttm2) + chge(ih1) * xh1 + chge(ih2) * xh2
     moldipy(i) = chge(mttm2) * yyy(mttm2) + chge(ih1) * yh1 + chge(ih2) * yh2
     moldipz(i) = chge(mttm2) * zzz(mttm2) + chge(ih1) * zh1 + chge(ih2) * zh2

  end do

!
! --- DMS if specified ---
!

  if( .not. ldms ) return

!
! --- initialize ESP ---
!
  
  vesp(:) = 0.0d0
  vesp_k(:) = 0.0d0;  vesp_r(:) = 0.0d0
  vesp_s(:) = 0.0d0;  vesp_c(:) = 0.0d0

  vesp_dc_k(:) = 0.0d0;  vesp_dc_r(:) = 0.0d0;  vesp_dc_c(:) = 0.0d0

!
! --- DMS setting ---
!  
! Some of this part were based on the reference code. 
! Modications were made so as to be consistent with the current code.
!

! gradients with respect to intra. coord  
  dp1dr1(:) = 0.0d0;  dp1dr2(:) = 0.0d0;  dp1dcos(:) = 0.0d0
  dp2dr1(:) = 0.0d0;  dp2dr2(:) = 0.0d0;  dp2dcos(:) = 0.0d0 
  
  mttm2 = listttm2(nttm2) 

  ttm2_mol: do i=1,n_water
     
     mttm2 = mttm2 + 1

     iox = listttm2(3*i-2) 
     ih1 = listttm2(3*i-1)
     ih2 = listttm2(3*i)  
     
     r1x = xxx(ih1) - xxx(iox)
     r1y = yyy(ih1) - yyy(iox)
     r1z = zzz(ih1) - zzz(iox)

     call images( imcon, 0, 1, 1, cell, r1x, r1y, r1z )

     r1 = sqrt( r1x**2 + r1y**2 + r1z**2 )

     r2x = xxx(ih2) - xxx(iox)
     r2y = yyy(ih2) - yyy(iox)
     r2z = zzz(ih2) - zzz(iox)
     
     call images( imcon, 0, 1, 1, cell, r2x, r2y, r2z )

     r2 = sqrt( r2x**2 + r2y**2 + r2z**2 )  ! angstrom
     
     cost = ( r1x * r2x + r1y * r2y + r1z * r2z ) / ( r1 * r2 )  
 
     x1 = ( r1 - reoh ) / reoh  ! ang/ang = dimensionless
     x2 = ( r2 - reoh ) / reoh
     x3 = cost - cose
!    theta = acos( cost )

     tmp_x = r1y*r2z - r1z*r2y
     tmp_y =-r1x*r2z + r1z*r2x
     tmp_z = r1x*r2y - r1y*r2x

     tmp = sqrt( tmp_x**2 + tmp_y**2 + tmp_z**2 )
     sint = tmp / (r1*r2)
     theta = atan2(sint,cost)

     fmat(:,:) = 0.0d0
     fmat(1,1) = 1.0d0
     fmat(1,2) = 1.0d0
     fmat(1,3) = 1.0d0

     do j=2,isump
     
        fmat(j,1) = fmat(j-1,1) * x1
        fmat(j,2) = fmat(j-1,2) * x2
        fmat(j,3)=  fmat(j-1,3) * x3

     end do

! q(r1,r2,theta), q(r2,r1,theta)

     p1 = 0.0d0
     p2 = 0.0d0

     do j=2,84

        p1 = p1 + coef(j) * fmat(idx(j,1),1) * fmat(idx(j,2),2) * fmat(idx(j,3),3)
        p2 = p2 + coef(j) * fmat(idx(j,2),1) * fmat(idx(j,1),2) * fmat(idx(j,3),3)

        dp1dr1(i) = dp1dr1(i) + coef(j) * ( idx(j,1) - 1 ) * &
             & fmat(idx(j,1)-1,1) * fmat(idx(j,2),2) * fmat(idx(j,3),3) / reoh
        dp1dr2(i) = dp1dr2(i) + coef(j) * ( idx(j,2) - 1 ) * &
             & fmat(idx(j,1),1) * fmat(idx(j,2)-1,2) * fmat(idx(j,3),3) / reoh

        dp2dr1(i) = dp2dr1(i) + coef(j) * ( idx(j,2) - 1 ) * &
             & fmat(idx(j,2)-1,1) * fmat(idx(j,1),2) * fmat(idx(j,3),3) / reoh
        dp2dr2(i) = dp2dr2(i) + coef(j) * ( idx(j,1) - 1 ) * &
             & fmat(idx(j,2),1) * fmat(idx(j,1)-1,2) * fmat(idx(j,3),3) / reoh

        dp1dcos(i) = dp1dcos(i) + coef(j) * ( idx(j,3) - 1 ) * &
             & fmat(idx(j,1),1) * fmat(idx(j,2),2) * fmat(idx(j,3)-1,3)
        dp2dcos(i) = dp2dcos(i) + coef(j) * ( idx(j,3) - 1 ) * &
             & fmat(idx(j,2),1) * fmat(idx(j,1),2) * fmat(idx(j,3)-1,3)
        
     end do

     pl1 = cost
     pl2 = 0.5d0 * ( 3.0d0 * pl1 * pl1 - 1.0d0 )

     damp = exp( -b1dms * ( ( r1 - reoh )**2 + ( r2 - reoh )**2 ) ) 

     ddampdr1 = -2.0d0 * b1dms * ( r1 - reoh ) * damp
     ddampdr2 = -2.0d0 * b1dms * ( r2 - reoh ) * damp

     pc0 = adms * ( ( r1**bdms ) + ( r2**bdms ) ) * ( c0 + pl1 * c1 + pl2 * c2 )

     dpc0dr1  = adms * ( bdms * r1**(bdms-1) ) * ( c0 + pl1 * c1 + pl2 * c2 )
     dpc0dr2  = adms * ( bdms * r2**(bdms-1) ) * ( c0 + pl1 * c1 + pl2 * c2 )
     dpc0dcos = adms * ( ( r1**bdms ) + ( r2**bdms ) ) * ( c1 + 3.0d0 * c2 * pl1 )

     dp1dr1(i)  = dp1dr1(i) * damp + p1 * ddampdr1 + dpc0dr1
     dp1dr2(i)  = dp1dr2(i) * damp + p1 * ddampdr2 + dpc0dr2

     dp2dr1(i)  = dp2dr1(i) * damp + p2 * ddampdr1 + dpc0dr1
     dp2dr2(i)  = dp2dr2(i) * damp + p2 * ddampdr2 + dpc0dr2

     dp1dcos(i) = dp1dcos(i) * damp + dpc0dcos
     dp2dcos(i) = dp2dcos(i) * damp + dpc0dcos

     p1 = coef(1) + p1 * damp + pc0  ! q(r1,r2,theta)
     p2 = coef(1) + p2 * damp + pc0  ! q(r2,r1,theta)

! For TTM3.
     if ( lttm3 ) then
        p1 = p1 + fac_r * ( r1-reoh_dip ) + fac_theta * ( theta-thetae_dip )
        p2 = p2 + fac_r * ( r2-reoh_dip ) + fac_theta * ( theta-thetae_dip )

        dp1dr1(i)  = dp1dr1(i) + fac_r
        dp2dr2(i)  = dp2dr2(i) + fac_r

        dp1dcos(i) = dp1dcos(i) - fac_theta / sint
        dp2dcos(i) = dp2dcos(i) - fac_theta / sint
     endif

     if( lttmfold ) then  ! original assign of charges

        chge(ih1)   = p1 / ( 1.0d0 - gamma )
        chge(ih2)   = p2 / ( 1.0d0 - gamma ) 
        chge(mttm2) = -( chge(ih1) + chge(ih2) )  ! M site

     else    ! new assign

        factor = gamma / 2.0d0 / ( 1.0d0 - gamma )

        chge(ih1)   = p1 + factor * ( p1 + p2 )
        chge(ih2)   = p2 + factor * ( p1 + p2 )
        chge(mttm2) = -( chge(ih1) + chge(ih2) ) ! M site

     end if

! PS dipole for check

     chkdipx(i) = p1 * r1x + p2 * r2x
     chkdipy(i) = p1 * r1y + p2 * r2y
     chkdipz(i) = p1 * r1z + p2 * r2z

  end do ttm2_mol

!
! --- new molecular dipole ---
!
  
  mttm2 = listttm2(nttm2) 

  do i=1,n_water

     mttm2 = mttm2 + 1

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
     
     moldipx(i) = chge(mttm2) * xxx(mttm2) + chge(ih1) * xh1 + chge(ih2) * xh2
     moldipy(i) = chge(mttm2) * yyy(mttm2) + chge(ih1) * yh1 + chge(ih2) * yh2
     moldipz(i) = chge(mttm2) * zzz(mttm2) + chge(ih1) * zh1 + chge(ih2) * zh2

     sumchge = chge(mttm2) + chge(ih1) + chge(ih2) + chge(iox)

     if( abs(sumchge) > 1.0d-10 ) then 

        write(6,*) 'WARN: sum of new charges is non-zero';  stop

     end if

  end do

  return

end subroutine qdetermine
