!
! === driver routine for partial charge derivative (with DMS) ===
!
! Last updated: 24 July 2006
!
! Ref. 1: C. J. Burnham and S. S. Xantheas, JCP 116, 5115 (2002)
! Ref. 2: G. S. Fanourgakis and S. S. Xantheas, JPCA 110, 4100 (2006)
! Ref. 3: H. Partridge and D. W. Schwenke, JCP 106, 4618 (1997)
!
! See README file for the details.
!

subroutine qdforce( listttm2, n_water, nttm2, natms, imcon, cell,&
     & gamma, xxx, yyy, zzz, fxx, fyy, fzz, stress, virqdf )

  use global_variables, only: mxatms
  use ps_type_dms,      only: lttmfold, vesp, &
                             dp1dr1, dp1dr2, dp2dr1, dp2dr2, dp1dcos, dp2dcos
  use multibead, only: bead_rank, bead_size
#ifdef HEAT_CURRENT
  use heatcurrent, only: update_stress_qd, update_forces
#endif

  implicit none

! arguments

  integer, intent(in) :: imcon
  integer, intent(in) :: listttm2(mxatms)
  integer, intent(in) :: nttm2, n_water
  integer, intent(in) :: natms

  real(8), intent(in) :: cell(9)
  real(8), intent(in) :: gamma
  real(8), intent(in) :: xxx(mxatms), yyy(mxatms), zzz(mxatms)

  real(8) :: fxx(mxatms), fyy(mxatms), fzz(mxatms)
  real(8) :: stress(9), virqdf

! local variables

  integer, parameter :: n_site = 4  ! 1: O, 2: H1, 3: H2, 4: M

  integer :: i, j, k, iox, ih1, ih2, imol
  integer :: mttm2
  integer :: label(n_site,n_site)

  real(8) :: cost, cott
  real(8) :: dcosdr12(n_water), dcosdr13(n_water), dcosdr23(n_water)
  real(8) :: dcostmp1, dcostmp2, dcostmp3
  real(8) :: dcostmp5, dcostmp6, dcostmp9
  real(8) :: dp1dsitex(n_site,n_water), dp2dsitex(n_site,n_water)
  real(8) :: dp1dsitey(n_site,n_water), dp2dsitey(n_site,n_water)
  real(8) :: dp1dsitez(n_site,n_water), dp2dsitez(n_site,n_water)
  real(8) :: dp1dr12, dp2dr12, dp1dr13, dp2dr13
  real(8) :: dqdx(n_site,n_site,n_water)
  real(8) :: dqdy(n_site,n_site,n_water)
  real(8) :: dqdz(n_site,n_site,n_water)
  real(8) :: dqdr12(n_site,n_water), dqdr13(n_site,n_water)
  real(8) :: dqdcos(n_site,n_water)
  real(8) :: dqdr_tmp
  real(8) :: drintx(3,4), drinty(3,4), drintz(3,4)
  real(8) :: factor
  real(8) :: r1x, r1y, r1z, r2x, r2y, r2z, rx, ry, rz, rjk
  real(8) :: rx1, ry1, rz1, rx2, ry2, rz2, rx3, ry3, rz3
  real(8) :: r1, r2, r3, rdot, r1sqr, r2sqr
  real(8) :: sint
  real(8) :: strs1, strs2, strs3, strs5, strs6, strs9
  real(8) :: tmp
  real(8) :: vtmp
#ifdef HEAT_CURRENT
  real(8), parameter :: third=1.d0/3.d0
  real(8) :: force_tmp(3)
  integer :: idx
#endif /* HEAT_CURRENT */

!
! --- initialization ---
!

  virqdf = 0.0d0

  label(:,:) = 0

  do j=1,n_site
     do k=1,n_site

        if( j == 1 .and. k == 2 ) label(j,k) = 1
        if( j == 1 .and. k == 3 ) label(j,k) = 2

     end do
  end do

!
! --- charge derivatives with respect to intra. coordinates ---
!

  dqdr12(:,:) = 0.0d0;  dqdr13(:,:) = 0.0d0
  dqdcos(:,:) = 0.0d0
  drintx(:,:) = 0.0d0;  drinty(:,:) = 0.0d0;  drintz(:,:) = 0.0d0

  mttm2 = listttm2(nttm2) + 1 + bead_rank

  ttm2_mol: do imol=bead_rank+1,n_water,bead_size

     iox = listttm2(3*imol-2)
     ih1 = listttm2(3*imol-1)
     ih2 = listttm2(3*imol)

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

     rdot = r1x * r2x + r1y * r2y + r1z * r2z
     cost = rdot / ( r1 * r2 )

     r3 = sqrt( r1**2 + r2**2 - 2.0d0 * r1 * r2 * cost )

     sint = sqrt( 1.0d0 - cost * cost )
     cott = cost / sint

     r1sqr = r1 * r1
     r2sqr = r2 * r2

     drintx(1,2) = r1x / r1
     drinty(1,2) = r1y / r1
     drintz(1,2) = r1z / r1

     drintx(1,1) = -drintx(1,2)
     drinty(1,1) = -drinty(1,2)
     drintz(1,1) = -drintz(1,2)

     drintx(2,3) = r2x / r2
     drinty(2,3) = r2y / r2
     drintz(2,3) = r2z / r2

     drintx(2,1) = -drintx(2,3)
     drinty(2,1) = -drinty(2,3)
     drintz(2,1) = -drintz(2,3)

     drintx(3,2) = cott * ( r1x / r1sqr - r2x / rdot )
     drinty(3,2) = cott * ( r1y / r1sqr - r2y / rdot )
     drintz(3,2) = cott * ( r1z / r1sqr - r2z / rdot )

     drintx(3,3) = cott * ( r2x / r2sqr - r1x / rdot )
     drinty(3,3) = cott * ( r2y / r2sqr - r1y / rdot )
     drintz(3,3) = cott * ( r2z / r2sqr - r1z / rdot )

     drintx(3,1) = -drintx(3,2) - drintx(3,3)
     drinty(3,1) = -drinty(3,2) - drinty(3,3)
     drintz(3,1) = -drintz(3,2) - drintz(3,3)

     do i=1,n_site

        dp1dsitex(i,imol) = dp1dr1(imol) * drintx(1,i) + dp1dr2(imol) * drintx(2,i) &
                        & - sint * dp1dcos(imol) * drintx(3,i)
        dp1dsitey(i,imol) = dp1dr1(imol) * drinty(1,i) + dp1dr2(imol) * drinty(2,i) &
                        & - sint * dp1dcos(imol) * drinty(3,i)
        dp1dsitez(i,imol) = dp1dr1(imol) * drintz(1,i) + dp1dr2(imol) * drintz(2,i) &
                        & - sint * dp1dcos(imol) * drintz(3,i)

        dp2dsitex(i,imol) = dp2dr1(imol) * drintx(1,i) + dp2dr2(imol) * drintx(2,i) &
                        & - sint * dp2dcos(imol) * drintx(3,i)
        dp2dsitey(i,imol) = dp2dr1(imol) * drinty(1,i) + dp2dr2(imol) * drinty(2,i) &
                        & - sint * dp2dcos(imol) * drinty(3,i)
        dp2dsitez(i,imol) = dp2dr1(imol) * drintz(1,i) + dp2dr2(imol) * drintz(2,i) &
                        & - sint * dp2dcos(imol) * drintz(3,i)

     end do

     dp1dr12 = dp1dr1(imol)
     dp2dr12 = dp2dr1(imol)

     dp1dr13 = dp1dr2(imol)
     dp2dr13 = dp2dr2(imol)

     dcosdr12(imol) = ( r1 - r2 * cost ) / r1 / r2
     dcosdr13(imol) = ( r2 - r1 * cost ) / r1 / r2
     dcosdr23(imol) = -r3 / r1 / r2

     if( lttmfold ) then  ! original assignment

        dqdr12(2,imol) =  dp1dr12 / ( 1.0d0 - gamma )
        dqdr13(2,imol) =  dp1dr13 / ( 1.0d0 - gamma )
        dqdcos(2,imol) =  dp1dcos(imol) / ( 1.0d0 - gamma )

        dqdr12(3,imol) =  dp2dr12 / ( 1.0d0 - gamma )
        dqdr13(3,imol) =  dp2dr13 / ( 1.0d0 - gamma )
        dqdcos(3,imol) =  dp2dcos(imol) / ( 1.0d0 - gamma )

        dqdr12(4,imol) = -( dp1dr12 + dp2dr12 ) / ( 1.0d0 - gamma )
        dqdr13(4,imol) = -( dp1dr13 + dp2dr13 ) / ( 1.0d0 - gamma )
        dqdcos(4,imol) = -( dp1dcos(imol) + dp2dcos(imol) ) / ( 1.0d0 - gamma )

     else  ! new assignment

        tmp = gamma / 2.0d0 / ( 1.0d0 - gamma )

        dqdr12(2,imol) = dp1dr12 + ( dp1dr12 + dp2dr12 ) * tmp
        dqdr13(2,imol) = dp1dr13 + ( dp1dr13 + dp2dr13 ) * tmp
        dqdcos(2,imol) = dp1dcos(imol) + ( dp1dcos(imol) + dp2dcos(imol) ) * tmp

        dqdr12(3,imol) = dp2dr12 + ( dp1dr12 + dp2dr12 ) * tmp
        dqdr13(3,imol) = dp2dr13 + ( dp1dr13 + dp2dr13 ) * tmp
        dqdcos(3,imol) = dp2dcos(imol) + ( dp1dcos(imol) + dp2dcos(imol) ) * tmp

        dqdr12(4,imol) = -( dp1dr12 + dp2dr12 ) / ( 1.0d0 - gamma )
        dqdr13(4,imol) = -( dp1dr13 + dp2dr13 ) / ( 1.0d0 - gamma )
        dqdcos(4,imol) = -( dp1dcos(imol) + dp2dcos(imol) ) / ( 1.0d0 - gamma )

     end if

     mttm2 = mttm2 + bead_size

  end do ttm2_mol

!
! --- charge derivatives with respect to site coordinates ---
!

  dqdx(:,:,:) = 0.0d0;  dqdy(:,:,:) = 0.0d0;  dqdz(:,:,:) = 0.0d0

  ttm2_mol2: do imol=bead_rank+1,n_water,bead_size

     do j=1,n_site

        if( lttmfold ) then  ! original assign of charges

! with respect to hydrogen 1 charge
           dqdx(2,j,imol) = dp1dsitex(j,imol) / ( 1.0d0 - gamma )
           dqdy(2,j,imol) = dp1dsitey(j,imol) / ( 1.0d0 - gamma )
           dqdz(2,j,imol) = dp1dsitez(j,imol) / ( 1.0d0 - gamma )

! with respect to hydrogen 2 charge
           dqdx(3,j,imol) = dp2dsitex(j,imol) / ( 1.0d0 - gamma )
           dqdy(3,j,imol) = dp2dsitey(j,imol) / ( 1.0d0 - gamma )
           dqdz(3,j,imol) = dp2dsitez(j,imol) / ( 1.0d0 - gamma )

! with respect to M site charge
           dqdx(4,j,imol) = -( dqdx(2,j,imol) + dqdx(3,j,imol) )
           dqdy(4,j,imol) = -( dqdy(2,j,imol) + dqdy(3,j,imol) )
           dqdz(4,j,imol) = -( dqdz(2,j,imol) + dqdz(3,j,imol) )

        else ! new assign of charges

           factor = gamma / 2.0d0 / ( 1.0d0 - gamma )

! with respect to hydrogen 1 charge
           dqdx(2,j,imol) = dp1dsitex(j,imol) * ( factor + 1.0d0 ) &
                        & + dp2dsitex(j,imol) * factor
           dqdy(2,j,imol) = dp1dsitey(j,imol) * ( factor + 1.0d0 ) &
                        & + dp2dsitey(j,imol) * factor
           dqdz(2,j,imol) = dp1dsitez(j,imol) * ( factor + 1.0d0 ) &
                        & + dp2dsitez(j,imol) * factor

! with respect to hydrogen 2 charge
           dqdx(3,j,imol) = dp2dsitex(j,imol) * ( factor + 1.0d0 ) &
                        & + dp1dsitex(j,imol) * factor
           dqdy(3,j,imol) = dp2dsitey(j,imol) * ( factor + 1.0d0 ) &
                        & + dp1dsitey(j,imol) * factor
           dqdz(3,j,imol) = dp2dsitez(j,imol) * ( factor + 1.0d0 ) &
                        & + dp1dsitez(j,imol) * factor

! with respect to M site charge
           dqdx(4,j,imol) = -( dqdx(2,j,imol) + dqdx(3,j,imol) )
           dqdy(4,j,imol) = -( dqdy(2,j,imol) + dqdy(3,j,imol) )
           dqdz(4,j,imol) = -( dqdz(2,j,imol) + dqdz(3,j,imol) )

        end if

     end do

  end do ttm2_mol2

!
! --- electrostatic potential --- (see forces.f)
!

!  vesp(:) = 0.0d0
!
!  do i=1,natms
!     vesp(i) = vesp_k(i) + vesp_r(i) + vesp_s(i) + vesp_c(i) &
!           & + vesp_dc_k(i) + vesp_dc_r(i) + vesp_dc_c(i)
!  end do

!
! --- add charge derivative contributions to total forces ---
!

  mttm2 = listttm2(nttm2) + 1 + bead_rank

  ttm2_mol3: do imol=bead_rank+1,n_water,bead_size

     iox = listttm2(3*imol-2)
     ih1 = listttm2(3*imol-1)
     ih2 = listttm2(3*imol)

     do j=1,n_site

        if( j == 1 ) vtmp = vesp(iox)    ! O
        if( j == 2 ) vtmp = vesp(ih1)    ! H1
        if( j == 3 ) vtmp = vesp(ih2)    ! H2
        if( j == 4 ) vtmp = vesp(mttm2)  ! M

        fxx(iox) = fxx(iox) - dqdx(j,1,imol) * vtmp
        fyy(iox) = fyy(iox) - dqdy(j,1,imol) * vtmp
        fzz(iox) = fzz(iox) - dqdz(j,1,imol) * vtmp

        fxx(ih1) = fxx(ih1) - dqdx(j,2,imol) * vtmp
        fyy(ih1) = fyy(ih1) - dqdy(j,2,imol) * vtmp
        fzz(ih1) = fzz(ih1) - dqdz(j,2,imol) * vtmp

        fxx(ih2) = fxx(ih2) - dqdx(j,3,imol) * vtmp
        fyy(ih2) = fyy(ih2) - dqdy(j,3,imol) * vtmp
        fzz(ih2) = fzz(ih2) - dqdz(j,3,imol) * vtmp

        fxx(mttm2) = fxx(mttm2) - dqdx(j,4,imol) * vtmp
        fyy(mttm2) = fyy(mttm2) - dqdy(j,4,imol) * vtmp
        fzz(mttm2) = fzz(mttm2) - dqdz(j,4,imol) * vtmp
#ifdef HEAT_CURRENT
        force_tmp=(/-dqdx(j,1,imol)*vtmp,-dqdy(j,1,imol)*vtmp,-dqdz(j,1,imol)*vtmp/)
        if (j == 4) then
          idx = n_water + imol
        else
          idx = 3*(imol-1) + j
        end if
        call update_forces(iox,idx,force_tmp)
        force_tmp=(/-dqdx(j,2,imol)*vtmp,-dqdy(j,2,imol)*vtmp,-dqdz(j,2,imol)*vtmp/)
        call update_forces(ih1,idx,force_tmp)
        force_tmp=(/-dqdx(j,3,imol)*vtmp,-dqdy(j,3,imol)*vtmp,-dqdz(j,3,imol)*vtmp/)
        call update_forces(ih2,idx,force_tmp)
        force_tmp=(/-dqdx(j,4,imol)*vtmp,-dqdy(j,4,imol)*vtmp,-dqdz(j,4,imol)*vtmp/)
        call update_forces(mttm2,idx,force_tmp)
#endif /*HEAT_CURRENT*/

     end do

     mttm2 = mttm2 + bead_size

  end do ttm2_mol3

!
! --- pressure tensor ---
!

  strs1 = 0.0d0;  strs2 = 0.0d0;  strs3 = 0.0d0
  strs5 = 0.0d0;  strs6 = 0.0d0
  strs9 = 0.0d0

  mttm2 = listttm2(nttm2) + 1 + bead_rank

  ttm2_mol4: do imol=bead_rank+1,n_water,bead_size

     iox = listttm2(3*imol-2)
     ih1 = listttm2(3*imol-1)
     ih2 = listttm2(3*imol)

     do_i: do i=1,n_site

        if( i == 1 ) vtmp = vesp(iox)    ! O
        if( i == 2 ) vtmp = vesp(ih1)    ! H1
        if( i == 3 ) vtmp = vesp(ih2)    ! H2
        if( i == 4 ) vtmp = vesp(mttm2)  ! M

        do_j: do j=1,n_site-1
           do_k: do k=j+1,n_site

              if( label(j,k) == 0 ) cycle

              if( label(j,k) == 1 ) then

                 dqdr_tmp = dqdr12(i,imol)

                 rx = xxx(iox) - xxx(ih1)
                 ry = yyy(iox) - yyy(ih1)
                 rz = zzz(iox) - zzz(ih1)

              else if( label(j,k) == 2 ) then

                 dqdr_tmp = dqdr13(i,imol)

                 rx = xxx(iox) - xxx(ih2)
                 ry = yyy(iox) - yyy(ih2)
                 rz = zzz(iox) - zzz(ih2)

              else

                 write(6,*) 'WRONG LABEL IN QDFORCE; STOP';  stop

              end if

              call images( imcon, 0, 1, 1, cell, rx, ry, rz )

              rjk = sqrt( rx * rx + ry * ry + rz * rz )

              strs1 = strs1 - vtmp * dqdr_tmp * rx * rx / rjk
              strs2 = strs2 - vtmp * dqdr_tmp * rx * ry / rjk
              strs3 = strs3 - vtmp * dqdr_tmp * rx * rz / rjk

              strs5 = strs5 - vtmp * dqdr_tmp * ry * ry / rjk
              strs6 = strs6 - vtmp * dqdr_tmp * ry * rz / rjk

              strs9 = strs9 - vtmp * dqdr_tmp * rz * rz / rjk

#ifdef HEAT_CURRENT
#ifdef HEAT_STRESS
          call update_stress_qd(iox,1,1,third*vtmp*dqdr_tmp*rx*rx/rjk)
          call update_stress_qd(iox,1,2,third*vtmp*dqdr_tmp*rx*ry/rjk)
          call update_stress_qd(iox,1,3,third*vtmp*dqdr_tmp*rx*rz/rjk)
          call update_stress_qd(iox,2,1,third*vtmp*dqdr_tmp*rx*ry/rjk)
          call update_stress_qd(iox,2,2,third*vtmp*dqdr_tmp*ry*ry/rjk)
          call update_stress_qd(iox,2,3,third*vtmp*dqdr_tmp*ry*rz/rjk)
          call update_stress_qd(iox,3,1,third*vtmp*dqdr_tmp*rx*rz/rjk)
          call update_stress_qd(iox,3,2,third*vtmp*dqdr_tmp*ry*rz/rjk)
          call update_stress_qd(iox,3,3,third*vtmp*dqdr_tmp*rz*rz/rjk)
          call update_stress_qd(ih1,1,1,third*vtmp*dqdr_tmp*rx*rx/rjk)
          call update_stress_qd(ih1,1,2,third*vtmp*dqdr_tmp*rx*ry/rjk)
          call update_stress_qd(ih1,1,3,third*vtmp*dqdr_tmp*rx*rz/rjk)
          call update_stress_qd(ih1,2,1,third*vtmp*dqdr_tmp*rx*ry/rjk)
          call update_stress_qd(ih1,2,2,third*vtmp*dqdr_tmp*ry*ry/rjk)
          call update_stress_qd(ih1,2,3,third*vtmp*dqdr_tmp*ry*rz/rjk)
          call update_stress_qd(ih1,3,1,third*vtmp*dqdr_tmp*rx*rz/rjk)
          call update_stress_qd(ih1,3,2,third*vtmp*dqdr_tmp*ry*rz/rjk)
          call update_stress_qd(ih1,3,3,third*vtmp*dqdr_tmp*rz*rz/rjk)
          call update_stress_qd(ih2,1,1,third*vtmp*dqdr_tmp*rx*rx/rjk)
          call update_stress_qd(ih2,1,2,third*vtmp*dqdr_tmp*rx*ry/rjk)
          call update_stress_qd(ih2,1,3,third*vtmp*dqdr_tmp*rx*rz/rjk)
          call update_stress_qd(ih2,2,1,third*vtmp*dqdr_tmp*rx*ry/rjk)
          call update_stress_qd(ih2,2,2,third*vtmp*dqdr_tmp*ry*ry/rjk)
          call update_stress_qd(ih2,2,3,third*vtmp*dqdr_tmp*ry*rz/rjk)
          call update_stress_qd(ih2,3,1,third*vtmp*dqdr_tmp*rx*rz/rjk)
          call update_stress_qd(ih2,3,2,third*vtmp*dqdr_tmp*ry*rz/rjk)
          call update_stress_qd(ih2,3,3,third*vtmp*dqdr_tmp*rz*rz/rjk)
#endif
#endif /* HEAT_CURRENT */

              virqdf = virqdf + vtmp * rjk * dqdr_tmp

           end do do_k
        end do do_j

        rx1 = xxx(iox) - xxx(ih1)
        ry1 = yyy(iox) - yyy(ih1)
        rz1 = zzz(iox) - zzz(ih1)

        rx2 = xxx(iox) - xxx(ih2)
        ry2 = yyy(iox) - yyy(ih2)
        rz2 = zzz(iox) - zzz(ih2)

        rx3 = xxx(ih1) - xxx(ih2)
        ry3 = yyy(ih1) - yyy(ih2)
        rz3 = zzz(ih1) - zzz(ih2)

        call images( imcon, 0, 1, 1, cell, rx1, ry1, rz1 )
        call images( imcon, 0, 1, 1, cell, rx2, ry2, rz2 )
        call images( imcon, 0, 1, 1, cell, rx3, ry3, rz3 )

        r1 = sqrt( rx1 * rx1 + ry1 * ry1 + rz1 * rz1 )
        r2 = sqrt( rx2 * rx2 + ry2 * ry2 + rz2 * rz2 )
        r3 = sqrt( rx3 * rx3 + ry3 * ry3 + rz3 * rz3 )

        dcostmp1 = dcosdr12(imol) * rx1 * rx1 / r1 &
               & + dcosdr13(imol) * rx2 * rx2 / r2 &
               & + dcosdr23(imol) * rx3 * rx3 / r3
        dcostmp2 = dcosdr12(imol) * rx1 * ry1 / r1 &
               & + dcosdr13(imol) * rx2 * ry2 / r2 &
               & + dcosdr23(imol) * rx3 * ry3 / r3
        dcostmp3 = dcosdr12(imol) * rx1 * rz1 / r1 &
               & + dcosdr13(imol) * rx2 * rz2 / r2 &
               & + dcosdr23(imol) * rx3 * rz3 / r3

        dcostmp5 = dcosdr12(imol) * ry1 * ry1 / r1 &
               & + dcosdr13(imol) * ry2 * ry2 / r2 &
               & + dcosdr23(imol) * ry3 * ry3 / r3
        dcostmp6 = dcosdr12(imol) * ry1 * rz1 / r1 &
               & + dcosdr13(imol) * ry2 * rz2 / r2 &
               & + dcosdr23(imol) * ry3 * rz3 / r3

        dcostmp9 = dcosdr12(imol) * rz1 * rz1 / r1 &
               & + dcosdr13(imol) * rz2 * rz2 / r2 &
               & + dcosdr23(imol) * rz3 * rz3 / r3

        strs1 = strs1 - vtmp * dqdcos(i,imol) * dcostmp1
        strs2 = strs2 - vtmp * dqdcos(i,imol) * dcostmp2
        strs3 = strs3 - vtmp * dqdcos(i,imol) * dcostmp3

        strs5 = strs5 - vtmp * dqdcos(i,imol) * dcostmp5
        strs6 = strs6 - vtmp * dqdcos(i,imol) * dcostmp6

        strs9 = strs9 - vtmp * dqdcos(i,imol) * dcostmp9

#ifdef HEAT_CURRENT
    call update_stress_qd(iox,1,1,third*vtmp*dqdcos(i,imol)*dcostmp1)
    call update_stress_qd(iox,1,2,third*vtmp*dqdcos(i,imol)*dcostmp2)
    call update_stress_qd(iox,1,3,third*vtmp*dqdcos(i,imol)*dcostmp3)
    call update_stress_qd(iox,2,1,third*vtmp*dqdcos(i,imol)*dcostmp2)
    call update_stress_qd(iox,2,2,third*vtmp*dqdcos(i,imol)*dcostmp5)
    call update_stress_qd(iox,2,3,third*vtmp*dqdcos(i,imol)*dcostmp6)
    call update_stress_qd(iox,3,1,third*vtmp*dqdcos(i,imol)*dcostmp3)
    call update_stress_qd(iox,3,2,third*vtmp*dqdcos(i,imol)*dcostmp6)
    call update_stress_qd(iox,3,3,third*vtmp*dqdcos(i,imol)*dcostmp9)
    call update_stress_qd(ih1,1,1,third*vtmp*dqdcos(i,imol)*dcostmp1)
    call update_stress_qd(ih1,1,2,third*vtmp*dqdcos(i,imol)*dcostmp2)
    call update_stress_qd(ih1,1,3,third*vtmp*dqdcos(i,imol)*dcostmp3)
    call update_stress_qd(ih1,2,1,third*vtmp*dqdcos(i,imol)*dcostmp2)
    call update_stress_qd(ih1,2,2,third*vtmp*dqdcos(i,imol)*dcostmp5)
    call update_stress_qd(ih1,2,3,third*vtmp*dqdcos(i,imol)*dcostmp6)
    call update_stress_qd(ih1,3,1,third*vtmp*dqdcos(i,imol)*dcostmp3)
    call update_stress_qd(ih1,3,2,third*vtmp*dqdcos(i,imol)*dcostmp6)
    call update_stress_qd(ih1,3,3,third*vtmp*dqdcos(i,imol)*dcostmp9)
    call update_stress_qd(ih2,1,1,third*vtmp*dqdcos(i,imol)*dcostmp1)
    call update_stress_qd(ih2,1,2,third*vtmp*dqdcos(i,imol)*dcostmp2)
    call update_stress_qd(ih2,1,3,third*vtmp*dqdcos(i,imol)*dcostmp3)
    call update_stress_qd(ih2,2,1,third*vtmp*dqdcos(i,imol)*dcostmp2)
    call update_stress_qd(ih2,2,2,third*vtmp*dqdcos(i,imol)*dcostmp5)
    call update_stress_qd(ih2,2,3,third*vtmp*dqdcos(i,imol)*dcostmp6)
    call update_stress_qd(ih2,3,1,third*vtmp*dqdcos(i,imol)*dcostmp3)
    call update_stress_qd(ih2,3,2,third*vtmp*dqdcos(i,imol)*dcostmp6)
    call update_stress_qd(ih2,3,3,third*vtmp*dqdcos(i,imol)*dcostmp9)
#endif /* HEAT_CURRENT */

        virqdf = virqdf + &
             & vtmp * ( dcostmp1 + dcostmp5 + dcostmp9 ) * dqdcos(i,imol)

     end do do_i

     mttm2 = mttm2 + bead_size

  end do ttm2_mol4

!  write(6,*) 'CHECK'
!  write(6,'(3f20.10)') strs1, strs2, strs3
!  write(6,'(3f20.10)') strs2, strs5, strs6
!  write(6,'(3f20.10)') strs3, strs6, strs9
!  write(6,*) 'CHECK2'
!  write(6,'(3f20.10)') stress(1), stress(2), stress(3)
!  write(6,'(3f20.10)') stress(2), stress(5), stress(6)
!  write(6,'(3f20.10)') stress(3), stress(6), stress(9)
!
!  write(6,*) 'CHECK', virqdf, ( strs1 + strs5 + strs9 ) / 3.0d0

  stress(1) = stress(1) + strs1
  stress(2) = stress(2) + strs2
  stress(3) = stress(3) + strs3
  stress(4) = stress(4) + strs2
  stress(5) = stress(5) + strs5
  stress(6) = stress(6) + strs6
  stress(7) = stress(7) + strs3
  stress(8) = stress(8) + strs6
  stress(9) = stress(9) + strs9

!VB:  reduced in forces.f
!  if (bead_size.gt.1) call gdsum(virqdf,1,vtmp)

  return

end subroutine qdforce
