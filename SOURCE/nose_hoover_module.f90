#include "assert.h"

!------------------------------------------------------------------------------
! nose_hoover.f90
!------------------------------------------------------------------------------
!
! integrator of Nose-Hoover chain thermostat based on a higher-order
! decomposition of the Liouville propagator
!
! Ref: Martyna, Klein, & Tuckerman, J. Chem. Phys. 97, 2635 (1992).
!      Martyna, Tuckerman, Tobias, & Klein, Mol. Phys. 87, 1117 (1996).
!      Tuckerman, Liu, Ciccotti, & Martyna, J. Chem. Phys. 115, 1678 (2001).
!
! Last updated: August 25, 2005
!
!------------------------------------------------------------------------------



!==============================================================================
module nose_hoover_module
!==============================================================================

  implicit none

  save

  character(*), parameter :: module_name = "nose_hoover_module:"

  logical :: module_init = .false.

  integer, parameter :: dp = 8   !! double precision

  !..................................................

  integer, parameter :: M = 9   !! chain length

  integer, parameter :: num_split  = 1   !! number of internal splitting
  integer, parameter :: num_suzuki = 1   !! number of Yoshida-Suzuki splitting [1,3,5]

  real(dp) :: coeff_suzuki(5,5)   !! coefficients of Yoshida-Suzuki splitting

  !..................................................

  type SystemCoordinate_type
     real(dp) :: mass
     real(dp), pointer :: vel
  endtype

  type Thermostat_type
     real(dp) :: eta(M)        !! position
     real(dp) :: v(M)          !! velocity
     real(dp) :: a(M)          !! acceleration (without friction term)
     real(dp) :: Q(M)          !! mass
     real(dp) :: Q_inv(M)      !! inverse mass
     real(dp) :: kT, Ndof_kT   !! kB * temperature
     type( SystemCoordinate_type ), pointer :: system_coords(:)
     integer :: num_system_coords
     integer :: system_coord_id = 0
     logical :: activated   !! flag to activate the thermostat
     logical :: initialized = .false.
  endtype

#ifdef AIX
  type(Thermostat_type), save :: thermo_lnv
#else
  type(Thermostat_type) :: thermo_lnv
#endif
  real(8) :: c2_lnv, mass_lnv, v_lnv, f_lnv, x_lnv, x_lnv_old

  !..................................................

  integer :: print_level = 0   !! output is suppressed
!!  integer :: print_level = 1

  !..................................................

!! Random number generator from Numerical Recipes.
  integer :: random_number_seed = 100
  integer :: idum_ran2

  !..................................................

contains


!------------------------------------------------------------------------------
subroutine nose_hoover_module_init
!------------------------------------------------------------------------------

  implicit none

  real(dp) :: tmp

  !..................................................

  if ( module_init ) return

  if ( print_level > 0 ) then
     write(6,*)
     write(6,*) "Initializing ", module_name
  endif

  module_init = .true.

!! Yoshida-Suzuki splitting.
  coeff_suzuki( 1, 1 ) = 1.0d0

  tmp = 1.0d0 / ( 2.0d0 - 2.0d0**( 1.0d0 / 3.0d0 ) )
  coeff_suzuki( 1, 3 ) = tmp
  coeff_suzuki( 2, 3 ) = 1.0d0 - 2.0d0 * tmp
  coeff_suzuki( 3, 3 ) = tmp

  tmp = 1.0d0 / ( 4.0d0 - 4.0d0**( 1.0d0 / 3.0d0 ) )
  coeff_suzuki( 1, 5 ) = tmp
  coeff_suzuki( 2, 5 ) = tmp
  coeff_suzuki( 3, 5 ) = 1.0d0 - 4.0d0 * tmp
  coeff_suzuki( 4, 5 ) = tmp
  coeff_suzuki( 5, 5 ) = tmp

!! Random number generator.
  idum_ran2 = - abs( random_number_seed )
  tmp = NR_ran2( idum_ran2 )

!! Output.
  if ( print_level > 0 ) then
     write(6,*)
     write(6,*) module_name
     write(6,*)
     write(6,*) "   max. chain length = ", M
     write(6,*) "   num_split    = ", num_split
     write(6,*) "   num_suzuki   = ", num_suzuki
  endif

end subroutine


!------------------------------------------------------------------------------
subroutine Thermostat_init  &
           (nchain,thermo,num_system_coords,num_constraints, &
            kT,tau,pos_init,vel_init,activate)

! initializes a Thermostat object.
!------------------------------------------------------------------------------

  implicit none

  type( Thermostat_type ) :: thermo

  integer, intent(in) :: nchain  !! number of oscillators in each chain
  integer, intent(in) :: num_system_coords   !! number of system coordinates
                                             !! to be thermostatted
  integer, intent(in) :: num_constraints   !! number of geometrical constraints
                                           !! in the system
  real(dp), intent(in) :: kT   !! external temperature
  real(dp), intent(in) :: tau   !! characteristic timescale of the system
  real(dp), intent(in), optional :: pos_init, vel_init   !! initial condition
                                                         !! of the thermostat
  logical, intent(in), optional :: activate   !! flag to activate the thermostat

  !..................................................

  integer :: Ndof, j

  !..................................................

  if ( .not.module_init ) call nose_hoover_module_init

  if ( nchain > M-1 ) then
     write(6,*) "Stop. Number of oscillators too large. nchain must be <= ", M-1
     stop
  end if
 
  thermo%initialized = .true.

  allocate( thermo%system_coords(num_system_coords) )

  thermo%num_system_coords = num_system_coords

  Ndof = num_system_coords - num_constraints

  thermo%Q(1) = kT * tau**2 * dble( Ndof )
  thermo%Q(2:nchain) = kT * tau**2

  thermo%Q_inv(1:nchain) = 1.0d0 / thermo%Q(1:nchain)

  thermo%Q(nchain+1:M) = 0.d0
  thermo%Q_inv(nchain+1:M) = 0.d0

  thermo%kT = kT
  thermo%Ndof_kT = dble( Ndof ) * kT

!! Initial position, velocity, and acceleration.
  thermo%eta(:) = 0.0d0

  do j = 1, nchain
     thermo%v(j) = sqrt( kT / thermo%Q(j) ) * gaussian()
  enddo

  if ( present( pos_init ) ) thermo%eta(:) = pos_init
  if ( present( vel_init ) ) thermo%v(:) = vel_init

  do j = 2, nchain
     thermo%a(j) = ( thermo%Q(j-1) * thermo%v(j-1)**2 - thermo%kT )  &
                     / thermo%Q(j)
     !! Note: thermo%a(1) is computed in (subr.)Thermostat_integrate.
  enddo

!! Flag to activate the thermostat.
  if ( present( activate ) ) then
     thermo%activated = activate
  else
     thermo%activated = .true.
  endif

end subroutine


!------------------------------------------------------------------------------
subroutine Thermostat_link  &
           (thermo,system_mass,system_velocity)

! establishes the link to system velocicities via pointer assignment.
!------------------------------------------------------------------------------

  implicit none

  type( Thermostat_type ) :: thermo

  real(dp), intent(in) :: system_mass

  real(dp), target :: system_velocity

  !..................................................

  integer :: id

  !..................................................

  if ( .not. thermo%initialized ) then
     write(6,*) module_name, " thermostat is not initialized. stop."
     stop
  endif

  thermo%system_coord_id = thermo%system_coord_id + 1

  id = thermo%system_coord_id

  if ( id <= thermo%num_system_coords ) then
     thermo%system_coords(id)%mass = system_mass
     thermo%system_coords(id)%vel => system_velocity   !! establish link
  else
     write(6,*) module_name, " thermo% num_system_coords is too small. stop."
     stop
  endif

end subroutine


!------------------------------------------------------------------------------
subroutine Thermostat_switch(thermo,activate)
!
! Switch of thermostats
!------------------------------------------------------------------------------

  implicit none

  type( Thermostat_type ) :: thermo

  logical, intent(in) :: activate

  !..................................................

  thermo%activated = activate

end subroutine


!------------------------------------------------------------------------------
subroutine Thermostat_reinit(thermo)
!
! Switch of thermostats
!------------------------------------------------------------------------------

  implicit none

  type( Thermostat_type ) :: thermo

  !..................................................

  thermo% system_coord_id = 0
  
end subroutine

! applies exp(iL_{T-baro}*dt)
subroutine Barostat_thermostat_integrate(nchain,timestep)

   implicit none

   integer, intent(in) :: nchain
   real(8), intent(in) :: timestep

   integer :: j, isplit, isuzuki
   real(8) :: dt, eps, eps2, eps4, kT, tmp

   dt = timestep/dble( num_split )

   do isplit  = 1, num_split
   do isuzuki = 1, num_suzuki

   eps = dt*coeff_suzuki(isuzuki, num_suzuki)

   eps2 = eps/2
   eps4 = eps/4

   kT = thermo_lnv%kT

   !! update force acting on the first barostat's thermostat oscillator
   thermo_lnv%a(1) = thermo_lnv%Q_inv(1)*(mass_lnv*v_lnv*v_lnv - kT)

   !! update thermostat's velocities
   thermo_lnv%v(nchain) = thermo_lnv%v(nchain) + eps2*thermo_lnv%a(nchain)
   do j = nchain - 1, 1, -1
      tmp = exp(-eps4*thermo_lnv%v(j+1))
      thermo_lnv%v(j) = tmp*(tmp*thermo_lnv%v(j) + eps2*thermo_lnv%a(j))
   end do

   !! advance thermostat's positions
   do j = 1, nchain
      thermo_lnv%eta(j) = thermo_lnv%eta(j) + eps*thermo_lnv%v(j)
   end do

   !! scale barostat's velocity
   v_lnv = v_lnv*exp(-eps*thermo_lnv%v(1))

   !! update thermostat's velocity
   thermo_lnv%a(1) = thermo_lnv%Q_inv(1)*(mass_lnv*v_lnv*v_lnv - kT)

   tmp = exp(-eps4*thermo_lnv%v(2))
   thermo_lnv%v(1) = tmp*(tmp*thermo_lnv%v(1) + eps2*thermo_lnv%a(1))

   !! update thermostat's velocities
   do j = 2, nchain-1
      thermo_lnv%a(j) = thermo_lnv%Q_inv(j) &
           *(thermo_lnv%Q(j-1)*thermo_lnv%v(j-1)**2 - kT)

      tmp = exp(-eps4*thermo_lnv%v(j+1))
      thermo_lnv%v(j) = tmp*(tmp*thermo_lnv%v(j) + eps2*thermo_lnv%a(j))
   enddo

   thermo_lnv%a(nchain) = thermo_lnv%Q_inv(nchain) &
        *(thermo_lnv%Q(nchain-1)*thermo_lnv%v(nchain-1)**2 - kT)

   thermo_lnv%v(nchain) = thermo_lnv%v(nchain) + eps2*thermo_lnv%a(nchain)

   end do ! isuzuki
   end do ! isplit

end subroutine

!------------------------------------------------------------------------------
! applies exp(iL_{T-part}*dt)
subroutine Thermostat_integrate  &
           (nchain,thermostats,num_thermostats,timestep)
!------------------------------------------------------------------------------

   implicit none

   integer, intent(in) :: nchain
   integer, intent(in) :: num_thermostats
   type( Thermostat_type ), target :: thermostats(num_thermostats)
   real(dp), intent(in) :: timestep

   !..................................................

   integer :: isplit, isuzuki, ithermo, j
   real(dp) :: dt, kT, Ndof_kT, eps, eps2, eps4, tmp, Ekin2, vel_scale_factor
   type( Thermostat_type ), pointer :: thermo
   type( SystemCoordinate_type ), pointer :: system(:)

   dt = timestep / dble( num_split )

   !------------------------------

   do ithermo = 1, num_thermostats

      thermo => thermostats(ithermo)
      if (.not.thermo%activated) cycle

      assert(thermo%initialized)
      assert(thermo%num_system_coords==thermo%system_coord_id)

      system => thermo%system_coords(:)

      Ekin2 = 0.d0
      do j = 1, thermo%num_system_coords
         Ekin2 = Ekin2 + system(j)%mass * system(j)%vel**2
      end do

      kT = thermo%kT
      Ndof_kT = thermo%Ndof_kT

      vel_scale_factor = 1.d0

      do isplit  = 1, num_split
      do isuzuki = 1, num_suzuki

      eps = dt * coeff_suzuki(isuzuki, num_suzuki)

      eps2 = eps * 0.5d0
      eps4 = eps * 0.25d0

      !! Update force acting on the first thermostat oscillator.
      thermo%a(1) = thermo%Q_inv(1) * ( Ekin2 - Ndof_kT )

      !! Update thermostat velocity.
      thermo%v(nchain) = thermo%v(nchain) + eps2 * thermo%a(nchain)
      do j = nchain-1, 1, -1
         tmp = exp( - eps4 * thermo%v(j+1) )
         thermo%v(j) = tmp * ( tmp * thermo%v(j) + eps2 * thermo%a(j) )
      enddo

      !! Update thermostat position.
      do j = 1, nchain
         thermo%eta(j) = thermo%eta(j) + eps * thermo%v(j)
      enddo

      !! Velocity scaling factor.
      tmp = exp( - eps * thermo%v(1) )
      vel_scale_factor = vel_scale_factor * tmp
      Ekin2 = Ekin2 * tmp**2

      !! Update thermostat velocity.
      thermo%a(1) = thermo%Q_inv(1) * ( Ekin2 - Ndof_kT )
      tmp = exp( - eps4 * thermo%v(2) )
      thermo%v(1) = tmp * ( tmp * thermo%v(1) + eps2 * thermo%a(1) )

      do j = 2, nchain-1
         thermo%a(j) = thermo%Q_inv(j)  &
            * ( thermo%Q(j-1) * thermo%v(j-1)**2 - kT )
         tmp = exp( - eps4 * thermo%v(j+1) )
         thermo%v(j) =  &
            tmp * ( tmp * thermo%v(j) + eps2 * thermo%a(j) )
      enddo

      thermo%a(nchain) = thermo%Q_inv(nchain)  &
         * ( thermo%Q(nchain-1) * thermo%v(nchain-1)**2 - kT )

      thermo%v(nchain) = thermo%v(nchain) + eps2 * thermo%a(nchain)

      end do ! isuzuki
      end do ! isplit

      !! Update system's velocity.
      do j = 1, thermo%num_system_coords
         system(j)%vel = system(j)%vel * vel_scale_factor
      enddo

   end do ! ithermo

end subroutine

!------------------------------------------------------------------------------
function Thermostat_hamiltonian  &
         (nchain,thermostats,num_thermostats,npt) result(E)

! computes the thermostat terms in the extended Hamiltonian.
!------------------------------------------------------------------------------

  implicit none

  logical, intent(in) :: npt
  integer, intent(in) :: nchain
  integer, intent(in) :: num_thermostats
  type( Thermostat_type ), intent(in), target :: thermostats(num_thermostats)
  real(dp) :: E

  !!..................................................

  integer :: j, ithermo
  type( Thermostat_type ), pointer :: thermo

  !!..................................................

  E = 0.0d0

  do ithermo = 1, num_thermostats
     thermo => thermostats(ithermo)
     E = E + 0.5d0 * thermo%Q(1) * thermo%v(1)**2  &
           + thermo%Ndof_kT * thermo%eta(1)
     do j = 2, nchain
        E = E + 0.5d0 * thermo%Q(j) * thermo%v(j)**2  &
              + thermo%kT * thermo%eta(j)
     enddo
  enddo

  if ( npt ) then
     E = E + 0.5d0 * mass_lnv * v_lnv**2 ! + P_ext * Volume
     do j = 1, nchain
        E = E + 0.5d0 * thermo_lnv%Q(j) * thermo_lnv%v(j)**2  &
              + thermo_lnv%kT*thermo_lnv%eta(j)
     end do
  end if

end function


!------------------------------------------------------------------------------
FUNCTION NR_ran2(idum,icontrol)

! from second edition of numerical recipes
! call with negative idum to initialize
! then do not alter idum from successive calls
!------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: idum
  INTEGER, OPTIONAL, INTENT(IN) :: icontrol

  !..................................................

  REAL(dp) :: NR_ran2

  INTEGER :: im1, im2, imm1, ia1, ia2, iq1, iq2, ir1, ir2, ntab, ndiv

  REAL(dp) :: am, eps, rnmx

  PARAMETER (im1=2147483563,im2=2147483399,am=1.0d0/im1,imm1=im1-1, &

       ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,ntab=32, &

       ndiv=1+imm1/ntab,eps=1.2d-7,rnmx=1.0d0-eps)

  INTEGER :: idum2, j, k, iv(ntab), iy, iounit

  SAVE :: iv, iy, idum2

  DATA idum2/123456789/, iv/ntab*0/, iy/0/, iounit/10/

  !..................................................

!
!       if icontrol=0, generate a random number normally
!       if icontrol=1, load old random number table, dont get a new #
!       if icontrol=2, save current random number table and dont get a
!                      new random number
!

  IF ( .NOT. PRESENT ( icontrol)) THEN
     IF (idum<=0) THEN
        idum = max(-idum,1)
        idum2 = idum
        DO j = ntab + 8, 1, -1
           k = idum/iq1
           idum = ia1*(idum-k*iq1) - k*ir1
           IF (idum<0) idum = idum + im1
           IF (j<=ntab) iv(j) = idum
        END DO
        iy = iv(1)
     END IF
     k = idum/iq1
     idum = ia1*(idum-k*iq1) - k*ir1
     IF (idum<0) idum = idum + im1
     k = idum2/iq2
     idum2 = ia2*(idum2-k*iq2) - k*ir2
     IF (idum2<0) idum2 = idum2 + im2
     j = 1 + iy/ndiv
     iy = iv(j) - idum2
     iv(j) = idum
     IF (iy<1) iy = iy + imm1
     NR_ran2 = min(am*iy,rnmx)
  ELSE IF (icontrol==1) THEN
     OPEN (unit=iounit,file='random.table')
     DO j = 1, ntab
        READ (iounit,*) iv(j)
     END DO
     READ (iounit,*) iy
     READ (iounit,*) idum
     READ (iounit,*) idum2
     CLOSE (iounit)
  ELSE IF (icontrol==2) THEN
     OPEN (unit=iounit,file='random.table')
     DO j = 1, ntab
        WRITE (iounit,*) iv(j)
     END DO
     WRITE (iounit,*) iy
     WRITE (iounit,*) idum
     WRITE (iounit,*) idum2
     CLOSE (iounit)
  END IF

END FUNCTION


!------------------------------------------------------------------------------
FUNCTION NR_gasdev(idum)
!------------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(dp) :: NR_gasdev

  INTEGER :: idum

  !..................................................

  REAL(dp) :: v1, v2, r, fac, gset
  INTEGER :: iset
  DATA iset/0/
  SAVE iset, gset

  !..................................................

  IF (iset==0) THEN
1    CONTINUE
     v1 = 2.0d0 * NR_ran2( idum ) - 1.0d0
     v2 = 2.0d0 * NR_ran2( idum ) - 1.0d0
     r = v1*v1 + v2*v2
     IF (r>=1.0d0) GO TO 1
     fac = sqrt(-2.0d0*log(r)/r)
     gset = v1*fac
     NR_gasdev = v2*fac
     iset = 1
  ELSE
     NR_gasdev = gset
     iset = 0
  END IF

END FUNCTION


!------------------------------------------------------------------------------
function gaussian () result ( ret )

! returns random Gaussian distribution ~ exp[ -0.5 * x**2 ]
!------------------------------------------------------------------------------

  implicit none

  real(dp) :: ret

  !..................................................

  ret = NR_gasdev ( idum_ran2 )

end function


end module


