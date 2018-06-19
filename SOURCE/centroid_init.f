c==============================================================================
c
      subroutine centroid_init_params 
     x   (mxatms,mxtmls,ntpmls,nummols,first_wat_ind,natms_r, 
     x    molmix_nskip,weight, temp)
c
c==============================================================================

      use centroid_module, only: 
     x                     pi, hbar,
     x                     file_out, file_pmf,
     x                     file_nh, file_npt, file_stress,
     x                     file_pimd, file_cmd,
     x                     file_history_pimd,
     x                     file_pos_cmd, file_vel_cmd, file_force_cmd,
     x                     file_dipole_pimd, file_dipole_cmd,
     x                     nbead, natom, nmol,
     x                     hod_in_h2o, hod_in_d2o,
     x                     d2o_in_h2o, h2o_in_d2o,
     x                     nmol_mix,
     x                     kB, kT, beta, NkT,
     x                     pimd_npt,
     x                     restart_run, cmd_run, trpmd_run,
     x                     adiab_param, Equal_part, Equal_part_prim,
     x                     cart_to_nmode, nmode_to_cart,  
     x                     tcart_to_tnmode,  
     x                     lambda_nmode, omega_nmode,  
     x                     mass_nmode, fict_mass_nmode,  
     x                     phys_mass,
     x                     mxnmbuff, nmbuff
      use Langevin_Thermostat,      only:
     x                     langevin,
     x                     freq_scale, force_scale


      use multibead, only: mb_rank
c     !..................................................

      implicit none

      integer, intent(in) :: mxatms, mxtmls, ntpmls

      integer, intent(in) :: first_wat_ind, natms_r, molmix_nskip

      integer, intent(in) :: nummols(mxtmls)

      real(8), intent(in) :: temp

      real(8) :: weight(mxatms)
      real(8) :: massH, massD, omega_P

      integer :: idim, ibead, iatom, imol, itpmls

c     !..................................................

      natom = natms_r

!     Total number of molecules (M site for TTM = last molecule).
      do itpmls = 1, ntpmls-1
         nmol = nmol + nummols(itpmls)
      enddo

      kT = kB * temp
      beta = 1.0d0 / kT
      NkT = dble( natom ) * kT

c     !---------------------------------------

      allocate( lambda_nmode(nbead) )
      allocate( omega_nmode(nbead) )

      allocate( mass_nmode(natom,nbead) )
      allocate( fict_mass_nmode(natom,nbead) )

      allocate( nmode_to_cart(nbead,nbead) )
      allocate( cart_to_nmode(nbead,nbead) )
      allocate( tcart_to_tnmode(nbead,nbead) )

c     !---------------------------------------

c     PIMD/CMD: allocate cartesian variables.
      allocate( phys_mass(natom) )

c     !---------------------------------------

c     Files for PIMD/CMD output.
      if (.not.mb_rank.eq.0) goto 100
!SR
#ifdef IPI
      goto 100    ! to skip writing output files for ipi 
#endif /* IPI */
      if (cmd_run) then

         open(file_out, file='OUTPUT_PIMD')
         open(file_nh, file='STATIS_NH')
         open(file_cmd, file='STATIS_CMD')
         open(file_pos_cmd, file='POSITION_CMD')
         open(file_vel_cmd, file='VELOCITY_CMD')
         open(file_force_cmd, file='FORCE_CMD')
         close(file_out)
         close(file_nh)
         close(file_cmd)
         close(file_pos_cmd)
         close(file_vel_cmd)
         close(file_force_cmd)
      
      elseif (.not.cmd_run) then
 
         open(file_out, file='OUTPUT_PIMD')
         open(file_nh, file='STATIS_NH')
         open(file_pimd, file='STATIS_PIMD')
         open(file_stress, file='STATIS_STRESS')
         close(file_out)
         close(file_nh)
         close(file_pimd)
         close(file_stress)
         
         if (pimd_npt) then
            open(file_npt, file='STATIS_NPT')
            close(file_npt)
         endif

      end if

 100  continue
c     !---------------------------------------

c     Compute constant factor for the quantum kinetic energy
c     using virial expression.
      Equal_part = 1.5d0 * dble( natom ) / beta
      Equal_part_prim = 1.5d0 * dble( nbead*natom ) / beta

c     !---------------------------------------

      if (hod_in_h2o) then
         write(78,*) 'first_wat_ind = ', first_wat_ind           
         massD = 2.0141d0
         do imol = 1, nmol_mix
            iatom = 1 + 3 * (imol - 1) * molmix_nskip + first_wat_ind
            write(78,*) 'orig weight(', iatom, ') = ', weight(iatom)
            weight(iatom) = massD
            write(78,*) 'new  weight(', iatom, ') = ', weight(iatom)
         enddo
         flush(78)
      else if (hod_in_d2o) then
         massH = 1.0079d0
         do imol = 1, nmol_mix
            iatom = 1 + 3 * (imol - 1) * molmix_nskip + first_wat_ind
            write(78,*) 'orig weight(', iatom, ') = ', weight(iatom)
            weight(iatom) = massH
            write(78,*) 'new  weight(', iatom, ') = ', weight(iatom)
         enddo
         flush(78)
      else if (h2o_in_d2o) then
         massH = 1.0079d0
         do imol = 1, nmol_mix
            iatom = 1 + 3 * (imol - 1) * molmix_nskip + first_wat_ind
            write(78,*) 'orig weight(', iatom, ') = ', weight(iatom)
            weight(iatom) = massH
            write(78,*) 'new  weight(', iatom, ') = ', weight(iatom)
            iatom = 2 + 3 * (imol - 1) * molmix_nskip + first_wat_ind
            write(78,*) 'orig weight(', iatom, ') = ', weight(iatom)
            weight(iatom) = massH
            write(78,*) 'new  weight(', iatom, ') = ', weight(iatom)
         enddo
         flush(78)
      else if (d2o_in_h2o) then
         massD = 2.0141d0
         do imol = 1, nmol_mix
            iatom = 1 + 3 * (imol - 1) * molmix_nskip + first_wat_ind
            write(78,*) 'orig weight(', iatom, ') = ', weight(iatom)
            weight(iatom) = massD
            write(78,*) 'new  weight(', iatom, ') = ', weight(iatom)
            iatom = 2 + 3 * (imol - 1) * molmix_nskip + first_wat_ind
            write(78,*) 'orig weight(', iatom, ') = ', weight(iatom)
            weight(iatom) = massD
            write(78,*) 'new  weight(', iatom, ') = ', weight(iatom)
         enddo
         flush(78)
      endif

c     Initialize physical masses.
      do iatom = 1, natom
         phys_mass(iatom) = weight(iatom)
      end do

c     !---------------------------------------

c     Get transformation matrix from cartesian to normal mode coordinates.
      call transformation_matrix

c     !---------------------------------------

c     Normal mode = 1 corresponds to the path-centroid. 
c     For the path-centroid use actual physical mass 
c     (as defined in transformation_matrix).
c     All other normal mode masses scaled by the adiabacity parameter.

      if ( trpmd_run ) then

         do ibead = 1, nbead
            do iatom = 1, natom
               fict_mass_nmode(iatom,ibead) = phys_mass(iatom)
            end do
         end do

      else

         do iatom = 1, natom
            fict_mass_nmode(iatom,1) = phys_mass(iatom)
         end do

         do ibead = 2, nbead
            do iatom = 1, natom
               fict_mass_nmode(iatom,ibead) =
     x            mass_nmode(iatom,ibead) * adiab_param**2 
            end do
         end do

      end if

c     !---------------------------------------

c     Normal-mode frequencies.
c     ibead = 1 -> path-centroid.

      omega_P = sqrt( dble(nbead) ) / hbar / beta

      omega_nmode(1) = omega_P

      do ibead = 2, nbead
         omega_nmode(ibead) = omega_P / adiab_param
      end do

!VB: buffer for nmode<->cart transforms

      mxnmbuff = 24*3*natom*nbead;
      allocate(nmbuff(mxnmbuff))

      end subroutine


c==============================================================================
c
      subroutine centroid_init_nose_hoover(iseed,nvx,nvy,nvz)
c
c     Initialize Nose-Hoover chain thermostats
c==============================================================================

       use multibead
       use nose_hoover_module, only:
     x                         nose_hoover_module_init,
     x                         Thermostat_init,
     x                         Thermostat_link,
     x                         gaussian

       use centroid_module, only: 
     x                      natom, nbead, 
     x                      kT, hbar, pi,  
     x                      nchain, nthermo, thermo, 
     x                      tau, omega_nmode,
     x                      fict_mass_nmode,
     x                      fac_nm

c     !..................................................

      implicit none

      integer, intent(in) :: iseed
      real(8) :: nvx(*),nvy(*),nvz(*)
      integer :: ibead, iatom, iatm1, iatm2, ithermo, j

      real(8), external :: NR_gauss
      real(8) :: dummy

c     !..................................................

      dummy = NR_gauss(iseed)

c     !---------------------------------------

      call nose_hoover_module_init

c     !---------------------------------------

      iatm1 = (bead_rank*natom)/bead_size+1
      iatm2 = ((bead_rank+1)*natom)/bead_size

      nthermo = 3 * (iatm2 - iatm1 + 1)

      allocate( thermo(3,iatm2 - iatm1 + 1) )

c     !---------------------------------------

      nbead = ring_size

      ithermo = 0

      do ibead = 1, nbead
         do iatom = 1, natom
            if (iatom.lt.iatm1.or.iatom.gt.iatm2
     x          .or.ibead.ne.(ring_rank+1)) then
               do j = 1, nchain
                  dummy = gaussian()
                  dummy = gaussian()
                  dummy = gaussian()
               end do
               cycle
            end if ! does not belong to this PE

            tau = 1.d0 / omega_nmode(ibead)
            tau = tau / fac_nm

            ithermo = ithermo + 1
            call Thermostat_init(nchain,thermo(1,ithermo),1,0,kT,tau)
            call Thermostat_link(thermo(1,ithermo),
     x                      fict_mass_nmode(iatom,ibead), nvx(iatom))

            call Thermostat_init(nchain,thermo(2,ithermo),1,0,kT,tau)
            call Thermostat_link(thermo(2,ithermo),
     x                      fict_mass_nmode(iatom,ibead), nvy(iatom))

            call Thermostat_init(nchain,thermo(3,ithermo),1,0,kT,tau)
            call Thermostat_link(thermo(3,ithermo),
     x                      fict_mass_nmode(iatom,ibead), nvz(iatom))
         enddo
      enddo

      end subroutine


c==============================================================================
c
      subroutine setup_vel_nmode(iseed,nvx,nvy,nvz)
c
c     Initialize normal mode velocities
c==============================================================================

      use multibead
      use centroid_module, only: natom, nbead, kB, kT, pi,  
     x                           fict_mass_nmode

c     !..................................................
  
      implicit none
 
      integer, intent(in) :: iseed

      real(8) :: nvx(*),nvy(*),nvz(*)
  
c   !..................................................
  
      integer :: iatom, ibead
  
      real(8), external :: NR_gauss
      real(8) :: sigma

c    !..................................................

      do ibead = 1, nbead
         do iatom = 1, natom
            sigma = sqrt( kT / fict_mass_nmode(iatom,ibead) )
            if (ibead.eq.(ring_rank+1)) then
               nvx(iatom) =  NR_gauss( iseed ) * sigma
               nvy(iatom) =  NR_gauss( iseed ) * sigma
               nvz(iatom) =  NR_gauss( iseed ) * sigma
            else
               sigma = NR_gauss( iseed )
               sigma = NR_gauss( iseed )
               sigma = NR_gauss( iseed )
            end if
         end do
      end do

      end subroutine


c==============================================================================
c
      subroutine fix_vel_nmode(nvx,nvy,nvz,temp)
c
c     Scale velocities to match the temperature
c==============================================================================

      use multibead
      use centroid_module, only: iseed, natom, nbead, kB,
     x                           fict_mass_nmode

c     !..................................................

      implicit none

      real(8), intent(in) :: temp
      real(8) :: nvx(*),nvy(*),nvz(*)

c     !..................................................

#include "mpif.h"

      integer :: ibead, iatom, iatm1, iatm2
      real(8) :: scale_factor
      real(8) :: E_kin, temp_ini
      real(8) :: tmp

c     !..................................................

      E_kin = 0.d0

      ibead = ring_rank + 1

      iatm1 = (bead_rank*natom)/bead_size+1
      iatm2 = ((bead_rank+1)*natom)/bead_size

      do iatom = iatm1, iatm2
         tmp = nvx(iatom)**2 + nvy(iatom)**2 + nvz(iatom)**2
         E_kin = E_kin + 0.5d0 * fict_mass_nmode(iatom,ibead) * tmp
      enddo

      call MPI_ALLREDUCE(E_kin, tmp, 1, MPI_DOUBLE_PRECISION,
     x                   MPI_SUM, comm_mb, iatom)

      E_kin = tmp

      temp_ini = 2.d0 * E_kin / 3.d0 / kB / dble( natom * nbead )
      scale_factor = sqrt( temp / temp_ini )

      do iatom = iatm1, iatm2
         nvx(iatom) = nvx(iatom) * scale_factor
         nvy(iatom) = nvy(iatom) * scale_factor
         nvz(iatom) = nvz(iatom) * scale_factor
      end do

      end subroutine


c==============================================================================
c
      subroutine scale_vel_centroid(nvx,nvy,nvz,temp)
c
c     Scale centroid velocities by fixing total momentum equal to zero
c==============================================================================

      use multibead
      use centroid_module, only: natom, nbead, kB,
     x                           fict_mass_nmode
  
c     !..................................................

      implicit none

      real(8), intent(in) :: temp
      real(8) :: nvx(*),nvy(*),nvz(*)

c     !..................................................

#include "mpif.h"

      integer :: ibead, iatom, iatm1, iatm2
      real(8) :: E_kin, temp_ini, tmp, scale_factor
      real(8) :: ptot(3), pmy(3)

      ibead = ring_rank+1
      if (ibead.ne.1) goto 111
  
c     !..................................................
  
c     Set total momentum equal to zero for path centroid.

      iatm1 = (bead_rank*natom)/bead_size+1
      iatm2 = ((bead_rank+1)*natom)/bead_size

      pmy(1:3) = 0.d0
      do iatom = iatm1, iatm2
         pmy(1) = pmy(1) + fict_mass_nmode(iatom,1) * nvx(iatom)
         pmy(2) = pmy(2) + fict_mass_nmode(iatom,1) * nvy(iatom)
         pmy(3) = pmy(3) + fict_mass_nmode(iatom,1) * nvz(iatom)
      enddo

      call MPI_ALLREDUCE(pmy, ptot, 3, MPI_DOUBLE_PRECISION,
     x                   MPI_SUM, comm_bead, iatom)

      ptot = ptot / dble( natom )

      E_kin = 0.d0

      do iatom = iatm1, iatm2
         nvx(iatom) =  
     x        nvx(iatom) - ptot(1) / fict_mass_nmode(iatom,1)
         nvy(iatom) =  
     x        nvy(iatom) - ptot(2) / fict_mass_nmode(iatom,1)
         nvz(iatom) =  
     x        nvz(iatom) - ptot(3) / fict_mass_nmode(iatom,1)

         tmp = nvx(iatom)**2 + nvy(iatom)**2 + nvz(iatom)**2
         E_kin = E_kin + 0.5d0 * fict_mass_nmode(iatom,1) * tmp
      enddo

      call MPI_ALLREDUCE(E_kin, tmp, 1, MPI_DOUBLE_PRECISION,
     x                   MPI_SUM, comm_bead, iatom)
      E_kin = tmp

      temp_ini = 2.d0 * E_kin / 3.d0 / kB / dble( natom )
      scale_factor = sqrt( temp / temp_ini )

c     !---------------------------------------

c     Scale velocity according to target temperature.
      E_kin = 0.d0

      do iatom = iatm1, iatm2
         nvx(iatom) = nvx(iatom) * scale_factor
         nvy(iatom) = nvy(iatom) * scale_factor
         nvz(iatom) = nvz(iatom) * scale_factor

         tmp = nvx(iatom)**2 + nvy(iatom)**2 + nvz(iatom)**2
         E_kin = E_kin + 0.5d0 * fict_mass_nmode(iatom,1) * tmp
      enddo

#if 0
      call MPI_ALLREDUCE(E_kin, tmp, 1, MPI_DOUBLE_PRECISION,
     x                   MPI_SUM, comm_bead, iatom)
      E_kin = tmp

      temp_ini = 2.d0 * E_kin / 3.d0 / kB / dble( natom )
      write(50+mb_rank,*) __FILE__,__LINE__,temp_ini
#endif

 111  continue

      call MPI_BARRIER(comm_mb, iatom)

      end subroutine


c==============================================================================
c
      subroutine setup_npt_pimd(iseed,tau_vol)
c
c==============================================================================

      use centroid_module, only: natom, kB, kT

      use nose_hoover_module, only : M, thermo_lnv, v_lnv,
     x                               f_lnv, c2_lnv, mass_lnv

      implicit none

      integer :: j

      integer, intent(in) :: iseed

      real(8), external :: NR_gauss

      real(8) :: tau_vol, sigma

      !...........................................

      c2_lnv = 1.d0 + 1.d0 / dble( natom )

      v_lnv = 0.d0

      f_lnv = 0.d0

      mass_lnv = 3.d0 * dble( natom+1 ) * kT * tau_vol * tau_vol

      thermo_lnv%Q(1:M) = kT * tau_vol * tau_vol

      thermo_lnv%Q_inv(1:M)= 1.0d0 / thermo_lnv%Q(1)

      thermo_lnv%kT = kT

      thermo_lnv%Ndof_kT = kT

      thermo_lnv%eta(:) = 0.0d0

      do j = 1, M
         sigma = sqrt( kT / thermo_lnv%Q(j) )                     
         thermo_lnv%v(j) = NR_gauss( iseed ) * sigma
      enddo

      do j = 2, M
         thermo_lnv%a(j) = thermo_lnv%Q_inv(j) *
     x        (thermo_lnv%Q(j-1)*thermo_lnv%v(j-1)**2-thermo_lnv%kT)
      enddo


      end subroutine


c==============================================================================
      function NR_gauss( idum )
c
c     Generate a gaussian distribution. From Numerical Recipes.
c==============================================================================

      implicit none

c     !! Return value
      real(8) :: NR_gauss

c     !! Arguments
      integer :: idum

c     !! Locals
      integer :: iset
      real(8) :: fac, gset, rsq, v1, v2, NR_ran1
      data iset/0/
      save iset, gset

c     !..................................................

      if ( iset==0 ) then
1     continue
       v1 = 2.d0 * NR_ran1( idum ) - 1.d0
       v2 = 2.d0 * NR_ran1( idum ) - 1.d0
       rsq = v1*v1 + v2*v2
      if ( rsq >= 1.d0 .or. rsq == 0.d0 ) go to 1
       fac = sqrt( -2.d0*log(rsq)/rsq )
       gset = v1*fac
       NR_gauss = v2*fac
       iset=1
      else
       NR_gauss = gset
       iset = 0
      endif

      end function


c==============================================================================
      function NR_ran1( idum )
c
c     Generate uniform distribution. From Numerical Recipes.
c==============================================================================

      implicit none

c     !! Arguments
      integer :: idum, idum1

c     !! Locals
      real(8) :: NR_ran1
      integer :: ia, im, iq, ir, ntab, ndiv
      real(8) :: am, eps, rnmx
      parameter(ia=16807, im=2147483647, am=1.d0/im, iq=127773, ir=2836)
      parameter(ntab=32, ndiv=1+(im-1)/ntab, eps=1.2d-7, rnmx=1.d0-eps)
      integer :: j, k, iv(ntab), iy
      data iv /ntab*0/, iy /0/
      save iv, iy

c     !..................................................

      if ( idum<=0 .or. iy==0 ) then
!     idum = max(-idum1,1)
      idum = max(-idum,1)
      do j = NTAB + 80 , 1, -1   !'80' is the num of warmup steps
        k = idum/iq
        idum = ia*(idum-k*iq) - ir*k
        if ( idum<0 ) idum = idum+im
        if ( j<=ntab ) iv(j) = idum
      end do
      iy = iv(1)
      endif
      k = idum/iq
      idum = ia*(idum-k*iq) - ir*k
      if ( idum<0 ) idum = idum + im
      j = 1 + iy/ndiv
      iy = iv(j)
      iv(j) = idum
      NR_ran1 = min( am*iy, rnmx )

      end function

