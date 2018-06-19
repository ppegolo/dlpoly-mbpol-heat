!------------------------------------------------------------------------------
module centroid_module
!
! internal variables used in CMD
!------------------------------------------------------------------------------

   use nose_hoover_module, only: Thermostat_type

   implicit none

   save

   !..................................................

   real(8), parameter :: conv_energy = 418.4d0

   real(8), parameter :: pi  = 3.141592653589793d0

   real(8), parameter :: hbar = 6.350780668d0
                                                                              
   real(8), parameter :: kB = 8.31451115d-1

   !..................................................

   !! Files for PIMD/CMD.

   integer :: file_out  = 101

   integer :: file_pimd = 1001

   integer :: file_cmd  = 1002

   integer :: file_nh   = 1003

   integer :: file_npt  = 1004

   integer :: file_pmf  = 1005

   integer :: file_ind  = 1006

   integer :: file_cv  = 1007

   integer :: file_stress  = 1008

   integer :: file_history_pimd = 2001

   integer :: file_dipole_pimd = 2003

   integer :: file_pos_cmd = 3001

   integer :: file_vel_cmd = 3002

   integer :: file_force_cmd = 3003

   integer :: file_dipole_cmd = 3004

   integer :: file_history_bead = 3101

   integer :: file_force_cv  = 3007

   integer :: file_initial_config_pimd = 4001

   integer :: file_final_config_pimd = 4002

   integer :: file_initial_dipole_pimd = 5001

   integer :: file_final_dipole_pimd = 5002

   !----------------------------------
  
   !! PIMD/CMD global parameters.

   logical :: restart_run = .false.

   logical :: restart_dipole = .false.

   logical :: restart_cmd = .false.

   logical :: pimd_npt = .false.

   logical :: cluster_run = .false.

   logical :: cmd_run = .true.

   logical :: cmd_nve = .true.

   logical :: trpmd_run = .false.

   logical :: w_bead = .false.

   logical :: set_vel = .false.

   logical :: activate = .true.

   integer :: iseed = -2209

   integer :: nbead = 1

   integer :: natom, nmol

   integer :: nequilib_cmd = 0
 
   integer :: nstep_cmd = 0

   real(8) :: tau_vol = 1.0d0

   real(8) :: adiab_param = 1.d0

   real(8) :: fac_nm = 1.d0

   real(8) :: kT, beta, NkT

   real(8) :: tau

   real(8) :: Epot_spring, Epot_deriv, Equal_part, Equal_part_prim

   real(8) :: Ekin_fict, Epot_fict, Etot_fict, temp_fict

   real(8) :: Ekin_pimd, Epot_pimd, Etot_pimd

   real(8) :: Ekin_prim, Etot_prim

   real(8) :: Ekin_cmd, Epot_cmd, Etot_cmd, temp_cmd

   real(8) :: Etot_nose

   real(8) :: vir_pimd, P_pimd

   !----------------------------------

   !! ISOTOPIC MIXTURES.

   logical :: hod_in_h2o = .false.
   logical :: hod_in_d2o = .false.
   logical :: h2o_in_d2o = .false.
   logical :: d2o_in_h2o = .false.
   integer :: nmol_mix = 0

   !----------------------------------

   !! PIMD/CMD: normal-mode variables.

   real(8), allocatable :: lambda_nmode(:)

   real(8), allocatable :: omega_nmode(:)

   real(8), allocatable :: mass_nmode(:,:)

   real(8), allocatable :: fict_mass_nmode(:,:)

   real(8), allocatable :: nmode_to_cart(:,:)

   real(8), allocatable :: cart_to_nmode(:,:)

   real(8), allocatable :: tcart_to_tnmode(:,:)

   !----------------------------------

   !! PIMD/CMD: cartesian variables.

   real(8), allocatable :: phys_mass(:)

   !----------------------------------

   !! PIMD/CMD: thermostat variables.

   integer :: nchain = 1

   integer :: nthermo

   type( Thermostat_type ), allocatable :: thermo(:,:)

   !----------------------------------

   real(8), allocatable :: nmbuff(:)
   integer              :: mxnmbuff

end module
