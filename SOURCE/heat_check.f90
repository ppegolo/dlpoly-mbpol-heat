!-----------------------------------------------------------------------------------------------------------------------
! This module writes the total energy and stress obtained as the sum of atomic contributions.
! Author: Paolo Pegolo
! Last update: 16-May-2018
!-----------------------------------------------------------------------------------------------------------------------
module heat_check
!
use heatcurrent
implicit none
!
contains
!
  subroutine write_check(time,engunit,strsunit)
    real(8), intent(in) :: time, engunit, strsunit
    integer :: check_file = 12, ew1p_file=13
    real(8) :: tot_pot, tot_kin, tot_strs(9), ew1p_strs(9)
    open(check_file, file='energy_stress_check.dat',position='append')
    open(ew1p_file, file='ew1p_stress.dat',position='append')
    call total_potential_energy(tot_pot)
    call total_kinetic_energy(tot_kin)
    call total_stress(tot_strs)

    call total_stress_ew1p(ew1p_strs)
    write(check_file,'(es12.5)') time
    write(check_file,'(2es15.7)') tot_kin/engunit, tot_pot/engunit
    write(check_file,'(9es15.7)') tot_strs*strsunit
    close(check_file)

    write(ew1p_file,'(9es15.7)') ew1p_strs
    close(ew1p_file)
  end subroutine write_check
!
end module
