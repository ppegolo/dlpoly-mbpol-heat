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
  subroutine write_check(time,engunit,strsunit,bead_suffix)
    real(8), intent(in) :: time, engunit, strsunit
    integer :: check_file = 12!,ew1p_file=13
    real(8) :: tot_pot, tot_kin, tot_strs(9), ew1p_strs(9), st_ew1p(1:3,1:3)
    character(8), intent(in) :: bead_suffix

    real(8) :: ecpe,e1,e2,e3,tot_ew1p
    open(check_file, file='CHECK_ENERGY_STRESS'//bead_suffix,position='append')
    open(check_file+1, file='CHECK_ENERGY_STRESS2'//bead_suffix,position='append')
    open(check_file+2, file='CHECK_ENERGY3'//bead_suffix,position='append')

    call total_potential_energy(tot_pot)
    call total_energy_ew1p(tot_ew1p)
    call total_kinetic_energy(tot_kin)
    call total_stress(tot_strs)

    call total_stress_ew1p(st_ew1p)
    ew1p_strs(1:3) = st_ew1p(1,1:3)
    ew1p_strs(4:6) = st_ew1p(2,1:3)
    ew1p_strs(7:9) = st_ew1p(3,1:3)

    write(check_file,'(es12.5)') time
    write(check_file,'(2es15.7)') tot_kin/engunit, tot_pot/engunit
    write(check_file,'(9es15.7)') (tot_strs-ew1p_strs)*strsunit
    close(check_file)

    write(check_file+1,'(es12.5)') time
    write(check_file+1,'(4es15.7)') tot_kin/engunit, (sum(atomic_potential_energy2)+sum(delta_potential_energy)&
    +tot_ew1p)/engunit, &
    sum(delta_potential_energy)/engunit, tot_pot/engunit
    close(check_file+1)

    write(check_file+2,'(2es15.7)') time, (tot_pot-(sum(atomic_potential_energy2)+tot_ew1p))/engunit
    close(check_file+2)


  end subroutine write_check
!
end module
