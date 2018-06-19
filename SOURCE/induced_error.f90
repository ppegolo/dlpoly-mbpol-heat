! computes L^2 and L^inf dipoles' errors

subroutine induced_error(natms, polr2, dipx, dipy, dipz, &
    & efieldkx, efieldky, efieldkz, emux, emuy, emuz, &
    & efddmurecx, efddmurecy, efddmurecz, efdcrecx, efdcrecy, efdcrecz, &
    & extx, exty, extz, err_L2, err_Linf)

   use multibead

   implicit none

   integer, intent(in) :: natms

   real(8), intent(in) :: polr2(*), dipx(*), dipy(*), dipz(*)
   real(8), intent(in) :: efieldkx(*), efieldky(*), efieldkz(*)
   real(8), intent(in) :: emux(*), emuy(*), emuz(*)
   real(8), intent(in) :: efddmurecx(*), efddmurecy(*), efddmurecz(*)
   real(8), intent(in) :: efdcrecx(*), efdcrecy(*), efdcrecz(*)
   real(8), intent(in) :: extx, exty, extz

   real(8), intent(out) :: err_L2, err_Linf

   real(8), parameter :: r4pie0 = 138935.4835d0
   real(8) :: dmux, dmuy, dmuz

   integer :: i, iatm1, iatm2

#ifdef MPI
#  include "mpif.h"
#endif


!  set up atoms numbers for this node

   iatm1 = (bead_rank*natms)/bead_size + 1
   iatm2 = ((bead_rank + 1)*natms)/bead_size

   err_L2 = 0.d0
   err_Linf = 0.d0

   do i = iatm1, iatm2
      if(polr2(i).le.1.d-6) cycle

      dmux = dipx(i) - polr2(i)*(efieldkx(i) + emux(i) &
          & + efddmurecx(i) + efdcrecx(i) + extx)/r4pie0
      dmuy = dipy(i) - polr2(i)*(efieldky(i) + emuy(i) &
          & + efddmurecy(i) + efdcrecy(i) + exty)/r4pie0
      dmuz = dipz(i) - polr2(i)*(efieldkz(i) + emuz(i) &
          & + efddmurecz(i) + efdcrecz(i) + extz)/r4pie0

      err_L2 = err_L2 + dmux**2 + dmuy**2 + dmuz**2

      err_Linf = max(err_Linf, abs(dmux))
      err_Linf = max(err_Linf, abs(dmuy))
      err_Linf = max(err_Linf, abs(dmuz))
   end do

   if (bead_size.gt.1) call gdsum(err_L2, 1, dmux)

   err_L2 = sqrt(err_L2)

#ifdef MPI
   call MPI_REDUCE(err_Linf, dmux, 1, MPI_DOUBLE_PRECISION, &
                   & MPI_MAX, 0, comm_bead, i)
   err_Linf = dmux
#endif /* MPI */

end subroutine induced_error
