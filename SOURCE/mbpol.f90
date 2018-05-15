#include "assert.h"

!#define VERBOSE0        ! this is replaced by DEBUG
#define ENABLE_X2B goahead
#define ENABLE_X3B aswell

!============================================================================

module mbpol

!----------------------------------------------------------------------------
#ifdef HEAT_CURRENT
use heatcurrent, only: update_stress_mbpol, update_energy_mbpol
#endif /* HEAT_CURRENT */


implicit none

real(8), private, save :: r2f = 0.d0
real(8), private, save :: r3f = 0.d0

!----------------------------------------------------------------------------

public :: mbpol_init
public :: mbpol_fini

public :: mbpol_forces ! does not MPI_REDUCE the forces

!----------------------------------------------------------------------------

real(8), private, allocatable, save :: monomers(:)

private :: do_dimer, do_trimer
private :: image_vectors, image_monomers

real(8), private, parameter :: engunit = 418.4d0

#ifdef HEAT_CURRENT
real(8), parameter :: third=1.d0/3.d0
#endif /* HEAT_CURRENT */

!============================================================================

contains

!============================================================================

subroutine mbpol_init(nrite, nttm2)

   use multibead, only: bead_rank, mb_abort

   implicit none

   integer, intent(in) :: nrite, nttm2

   integer :: ierr, nw

   assert(.not.allocated(monomers))

   call mbpol_2b_cutoff(r2f)
   call mbpol_3b_cutoff(r3f)

   if (bead_rank.eq.0) then
      write(nrite,*) ' ** Using MB-pol for water-water interactions **'
#     ifdef ENABLE_X2B
      write(nrite,*) ' ** 2B cutoff = ', r2f
#     endif
#     ifdef ENABLE_X3B
      write(nrite,*) ' ** 3B cutoff = ', r3f
#     endif
   end if

   nw = nttm2/3

   allocate(monomers(9*nw), stat=ierr)
   if (ierr.ne.0) then
      write(nrite,*) ' ** Error ** : could not allocate memory'
      call mb_abort()
   end if

end subroutine mbpol_init

!----------------------------------------------------------------------------

subroutine mbpol_fini()

   implicit none
   deallocate(monomers)

end subroutine mbpol_fini

!----------------------------------------------------------------------------

subroutine mbpol_forces(imcon,cell,nttm2,listttm2,wat_monomer_ener, &
   & xxx,yyy,zzz,fxx,fyy,fzz,engsrp,virsrp,stress) !FP: added

   use multibead, only: bead_rank, bead_size

   implicit none

   integer, intent(in) :: imcon
   real(8), intent(in) :: cell(*)
   integer, intent(in) :: nttm2, listttm2(:)
   real(8), intent(in) :: xxx(*), yyy(*), zzz(*)
   real(8), intent(in) :: wat_monomer_ener(*)  ! SR

   real(8), intent(inout) :: fxx(*), fyy(*), fzz(*)
   real(8), intent(inout) :: engsrp, virsrp, stress(9) !FP: added

   real(8) :: eng2, vir2
   real(8) :: eng3, vir3
   real(8) :: strs(9) !FP: added
   real(8) :: wat_monomer_ceiling  !SR
   integer :: t1,t2,l           ! SR

   real(8) :: buffer(8)

   integer :: i, j, k, o1, o2, o3,tmp_nw
   integer(kind=8) :: nw, nw2, nw3,t,tmp_bead_size
   integer(kind=8) :: istart, iend, istep


   real(8) :: timeint,timetot


   !SR    ceiling value --  set to 60 kcal/mol

   wat_monomer_ceiling=60          ! kcal/mol
   wat_monomer_ceiling=wat_monomer_ceiling*engunit          ! internal units

   !
   ! copy waters into the monomers(:) and image every molecule
   !

!SR adding ceiling for mbpol potential mono energ
! but it is not helpful for optimization to skip monomers if energy is
! greater than ceiling value
! so fixing it temporarily


   nw = nttm2/3

#ifdef OPT
   l=0

! SR: Water monomer ceiling is removed ...
! for geometry optimization with MBpol,  the following criteria
! is enabled... If the engsrp is <-1.d7, it won't evalaute loops
      if(engsrp<-1.d7) then
         engsrp=abs(engsrp)
         goto 24050
      endif
!  do t1 = 1, nttm2/3

!   if(wat_monomer_ener(t1)>wat_monomer_ceiling) then
!      engsrp=abs(engsrp)
!      goto 24050
!     nw=nw-1
!     cycle
!   endif
! end do

!  do t2=1,3              ! 3 atoms
!   t=(t1-1)*3+t2
!   i = listttm2(t)
!   l=l+1
!   monomers(3*l - 2) = xxx(i)
!   monomers(3*l - 1) = yyy(i)
!   monomers(3*l - 0) = zzz(i)
!  end do
!end do
!
!  tmp_nw=nw          ! to pass same integer kind
#endif /* OPT */


 do t = 1, nttm2
    i = listttm2(t)
    monomers(3*t - 2) = xxx(i)
    monomers(3*t - 1) = yyy(i)
    monomers(3*t - 0) = zzz(i)
 end do



  tmp_nw=nw          ! to pass same integer kind


  call image_monomers(imcon, cell, tmp_nw, monomers)

   !
   ! loop over the pairs
   !

   nw2 = (nw*(nw - 1))/2

   eng2 = 0.d0
   vir2 = 0.d0

#  ifdef ENABLE_X2B

#ifdef DEBUG
      call timchk(0,timeint)
      timetot=timeint
#endif

   i = 1
   j = 2 + bead_rank

   tmp_bead_size=bead_size
   istart = bead_rank
   iend   = nw2 - 1
   istep  = bead_size
   do t = istart, iend, istep
!   do t = 1, nw2
!      if (mod(t,tmp_bead_size).eq.bead_rank) then
         o1 = listttm2(3*i - 2)
         o2 = listttm2(3*j - 2)
#ifdef HEAT_CURRENT
         call do_dimer(imcon,cell, &
              & monomers(9*(i-1)+1),monomers(9*(j-1)+1), &
              & fxx(o1),fyy(o1),fzz(o1),fxx(o2),fyy(o2),fzz(o2), &
              & eng2,vir2,strs,o1,o2) !FP: added.
#else
          call do_dimer(imcon,cell, &
               & monomers(9*(i-1)+1),monomers(9*(j-1)+1), &
               & fxx(o1),fyy(o1),fzz(o1),fxx(o2),fyy(o2),fzz(o2), &
               & eng2,vir2,strs)
#endif /* HEAT_CURRENT */

         stress(1) = stress(1) + strs(1) !FP: added
         stress(2) = stress(2) + strs(2) !FP: added
         stress(3) = stress(3) + strs(3) !FP: added
         stress(4) = stress(4) + strs(4) !FP: added
         stress(5) = stress(5) + strs(5) !FP: added
         stress(6) = stress(6) + strs(6) !FP: added
         stress(7) = stress(7) + strs(7) !FP: added
         stress(8) = stress(8) + strs(8) !FP: added
         stress(9) = stress(9) + strs(9) !FP: added
!      end if ! mod(t,bead_size).eq.bead_rank

!      j = j+1
!      if (j.gt.nw) then
!         i = i+1
!         j = i+1
!      end if

      j = j + bead_size
      if (j > nw) then
         do while (j > nw) ! -- only needed near end of pair list
            i = i + 1
            j = i + j - nw
         end do
      end if

   end do ! t = 1, nw2

#ifdef DEBUG
      call timchk(0,timeint)
      write(4444,"(a,F10.3,a)") 'Time taken for MB-pol 2B       = ',timeint-timetot,' sec'
      write(4444,*)
#endif


#  endif /* ENABLE_X2B */

   !
   ! loop over the triples
   !

   nw3 = (nw*(nw - 1)*(nw - 2))/6

   eng3 = 0.d0
   vir3 = 0.d0

#  ifdef ENABLE_X3B

#ifdef DEBUG
      call timchk(0,timeint)
      timetot=timeint
#endif


   i = 1
   j = 2
   k = 3 + bead_rank

!   tmp_bead_size=bead_size
   istart = bead_rank
   iend   = nw3 - 1
   istep  = bead_size
   do t = istart, iend, istep
!   do t = 1, nw3
!      if (mod(t, tmp_bead_size).eq.bead_rank) then
         o1 = listttm2(3*i - 2)
         o2 = listttm2(3*j - 2)
         o3 = listttm2(3*k - 2)
!         write(*,*) 'wat-wat',o1,o2,o3
#ifdef HEAT_CURRENT
         call do_trimer(imcon,cell, &
          & monomers(9*(i-1)+1),monomers(9*(j-1)+1),monomers(9*(k-1)+1), &
          & fxx(o1),fyy(o1),fzz(o1), &
          & fxx(o2),fyy(o2),fzz(o2), &
          & fxx(o3),fyy(o3),fzz(o3), &
          & eng3,vir3,strs,o1,o2,o3)
#else
          call do_trimer(imcon,cell, &
           & monomers(9*(i-1)+1),monomers(9*(j-1)+1),monomers(9*(k-1)+1), &
           & fxx(o1),fyy(o1),fzz(o1), &
           & fxx(o2),fyy(o2),fzz(o2), &
           & fxx(o3),fyy(o3),fzz(o3), &
           & eng3,vir3,strs)
#endif /* HEAT_CURRENT */

         stress(1) = stress(1) + strs(1) !FP: added
         stress(2) = stress(2) + strs(2) !FP: added
         stress(3) = stress(3) + strs(3) !FP: added
         stress(4) = stress(4) + strs(4) !FP: added
         stress(5) = stress(5) + strs(5) !FP: added
         stress(6) = stress(6) + strs(6) !FP: added
         stress(7) = stress(7) + strs(7) !FP: added
         stress(8) = stress(8) + strs(8) !FP: added
         stress(9) = stress(9) + strs(9) !FP: added
!      end if ! mod(t, bead_size).eq.bead_rank

!      k = k+1
!      if (k.gt.nw) then
!         j = j+1
!         if ((j+1).gt.nw) then
!            i = i+1
!            j = i+1
!         end if
!         k = j+1
!      end if ! k.gt.nw

      k = k + bead_size
      if (k > nw) then
         do while (k > nw) ! do while() only needed near end of triplet list
            j = j + 1
            if ((j+1) > nw) then
               i = i+1
               j = i+1
            end if
            k = j + (k-nw)
         end do
      end if ! k.gt.nw

   end do ! t = 1, nw3


#ifdef DEBUG
      call timchk(0,timeint)
      write(4444,"(a,F10.3,a)") 'Time taken for MB-pol 3B       = ',timeint-timetot,' sec'
      write(4444,*)
#endif

#  endif /* ENABLE_X3B */

   !
   ! reduce energy/virial and return
   !

   if (bead_size.gt.1) then
      buffer(5) = eng2
      buffer(6) = vir2
      buffer(7) = eng3
      buffer(8) = vir3

      !buffer(14) = eng2
      !buffer(15) = vir2
      !buffer(16) = eng3
      !buffer(17) = vir3
      !buffer(18) = stress(1)
      !buffer(19) = stress(2)
      !buffer(20) = stress(3)
      !buffer(21) = stress(4)
      !buffer(22) = stress(5)
      !buffer(23) = stress(6)
      !buffer(24) = stress(7)
      !buffer(25) = stress(8)
      !buffer(26) = stress(9)

      call gdsum(buffer(5), 4, buffer(1))
      !call gdsum(buffer(14), 13, buffer(1))

      eng2 = buffer(1)
      vir2 = buffer(2)
      eng3 = buffer(3)
      vir3 = buffer(4)
      !stress(1) = buffer(5)
      !stress(2) = buffer(6)
      !stress(3) = buffer(7)
      !stress(4) = buffer(8)
      !stress(5) = buffer(9)
      !stress(6) = buffer(10)
      !stress(7) = buffer(11)
      !stress(8) = buffer(12)
      !stress(9) = buffer(13)
   endif


   engsrp = engsrp+eng2+eng3
   virsrp = virsrp+vir2+vir3

#ifdef DEBUG
 if(bead_rank==0) then
   write(4444,*) 'MB-Pol polynomial contributions to the potential energy'
   write(4444,*) 'water-water 2B',eng2/engunit
   write(4444,*) 'water-water 3B',eng3/engunit
 endif
#endif

24050 continue

!   if(bead_rank.eq.0)write(200,*)eng2,eng3

end subroutine mbpol_forces

!----------------------------------------------------------------------------
#ifdef HEAT_CURRENT
subroutine do_dimer(imcon,cell,w1,w2,fx1,fy1,fz1,fx2,fy2,fz2,engacc,viracc,strs,o1,o2)
#else
subroutine do_dimer(imcon,cell,w1,w2,fx1,fy1,fz1,fx2,fy2,fz2,engacc,viracc,strs)
#endif /*HEAT_CURRENT*/


   implicit none

   integer, intent(in) :: imcon

#ifdef HEAT_CURRENT
   integer, intent(in) :: o1, o2
#endif /* HEAT_CURRENT */

   real(8), intent(in) :: cell(*)

   real(8), intent(inout) :: w1(*)
   real(8), intent(inout) :: w2(*)

   real(8), intent(inout) :: fx1(*),fy1(*),fz1(*)
   real(8), intent(inout) :: fx2(*),fy2(*),fz2(*)

   real(8), intent(inout) :: engacc,viracc,strs(9) !FP: added

   integer :: k

   real(8) :: r12(3),g1(9),g2(9),e2

!FP: added
#ifdef STRESS
   ! initialize stress tensor accumulators
   strs(1:9) = 0.d0
#endif

   r12(1) = w1(1) - w2(1)
   r12(2) = w1(2) - w2(2)
   r12(3) = w1(3) - w2(3)

   call image_vectors(imcon,cell,1,r12)

   if ((r12(1)**2 + r12(2)**2 + r12(3)**2).gt.r2f**2) &
   &   return

   do k = 1, 3
      w1(3+k) = (w1(3+k) - w1(k)) + (w2(k) + r12(k))
      w1(6+k) = (w1(6+k) - w1(k)) + (w2(k) + r12(k))
      w1(k) = w2(k) + r12(k)
   end do

   call mbpol_2b_poly(w1,w2,e2,g1,g2)

   engacc = engacc + e2*engunit

#ifdef HEAT_CURRENT
   call update_energy_mbpol(o1,0.5*e2*engunit)
   call update_energy_mbpol(o2,0.5*e2*engunit)
#endif /* HEAT_CURRENT */

   do k = 1, 9
      g1(k) = g1(k)*engunit
      viracc = viracc + w1(k)*g1(k)
      g2(k) = g2(k)*engunit
      viracc = viracc + w2(k)*g2(k)
   end do

   ! forces

   do k = 1, 3
      fx1(k) = fx1(k) - g1(3*k - 2)
      fy1(k) = fy1(k) - g1(3*k - 1)
      fz1(k) = fz1(k) - g1(3*k - 0)

      fx2(k) = fx2(k) - g2(3*k - 2)
      fy2(k) = fy2(k) - g2(3*k - 1)
      fz2(k) = fz2(k) - g2(3*k - 0)
   end do

#ifdef STRESS

   ! calculate stress tensor

   strs(1) = w1(1) * g1(1) &
           + w1(4) * g1(4) &
           + w1(7) * g1(7) &
           + w2(1) * g2(1) &
           + w2(4) * g2(4) &
           + w2(7) * g2(7)

   strs(2) = w1(1) * g1(2) &
           + w1(4) * g1(5) &
           + w1(7) * g1(8) &
           + w2(1) * g2(2) &
           + w2(4) * g2(5) &
           + w2(7) * g2(8)

   strs(3) = w1(1) * g1(3) &
           + w1(4) * g1(6) &
           + w1(7) * g1(9) &
           + w2(1) * g2(3) &
           + w2(4) * g2(6) &
           + w2(7) * g2(9)

   strs(4) = strs(2)

   strs(5) = w1(2) * g1(2) &
           + w1(5) * g1(5) &
           + w1(8) * g1(8) &
           + w2(2) * g2(2) &
           + w2(5) * g2(5) &
           + w2(8) * g2(8)

   strs(6) = w1(2) * g1(3) &
           + w1(5) * g1(6) &
           + w1(8) * g1(9) &
           + w2(2) * g2(3) &
           + w2(5) * g2(6) &
           + w2(8) * g2(9)

   strs(7) = strs(3)

   strs(8) = strs(6)

   strs(9) = w1(3) * g1(3) &
           + w1(6) * g1(6) &
           + w1(9) * g1(9) &
           + w2(3) * g2(3) &
           + w2(6) * g2(6) &
           + w2(9) * g2(9)

   strs(:) = -strs(:)
#endif

#ifdef HEAT_CURRENT
  call update_stress_mbpol(o1,1,1,0.5d0*strs(1))
  call update_stress_mbpol(o1,1,2,0.5d0*strs(2))
  call update_stress_mbpol(o1,1,3,0.5d0*strs(3))
  call update_stress_mbpol(o1,2,1,0.5d0*strs(4))
  call update_stress_mbpol(o1,2,2,0.5d0*strs(5))
  call update_stress_mbpol(o1,2,3,0.5d0*strs(6))
  call update_stress_mbpol(o1,3,1,0.5d0*strs(7))
  call update_stress_mbpol(o1,3,2,0.5d0*strs(8))
  call update_stress_mbpol(o1,3,3,0.5d0*strs(9))
  call update_stress_mbpol(o2,1,1,0.5d0*strs(1))
  call update_stress_mbpol(o2,1,2,0.5d0*strs(2))
  call update_stress_mbpol(o2,1,3,0.5d0*strs(3))
  call update_stress_mbpol(o2,2,1,0.5d0*strs(4))
  call update_stress_mbpol(o2,2,2,0.5d0*strs(5))
  call update_stress_mbpol(o2,2,3,0.5d0*strs(6))
  call update_stress_mbpol(o2,3,1,0.5d0*strs(7))
  call update_stress_mbpol(o2,3,2,0.5d0*strs(8))
  call update_stress_mbpol(o2,3,3,0.5d0*strs(9))
#endif /* HEAT_CURRENT */

end subroutine do_dimer

!----------------------------------------------------------------------------

#ifdef HEAT_CURRENT
subroutine do_trimer(imcon,cell,w1,w2,w3, &
 & fx1,fy1,fz1,fx2,fy2,fz2,fx3,fy3,fz3,engacc,viracc,strs,o1,o2,o3)
#else
subroutine do_trimer(imcon,cell,w1,w2,w3, &
 & fx1,fy1,fz1,fx2,fy2,fz2,fx3,fy3,fz3,engacc,viracc,strs)
#endif /* HEAT_CURRENT */

   implicit none

#ifdef HEAT_CURRENT
integer, intent(in) :: o1,o2,o3
#endif
   integer, intent(in) :: imcon
   real(8), intent(in) :: cell(*)

   real(8), intent(inout) :: w1(*)
   real(8), intent(inout) :: w2(*)
   real(8), intent(inout) :: w3(*)

   real(8), intent(inout) :: fx1(*),fy1(*),fz1(*)
   real(8), intent(inout) :: fx2(*),fy2(*),fz2(*)
   real(8), intent(inout) :: fx3(*),fy3(*),fz3(*)

   real(8), intent(inout) :: engacc,viracc,strs(9)

   real(8) :: rab(9),r12,r23,r13
   real(8) :: e3b,g1(9),g2(9),g3(9)

   logical :: istoobig1, istoobig2, istoobig3

   integer :: k

!FP: added
#ifdef STRESS
   ! initialize stress tensor accumulators
   strs(1:9) = 0.d0
#endif

   ! O-O distances

   do k = 1, 3
     rab(k + 0) = w1(k) - w2(k)
     rab(k + 3) = w2(k) - w3(k)
     rab(k + 6) = w1(k) - w3(k)
   end do

   call image_vectors(imcon,cell,3,rab)

   r12 = 0.d0
   r23 = 0.d0
   r13 = 0.d0

   do k = 1, 3
      r12 = r12 + rab(k + 0)**2
      r23 = r23 + rab(k + 3)**2
      r13 = r13 + rab(k + 6)**2
   end do

   r12 = sqrt(r12)
   r23 = sqrt(r23)
   r13 = sqrt(r13)

   istoobig1 = r12.gt.r3f.or.r13.gt.r3f
   istoobig2 = r12.gt.r3f.or.r23.gt.r3f
   istoobig3 = r13.gt.r3f.or.r23.gt.r3f

   if (istoobig1.and.istoobig2.and.istoobig3) return

   ! image the trimer [monomers are already imaged]
   ! w1 = w2 + rab(1:3)
   ! w3 = w2 - rab(4:6)

   do k = 1, 3
      w1(3+k) = (w1(3+k) - w1(k)) + (w2(k) + rab(k))
      w1(6+k) = (w1(6+k) - w1(k)) + (w2(k) + rab(k))
      w1(k) = w2(k) + rab(k)
      w3(3+k) = (w3(3+k) - w3(k)) + (w2(k) - rab(k+3))
      w3(6+k) = (w3(6+k) - w3(k)) + (w2(k) - rab(k+3))
      w3(k) = w2(k) - rab(k+3)
   end do

   g1(1:9) = 0.d0
   g2(1:9) = 0.d0
   g3(1:9) = 0.d0

   call mbpol_3b_poly(w1,w2,w3,e3b,g1,g2,g3)

   engacc = engacc + e3b*engunit
#ifdef HEAT_CURRENT
    call update_energy_mbpol(o1,third*e3b*engunit)
    call update_energy_mbpol(o2,third*e3b*engunit)
    call update_energy_mbpol(o3,third*e3b*engunit)
#endif

!   if (e3b.lt.-3.d0) &
!       & call drop_trimer(e3b,w1,w2,w3)

   do k = 1, 9
      g1(k) = g1(k)*engunit
      viracc = viracc + w1(k)*g1(k)
      g2(k) = g2(k)*engunit
      viracc = viracc + w2(k)*g2(k)
      g3(k) = g3(k)*engunit
      viracc = viracc + w3(k)*g3(k)
   end do

   ! forces

   do k = 1, 3
      fx1(k) = fx1(k) - g1(3*k - 2)
      fy1(k) = fy1(k) - g1(3*k - 1)
      fz1(k) = fz1(k) - g1(3*k - 0)

      fx2(k) = fx2(k) - g2(3*k - 2)
      fy2(k) = fy2(k) - g2(3*k - 1)
      fz2(k) = fz2(k) - g2(3*k - 0)

      fx3(k) = fx3(k) - g3(3*k - 2)
      fy3(k) = fy3(k) - g3(3*k - 1)
      fz3(k) = fz3(k) - g3(3*k - 0)
   end do

#ifdef STRESS

   ! calculate stress tensor

   strs(1) = w1(1) * g1(1) &
           + w1(4) * g1(4) &
           + w1(7) * g1(7) &
           + w2(1) * g2(1) &
           + w2(4) * g2(4) &
           + w2(7) * g2(7) &
           + w3(1) * g3(1) &
           + w3(4) * g3(4) &
           + w3(7) * g3(7)

   strs(2) = w1(1) * g1(2) &
           + w1(4) * g1(5) &
           + w1(7) * g1(8) &
           + w2(1) * g2(2) &
           + w2(4) * g2(5) &
           + w2(7) * g2(8) &
           + w3(1) * g3(2) &
           + w3(4) * g3(5) &
           + w3(7) * g3(8)

   strs(3) = w1(1) * g1(3) &
           + w1(4) * g1(6) &
           + w1(7) * g1(9) &
           + w2(1) * g2(3) &
           + w2(4) * g2(6) &
           + w2(7) * g2(9) &
           + w3(1) * g3(3) &
           + w3(4) * g3(6) &
           + w3(7) * g3(9)

   strs(4) = strs(2)

   strs(5) = w1(2) * g1(2) &
           + w1(5) * g1(5) &
           + w1(8) * g1(8) &
           + w2(2) * g2(2) &
           + w2(5) * g2(5) &
           + w2(8) * g2(8) &
           + w3(2) * g3(2) &
           + w3(5) * g3(5) &
           + w3(8) * g3(8)

   strs(6) = w1(2) * g1(3) &
           + w1(5) * g1(6) &
           + w1(8) * g1(9) &
           + w2(2) * g2(3) &
           + w2(5) * g2(6) &
           + w2(8) * g2(9) &
           + w3(2) * g3(3) &
           + w3(5) * g3(6) &
           + w3(8) * g3(9)

   strs(7) = strs(3)

   strs(8) = strs(6)

   strs(9) = w1(3) * g1(3) &
           + w1(6) * g1(6) &
           + w1(9) * g1(9) &
           + w2(3) * g2(3) &
           + w2(6) * g2(6) &
           + w2(9) * g2(9) &
           + w3(3) * g3(3) &
           + w3(6) * g3(6) &
           + w3(9) * g3(9)

   strs(:) = -strs(:)

#endif

#ifdef HEAT_CURRENT
call update_stress_mbpol(o1,1,1,third*strs(1))
call update_stress_mbpol(o1,1,2,third*strs(2))
call update_stress_mbpol(o1,1,3,third*strs(3))
call update_stress_mbpol(o1,2,1,third*strs(4))
call update_stress_mbpol(o1,2,2,third*strs(5))
call update_stress_mbpol(o1,2,3,third*strs(6))
call update_stress_mbpol(o1,3,1,third*strs(7))
call update_stress_mbpol(o1,3,2,third*strs(8))
call update_stress_mbpol(o1,3,3,third*strs(9))
call update_stress_mbpol(o2,1,1,third*strs(1))
call update_stress_mbpol(o2,1,2,third*strs(2))
call update_stress_mbpol(o2,1,3,third*strs(3))
call update_stress_mbpol(o2,2,1,third*strs(4))
call update_stress_mbpol(o2,2,2,third*strs(5))
call update_stress_mbpol(o2,2,3,third*strs(6))
call update_stress_mbpol(o2,3,1,third*strs(7))
call update_stress_mbpol(o2,3,2,third*strs(8))
call update_stress_mbpol(o2,3,3,third*strs(9))
call update_stress_mbpol(o3,1,1,third*strs(1))
call update_stress_mbpol(o3,1,2,third*strs(2))
call update_stress_mbpol(o3,1,3,third*strs(3))
call update_stress_mbpol(o3,2,1,third*strs(4))
call update_stress_mbpol(o3,2,2,third*strs(5))
call update_stress_mbpol(o3,2,3,third*strs(6))
call update_stress_mbpol(o3,3,1,third*strs(7))
call update_stress_mbpol(o3,3,2,third*strs(8))
call update_stress_mbpol(o3,3,3,third*strs(9))
#endif /* HEAT_CURRENT */

end subroutine do_trimer

!----------------------------------------------------------------------------

subroutine image_monomers(imcon, cell, nw, xyz)

   implicit none

   integer, intent(in) :: imcon
   real(8), intent(in) :: cell(*)
   integer, intent(in) :: nw

   real(8), intent(inout) :: xyz(*)

   integer :: i, i9, k
   real(8) :: oh(6)

   do i = 0, nw - 1
      i9 = 9*i
      do k = 1, 3
         oh(k)   = xyz(i9+3+k) - xyz(i9+k)
         oh(k+3) = xyz(i9+6+k) - xyz(i9+k)
      end do
      call image_vectors(imcon,cell,2,oh)
      do k = 1, 3
         xyz(i9+3+k) = xyz(i9+k) + oh(k)
         xyz(i9+6+k) = xyz(i9+k) + oh(k+3)
      end do
   end do

end subroutine image_monomers

!----------------------------------------------------------------------------

subroutine image_vectors(imcon,cell,n,dx)

   implicit none

   integer, intent(in) :: imcon
   real(8), intent(in) :: cell(*)
   integer, intent(in) :: n

   real(8), intent(inout) :: dx(*)

   integer :: i, i3
   real(8) :: a1, b1, c1
   real(8) :: ssx, ssy, ssz, xss, yss, zss, det
   real(8) :: rcell(9)

   if (imcon.eq.1) then
      a1 = 1.d0/cell(1)
      do i = 1, n
         i3 = 3*i-2
         dx(i3)   = dx(i3)   - cell(1)*nint(dx(i3)*a1)
         dx(i3+1) = dx(i3+1) - cell(1)*nint(dx(i3+1)*a1)
         dx(i3+2) = dx(i3+2) - cell(1)*nint(dx(i3+2)*a1)
      end do
   else if (imcon.eq.2) then
      a1 = 1.d0/cell(1)
      b1 = 1.d0/cell(5)
      c1 = 1.d0/cell(9)
      do i = 1, n
         i3 = 3*i-2
         dx(i3)   = dx(i3)   - cell(1)*nint(dx(i3)*a1)
         dx(i3+1) = dx(i3+1) - cell(5)*nint(dx(i3+1)*b1)
         dx(i3+2) = dx(i3+2) - cell(9)*nint(dx(i3+2)*c1)
      end do
   else if (imcon.eq.3) then

      call invert(cell,rcell,det)
      do i = 1, n
         i3 = 3*i-2

         ssx=(rcell(1)*dx(i3)+rcell(4)*dx(i3+1)+rcell(7)*dx(i3+2))
         ssy=(rcell(2)*dx(i3)+rcell(5)*dx(i3+1)+rcell(8)*dx(i3+2))
         ssz=(rcell(3)*dx(i3)+rcell(6)*dx(i3+1)+rcell(9)*dx(i3+2))

         xss=ssx-nint(ssx)
         yss=ssy-nint(ssy)
         zss=ssz-nint(ssz)

         dx(i3  )=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
         dx(i3+1)=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
         dx(i3+2)=(cell(3)*xss+cell(6)*yss+cell(9)*zss)

       end do
   else
       call fixme(__FILE__,__LINE__)
   end if

end subroutine image_vectors

subroutine drop_trimer(e3b,w1,w2,w3)

   use multibead, only: bead_rank, bead_size

   implicit none

   real(8), intent(in) :: e3b

   real(8), intent(in) :: w1(*)
   real(8), intent(in) :: w2(*)
   real(8), intent(in) :: w3(*)

#  if 0
   write(500+bead_rank,*) '9'
   write(500+bead_rank,*) e3b
   write(500+bead_rank,*) 'O', w1(1:3)
   write(500+bead_rank,*) 'H', w1(4:6)
   write(500+bead_rank,*) 'H', w1(7:9)
   write(500+bead_rank,*) 'O', w2(1:3)
   write(500+bead_rank,*) 'H', w2(4:6)
   write(500+bead_rank,*) 'H', w2(7:9)
   write(500+bead_rank,*) 'O', w3(1:3)
   write(500+bead_rank,*) 'H', w3(4:6)
   write(500+bead_rank,*) 'H', w3(7:9)
#  else
   write(500+bead_rank,*) e3b
#  endif

end subroutine drop_trimer

!----------------------------------------------------------------------------

end module mbpol

!============================================================================
