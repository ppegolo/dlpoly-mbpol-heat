#include "assert.h"

#define ENABLE_X2B goahead
!#define ENABLE_X3B aswell

!============================================================================

module mbnrg

!----------------------------------------------------------------------------

implicit none

real(8), private, save :: r2f = 0.d0
real(8), private, save :: r3f = 0.d0

!----------------------------------------------------------------------------

public :: mbnrg_init
public :: mbnrg_fini

public :: mbnrg_forces ! does not MPI_REDUCE the forces

!----------------------------------------------------------------------------

real(8), private, allocatable, save :: water_monomers(:),ions(:)
integer, private, allocatable, save :: polynomol(:),selectpoly(:)

private :: do_dimer, do_trimer
private :: image_vectors, image_water_monomers

real(8), private, parameter :: engunit = 418.4d0

!============================================================================

contains

!============================================================================

subroutine mbnrg_init(nrite, nttm2, nions, mxtmls, key, nummols)

   use multibead, only: bead_rank, mb_abort

   implicit none

   integer, intent(in) :: nrite, nttm2, nions, mxtmls

   integer,intent(in)  :: key(mxtmls),nummols(mxtmls)

   integer :: ierr, nw,itmp,i,j

   assert(.not.allocated(water_monomers))
   assert(.not.allocated(ions))
 
   allocate(polynomol(mxtmls),selectpoly(nions))

   polynomol=nummols
   
  itmp=0 
  do i=1,mxtmls     ! types 

     do j=1,nummols(i)   ! no of ions in each type 

        itmp=itmp+1 

        if(itmp>nions) goto 4500

        selectpoly(itmp)=key(i)
 
      enddo

  enddo   

  4500 continue 

! SR: activate r3f line when 3B polynomials are ready 

! MRR & DZ - Add multion mbnrg
! Set r2f to 10.0, so we avoid problems 
!   select case(polykey)
!
!      case(1) 
!
!         call  mbnrg_2b_h2o_f_cutoff(r2f)
!!         call  mbnrg_3b_h2o_f_cutoff(r3f)
!
!      case(2) 
!
!         call  mbnrg_2b_h2o_cl_cutoff(r2f)
!!         call  mbnrg_3b_h2o_cl_cutoff(r3f)
!      
!      case(3) 
!
!         call  mbnrg_2b_h2o_br_cutoff(r2f)
!!         call  mbnrg_3b_h2o_br_cutoff(r3f)
!     
!      case(4) 
!
!         call  mbnrg_2b_h2o_i_cutoff(r2f)
!!         call  mbnrg_3b_h2o_ion_cutoff(r3f)
!!         call  mbnrg_3b_h2o_i_cutoff(r3f)
!    
!      case(5) 
!
!         call  mbnrg_2b_h2o_li_cutoff(r2f)
!!         call  mbnrg_3b_h2o_f_cutoff(r3f)
!
!      case(6) 
!
!         call  mbnrg_2b_h2o_na_cutoff(r2f)
!!         call  mbnrg_3b_h2o_f_cutoff(r3f)
!
!      case(7) 
!
!         call  mbnrg_2b_h2o_k_cutoff(r2f)
!!         call  mbnrg_3b_h2o_f_cutoff(r3f)
!
!      case(8) 
!
!         call  mbnrg_2b_h2o_rb_cutoff(r2f)
!!         call  mbnrg_3b_h2o_f_cutoff(r3f)
!
!      case(9) 
!
!         call  mbnrg_2b_h2o_cs_cutoff(r2f)
!!         call  mbnrg_3b_h2o_f_cutoff(r3f)
!
!
!    end select 
   r2f = 10.0

! END MRR & FP
!DZ: set r3f to 7 for iodide

   r3f = 7.0
   

   if (bead_rank.eq.0) then
      write(nrite,*) ' ** Using MB-nrg for ion-water interactions **'
#     ifdef ENABLE_X2B
      write(nrite,*) ' ** 2B cutoff = ', r2f
#     endif
#     ifdef ENABLE_X3B
      write(nrite,*) ' ** 3B cutoff = ', r3f
#     endif
   end if

   nw = nttm2/3

   allocate(water_monomers(9*nw), stat=ierr)
   allocate(ions(3*nions), stat=ierr)
   if (ierr.ne.0) then
      write(nrite,*) ' ** Error ** : could not allocate memory'
      call mb_abort()
   end if

end subroutine mbnrg_init

!----------------------------------------------------------------------------

subroutine mbnrg_fini()

   implicit none
   deallocate(water_monomers,ions)

end subroutine mbnrg_fini

!----------------------------------------------------------------------------

subroutine mbnrg_forces(imcon,cell,nttm2,listttm2,nions,listions, &
   & xxx,yyy,zzz,fxx,fyy,fzz,engsrp,virsrp,stress) !FP: added

   use multibead, only: bead_rank, bead_size

   implicit none

   integer, intent(in) :: imcon
   real(8), intent(in) :: cell(*)
   integer, intent(in) :: nttm2, listttm2(:)
   integer, intent(in) :: nions, listions(:)
   real(8), intent(in) :: xxx(*), yyy(*), zzz(*)

   real(8), intent(inout) :: fxx(*), fyy(*), fzz(*)
   real(8), intent(inout) :: engsrp, virsrp, stress(9) !FP: added

   real(8) :: eng2, vir2
   real(8) :: eng3, vir3
   real(8) :: strs(9) !FP: added
   real(8) :: wat_monomer_ceiling  !SR
   integer :: t1,t2,l           ! SR 

   real(8) :: buffer(8)

   integer :: i, j, k, o1, o2, o3,tmp_nw,tmppoly
   integer(kind=8) :: nw, nw2, nw3,t,tmp_bead_size
   integer(kind=8) :: istart, iend, istep

   
   !SR    ceiling value --  set to 60 kcal/mol 

!   wat_monomer_ceiling=60          ! kcal/mol
!   wat_monomer_ceiling=wat_monomer_ceiling*engunit          ! internal units

   !
   ! copy waters into the water_monomers(:) and image every molecule
   !

!SR adding ceiling for mbnrg potential mono energ
! but it is not helpful for optimization to skip water_monomers if energy is
! greater than ceiling value
! so fixing it temporarily


   nw = nttm2/3

 do t = 1, nttm2
    i = listttm2(t)
    water_monomers(3*t - 2) = xxx(i)
    water_monomers(3*t - 1) = yyy(i)
    water_monomers(3*t - 0) = zzz(i)
 end do


!!DZ: ions are always listed first, so we can loop through the ions in the
!beginning
 do t = 1, nions
!    i = listions(t)
!    ions(3*t - 2) = xxx(i)
!    ions(3*t - 1) = yyy(i)
!    ions(3*t - 0) = zzz(i)
    ions(3*t - 2) = xxx(t)
    ions(3*t - 1) = yyy(t)
    ions(3*t - 0) = zzz(t)
 end do



  tmp_nw=nw          ! to pass same integer kind


  call image_water_monomers(imcon, cell, tmp_nw, water_monomers)

   !
   ! loop over the pairs
   !

   nw2 = nw*nions

   eng2 = 0.d0
   vir2 = 0.d0

#  ifdef ENABLE_X2B

!   i = 1
!   j = 1 + bead_rank
!   write(*,*) 'bead',bead_rank 

  tmp_bead_size=bead_size
  istart=bead_rank*(nw/bead_size)+1+bead_rank
  iend=istart+(nw/bead_size)
  if(bead_rank==bead_size-1) then
      iend=nw
  endif

  do j=1,nions

!    tmppoly=selectpoly(j)

!   istart = bead_rank*(bead_rank/bead_size)
!   iend   = nw*(bead_rank+1)/bead_size
!   istep  = bead_size
!    write(*,*) 'test',istart,iend,j,bead_rank
! Debbie: debug
   do i = istart, iend
!   do t = 1, nw2
!      if (mod(t,tmp_bead_size).eq.bead_rank) then
         o1 = listttm2(3*i - 2)
         o2 = listions(j)
!         o2 = listttm2(3*j - 2)
!         call do_dimer(tmppoly,imcon,cell, &
!              & water_monomers(9*(i-1)+1),ions(3*(j-1)+1), &
!              & fxx(o1),fyy(o1),fzz(o1),fxx(o2),fyy(o2),fzz(o2), &
!              & eng2,vir2,strs) !FP: added
          call do_dimer(selectpoly(j),0,imcon,cell, &
              & water_monomers(9*(i-1)+1),ions(3*(j-1)+1), &
              & fxx(o1),fyy(o1),fzz(o1),fxx(o2),fyy(o2),fzz(o2), &
              & eng2,vir2,strs) !DZ: modified for ion-ion interactions
!         call do_dimer(imcon,cell, &
!              & water_monomers(9*(i-1)+1),water_monomers(9*(j-1)+1), &
!              & fxx(o1),fyy(o1),fzz(o1),fxx(o2),fyy(o2),fzz(o2), &
!              & eng2,vir2,strs) !FP: added
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

!     j = j + bead_size 
!     if (j > nw) then
!        do while (j > nw) ! -- only needed near end of pair list
!           i = i + 1
!           j = i + j - nw
!        end do
!     end if
!Debbie: debug
     end do ! t = 1, nw2
     
     do k = j+1,nions
!!add other stuff
       o1 = listions(j)
       o2 = listions(k)
       call do_dimer(selectpoly(j),selectpoly(k),imcon,cell, &
              & ions(3*(j-1)+1),ions(3*(k-1)+1), &
              & fxx(o1),fyy(o1),fzz(o1),fxx(o2),fyy(o2),fzz(o2), &
              & eng2,vir2,strs)
       stress(1) = stress(1) + strs(1) !FP: added
       stress(2) = stress(2) + strs(2) !FP: added
       stress(3) = stress(3) + strs(3) !FP: added
       stress(4) = stress(4) + strs(4) !FP: added
       stress(5) = stress(5) + strs(5) !FP: added
       stress(6) = stress(6) + strs(6) !FP: added
       stress(7) = stress(7) + strs(7) !FP: added
       stress(8) = stress(8) + strs(8) !FP: added
       stress(9) = stress(9) + strs(9) !FP: added

     end do 
    !! add loop that goes through the remaining (nions-j) ions to find the
    !potential, do i write another function or do i modify the original
    !function?
   end do ! t = 1, nw2

#  endif /* ENABLE_X2B */

   !
   ! loop over the triples
   !

   nw3 = (nw*(nw - 1)*nions)/2

   eng3 = 0.d0
   vir3 = 0.d0

#  ifdef ENABLE_X3B

!   i = 1
!   j = 2
!   k = 1 + bead_rank         ! for ion 

!   tmp_bead_size=bead_size

  tmp_bead_size=bead_size
  istart=bead_rank*(nw/bead_size)+1+bead_rank
  iend=istart+(nw/bead_size)
  if(bead_rank==bead_size-1) then
    iend=nw
  endif

  do k=1,nions

   do i = istart, iend

     do j = i+1, nw
!   do t = 1, nw3
!      if (mod(t, tmp_bead_size).eq.bead_rank) then
!        write(*,*) 'test',i,j,k
         o1 = listttm2(3*i - 2)
         o2 = listttm2(3*j - 2)
         o3 = listions(k)
!         o3 = listttm2(3*k - 2)
! write(*,*) 'index',o1,o2,o3
! DZ (11/30/17):
         call do_trimer(selectpoly(k), imcon,cell, &
          & water_monomers(9*(i-1)+1),water_monomers(9*(j-1)+1),ions(3*(k-1)+1), &
          & fxx(o1),fyy(o1),fzz(o1), &
          & fxx(o2),fyy(o2),fzz(o2), &
          & fxx(o3),fyy(o3),fzz(o3), &
          & eng3,vir3,strs)
!  write(24000,*) 'water-ion ',eng3/engunit 
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

!      j = j + bead_size 
!      k = k + bead_size 
!      if (k > nw) then
!         do while (j > nw) ! do while() only needed near end of triplet list
!            i = i + 1
!            if ((j+1) > nw) then
!               i = i+1
!               j = i+j-nw
!            end if
!            k = j + (k-nw)
         end do
      end do
!      end if ! k.gt.nw

   end do ! t = 1, nw3

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
   write(4444,*) 'MB-nrg PES polynomial contributions to the potential energy'
   write(4444,*) 'wat-ion 2B',eng2/engunit
   write(4444,*) 'wat-ion 3B',eng3/engunit
 endif 
#endif 
24050 continue 

!   if(bead_rank.eq.0)write(200,*)eng2,eng3

end subroutine mbnrg_forces

!----------------------------------------------------------------------------

subroutine do_dimer(tmppoly1,tmppoly2,imcon,cell,w1,w2,fx1,fy1,fz1,fx2,fy2,fz2,engacc,viracc,strs)

   implicit none

   integer, intent(in) :: imcon,tmppoly1,tmppoly2
   real(8), intent(in) :: cell(*)

   real(8), intent(inout) :: w1(*)
   real(8), intent(inout) :: w2(*)

   real(8), intent(inout) :: fx1(*),fy1(*),fz1(*)
   real(8), intent(inout) :: fx2(*),fy2(*),fz2(*)

   real(8), intent(inout) :: engacc,viracc,strs(9) !FP: added

   integer :: k

   real(8) :: r12(3),g1(9),g2(3),e2



!SR: for some reason, though g2 variable not initialized to zero, it worked on MB-pol.
! but the same is not working for MB-nrg. So, added g2 below

   g2=0.d0 ; g1=0.d0 

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
      if (tmppoly2 == 0) then
          w1(3+k) = (w1(3+k) - w1(k)) + (w2(k) + r12(k))
!          write(*,*) "w1(3+k) = ", w1(3+k)
          w1(6+k) = (w1(6+k) - w1(k)) + (w2(k) + r12(k))
!          write(*,*) "w1(6+k) = ", w1(6+k)
      endif
      w1(k) = w2(k) + r12(k)
!      write(*,*) "w1(k) = ", w1(k)
   end do
   
  !Debbie:

   

   if (tmppoly2 == 0) then
      select case(tmppoly1)

         case(1) 
            call mbnrg_2b_h2o_f_poly(w1,w2,e2,g1,g2)
         case(2)
            call mbnrg_2b_h2o_cl_poly(w1,w2,e2,g1,g2)
         case(3) 
            call mbnrg_2b_h2o_br_poly(w1,w2,e2,g1,g2)
         case(4) 
            call mbnrg_2b_h2o_i_poly(w1,w2,e2,g1,g2)
         case(5) 
            call mbnrg_2b_h2o_li_poly(w1,w2,e2,g1,g2)
         case(6)
            call mbnrg_2b_h2o_na_poly(w1,w2,e2,g1,g2)
         case(7) 
            call mbnrg_2b_h2o_k_poly(w1,w2,e2,g1,g2)
         case(8) 
            call mbnrg_2b_h2o_rb_poly(w1,w2,e2,g1,g2)
         case(9) 
            call mbnrg_2b_h2o_cs_poly(w1,w2,e2,g1,g2)

       end select 
   else if (tmppoly2 == 1) then
      select case(tmppoly1)

         case(1)
            call mbnrg_2b_f_f_poly(w1,w2,e2,g1,g2)
         case(5)
            call mbnrg_2b_li_f_poly(w1,w2,e2,g1,g2)
         case(6)
            call mbnrg_2b_na_f_poly(w1,w2,e2,g1,g2)
         case(7)
            call mbnrg_2b_k_f_poly(w1,w2,e2,g1,g2)
         case(8)
            call mbnrg_2b_rb_f_poly(w1,w2,e2,g1,g2)
         case(9)
            call mbnrg_2b_cs_f_poly(w1,w2,e2,g1,g2)

       end select
   elseif (tmppoly2 == 2) then
      select case(tmppoly1)

         case(2)
            call mbnrg_2b_cl_cl_poly(w1,w2,e2,g1,g2)
         case(5)
            call mbnrg_2b_li_cl_poly(w1,w2,e2,g1,g2)
         case(6)
            call mbnrg_2b_na_cl_poly(w1,w2,e2,g1,g2)
         case(7)
            call mbnrg_2b_k_cl_poly(w1,w2,e2,g1,g2)
         case(8)
            call mbnrg_2b_rb_cl_poly(w1,w2,e2,g1,g2)
         case(9)
            call mbnrg_2b_cs_cl_poly(w1,w2,e2,g1,g2)

       end select
   elseif (tmppoly2 == 3) then
      select case(tmppoly1)
         case(3)
            call mbnrg_2b_br_br_poly(w1,w2,e2,g1,g2)
         case(5)
            call mbnrg_2b_li_br_poly(w1,w2,e2,g1,g2)
         case(6)
            call mbnrg_2b_na_br_poly(w1,w2,e2,g1,g2)
         case(7)
            call mbnrg_2b_k_br_poly(w1,w2,e2,g1,g2)
         case(8)
            call mbnrg_2b_rb_br_poly(w1,w2,e2,g1,g2)
         case(9)
            call mbnrg_2b_cs_br_poly(w1,w2,e2,g1,g2)

       end select
   elseif (tmppoly2 == 4) then
      select case(tmppoly1)

         case(4)
            call mbnrg_2b_i_i_poly(w1,w2,e2,g1,g2)
         case(5)
            call mbnrg_2b_li_i_poly(w1,w2,e2,g1,g2)
         case(6)
            call mbnrg_2b_na_i_poly(w1,w2,e2,g1,g2)
         case(7)
            call mbnrg_2b_k_i_poly(w1,w2,e2,g1,g2)
         case(8)
            call mbnrg_2b_rb_i_poly(w1,w2,e2,g1,g2)
         case(9)
            call mbnrg_2b_cs_i_poly(w1,w2,e2,g1,g2)
       end select
   elseif (tmppoly2 == 5) then
      select case(tmppoly1)

         case(1)
            call mbnrg_2b_li_f_poly(w1,w2,e2,g1,g2)
         case(2)
            call mbnrg_2b_li_cl_poly(w1,w2,e2,g1,g2)
         case(3)
            call mbnrg_2b_li_br_poly(w1,w2,e2,g1,g2)
         case(4)
            call mbnrg_2b_li_i_poly(w1,w2,e2,g1,g2)
         case(5)
            call mbnrg_2b_li_li_poly(w1,w2,e2,g1,g2)
       end select
   elseif (tmppoly2 == 6) then
      select case(tmppoly1)

         case(1)
            call mbnrg_2b_na_f_poly(w1,w2,e2,g1,g2)
         case(2)
            call mbnrg_2b_na_cl_poly(w1,w2,e2,g1,g2)
         case(3)
            call mbnrg_2b_na_br_poly(w1,w2,e2,g1,g2)
         case(4)
            call mbnrg_2b_na_i_poly(w1,w2,e2,g1,g2)
         case(6)
            call mbnrg_2b_na_na_poly(w1,w2,e2,g1,g2)
       end select
   elseif (tmppoly2 == 7) then
      select case(tmppoly1)

         case(1)
            call mbnrg_2b_k_f_poly(w1,w2,e2,g1,g2)
         case(2)
            call mbnrg_2b_k_cl_poly(w1,w2,e2,g1,g2)
         case(3)
            call mbnrg_2b_k_br_poly(w1,w2,e2,g1,g2)
         case(4)
            call mbnrg_2b_k_i_poly(w1,w2,e2,g1,g2)
         case(7)
            call mbnrg_2b_k_k_poly(w1,w2,e2,g1,g2)
       end select
   elseif (tmppoly2 == 8) then
      select case(tmppoly1)

         case(1)
            call mbnrg_2b_rb_f_poly(w1,w2,e2,g1,g2)
         case(2)
            call mbnrg_2b_rb_cl_poly(w1,w2,e2,g1,g2)
         case(3)
            call mbnrg_2b_rb_br_poly(w1,w2,e2,g1,g2)
         case(4)
            call mbnrg_2b_rb_i_poly(w1,w2,e2,g1,g2)
         case(8)
            call mbnrg_2b_rb_rb_poly(w1,w2,e2,g1,g2)
       end select
   else
      select case(tmppoly1)

         case(1)
            call mbnrg_2b_cs_f_poly(w1,w2,e2,g1,g2)
         case(2)
            call mbnrg_2b_cs_cl_poly(w1,w2,e2,g1,g2)
         case(3)
            call mbnrg_2b_cs_br_poly(w1,w2,e2,g1,g2)
         case(4)
            call mbnrg_2b_cs_i_poly(w1,w2,e2,g1,g2)
         case(9)
            call mbnrg_2b_cs_cs_poly(w1,w2,e2,g1,g2)
       end select
   end if
        
      !!!select other cases

!   call mbnrg_2b_h2o_ion_poly(w1,w2,e2,g1,g2)

   engacc = engacc + e2*engunit

   do k = 1, 9
      g1(k) = g1(k)*engunit
      viracc = viracc + w1(k)*g1(k)
   enddo
   
   do k=1,3 
      g2(k) = g2(k)*engunit
      viracc = viracc + w2(k)*g2(k)
   end do

   ! forces

   do k = 1, 3
      fx1(k) = fx1(k) - g1(3*k - 2)
      fy1(k) = fy1(k) - g1(3*k - 1)
      fz1(k) = fz1(k) - g1(3*k - 0)
   enddo

   fx2(1) = fx2(1) - g2(1)
   fy2(1) = fy2(1) - g2(2)
   fz2(1) = fz2(1) - g2(3)

#ifdef STRESS
      
   ! calculate stress tensor

   strs(1) = w1(1) * g1(1) &
           + w1(4) * g1(4) &
           + w1(7) * g1(7) &
           + w2(1) * g2(1) 

   strs(2) = w1(1) * g1(2) &
           + w1(4) * g1(5) &
           + w1(7) * g1(8) &
           + w2(1) * g2(2) 

   strs(3) = w1(1) * g1(3) &
           + w1(4) * g1(6) &
           + w1(7) * g1(9) &
           + w2(1) * g2(3) 

   strs(4) = strs(2)
 
   strs(5) = w1(2) * g1(2) &
           + w1(5) * g1(5) &
           + w1(8) * g1(8) &
           + w2(2) * g2(2) 

   strs(6) = w1(2) * g1(3) &
           + w1(5) * g1(6) &
           + w1(8) * g1(9) &
           + w2(2) * g2(3) 

   strs(7) = strs(3)

   strs(8) = strs(6)

   strs(9) = w1(3) * g1(3) &
           + w1(6) * g1(6) &
           + w1(9) * g1(9) &
           + w2(3) * g2(3) 

   strs(:) = -strs(:)

#endif

end subroutine do_dimer

!----------------------------------------------------------------------------
!DZ (11/30/17): modifying do_trimer to take in multiple ions
subroutine do_trimer(tmppoly, imcon,cell,w1,w2,w3, &
 & fx1,fy1,fz1,fx2,fy2,fz2,fx3,fy3,fz3,engacc,viracc,strs)

   implicit none

   integer, intent(in) :: imcon, tmppoly
   real(8), intent(in) :: cell(*)

   real(8), intent(inout) :: w1(*)
   real(8), intent(inout) :: w2(*)
   real(8), intent(inout) :: w3(*)

   real(8), intent(inout) :: fx1(*),fy1(*),fz1(*)
   real(8), intent(inout) :: fx2(*),fy2(*),fz2(*)
   real(8), intent(inout) :: fx3(*),fy3(*),fz3(*)

   real(8), intent(inout) :: engacc,viracc,strs(9)

   real(8) :: rab(9),r12,r23,r13,eng_old
   real(8) :: e3b,g1(9),g2(9),g3(3)

   logical :: istoobig1, istoobig2, istoobig3

   integer :: k

!FP: added
#ifdef STRESS
   ! initialize stress tensor accumulators
   strs(1:9) = 0.d0
#endif

   ! O-O-ion distances

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

!   istoobig1 = r12.gt.r3f.or.r13.gt.r3f
!   istoobig2 = r12.gt.r3f.or.r23.gt.r3f
!   istoobig3 = r13.gt.r3f.or.r23.gt.r3f

   istoobig3 = r13.lt.r3f.and.r23.lt.r3f
!   istoobig2 = r12.gt.r3f.or.r23.gt.r3f
   if (.not.istoobig3) return

   ! image the trimer [water_monomers are already imaged]
   ! w1 = w2 + rab(1:3)
   ! w3 = w2 - rab(4:6)

!  write(245532,*) '3b',engacc/418.4 ,'before'
!write(245532,*) w1(1:9)
!write(245532,*) w2(1:9)
!write(245532,*) w3(1:3)

   do k = 1, 3
      w1(3+k) = (w1(3+k) - w1(k)) + (w2(k) + rab(k))
      w1(6+k) = (w1(6+k) - w1(k)) + (w2(k) + rab(k))

      w1(k) = w2(k) + rab(k)

!      w3(3+k) = (w3(3+k) - w3(k)) + (w2(k) - rab(k+3))
!      w3(6+k) = (w3(6+k) - w3(k)) + (w2(k) - rab(k+3))
      w3(k) = w2(k) - rab(k+3)
   end do

   g1(1:9) = 0.d0
   g2(1:9) = 0.d0
   g3(1:3) = 0.d0
!   write(245532,*) '3b',engacc/418.4 , 'after'
!write(245532,*) w1(1:9)
!write(245532,*) w2(1:9)
!write(245532,*) w3(1:3)
!

!DZ (11/30/17): adding different cases for different ions

   select case(tmppoly)
      case(1)
         call mbnrg_3b_h2o_h2o_f_poly(w1,w2,w3,e3b,g1,g2,g3)
      case(2)
         call mbnrg_3b_h2o_h2o_cl_poly(w1,w2,w3,e3b,g1,g2,g3)
      case(3)
         call mbnrg_3b_h2o_h2o_br_poly(w1,w2,w3,e3b,g1,g2,g3)
      case(4)
         call mbnrg_3b_h2o_h2o_i_poly(w1,w2,w3,e3b,g1,g2,g3)
      case(5)
         call mbnrg_3b_h2o_h2o_li_poly(w1,w2,w3,e3b,g1,g2,g3)
      case(6)
         call mbnrg_3b_h2o_h2o_na_poly(w1,w2,w3,e3b,g1,g2,g3)
      case(7)
         call mbnrg_3b_h2o_h2o_k_poly(w1,w2,w3,e3b,g1,g2,g3)
      case(8)
         call mbnrg_3b_h2o_h2o_rb_poly(w1,w2,w3,e3b,g1,g2,g3)
      case(9)
         call mbnrg_3b_h2o_h2o_cs_poly(w1,w2,w3,e3b,g1,g2,g3)
    end select

!   call mbnrg_3b_h2o_ion_poly(w1,w2,w3,e3b,g1,g2,g3)
!write(24000,*) 'e3b',e3b 
   engacc = engacc + e3b*engunit

!   if (e3b.lt.-3.d0) &
!       & call drop_trimer(e3b,w1,w2,w3)

   do k = 1, 9
      g1(k) = g1(k)*engunit
      viracc = viracc + w1(k)*g1(k)
      g2(k) = g2(k)*engunit
      viracc = viracc + w2(k)*g2(k)
   enddo

   do k=1,3 
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
  enddo

  fx3(1) = fx3(1) - g3(1)
  fy3(1) = fy3(1) - g3(2)
  fz3(1) = fz3(1) - g3(3)

#ifdef STRESS
      
   ! calculate stress tensor

   strs(1) = w1(1) * g1(1) &
           + w1(4) * g1(4) &
           + w1(7) * g1(7) &
           + w2(1) * g2(1) &
           + w2(4) * g2(4) &
           + w2(7) * g2(7) &
           + w3(1) * g3(1) 

   strs(2) = w1(1) * g1(2) &
           + w1(4) * g1(5) &
           + w1(7) * g1(8) &
           + w2(1) * g2(2) &
           + w2(4) * g2(5) &
           + w2(7) * g2(8) &
           + w3(1) * g3(2) 

   strs(3) = w1(1) * g1(3) &
           + w1(4) * g1(6) &
           + w1(7) * g1(9) &
           + w2(1) * g2(3) &
           + w2(4) * g2(6) &
           + w2(7) * g2(9) &
           + w3(1) * g3(3) 

   strs(4) = strs(2)
 
   strs(5) = w1(2) * g1(2) &
           + w1(5) * g1(5) &
           + w1(8) * g1(8) &
           + w2(2) * g2(2) &
           + w2(5) * g2(5) &
           + w2(8) * g2(8) &
           + w3(2) * g3(2) 

   strs(6) = w1(2) * g1(3) &
           + w1(5) * g1(6) &
           + w1(8) * g1(9) &
           + w2(2) * g2(3) &
           + w2(5) * g2(6) &
           + w2(8) * g2(9) &
           + w3(2) * g3(3) 

   strs(7) = strs(3)

   strs(8) = strs(6)

   strs(9) = w1(3) * g1(3) &
           + w1(6) * g1(6) &
           + w1(9) * g1(9) &
           + w2(3) * g2(3) &
           + w2(6) * g2(6) &
           + w2(9) * g2(9) &
           + w3(3) * g3(3) 

   strs(:) = -strs(:)

#endif

end subroutine do_trimer

!----------------------------------------------------------------------------

subroutine image_water_monomers(imcon, cell, nw, xyz)

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

end subroutine image_water_monomers

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

end module mbnrg

!============================================================================
