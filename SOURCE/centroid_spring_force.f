c==============================================================================
c
      subroutine centroid_spring_force
     x   (xxx,yyy,zzz,npx,npy,npz,fxx,fyy,fzz,Epot_spring,Epot_deriv)
c
c==============================================================================

       use multibead
       use centroid_module, only: nmbuff,
     x                      nbead, natom, hbar, beta,
     x                      phys_mass
       use Langevin_Thermostat , only: langevin,add_spring_to_forces,
     x     force_scale,freq_scale 

c     !................................................................

      implicit none

      real(8), intent(in) :: xxx(*), yyy(*), zzz(*)
      real(8), intent(in) :: npx(*), npy(*), npz(*)

      real(8) :: fxx(*), fyy(*), fzz(*)
      real(8) :: Epot_spring, Epot_deriv
      real(8) :: fac1,fac3

c     !................................................................

#include "mpif.h"

      integer :: status(MPI_STATUS_SIZE)
      integer :: iatm1, iatm2, iatm
      integer :: j, ierr, next_rank, prev_rank
      integer :: inext, iprev, iatom

      integer, parameter ::  cw_tag = 27606
      integer, parameter :: ccw_tag = 27607

      real(8) :: coeff, dx(3), dy(3), dz(3)
      real(8) :: sum_e, sum_d, sum_v, tmp

c     !................................................................

      coeff = dble( nbead ) / ( hbar*beta )**2

c     !----------------------------------------------------------------

      iatm1 = (bead_rank*natom)/bead_size+1
      iatm2 = ((bead_rank+1)*natom)/bead_size

      iatm = iatm2 - iatm1 + 1

!     1) broadcast centroid to other beads

      fac1 = 1.d0
      fac3 = force_scale 
      if(langevin) then
        ! normalization of NM transformation is 1/sqrt(nbead)
        fac1 = 1.0d0 / sqrt(dble(nbead)) 
        if(.not.add_spring_to_forces) then
          fac3 = 0.d0
        endif
      endif
      j = 0
      do iatom = iatm1, iatm2
         nmbuff(j + 1) = npx(iatom)*fac1
         nmbuff(j + 2) = npy(iatom)*fac1
         nmbuff(j + 3) = npz(iatom)*fac1
         j = j + 3
      end do

      call MPI_BCAST(nmbuff, 3*iatm, MPI_DOUBLE_PRECISION,
     x               0, comm_ring, ierr)

!     2) get coordinates of the prev/next beads

      next_rank = ring_rank + 1
      prev_rank = ring_rank - 1

      if (next_rank.eq.ring_size) next_rank = 0
      if (prev_rank.lt.0) prev_rank = ring_size - 1

      iprev = 3*iatm
      inext = 2*3*iatm

      j = 3*3*iatm
      do iatom = iatm1, iatm2
         nmbuff(j + 1) = xxx(iatom)
         nmbuff(j + 2) = yyy(iatom)
         nmbuff(j + 3) = zzz(iatom)
         j = j + 3
      end do

      call MPI_SENDRECV(nmbuff(3*3*iatm + 1), 3*iatm,
     x                  MPI_DOUBLE_PRECISION, next_rank, cw_tag,
     x                  nmbuff(iprev + 1), 3*iatm,
     x                  MPI_DOUBLE_PRECISION, prev_rank, cw_tag,
     x                  comm_ring, status, ierr)

      call MPI_SENDRECV(nmbuff(3*3*iatm + 1), 3*iatm,
     x                  MPI_DOUBLE_PRECISION, prev_rank, ccw_tag,
     x                  nmbuff(inext + 1), 3*iatm,
     x                  MPI_DOUBLE_PRECISION, next_rank, ccw_tag,
     x                  comm_ring, status, ierr)

!     3) compute the thing

      sum_e = 0.d0
      sum_d = 0.d0

      j = 0

      do iatom = iatm1, iatm2

         dx(1) = xxx(iatom) - nmbuff(j + 1)
         dy(1) = yyy(iatom) - nmbuff(j + 2)
         dz(1) = zzz(iatom) - nmbuff(j + 3)

         dx(2) = xxx(iatom) - nmbuff(iprev + j + 1)
         dy(2) = yyy(iatom) - nmbuff(iprev + j + 2)
         dz(2) = zzz(iatom) - nmbuff(iprev + j + 3)

         dx(3) = xxx(iatom) - nmbuff(inext + j + 1)
         dy(3) = yyy(iatom) - nmbuff(inext + j + 2)
         dz(3) = zzz(iatom) - nmbuff(inext + j + 3)

         sum_d = sum_d - dx(1) * fxx(iatom)
     x                 - dy(1) * fyy(iatom)
     x                 - dz(1) * fzz(iatom)

         sum_e = sum_e
     x         + phys_mass(iatom) * (dx(2)**2 + dy(2)**2 + dz(2)**2)

         fxx(iatom) = fxx(iatom)
     x              - coeff * phys_mass(iatom) * (dx(2) + dx(3)) * fac3
         fyy(iatom) = fyy(iatom)
     x              - coeff * phys_mass(iatom) * (dy(2) + dy(3)) * fac3
         fzz(iatom) = fzz(iatom)
     x              - coeff * phys_mass(iatom) * (dz(2) + dz(3)) * fac3 

         j = j + 3

      enddo

      Epot_spring = 0.5d0 * coeff * sum_e 

      Epot_deriv = 0.5d0 * sum_d / force_scale

      nmbuff(1) = Epot_spring
      nmbuff(2) = Epot_deriv

      call MPI_ALLREDUCE(nmbuff(1),nmbuff(3),2,MPI_DOUBLE_PRECISION,
     x                   MPI_SUM, comm_mb, ierr)

      Epot_spring = nmbuff(3)
      Epot_deriv = nmbuff(4)

      end subroutine
