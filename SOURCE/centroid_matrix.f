c==============================================================================
c
      subroutine transformation_matrix
c
c==============================================================================

      use centroid_module, only: nbead, natom, 
     x                           kT, beta, hbar, pi,  
     x                           phys_mass,
     x                           cart_to_nmode, nmode_to_cart,  
     x                           tcart_to_tnmode,  
     x                           lambda_nmode, omega_nmode,  
     x                           mass_nmode, fict_mass_nmode
      use Langevin_Thermostat,   only:
     x                           langevin,freq_scale, force_scale  

c     !..................................................

      implicit none

      integer :: n, n1, n2, iatom, ibead

      real(8) :: omega_P, tmp, tmp1

c     !..................................................

      omega_P = sqrt(dble(nbead)) / hbar / beta * freq_scale

c     !---------------------------------------

c     Transformation matrix from normal modes to cartesian coordinates.
c     Check: it should come from one of the first papers by Cao and Voth.

      nmode_to_cart(:, 1) = 1.d0

      do n1 = 1, nbead/2
         nmode_to_cart(2*n1,nbead) = -1.0d0
         nmode_to_cart(2*n1-1,nbead) =  1.0d0
      end do


      do n1 = 1, nbead

         tmp = 2.0d0 * pi * dble(n1-1) / dble(nbead)

         do n2 = 1, nbead/2 - 1
            nmode_to_cart(n1,2*n2) =
     x         sqrt(2.d0) * cos(tmp * dble(n2))
            nmode_to_cart(n1,2*n2+1) =
     x         - sqrt(2.d0) * sin(tmp * dble(n2))
         enddo

      enddo

c     !---------------------------------------

c     ! Transformation matrix from cartesian coordinates to normal modes.

      do n2 = 1, nbead
         do n1 = 1, nbead
            cart_to_nmode(n2,n1) = 
     x          nmode_to_cart(n1,n2) / dble(nbead)
            tcart_to_tnmode(n2,n1) = nmode_to_cart(n1,n2)
         enddo
      enddo

      ! change ordering of modes for consistency (specially in PIGLET) 
      !if(langevin.and..false.) then
      if(langevin) then

        ! nmode_to_cart(n1,n2) = C(n1,n2), --> n1=j, n2=k 
        do n1 = 1, nbead
           tmp1 = dsqrt(1.d0/nbead)
           !tmp1 = 1.d0
           nmode_to_cart(n1,1) = tmp1 
           nmode_to_cart(n1,nbead/2+1) = tmp1*(-1.0d0)**(n1)
           tmp1 = dsqrt(2.d0/nbead)
           !tmp1 = dsqrt(2.d0) 
           tmp = 2.0d0 * pi * dble( n1 ) / dble( nbead )
           do n2 = 2, nbead/2
            nmode_to_cart(n1,n2) =  
     x         tmp1 * cos( tmp * dble( n2-1 ) )
            nmode_to_cart(n1,n2+nbead/2) =  
     x         tmp1 * sin( tmp * dble( n2+nbead/2.0d0-1.0d0 ) )
           enddo
        enddo  

c       ! Transformation matrix from cartesian coordinates to normal modes.
        do n2 = 1, nbead
           do n1 = 1, nbead
              cart_to_nmode(n2,n1) = nmode_to_cart(n1,n2) 
              tcart_to_tnmode(n2,n1) = nmode_to_cart(n1,n2) 
           enddo
        enddo

      endif

c     !---------------------------------------

c     Normal-mode eigenvalues. 

      if (langevin) then
         ! these are not actual frequencies, but they produce the
         ! appropriate masses. These should not be used anywhere anyway 
         do n = 1, nbead
            lambda_nmode(n) = 1.0d0
         enddo

      else  

         lambda_nmode(1) = 0.0d0
         do n = 1, nbead/2 - 1
            tmp = 2.0d0 * pi * dble(n) / dble(nbead)
            tmp = 2.0d0 * (1.0d0 - cos(tmp))
            lambda_nmode(2*n) = tmp * dble(nbead)
            lambda_nmode(2*n+1) = tmp * dble(nbead)
         enddo
         lambda_nmode(nbead) = 4.d0 * dble(nbead)

      end if

c     !---------------------------------------

c     Normal-mode masses.
c     ibead = 1 --> path-centroid.

      do iatom = 1, natom
         mass_nmode(iatom,1) = phys_mass(iatom)*lambda_nmode(1)
         do ibead = 2, nbead
            mass_nmode(iatom,ibead) =  
     x           phys_mass(iatom) * lambda_nmode(ibead)
         enddo
      enddo

      end subroutine


c==============================================================================
c
      subroutine transform_position_from_cart_to_nmode
     x (xxx,yyy,zzz,npx,npy,npz)
c
c     Transform cartesian positions into normal mode positions.
c==============================================================================

      use multibead
      use centroid_module, only: natom, nbead, cart_to_nmode,
     x                           nmbuff, mxnmbuff

c     !..................................................

      implicit none

      real(8), intent(in) :: xxx(natom),yyy(natom),zzz(natom)
      real(8)             :: npx(natom),npy(natom),npz(natom)

#include "mpif.h"

c     !..................................................

      integer :: ierr, iatm1, iatm2, iatm
      integer :: i, j, b1, ibead

c     !..................................................

      iatm1 = (bead_rank*natom)/bead_size+1
      iatm2 = ((bead_rank+1)*natom)/bead_size

      iatm = iatm2-iatm1+1

      ! gather coordinates
      j=0
      do i=iatm1,iatm2
         nmbuff(j+1)=xxx(i)
         nmbuff(j+2)=yyy(i)
         nmbuff(j+3)=zzz(i)
         j=j+3
      end do

      call MPI_ALLGATHER(nmbuff(1),3*iatm,MPI_DOUBLE_PRECISION,
     x            nmbuff(3*iatm+1),3*iatm,MPI_DOUBLE_PRECISION,
     x                   comm_ring,ierr)

      npx(iatm1:iatm2)=0.d0
      npy(iatm1:iatm2)=0.d0
      npz(iatm1:iatm2)=0.d0

      ! transform
      do ibead=1,nbead
          b1=3*iatm*ibead
          j=0
          do i=iatm1,iatm2
             npx(i)=npx(i)
     x             +cart_to_nmode(ring_rank+1,ibead)*nmbuff(b1+j+1)
             npy(i)=npy(i)
     x             +cart_to_nmode(ring_rank+1,ibead)*nmbuff(b1+j+2)
             npz(i)=npz(i)
     x             +cart_to_nmode(ring_rank+1,ibead)*nmbuff(b1+j+3)
             j=j+3
          enddo
      enddo

      end subroutine


c==============================================================================
c
      subroutine transform_position_from_nmode_to_cart
     x (npx,npy,npz,xxx,yyy,zzz)
c
c     Transform normal mode positions into cartesian positions.
c==============================================================================

      use multibead
      use centroid_module, only: natom, nbead, nmode_to_cart,
     x                           nmbuff, mxnmbuff

c     !..................................................

      implicit none

      real(8), intent(in) :: npx(natom),npy(natom),npz(natom)
      real(8)             :: xxx(natom),yyy(natom),zzz(natom)

#include "mpif.h"

c     !..................................................

      integer :: ierr, iatm1, iatm2, iatm
      integer :: i, j, b, ibead

c     !..................................................

      iatm1 = (bead_rank*natom)/bead_size+1
      iatm2 = ((bead_rank+1)*natom)/bead_size

      iatm = iatm2-iatm1+1

      j=0
      do i=iatm1,iatm2
         nmbuff(j+1)=npx(i)
         nmbuff(j+2)=npy(i)
         nmbuff(j+3)=npz(i)
         j=j+3
      end do

      call MPI_ALLGATHER(nmbuff(1),3*iatm,MPI_DOUBLE_PRECISION,
     x            nmbuff(3*iatm+1),3*iatm,MPI_DOUBLE_PRECISION,
     x                   comm_ring,ierr)

      xxx(iatm1:iatm2)=0.d0
      yyy(iatm1:iatm2)=0.d0
      zzz(iatm1:iatm2)=0.d0

      ! transform
      do ibead=1,nbead
          b=3*iatm*ibead
          j=0
          do i=iatm1,iatm2
             xxx(i)=xxx(i)
     x             +nmode_to_cart(ring_rank+1,ibead)*nmbuff(b+j+1)
             yyy(i)=yyy(i)
     x             +nmode_to_cart(ring_rank+1,ibead)*nmbuff(b+j+2)
             zzz(i)=zzz(i)
     x             +nmode_to_cart(ring_rank+1,ibead)*nmbuff(b+j+3)
             j=j+3
          enddo
      enddo

      end subroutine
      

c==============================================================================
c
      subroutine transform_force_from_cart_to_nmode
     x (fxx,fyy,fzz,nfx,nfy,nfz)
c
c     Transform cartesian forces into normal mode forces
c==============================================================================

      use multibead
      use centroid_module, only: natom, nbead, tcart_to_tnmode,
     x                           nmbuff, mxnmbuff

c     !..................................................

      implicit none

      real(8), intent(in) :: fxx(natom),fyy(natom),fzz(natom)
      real(8)             :: nfx(natom),nfy(natom),nfz(natom)

#include "mpif.h"

c     !..................................................

      integer :: ierr, iatm1, iatm2, iatm
      integer :: i, j, b, ibead

c     !..................................................

      iatm1 = (bead_rank*natom)/bead_size+1
      iatm2 = ((bead_rank+1)*natom)/bead_size

      iatm = iatm2-iatm1+1

      ! gather the forces
      j=0
      do i=iatm1,iatm2
         nmbuff(j+1)=fxx(i)
         nmbuff(j+2)=fyy(i)
         nmbuff(j+3)=fzz(i)
         j=j+3
      end do

      call MPI_ALLGATHER(nmbuff(1),3*iatm,MPI_DOUBLE_PRECISION,
     x            nmbuff(3*iatm+1),3*iatm,MPI_DOUBLE_PRECISION,
     x                   comm_ring,ierr)

      nfx(iatm1:iatm2)=0.d0
      nfy(iatm1:iatm2)=0.d0
      nfz(iatm1:iatm2)=0.d0

      ! transform
      do ibead=1,nbead
          b=3*iatm*ibead
          j=0
          do i=iatm1,iatm2
             nfx(i)=nfx(i)
     x             +tcart_to_tnmode(ring_rank+1,ibead)*nmbuff(b+j+1)
             nfy(i)=nfy(i)
     x             +tcart_to_tnmode(ring_rank+1,ibead)*nmbuff(b+j+2)
             nfz(i)=nfz(i)
     x             +tcart_to_tnmode(ring_rank+1,ibead)*nmbuff(b+j+3)
             j=j+3
          enddo
      enddo

      end subroutine


c==============================================================================
c
      subroutine transform_velocity_from_nmode_to_cart
     x (nvx,nvy,nvz,vxx,vyy,vzz)
c
c     Transform normal mode positions into cartesian positions.
c==============================================================================

      use multibead
      use centroid_module, only: natom, nbead, nmode_to_cart,
     x                           nmbuff, mxnmbuff

c     !..................................................

      implicit none

      real(8), intent(in) :: nvx(natom),nvy(natom),nvz(natom)
      real(8)             :: vxx(natom),vyy(natom),vzz(natom)

#include "mpif.h"

c     !..................................................

      integer :: ierr, iatm1, iatm2, iatm
      integer :: i, j, b, ibead

c     !..................................................

      iatm1 = (bead_rank*natom)/bead_size+1
      iatm2 = ((bead_rank+1)*natom)/bead_size

      iatm = iatm2-iatm1+1

      ! gather the forces
      j=0
      do i=iatm1,iatm2
         nmbuff(j+1)=nvx(i)
         nmbuff(j+2)=nvy(i)
         nmbuff(j+3)=nvz(i)
         j=j+3
      end do

      call MPI_ALLGATHER(nmbuff(1),3*iatm,MPI_DOUBLE_PRECISION,
     x            nmbuff(3*iatm+1),3*iatm,MPI_DOUBLE_PRECISION,
     x                   comm_ring,ierr)

      vxx(iatm1:iatm2)=0.d0
      vyy(iatm1:iatm2)=0.d0
      vzz(iatm1:iatm2)=0.d0

      ! transform
      do ibead=1,nbead
          b=3*iatm*ibead
          j=0
          do i=iatm1,iatm2
             vxx(i)=vxx(i)
     x             +nmode_to_cart(ring_rank+1,ibead)*nmbuff(b+j+1)
             vyy(i)=vyy(i)
     x             +nmode_to_cart(ring_rank+1,ibead)*nmbuff(b+j+2)
             vzz(i)=vzz(i)
     x             +nmode_to_cart(ring_rank+1,ibead)*nmbuff(b+j+3)
             j=j+3
          enddo
      enddo

      end subroutine


c==============================================================================
c
      subroutine transform_velocity_from_cart_to_nmode
     x (vxx,vyy,vzz,nvx,nvy,nvz)
c
c     Transform normal mode positions into cartesian positions.
c==============================================================================

      use multibead
      use centroid_module, only: natom, nbead, cart_to_nmode,
     x                           nmbuff, mxnmbuff

c     !..................................................

      implicit none

      real(8), intent(in) :: vxx(natom),vyy(natom),vzz(natom)
      real(8)             :: nvx(natom),nvy(natom),nvz(natom)

#include "mpif.h"

c     !..................................................

      integer :: ierr, iatm1, iatm2, iatm
      integer :: i, j, b, ibead

c     !..................................................

      iatm1 = (bead_rank*natom)/bead_size+1
      iatm2 = ((bead_rank+1)*natom)/bead_size

      iatm = iatm2-iatm1+1

      ! gather the velocities
      j=0
      do i=iatm1,iatm2
         nmbuff(j+1)=vxx(i)
         nmbuff(j+2)=vyy(i)
         nmbuff(j+3)=vzz(i)
         j=j+3
      end do

      call MPI_ALLGATHER(nmbuff(1),3*iatm,MPI_DOUBLE_PRECISION,
     x            nmbuff(3*iatm+1),3*iatm,MPI_DOUBLE_PRECISION,
     x                   comm_ring,ierr)

      nvx(iatm1:iatm2)=0.d0
      nvy(iatm1:iatm2)=0.d0
      nvz(iatm1:iatm2)=0.d0

      ! transform
      do ibead=1,nbead
          b=3*iatm*ibead
          j=0
          do i=iatm1,iatm2
             nvx(i)=nvx(i)
     x             +cart_to_nmode(ring_rank+1,ibead)*nmbuff(b+j+1)
             nvy(i)=nvy(i)
     x             +cart_to_nmode(ring_rank+1,ibead)*nmbuff(b+j+2)
             nvz(i)=nvz(i)
     x             +cart_to_nmode(ring_rank+1,ibead)*nmbuff(b+j+3)
             j=j+3
          enddo
      enddo

      end subroutine

c==============================================================================
