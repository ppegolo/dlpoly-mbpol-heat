!
! === output data (TTM2 model) ===
!
! Last updated: 23 Feb 2006 by S. Iuchi
!

subroutine output_ttm2_dip( listttm2, keybin, n_water, nttm2, dipx, dipy, dipz, lcmd)
  
  use multibead, only: bead_suffix, comm_ring, ring_rank, ring_size
  use dipole_moments,   only: moldipx, moldipy, moldipz, chkdipx, chkdipy, chkdipz
  use global_variables, only: ndip1, ndip2, ndip3, mxatms
!  use ps_type_dms,      only: ldms
  use unit_parameters,  only: eatd

  implicit none

! arguments
  
  integer, intent(in) :: n_water, listttm2(mxatms), nttm2, keybin 

  real(8), intent(in) :: dipx(mxatms), dipy(mxatms), dipz(mxatms)

  logical, intent(in) :: lcmd

#include "mpif.h"

! local variables

  integer :: i, mttm2, iox, ih1, ih2, imm

! for form in open statement 

  character(len=50) :: fileformat
  
!  real(8) :: chkdip(n_water)
  real(8) :: indx(n_water), indy(n_water), indz(n_water), ind(n_water)
  real(8) :: moldip(n_water), scratch(n_water)

!
! --- initialize induced dipole moments ---
!

  indx(:) = 0.0d0;  indy(:) = 0.0d0;  indz(:) = 0.0d0

  ind(:) = 0.0d0

!
! --- output water dipole moments ---
!

  mttm2 = listttm2(nttm2)

  do i=1,n_water
     
     mttm2 = mttm2 + 1
     
     iox = listttm2(3*i-2)
     ih1 = listttm2(3*i-1)
     ih2 = listttm2(3*i)
!    imm = n_water*3 + i
     imm = listttm2(3*n_water) + i

     indx(i) = indx(i) + dipx(iox) + dipx(ih1) + dipx(ih2) + dipx(imm)
     indy(i) = indy(i) + dipy(iox) + dipy(ih1) + dipy(ih2) + dipy(imm)
     indz(i) = indz(i) + dipz(iox) + dipz(ih1) + dipz(ih2) + dipz(imm)

  end do

!
! --- open files ---
!

  if(keybin==0) fileformat='formatted'
  if(keybin==1) fileformat='unformatted'

  if (lcmd) then

    ! reduce the dipoles at bead #1 (ring_rank == 0)

    call MPI_REDUCE(moldipx, scratch, n_water, &
        & MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm_ring, i)
    do i=1,n_water
      moldipx(i) = scratch(i)/ring_size
    end do

    call MPI_REDUCE(moldipy, scratch, n_water, &
        & MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm_ring, i)
    do i=1,n_water
      moldipy(i) = scratch(i)/ring_size
    end do

    call MPI_REDUCE(moldipz, scratch, n_water, &
        & MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm_ring, i)
    do i=1,n_water
      moldipz(i) = scratch(i)/ring_size
    end do

    call MPI_REDUCE(indx, scratch, n_water, &
        & MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm_ring, i)
    do i=1,n_water
      indx(i) = scratch(i)/ring_size
    end do

    call MPI_REDUCE(indy, scratch, n_water, &
        & MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm_ring, i)
    do i=1,n_water
      indy(i) = scratch(i)/ring_size
    end do

    call MPI_REDUCE(indz, scratch, n_water, &
        & MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm_ring, i)
    do i=1,n_water
      indz(i) = scratch(i)/ring_size
    end do

    if (ring_rank.eq.0) then ! this is centroid bead and lhead is true

      open(ndip1,file='DIPIND_CMD',position='append',form=fileformat)
      open(ndip2,file='DIPMOL_CMD',position='append',form=fileformat)

      do i=1,n_water

        ind(i)    = sqrt( indx(i)**2    + indy(i)**2    + indz(i)**2    )
        moldip(i) = sqrt( moldipx(i)**2 + moldipy(i)**2 + moldipz(i)**2 )
      end do

        if(keybin==0) then
         do i=1,n_water
           write(ndip1,'(i5,4f20.12)') i, indx(i) * eatd, indy(i) * eatd,&
             & indz(i) * eatd, ind(i) * eatd
           write(ndip2,'(i5,4f20.12)') i, moldipx(i) * eatd, moldipy(i) * eatd,&
             & moldipz(i) * eatd, moldip(i) * eatd
         end do
        elseif(keybin==1) then 
           write(ndip1) (i,i=1,n_water), indx(1:n_water) * eatd, indy(1:n_water) * eatd,&
             & indz(1:n_water) * eatd, ind(1:n_water) * eatd
           write(ndip2) (i,i=1,n_water), moldipx(1:n_water) * eatd, moldipy(1:n_water) * eatd,&
             & moldipz(1:n_water) * eatd, moldip(1:n_water) * eatd
        endif

      close(ndip1)
      close(ndip2)

    end if ! ring_rank.eq.0

  else ! not CMD

    open(ndip1,file='DIPIND'//bead_suffix,position='append',form=fileformat)
    open(ndip2,file='DIPMOL'//bead_suffix,position='append',form=fileformat)

!  if( ldms ) open(ndip3,file='DIPCHK',position='append')
!
!
! --- output molecular and induced dipole moments ---
!

    do i=1,n_water
       ind(i)    = sqrt( indx(i)**2    + indy(i)**2    + indz(i)**2    )
       moldip(i) = sqrt( moldipx(i)**2 + moldipy(i)**2 + moldipz(i)**2 )
    enddo

!     if( ldms ) chkdip(i) = sqrt( chkdipx(i)**2 + chkdipy(i)**2 + chkdipz(i)**2 )

       if(keybin==0) then
         do i=1,n_water
          write(ndip1,'(i5,4f20.12)') i, indx(i) * eatd, indy(i) * eatd,&
             & indz(i) * eatd, ind(i) * eatd
          write(ndip2,'(i5,4f20.12)') i, moldipx(i) * eatd, moldipy(i) * eatd,&
             & moldipz(i) * eatd, moldip(i) * eatd
         end do
       elseif(keybin==1) then
!         do i=1,n_water
!         write(ndip1) i, indx(i) * eatd, indy(i) * eatd,&
!            & indz(i) * eatd, ind(i) * eatd
!         write(ndip2) i, moldipx(i) * eatd, moldipy(i) * eatd,&
!            & moldipz(i) * eatd, moldip(i) * eatd
!        end do
         write(ndip1) (i,i=1,n_water), indx(1:n_water) * eatd, indy(1:n_water) * eatd,&
            & indz(1:n_water) * eatd, ind(1:n_water) * eatd
         write(ndip2) (i,i=1,n_water), moldipx(1:n_water) * eatd, moldipy(1:n_water) * eatd,&
            & moldipz(1:n_water) * eatd, moldip(1:n_water) * eatd
       endif
!     if( ldms ) write(ndip3,'(i5,4f20.12)') i, chkdipx(i) * eatd, chkdipy(i) * eatd,&
!          & chkdipz(i) * eatd, chkdip(i) * eatd
     

!
! --- close files ---
!

    close(ndip1);  close(ndip2)

!  if( ldms ) close(ndip3)

  end if ! lcmd

  return

end subroutine output_ttm2_dip
