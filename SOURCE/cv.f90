!**********************************************************************
!
!**********************************************************************
!----------------------------------------------------------------------------
MODULE CV
  !----------------------------------------------------------------------------
  SAVE
  !
  LOGICAL :: debug 
  LOGICAL :: lhead,pimd_head 
  INTEGER :: nions,nColVar
  !
  LOGICAL :: CV_Ql = .false. 
  !
  ! pre-allocating makes things simpler 
  REAL*8 :: ev_pos(10), ev_force(10),ev_kappa(10),cv_pos(10),Ql_En(10) 
  REAL*8,ALLOCATABLE :: CV_FORCE(:,:),buffer_in(:),buffer_out(:) 
  REAL*8 :: CV_STRESS(9)
  REAL*8 :: rmin_ql,rmax_ql,nlist_cutoff 
  !
  PUBLIC :: CV_init,CV_eval,CV_STRESS,CV_Ql,ev_pos,ev_force,ev_kappa
  PUBLIC :: nColVar,CV_FORCE,rmin_ql,rmax_ql,nlist_cutoff,cv_pos,Ql_En 
  !
  CONTAINS
    !
    ! ... public methods
    !
    !------------------------------------------------------------------------
    SUBROUTINE CV_init(lhead_,pimd_head_,nat_)
    !------------------------------------------------------------------------
      !
      !
      IMPLICIT NONE
      !
      logical, intent(in) :: lhead_,pimd_head_
      integer, intent(in) :: nat_
      !
      INTEGER i,j,IERR
      !

      if(.not.CV_Ql) return

      lhead = lhead_
      pimd_head = pimd_head_
      nions = nat_

      allocate(CV_FORCE(3,nions),buffer_in(3*nions+60),buffer_out(3*nions+60))
     
    !------------------------------------------------------------------------
    END SUBROUTINE CV_init
    !------------------------------------------------------------------------
    SUBROUTINE CV_eval(cell,xxx,yyy,zzz,cv_energy)
    !------------------------------------------------------------------------
      !
      !
      use multibead,   only: comm_bead,bead_rank,bead_size
      !
      IMPLICIT NONE
      !
#ifdef MPI
#  include "mpif.h"
#endif
      !
      real*8, intent(in) :: xxx(:),yyy(:),zzz(:) 
      real*8, intent(in) :: cell(9) 
      real*8, intent(out) :: cv_energy 
      ! 
      INTEGER :: i,j,ierr,nsend 

      IF(.NOT.CV_Ql) return

      Ql_En(:) = 0.d0
      cv_pos(:)=0.d0
      ev_force(:)=0.d0
      CV_FORCE(:,:) = 0.d0
      CV_STRESS(:)=0.d0     

      ! parallelize code
!      CALL Ql_wrapper(comm_bead,bead_rank,bead_size,nions, nlist_cutoff, cell, & 
!                       xxx, yyy, zzz, ev_pos, ev_force,   &
!                       ev_kappa,cv_pos,rmin_ql,rmax_ql,CV_FORCE,CV_STRESS,Ql_En) 
#ifdef MPI
      nsend=1
      do i=1,nions 
        do j=1,3
          buffer_out(nsend) = CV_FORCE(j,i)
          nsend = nsend + 1 
        enddo
      enddo
      do i=1,9
          buffer_out(nsend) = CV_STRESS(i)
          nsend = nsend + 1
      enddo
      do i=1,nColVar
        buffer_out(nsend) = cv_pos(i)
        buffer_out(nsend+1) = ev_force(i)
        buffer_out(nsend+2) = Ql_En(i)
        nsend = nsend + 3 
      enddo
      call MPI_ALLREDUCE(buffer_out,buffer_in,nsend,MPI_DOUBLE_PRECISION,   &
                      MPI_SUM, comm_bead ,ierr)
      nsend=1
      do i=1,nions  
        do j=1,3
          CV_FORCE(j,i) = buffer_in(nsend)
          nsend = nsend + 1
        enddo
      enddo
      do i=1,9
          CV_STRESS(i) = buffer_in(nsend)
          nsend = nsend + 1
      enddo
      do i=1,nColVar
        cv_pos(i) = buffer_in(nsend)
        ev_force(i) = buffer_in(nsend+1)
        Ql_En(i) = buffer_in(nsend+2)
        nsend = nsend + 3
      enddo
#endif

      ! add energy contributions
      cv_energy = 0.d0
      do i=1,nColVar 
        cv_energy = cv_energy + Ql_En(i)
      enddo

      RETURN
    !------------------------------------------------------------------------
    END SUBROUTINE CV_EVAL 
    !------------------------------------------------------------------------
  END MODULE CV
    !------------------------------------------------------------------------
