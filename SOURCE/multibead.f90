#include "assert.h"

module multibead

implicit none

integer, save :: comm_mb, & ! mb = multibead
                 comm_bead, & ! bead's internal communicator
                 comm_ring ! inter-bead atom ring iatm1:iatm2

integer, save :: mb_rank,mb_size, & ! comm_mb
                 bead_rank,bead_size, & ! comm_bead
                 ring_rank,ring_size ! comm_ring

character(8), save :: bead_suffix = ''

public :: mb_init
public :: mb_finalize

public :: mb_abort
public :: is_bead_head

contains

subroutine mb_init(nbeads)

   implicit none

   integer, optional, intent(in) :: nbeads

#ifdef MPI
#  include "mpif.h"

   integer :: ierr,ibead,nbead
   character(8) :: str

   call MPI_INIT(ierr)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)

   call MPI_COMM_DUP(MPI_COMM_WORLD,comm_mb,ierr)

   call MPI_COMM_RANK(comm_mb,mb_rank,ierr)
   call MPI_COMM_SIZE(comm_mb,mb_size,ierr)

   if (present(nbeads)) then
      assert(nbeads.gt.0)
      nbead = nbeads
   else
      str = ''
      call getenv('PIMD_CMD_NUM_BEADS',str)
      read(str,'(i3)') nbead
      nbead = max(1,nbead)
      call MPI_BCAST(nbead,1,MPI_INTEGER,0,comm_mb,ierr)
   end if ! present(nbeads)

   if (mb_size.lt.nbead) then
      if (mb_rank.eq.0) &
        print '(/a,i2,a/)', ' ** Error ** : too few threads for ',nbead,' beads'
      call MPI_ABORT(comm_mb,0,ierr)
   end if

   if (mod(mb_size,nbead).ne.0) then
      if (mb_rank.eq.0) &
        print '(/a/)', ' ** Error ** : incommensurate number of beads/threads'
      call MPI_ABORT(comm_mb,0,ierr)
   end if

   ibead = mb_rank/(mb_size/nbead) + 1
   if (.not.present(nbeads)) then
      write (bead_suffix,"('.',i2.2)") ibead
   else
      bead_suffix = ''
   end if

   call MPI_COMM_SPLIT(comm_mb,ibead,mb_rank,comm_bead,ierr)

   call MPI_COMM_RANK(comm_bead,bead_rank,ierr)
   call MPI_COMM_SIZE(comm_bead,bead_size,ierr)

   call MPI_COMM_SPLIT(comm_mb,mod(mb_rank,mb_size/nbead), &
                       mb_rank,comm_ring,ierr)

   call MPI_COMM_RANK(comm_ring,ring_rank,ierr)
   call MPI_COMM_SIZE(comm_ring,ring_size,ierr)
#else
   mb_size = 1
   mb_rank = 0
   bead_size = 1
   bead_rank = 0
   ring_size = 1
   ring_rank = 0
#endif
end subroutine mb_init

subroutine mb_finalize()

   implicit none

#ifdef MPI
   integer :: ierr

   call MPI_COMM_FREE(comm_bead,ierr)
   call MPI_COMM_FREE(comm_ring,ierr)
   call MPI_COMM_FREE(comm_mb,ierr)

   call MPI_FINALIZE(ierr)
#endif

end subroutine mb_finalize

subroutine mb_abort()

   implicit none

#ifdef MPI
#  include "mpif.h"

   integer :: ierr

   call MPI_ABORT(comm_mb,0,ierr)
#endif

end subroutine mb_abort

logical function is_bead_head()

   implicit none

   is_bead_head = bead_rank.eq.0

end function is_bead_head

end module multibead
