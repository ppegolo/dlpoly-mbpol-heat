      subroutine afailed(filename, linenum)

        implicit none

        character(*), intent(in) :: filename
        integer, intent(in) :: linenum
#ifdef MPI
        integer :: ierr
#       include "mpif.h"
#endif

        write (6, fmt='(a,a,a,i0,a)') ' ** Error ** : ', filename,
     x     ':', linenum, ': assert() failed'
#if defined(AIX)
        call flush_(6)
#else
        call flush(6)
#endif
#ifdef MPI
        call MPI_ABORT(MPI_COMM_WORLD,0,ierr)
#else
        stop
#endif

      end subroutine afailed
