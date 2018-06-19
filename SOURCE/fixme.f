      subroutine fixme(filename,line)
      use multibead, only: mb_rank, mb_abort
      implicit none
      character(*) :: filename
      integer :: line
      write(5000+mb_rank,*) ' ** FIXME ** at ', filename, line
#if defined(AIX)
      call flush_(5000+mb_rank)
#else
      call flush(5000+mb_rank)
#endif
      stop
      end
