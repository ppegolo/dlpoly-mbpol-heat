      subroutine timchk(ktim,time)
      use multibead, only: is_bead_head
c     
c***********************************************************************
c     
c     dlpoly timing routine.
c     time elapsed in seconds
c     makes use of intrinsic function 'etime' for non MPI version
c     
c     copyright daresbury laboratory sept 1994
c     author t.forester sept 1994
c
c     wl
c     2001/08/31 11:13:53
c     1.8
c     Exp
c     
c***********************************************************************
c     
#include "dl_params.inc"
#include "comms.inc"
      real*4 second(2),etime

      save init,tzero

      data init/0/

   10 format(/,' time elapsed since job start = ',f15.3,' seconds',/)

c
c     determine clock time

#ifdef SGISHMEM
      if(init.eq.0) then
        init = 1
        tzero= MPI_wtime()
      endif
      time = MPI_wtime()-tzero
#endif

#ifdef MPI
#ifdef MPIU
#define MPI_wtime MPI_wtime_
#endif
      if(init.eq.0) then
        init = 1
        tzero= MPI_wtime()
      endif
      time = MPI_wtime()-tzero
#endif
#ifdef INTEL
      if(init.eq.0)init=mclock()
      itime=mclock()-init
      time=dble(itime)/1000.d0
#endif
#ifdef CRAYT3D
      if(init.eq.0)then
c       call rtc(time)
c       tzero=time*1.6666666666666d-9
        init=1
        call timef(tzero)
        time = 0.0d0
      endif
c     call rtc(time)
c     time=time*1.6666666666666d-9-tzero
      call timef(tst)
      tst = tst * 1.0d-3
      time = tst - tzero
#endif
#ifdef SERIAL
      time=etime(second)
#endif
      if(ktim.gt.0.and.is_bead_head()) write(nrite,10)time

      return
      end
