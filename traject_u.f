        subroutine traject_u
     x  (ltraj,cfgname,atmnam,idnode,imcon,istraj,keytrj,natms,
     x  nstraj,nstep,tstep,cell,chge,weight,xxx,yyy,zzz,vxx,vyy,
     x  vzz,fxx,fyy,fzz)

c       
c***********************************************************************
c       
c       dl_poly subroutine for writing history file at selected
c       intervals in simulation
c
c       Unformatted, double precision version
c       
c       copyright - daresbury laboratory 1992
c       author    - w. smith dec 1992.
c       
c       wl
c       2001/05/30 12:40:27
c       1.5
c       Exp
c       
c***********************************************************************
c       

#include "dl_params.inc"

        logical newjob,ltraj

        character*80 cfgname
        character*8 atmnam(mxatms)

        dimension cell(9)
        dimension chge(mxatms),weight(mxatms)
        dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
        dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
        dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)

        save newjob
        data newjob/.true./

#ifdef VAMPIR
      call VTBEGIN(163, ierr)
#endif
        if(ltraj.and.idnode.eq.0)then

          if(newjob)  then
            newjob = .false.
c       
c       open the history file if new job or file closed

            open(nhist,file='HISTORY',form='unformatted',
     x       position='append')

          endif
          if(nstep.eq.nstraj)then
          
            write(nhist) cfgname
            write(nhist) dble(natms)
            write(nhist) (atmnam(i),i=1,natms)
            write(nhist) (weight(i),i=1,natms)
            write(nhist) (chge(i),i=1,natms)
          
          endif

          if(mod(nstep-nstraj,istraj).eq.0)then

            write(nhist)dble(nstep),dble(natms),dble(keytrj),
     x           dble(imcon),tstep

            if(imcon.gt.0) write(nhist) cell

            write(nhist) (xxx(i),i = 1,natms)
            write(nhist) (yyy(i),i = 1,natms)
            write(nhist) (zzz(i),i = 1,natms)

            if(keytrj.ge.1)then
              write(nhist) (vxx(i),i = 1,natms)
              write(nhist) (vyy(i),i = 1,natms)
              write(nhist) (vzz(i),i = 1,natms)
            endif
            if(keytrj.ge.2)then
              write(nhist) (fxx(i),i = 1,natms)
              write(nhist) (fyy(i),i = 1,natms)
              write(nhist) (fzz(i),i = 1,natms)
            endif

          endif
c       
c       close history file at regular intervals

          if(.not.newjob.and.mod(nstep,ndump).eq.0)then

            close (nhist)
            newjob=.true.

          endif

        endif

#ifdef VAMPIR
      call VTEND(163, ierr)
#endif
        return
        end
