! 13 Sep 06 - IUCHI - EMPLOY HIGHER PRECISION
! 20 OCT 05 - IUCHI - ADD MOVIE OUTPUT PART
!
      subroutine traject
     x  (ltraj,cfgname,atmnam,idnode,imcon,istraj,keytrj,natms,
     x  nstraj,nstep,tstep,cell,chge,weight,xxx,yyy,zzz,vxx,vyy,
     x  vzz,fxx,fyy,fzz)
      use multibead, only: bead_suffix
c     
c***********************************************************************
c     
c     dl_poly subroutine for writing history file at selected
c     intervals in simulation
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith dec 1992.
c
c     wl
c     2001/05/30 12:40:27
c     1.7
c     Exp
c
!     Last updated: Oct 20, 2005 by S. Iuchi
!     
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
      call VTBEGIN(70, ierr)
#endif
      if(ltraj.and.idnode.eq.0)then
        
c     
c     open the history file if new job or file closed
        
        if(newjob)then
          
          newjob = .false.

          open(nhist,file='HISTORY'//bead_suffix,position='append')
! add MOVIE file
          open(nmovi,file='MOVIE'//bead_suffix,position='append')

        endif
        
        if(nstep.eq.nstraj)then
          
          write(nhist,'(a80)') cfgname
          write(nhist,'(3i10)') keytrj,imcon,natms
          
        endif
        
        if(mod(nstep-nstraj,istraj).eq.0)then

          write(nhist,'(a8,4i10,f12.6)') 'timestep',
     x         nstep,natms,keytrj,imcon,tstep

! add MOVIE file
          write(nmovi,'(i10)') natms
          write(nmovi,'(a9,i12)') 'timestep:',nstep

          if(imcon.gt.0) write(nhist,'(3g12.4)') cell

          do i = 1,natms
            write(nhist,'(a8,i10,2f12.6)')
     x        atmnam(i),i,weight(i),chge(i)
            write(nhist,'(1p,3e20.10)') xxx(i),yyy(i),zzz(i)
            if(keytrj.ge.1)then
              write(nhist,'(1p,3e20.10)') vxx(i),vyy(i),vzz(i)
            endif
            if(keytrj.ge.2)then
              write(nhist,'(1p,3e20.10)') fxx(i),fyy(i),fzz(i)
            endif

! add MOVIE file
            write(nmovi,'(a2,2x,1p,3e12.4)') atmnam(i)(1:2),
     $           xxx(i),yyy(i),zzz(i)    

          enddo

        endif

c     
c     close history file at regular intervals
        
        if(.not.newjob.and.mod(nstep,ndump).eq.0)then
          
          close (nhist)
          close (nmovi)  ! for MOVIE file
          newjob = .true.
          
        endif
        
      endif

#ifdef VAMPIR
      call VTEND(70, ierr)
#endif
      return
      end
