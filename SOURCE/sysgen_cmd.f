! 07 NOV 05 - IUCHI - SET KEYRES=3 IF BELOW CASE
! 20 OCT 05 - IUCHI - SET KEYRES=2 IF LEVCFG>=1 AND KEYRES=0
! 
      subroutine sysgen_cmd
     x  (lhead,loglnk,lneut,lopt2,atmnam,cfgname,sitnam,idnode,imcon,
     x  mbnrg_list,mbnrg_nmol,mbnrg_index,mbnrg_key,
     x  keyens,keyfce,keyres,levcfg,multt,mxnode,ntpmls,delr,rcut,
     x  volm,lfzsit,lstfrz,ltpsit,ltype,lpolar,nummols,numsit,nugrp,
     x  lstneu,buffer,cell,chge,chgsit,polr,polr2,polarsit,polarsit2,
     x  fxx,fyy,fzz,vxx,vyy,vzz,
     x  weight,wgtsit,xxx,yyy,zzz)
      use multibead, only: bead_suffix
c     
c***********************************************************************
c     
c     dl_poly subroutine for reading the configuration data file
c     
c     copyright - daresbury laboratory 1997
c     author    - w. smith         july 1997
c     
c     $Author: souda $
c     $Date: 2002/06/27 21:45:59 $
c     $Revision: 1.8 $
c     $State: Exp $
!
!     Last updated: 21 Oct 2005 by S. Iuchi
c     
c***********************************************************************
c     

#include "dl_params.inc"

      character*256 record
      character*80 cfgname
      character*8 tmp_atname 
      character*8  atmnam(mxatms)
      character*8  atname,sitnam(mxsite)

      logical lhead
      logical loglnk,safe,lneut,lopt2,lpolar

      dimension lstneu(mxatms),nugrp(mxsite),lfzsit(mxsite)
      dimension cell(9),celprp(10),lstfrz(mxatms)
      dimension nummols(mxtmls),numsit(mxtmls)
      dimension ltype(mxatms),ltpsit(mxsite),wgtsit(mxsite)
      dimension chge(mxatms),weight(mxatms),chgsit(mxsite)
      dimension polr(mxatms),polarsit(mxsite)
      dimension polr2(mxatms),polarsit2(mxsite)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension buffer(mxbuff)

      character*8, allocatable :: carr(:)
      integer,     allocatable :: iarr(:)
      real*8,      allocatable :: rarr(:)

c mbnrg
c    second index is equal to 1 -- temporary change it if molecule has
c    more than one atom
      integer :: mbnrg_index(mxtmls)
      integer :: mbnrg_list(mxtmls,1),mbnrg_nmol(mxtmls)

      ! sandeep 

      integer,parameter :: nofields=100
      integer           :: countfields,fieldlength,mbnrg_key
      integer           :: starts(nofields),ends(nofields)
      character(len=256):: tmprecord

      allocate (carr(20), stat=memcha)
      allocate (iarr(20), stat=memint)
      allocate (rarr(20), stat=memrea)
      if (memcha.ne.0.or.memint.ne.0.or.memrea.ne.0) then
         if (lhead) write(*,*) 'SYSGEN: memory allocation error..'
         stop
      endif

#ifdef VAMPIR
      call VTBEGIN(4, ierr)
#endif

c
c     open the system input file

      if (.true.) open (nconf,file='CONFIG'//bead_suffix)
c     
c     read the CONFIG file header

c$$$  read(nconf,'(a80)',end=100) cfgname
      call getrec(safe,idnode,mxnode,nconf,record)
      if (.not.safe) goto 100

      cfgname=record(1:80)
      if(lhead) write(nrite,
     x  "(/,1x,'configuration file name: ',/,/,10x,a80)") cfgname

c$$$  read(nconf,'(2i10)',end=100) levcfg,imcon
      call getrec(safe,idnode,mxnode,nconf,record)
      if (.not.safe) goto 100
      call cal_field(record,nofields,countfields,starts,ends)

      tmprecord=record(starts(1):ends(1))
      read(tmprecord,*) levcfg
      tmprecord=record(starts(2):ends(2))
      read(tmprecord,*) imcon

c      call readrec(record,'ii ',carr,iarr,rarr)
c      levcfg =iarr(1)
c      imcon  =iarr(2)
      if(lhead) write(nrite,
     x  "(/,/,1x,'selected image convention',6x,i10)") imcon

      if((imcon.eq.0.or.imcon.eq.6).and.
     x  (keyfce/2.eq.1.or.keyfce/2.eq.6))
     x  call error(idnode,180)

      if(imcon.eq.0.and.(.not.lneut).and.(keyfce.gt.1)
     x  .and.(multt.eq.1)) call warning(idnode,30,0.d0,0.d0,0.d0)

      if(imcon.eq.0.and.(keyens.ge.4.and.keyens.le.7))
     x  call error(idnode,390)
      if(imcon.le.2.and.(keyens.eq.6.or.keyens.eq.7)) imcon = 3

c     
c     specify molecular dynamics simulation cell

      if(imcon.eq.0)then

c     if no periodic boundaries - set zero values for cell 
c     vectors and cell volume
        
        do i=1,9
          cell(i)=0.d0
        enddo

        volm=0.d0

      else
c     
c     read and write cell vectors

c$$$  read(nconf,'(3f20.0)',end=100) cell

        do i=1,9,3

           call getrec(safe,idnode,mxnode,nconf,record)
           if (.not.safe) goto 100
           call cal_field(record,nofields,countfields,starts,ends)

           tmprecord=record(starts(1):ends(1))
           read(tmprecord,*) cell(i)
           tmprecord=record(starts(2):ends(2))
           read(tmprecord,*) cell(i+1)
           tmprecord=record(starts(3):ends(3))
           read(tmprecord,*) cell(i+2)


c           call readrec(record,'rrr ',carr,iarr,rarr)

c           do k=1,3
c              cell(i+k-1) = rarr(k)
c           enddo

        enddo

        if(lhead)then
          write(nrite,"(/,/,1x,'simulation cell vectors'/,/)")
          write(nrite,"(21x,3f12.6)") cell
        endif
c     
c     check integrity of cell vectors : for cubic, TO and RD cases
c     ie. cell(1)=cell(5)=cell(9) (or cell(9)/sqrt(2) for RD)

        if((imcon.eq.1).or.(imcon.eq.4).or.(imcon.eq.5)) then

          axx = (abs(cell(1))+abs(cell(5)))/2.d0
          test = 1.d-8*axx
          if(abs(cell(1)-axx).gt.test) call error(idnode,410)
          if(abs(cell(5)-axx).gt.test) call error(idnode,410)
          if(imcon.eq.5)then
            if(abs(cell(9)-axx*sqrt(2.d0)).gt.test) 
     x        call error(idnode,410)
          else
            if(abs(cell(9)-axx).gt.test) call error(idnode,410)
          endif

        endif
c
c     check integrity of hexagonal prism cell vectors

        if(imcon.eq.7)then

          rt3=sqrt(3.d0)
          if(abs(cell(1)-rt3*cell(5)).ge.1.d-6)
     x    call error(idnode,410)

        endif
c     
c     check for diagonal cell matrix if appropriate

        if((imcon.eq.1).or.(imcon.eq.2).or.(imcon.eq.4).or.
     x    (imcon.eq.5).or.(imcon.eq.7)) then

          if(cell(2).ne.0.d0) call error(idnode,410)
          if(cell(3).ne.0.d0) call error(idnode,410)
          if(cell(4).ne.0.d0) call error(idnode,410)
          if(cell(6).ne.0.d0) call error(idnode,410)
          if(cell(7).ne.0.d0) call error(idnode,410)
          if(cell(8).ne.0.d0) call error(idnode,410)

        endif
c     
c     calculate dimensional properties of simulation cell

        call dcell(cell,celprp)

        volm=celprp(10)

        if(imcon.eq.4)then

          volm=0.5d0*celprp(10)

        elseif(imcon.eq.5)then

          volm=0.5d0*celprp(10)

        elseif(imcon.eq.7)then

          volm=0.5d0*celprp(10)

        endif

      endif

      if(lhead) write(nrite,
     x  "(/,/,1x,'system volume     ',2x,1p,g22.12)") volm
c     
c     check value of cutoff and reset if necessary

      if(imcon.gt.0)then

        width=min(celprp(7),celprp(8),celprp(9))/2.d0
        if(imcon.eq.4)width=sqrt(3.d0)*cell(1)/4.d0
        if(imcon.eq.5)width=cell(1)/2.d0
        if(imcon.eq.6)width=min(celprp(7),celprp(8))/2.d0
c     
c     halt program if potential cutoff exceeds cell width

        if(rcut.gt.width) call error(idnode,95)

      endif

c     
c     check on validity of config file contents

      if(keyres.gt.0.and.levcfg.lt.1)call error(idnode,85)
!      if(levcfg.ge.1.and.keyres.eq.0) keyres = 2   
      if(levcfg.ge.1.and.keyres.eq.0) keyres = 3   

      indatm=0
      indnam=0
      indneu=0
      tmp_index=0
      safe=.true.

      do k=1,ntpmls
c mbnrg 
c    right now, i am not using mbnrg_nmol variable later on, but it
c    should be taken care when it has more than one atom in the molecule
        if(mbnrg_index(k)>0) then
           mbnrg_nmol(mbnrg_index(k))=nummols(k)
        endif 

        do l=1,nummols(k)

          do m=1,numsit(k)

            indatm=indatm+1

            if (indatm.gt.mxatms) call error(idnode,45)

            xxx(indatm)  = 0.d0
            yyy(indatm)  = 0.d0
            zzz(indatm)  = 0.d0
            vxx(indatm)  = 0.d0
            vyy(indatm)  = 0.d0
            vzz(indatm)  = 0.d0
            fxx(indatm)  = 0.d0
            fyy(indatm)  = 0.d0
            fzz(indatm)  = 0.d0
            xcoord = 0.d0
            ycoord = 0.d0
            zcoord = 0.d0
            xveloc = 0.d0
            yveloc = 0.d0
            zveloc = 0.d0
            xforce = 0.d0
            yforce = 0.d0
            zforce = 0.d0

c           Read atom name
c           --------------
            call getrec(safe,idnode,mxnode,nconf,record)
            if (.not.safe) goto 100
            call cal_field(record,nofields,countfields,starts,ends)
            
            atname = record(starts(1):ends(1))

c           Read cartesian coordinates
c           --------------------------
            call getrec(safe,idnode,mxnode,nconf,record)
            if (.not.safe) goto 100
            if (.true.) then
           
              call cal_field(record,nofields,countfields,starts,ends)

              tmprecord=record(starts(1):ends(1))
              read(tmprecord,*) xcoord
              tmprecord=record(starts(2):ends(2))
              read(tmprecord,*) ycoord
              tmprecord=record(starts(3):ends(3))
              read(tmprecord,*) zcoord

            endif


            if (levcfg.gt.0) then

c              Read components of velocities
c              -----------------------------
               call getrec(safe,idnode,mxnode,nconf,record)
               if (.not.safe) goto 100

               call cal_field(record,nofields,countfields,starts,ends)
               if (.true.) then
c                  call readrec(record,'rrr ',carr,iarr,rarr)

                  tmprecord=record(starts(1):ends(1))
                  read(tmprecord,*) xveloc
                  tmprecord=record(starts(2):ends(2))
                  read(tmprecord,*) yveloc
                  tmprecord=record(starts(3):ends(3))
                  read(tmprecord,*) zveloc
 
               endif

            endif

            if (levcfg.gt.1) then

c              Read components of forces
c              -------------------------
               call getrec(safe,idnode,mxnode,nconf,record)
               if (.not.safe) goto 100
               call cal_field(record,nofields,countfields,starts,ends)

               if (.true.) then
c                  call readrec(record,'rrr ',carr,iarr,rarr)

                  tmprecord=record(starts(1):ends(1))
                  read(tmprecord,*) xforce
                  tmprecord=record(starts(2):ends(2))
                  read(tmprecord,*) yforce
                  tmprecord=record(starts(3):ends(3))
                  read(tmprecord,*) zforce
 
                endif

            endif

            if (atname.eq.sitnam(indnam+m)) then

               xxx(indatm) = xcoord
               yyy(indatm) = ycoord
               zzz(indatm) = zcoord

c SR: mbnrg  collect molecule or ion coordinates in a separate array 
c WARNING: right now, it assumes only one atom and 1 molecule If it
c changes, this part must be changed 
               if(mbnrg_index(k)>0) then
!                 tmp_index=tmp_index+1 
!                    write(*,*) 'TEST',indatm,mbnrg_index(1)
!                 mbnrg_list(mbnrg_index(k),1)=indatm
                  mbnrg_list(1,1)=indatm
                  tmp_atname=trim(adjustl(atname))
                  call lowcase(tmp_atname,80)
! MRR & DZ - Adding multiple mbnrg ions
                  if(tmp_atname=="f") then

!                      mbnrg_key=1
                      mbnrg_index(k)=1

                    elseif(tmp_atname=="cl") then
!                      mbnrg_key=2
                      mbnrg_index(k)=2

                    elseif(tmp_atname=="br") then

!                      mbnrg_key=3
                      mbnrg_index(k)=3

                    elseif(tmp_atname=="i") then  

!                     mbnrg_key=4
                      mbnrg_index(k)=4

                    elseif(tmp_atname=="li") then  

!                      mbnrg_key=5
                      mbnrg_index(k)=5

                    elseif(tmp_atname=="na") then  
!                      mbnrg_key=6
                      mbnrg_index(k)=6

                    elseif(tmp_atname=="k") then  

!                      mbnrg_key=7
                     mbnrg_index(k)=7

                    elseif(tmp_atname=="rb") then  

!                      mbnrg_key=8
                     mbnrg_index(k)=8

                    elseif(tmp_atname=="cs") then  

!                      mbnrg_key=9
                      mbnrg_index(k)=9

! END MRR & DZ

                    else

                     write(*,*) 'Polynomial for MB-nrg is required, but
     x               it is not available for requested atom type'
                     call exit(-1) 

                  endif
                endif 

               if (levcfg.gt.0) then
                  vxx(indatm) = xveloc
                  vyy(indatm) = yveloc
                  vzz(indatm) = zveloc
               endif

               if (levcfg.gt.1) then
                  fxx(indatm) = xforce
                  fyy(indatm) = yforce
                  fzz(indatm) = zforce
               endif

            else

               write(nrite,"(/,/,'unidentified atom label :',a8,
     x                     ': atom number ',i5)") atname,indatm
               safe=.false.

            endif

            call gstate(safe)
            if(.not.safe) call error(idnode,25)

            ltype(indatm)  = ltpsit(indnam+m)
            weight(indatm) = wgtsit(indnam+m)
            chge(indatm)   = chgsit(indnam+m)
            if (lpolar) then
               polr(indatm)=polarsit(indnam+m)
               polr2(indatm)=polarsit2(indnam+m)
            endif
            atmnam(indatm) = sitnam(indnam+m)
            lstfrz(indatm) = lfzsit(indnam+m)
            if(lneut) lstneu(indatm) = nugrp(indnam+m)+indneu

          enddo

          indneu = indneu + nugrp(indnam+numsit(k))

        enddo

        indnam = indnam + numsit(k)
        tmp_index=0 

      enddo

      if (.false.) then !VB getrec() broadcasts this

         call gdsum(xxx,indatm,buffer)
         call gdsum(yyy,indatm,buffer)
         call gdsum(zzz,indatm,buffer)

         if (levcfg.gt.0) then
            call gdsum(vxx,indatm,buffer)
            call gdsum(vyy,indatm,buffer)
            call gdsum(vzz,indatm,buffer)
         endif

         if (levcfg.gt.1) then
            call gdsum(fxx,indatm,buffer)
            call gdsum(fyy,indatm,buffer)
            call gdsum(fzz,indatm,buffer)
         endif

      endif

c     
c     decide on whether to use link cells for verlet list constructor

c
c     set widths if unset - needed for check on link cells below

      if(imcon.eq.0.or.imcon.eq.6) then

        xhi=dabs(xxx(1))
        yhi=dabs(yyy(1))
        zhi=dabs(zzz(1))

        do i=2,indatm
          xhi=max(xhi,dabs(xxx(i)))
          yhi=max(yhi,dabs(yyy(i)))
          zhi=max(zhi,dabs(zzz(i)))
        enddo

        if (imcon.eq.0) then
           cell(1) = max(2.d0*xhi+rcut+delr,3.d0*(rcut+delr))
           cell(5) = max(2.d0*yhi+rcut+delr,3.d0*(rcut+delr))
           cell(9) = max(2.d0*zhi+rcut+delr,3.d0*(rcut+delr))
        endif

        if (imcon.eq.6)then
           cell(9) = max(2.d0*zhi+rcut+delr,3.d0*(rcut+delr),cell(9))
        endif

        call dcell(cell,celprp)
        volm=celprp(10)

      endif

c      loglnk = .true.
      loglnk = .false. 

      ilx = max(3,int(celprp(7)/(rcut+delr)))
      ily = max(3,int(celprp(8)/(rcut+delr)))
      ilz = max(3,int(celprp(9)/(rcut+delr)))
      ncells = ilx*ily*ilz
c      if (ncells.eq.27.and.mxnode.eq.1) loglnk=.false.
c      if (lneut.and.ncells.le.36.and.mxnode.eq.1) loglnk=.false.
c      if (imcon.eq.4.or.imcon.eq.5.or.imcon.eq.7) loglnk = .false.
      if(lopt2) loglnk=.false.
      if (loglnk.and.ncells.gt.mxcell) then
         dum1=dble(ncells)
         dum2=dble(mxcell)
         call warning(idnode,90,dum1,dum2,dum2)
         loglnk=.false.
      endif

      if(loglnk.and.lhead) 
     x  write(nrite,"(/,/,' link cell algorithm in use')")

      if (lhead) close (nconf)

#ifdef VAMPIR
      call VTEND(4, ierr)
#endif

      deallocate (carr, stat=memcha)
      deallocate (iarr, stat=memint)
      deallocate (rarr, stat=memrea)
      if (memcha.ne.0.or.memint.ne.0.or.memrea.ne.0) then
         if (lhead)
     x      write(*,*) 'SYSDEF: memory deallocation error...'
         stop
      endif

      return
c     
c     error exit for config file read

  100 continue
      deallocate (carr, stat=memcha)
      deallocate (iarr, stat=memint)
      deallocate (rarr, stat=memrea)
      call error(idnode,55)

      end
