      subroutine cfgscan
     x  (lhead,idnode,mxnode,nconf,imcon,volm,xhi,yhi,zhi,cell,buffer)
      use multibead, only: bead_suffix
c***********************************************************************
c     
c     dl_poly subroutine for scanning the initial configuration
c     file to determine the number of atoms present
c     
c     copyright - daresbury laboratory 1997
c     author    - w. smith  june       1997
c     
c     note: volm is volume containing all particles, not system volume
c
c     $Author: souda $
c     $Date: 2002/06/27 21:45:49 $
c     $Revision: 1.6 $
c     $State: Exp $
c
c***********************************************************************
c     
      
      implicit real*8(a-h,o-z)

      parameter (mega=1000000)
      
      character*256 record
      character*80 header
      character*8 name
      logical lvolm,safe,lhead
      dimension cell(9),celprp(10),buffer(10),extra(4)

      character*8, allocatable :: carr(:)
      integer,     allocatable :: iarr(:)
      real*8,      allocatable :: rarr(:)

      ! sandeep 

      integer,parameter :: nofields=100
      integer           :: countfields,fieldlength
      integer           :: starts(nofields),ends(nofields)
      character(len=256):: tmprecord
      character(len=1):: a


      allocate (carr(20), stat=memcha)
      allocate (iarr(20), stat=memint)
      allocate (rarr(20), stat=memrea)

      if (memcha.ne.0.or.memint.ne.0.or.memrea.ne.0) then
         if (lhead)
     x      write(*,*) 'CFGSCAN: memory allocation error...'
         stop
      endif

      imcon=0
      xhi=0.d0
      yhi=0.d0
      zhi=0.d0
      volm=0.d0
      do i=1,9
        cell(i)=0.d0
      enddo

ccc      if (lhead) then

         if (lhead) open (nconf,file='CONFIG'//bead_suffix)
c     
c        read the CONFIG file header
        
ccc         if (lhead) read(nconf,'(a80)',end=100) header
         call getrec(safe,idnode,mxnode,nconf,record)
         if (.not.safe) goto 999
         header = record(1:80)

         call getrec(safe,idnode,mxnode,nconf,record)
         if (.not.safe) goto 999
         call cal_field(record,nofields,countfields,starts,ends)
    
         tmprecord=record(starts(1):ends(1))
         read(tmprecord,*) levcfg
         tmprecord=record(starts(2):ends(2))
         read(tmprecord,*) imcon

         lvolm=(imcon.eq.0.or.imcon.eq.6)
        
c     
c        specify molecular dynamics simulation cell
        
         if (imcon.gt.0) then

            do i=1,9,3
               call getrec(safe,idnode,mxnode,nconf,record)
               if (.not.safe) goto 999

               call cal_field(record,nofields,countfields,starts,ends)

               tmprecord=record(starts(1):ends(1))
               read(tmprecord,*) cell(i)
               tmprecord=record(starts(2):ends(2))
               read(tmprecord,*) cell(i+1)
               tmprecord=record(starts(3):ends(3))
               read(tmprecord,*) cell(i+2)

            enddo
            call dcell(cell,celprp)
    
         endif
        
         if (.not.lvolm) then
            volm=celprp(10)
            if (imcon.eq.4) then
               volm=0.5d0*celprp(10)
            elseif (imcon.eq.5) then
               volm=0.5d0*celprp(10)
            elseif (imcon.eq.7) then
               volm=0.5d0*celprp(10)
            endif
         endif

         do i=1,mega

c           Read atom name
c           --------------
            call getrec(safe,idnode,mxnode,nconf,record)
            if (.not.safe) goto 100

            call cal_field(record,nofields,countfields,starts,ends)
            
            name = record(starts(1):ends(1))

c           Read cartesian coordinates
c           --------------------------
            call getrec(safe,idnode,mxnode,nconf,record)
            if (.not.safe) goto 100

             call cal_field(record,nofields,countfields,starts,ends)

             tmprecord=record(starts(1):ends(1))
             read(tmprecord,*) xxx
             tmprecord=record(starts(2):ends(2))
             read(tmprecord,*) yyy
             tmprecord=record(starts(3):ends(3))
             read(tmprecord,*) zzz


            if (levcfg.gt.0) then

c              Read components of velocities
c              -----------------------------
               call getrec(safe,idnode,mxnode,nconf,record)
               if (.not.safe) goto 100

               call cal_field(record,nofields,countfields,starts,ends)

               tmprecord=record(starts(1):ends(1))
               read(tmprecord,*) uuu
               tmprecord=record(starts(2):ends(2))
               read(tmprecord,*) vvv
               tmprecord=record(starts(3):ends(3))
               read(tmprecord,*) www

            endif

            if (levcfg.gt.1) then

c              Read components of forces
c              -------------------------
               call getrec(safe,idnode,mxnode,nconf,record)
               if (.not.safe) goto 100
          
               call cal_field(record,nofields,countfields,starts,ends)

               tmprecord=record(starts(1):ends(1))
               read(tmprecord,*) uuu
               tmprecord=record(starts(2):ends(2))
               read(tmprecord,*) vvv
               tmprecord=record(starts(3):ends(3))
               read(tmprecord,*) www

            endif

            if (lvolm) then
               if (i.eq.1) then
                  xhi = dabs(xxx)
                  yhi = dabs(yyy)
                  zhi = dabs(zzz)
               else
                 xhi = max(xhi,dabs(xxx))
                 yhi = max(yhi,dabs(yyy))
                 zhi = max(zhi,dabs(zzz))
               endif
            endif
          
         enddo
        
  100    continue

         if (imcon.eq.0) then

            volm=8.d0*xhi*yhi*zhi

         elseif (imcon.eq.6) then

            coz=(cell(1)*cell(4)+cell(2)*cell(5)+cell(3)*cell(6))/
     x      (celprp(1)*celprp(2))
            volm=2.d0*zhi*celprp(1)*celprp(2)*sqrt(1.d0-coz**2)

         endif

         if(lhead)close (nconf)

ccc      endif

c      extra(1)=dble(imcon)
c      extra(2)=xhi
c      extra(3)=yhi
c      extra(4)=zhi
c      call gdsum(extra,4,buffer)
c      call gdsum(volm,1,buffer)
c      call gdsum(cell,9,buffer)
c      imcon=nint(extra(1))
c      xhi=extra(2)
c      yhi=extra(3)
c      zhi=extra(4)

c      if (lhead) then
c         write(*,*) 'From cfgscan:'
c	 write(*,*) 'cell:',cell
c	 write(*,*) 'imcon:',imcon
c	 write(*,*) 'volm:',volm
c	 write(*,*) 'xhi,yhi,zhi:',xhi,yhi,zhi
c      endif

      deallocate (carr, stat=memcha)
      deallocate (iarr, stat=memint)
      deallocate (rarr, stat=memrea)
      if (memcha.ne.0.or.memint.ne.0.or.memrea.ne.0) then
         if (lhead)
     x      write(*,*) 'CFGSCAN: memory deallocation error...'
         stop
      endif
      return

  999 continue
      deallocate (carr, stat=memcha)
      deallocate (iarr, stat=memint)
      deallocate (rarr, stat=memrea)
      call error(idnode,55)

      end
