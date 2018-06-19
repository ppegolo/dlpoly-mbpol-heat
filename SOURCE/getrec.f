      subroutine getrec(safe,idnode,mxnode,ifile,record)
c     
c*********************************************************************
c     
c     dl_poly subroutine to read a character string on one node
c     and broadcast it to all other nodes
c     
c     copyright daresbury laboratory 1994
c     author w.smith december 1994
c     
c     wl
c     1996/02/15 14:32:56
c     1.1.1.1
c     Exp
c     
c*********************************************************************
c     
      
      character*256 record
      common/recget/export(256),import(256)
      integer export,import,idnode,mxnode,ifile
      logical safe
      
      safe=.true.
      
      call gsync()
      
      if(idnode.eq.0)then
        
        read(ifile,'(a256)',end=100)record
        
        do i=1,256
          
          export(i)=ichar(record(i:i))
          
        enddo
        
        call gstate(safe)
        call gisum(export,256,import)
        
        return
        
  100   safe=.false.
        
        call gstate(safe)
        
      else
        
        call gstate(safe)
        if(.not.safe)return

        do i=1,256

          export(i)=0

        enddo

        call gisum(export,256,import)
        
        do i=1,256
          
          record(i:i)=char(export(i))
          
        enddo
        
        return
        
      endif
      
      end






