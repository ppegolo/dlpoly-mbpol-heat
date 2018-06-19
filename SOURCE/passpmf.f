      subroutine passpmf
     x  (idnode,mxnode,natms,nspmf,listpm,listin,lstpmt,lstpmf,npmf)
c     
c*********************************************************************
c     
c     dl_poly subroutine for passing information about PMF
c     constraints between nodes
c     
c     parallel replicated data version assuming direct node-node
c     connection (i.e. this version may be intel specific)
c     
c     copyright - daresbury laboratory 1995
c     author    - t.forester august 1995.
c     
c     wl
c     2001/05/30 12:40:21
c     1.5
c     $Sate: Exp $
c     
c***********************************************************************
c     

#include "dl_params.inc"

      dimension listpm(mxpmf),listin(mxatms),lstpmt(mxpmf)
      dimension lstpmf(mxspmf,mspmf),npmf(2)
#ifdef VAMPIR
      call VTBEGIN(154, ierr)
#endif
      if(mxproc.lt.mxnode)call error(idnode,102)
      if(mxpmf.lt.natms) call error(idnode,490)

      do i=1,natms
        
        listpm(i)=0
        
      enddo
      
      do k=1,nspmf
        
        do j = 1,npmf(1)+npmf(2)

          i=lstpmf(j,k)
          listpm(i)= 1
          
        enddo

      enddo
c     
c     keep record of all atoms subject to pmf constraints
      
      do i=1,natms
        
        if(listpm(i).gt.0)then
          
          lstpmt(i)=1
          
        else
          
          lstpmt(i)=0
          
        endif
        
      enddo
      
      if(mxnode.gt.1)call gisum(lstpmt,natms,listin)

#ifdef VAMPIR
      call VTEND(154, ierr)
#endif
      return
      end
