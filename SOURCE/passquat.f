      subroutine passquat
     x  (lcnb,idnode,mxnode,natms,ngrp,nscons,ntpmls,listin,
     x  listcon,lstrgd,lstout,lstcsit,lstgtp,nummols,numgrp,numgsit)
c     
c*********************************************************************
c     
c     dl_poly subroutine for passing information about rigid body 
c     atoms involved in bond constraints between nodes
c     
c     parallel replicated data version assuming direct node-node
c     connection
c     
c     copyright - daresbury laboratory 1995
c     author    - t. forester december 1995.
c     
c     wl
c     2001/05/30 12:40:21
c     1.3
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
#include "comms.inc"

      logical lcnb
      
      dimension listin(mxatms)
      dimension listcon(mxcons,3),lstcsit(2*mxcons)
      dimension lstout(mxatms),lstrgd(mxgatm)
      dimension nummols(mxtmls),numgrp(mxtmls),numgsit(mxungp)
      dimension lstgtp(mxgrp)

#ifdef MPI
      dimension status(MPI_STATUS_SIZE)
#endif
      
#ifdef VAMPIR
      call VTBEGIN(155, ierr)
#endif
      if(mxproc.lt.mxnode)call error(idnode,102)
      
c
c     block indices for groups
      
      igrp1 = (idnode*ngrp)/mxnode + 1
      igrp2 = ((idnode+1)*ngrp)/mxnode
      
c
c     locate site indices of atoms in constraints

      do i = 1,natms
        listin(i) = 0
      enddo
c     
c     loop over molecule types

      jr = 0 
      igrp = 0
      do itmols=1,ntpmls
c     
c     loop over molecules in system
        
        do imols=1,nummols(itmols)
c     
c     construct rigid body site list: each processor has a different copy
          
          do lgrp=1,numgrp(itmols)
            
            igrp=igrp+1
            
            if((igrp.ge.igrp1).and.(igrp.le.igrp2)) then
                
              id = lstgtp(igrp)
              do jj = 1,numgsit(id)
                  
                jr = jr +1
                i = lstrgd(jr)
                listin(i) = jj

              enddo
            endif
          enddo
        enddo
      enddo

      if(mxnode.gt.1) call gisum(listin,natms,lstout)

      lcnb = .true.
      ik = 0
      do k=1,nscons
        
        i=listcon(k,2)
        j=listcon(k,3)

        if(listin(i).ne.0) then
          ik = ik + 1
          lstcsit(ik) = listin(i)
          lcnb = .false.
        endif

        if(listin(j).ne.0) then
          ik = ik + 1
          lstcsit(ik) = listin(j)
          lcnb = .false.
        endif

      enddo
c
c     lcnb flags bodies connected by constraints

      if(mxnode.gt.1) call gstate(lcnb)
      lcnb = (.not.lcnb)

#ifdef VAMPIR
      call VTEND(155, ierr)
#endif
      return
      end



