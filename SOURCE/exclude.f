      subroutine exclude
     x    (idnode,mxnode,natms,ntpmls,keybnd,keyang,lexatm,lexsit,
     x    lstang,lstbnd,lstcon,lstdih,lstinv,nexatm,nexsit,numang,
     x    numbonds,numcon,numdih,numinv,numsit,numgrp,listyp,lstgst,
     x    numgsit,numshl,lstshl,prmdih)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for constructing the excluded pair
c     interaction list of the system to be simulated
c     
c     keybnd < 0 distance restraint so not excluded
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith        june 1992
c     
c     rigid body exclusions added : t.forester nov 1993
c     check on 1..4 scale factors : t.forester feb 1994
c     inversion terms added       : w.smith    jul 1996
c     
c     wl
c     2001/05/30 12:40:04
c     1.6
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical newjob,check,safe,lchk
      
      dimension numsit(mxtmls),numcon(mxtmls),numang(mxtmls)
      dimension numbonds(mxtmls),lstcon(mxtcon,2)
      dimension lstang(mxtang,3),lexatm(msatms,mxexcl),nexatm(msatms)
      dimension lstbnd(mxtbnd,2),lexsit(mxsite,mxexcl),nexsit(mxsite)
      dimension keybnd(mxtbnd),numgrp(mxtmls),numshl(mxtmls)
      dimension numdih(mxtmls),lstdih(mxtdih,4),lstshl(mxtshl,2)
      dimension lstgst(mxungp,mxngp),listyp(mxungp),numgsit(mxungp)
      dimension prmdih(mxtdih,mxpdih)
      dimension keyang(mxtang),numinv(mxtmls),lstinv(mxtinv,4)
      
      save newjob
      
      data newjob/.true./
      
#ifdef VAMPIR
      call VTBEGIN(77, ierr)
#endif
      if(newjob)then
c     
c     check on array allocations
        
        nsatms=(natms+mxnode-1)/mxnode
        if(nsatms.gt.msatms)then
          
          if(idnode.eq.0) write(nrite,*) 'make msatms >= ',nsatms
          call error(idnode,100)
          
        endif
        
        newjob=.false.
        
      endif
c
c     variables for array bound checking

      ibig = 0
      safe = .true.
c     
c     initialise excluded atom arrays
      
      do i=1,mxsite
        
        nexsit(i)=0
        
      enddo
      
      do i=1,msatms
        
        nexatm(i)=0
        
      enddo

      do j=1,mxexcl
        
        do i=1,mxsite
          
          lexsit(i,j)=0
          
        enddo
        
        do i=1,msatms
          
          lexatm(i,j)=0
          
        enddo
        
      enddo
      
c     
c     loop over molecules in system
      
      ibonds=0
      iangle=0
      iconst=0
      idihdr=0
      invers=0
      igrp  =0
      isite =0
      ishels=0
      
      do itmols=1,ntpmls
        
c     
c     exclude sites on basis of chemical bonds
        
        do i=1,numbonds(itmols)
          
          ibonds=ibonds+1
          
          if(keybnd(ibonds).gt.0) then
            
            ia=lstbnd(ibonds,1)+isite
            ib=lstbnd(ibonds,2)+isite
c
c     check interaction not already included

            lchk = .true.
            do jz = 1,min(nexsit(ia),mxexcl)
              if(lexsit(ia,jz).eq.ib-isite) lchk=.false.
            enddo
            if(lchk) then
              nexsit(ia)=nexsit(ia)+1
              nexsit(ib)=nexsit(ib)+1
              if(max(nexsit(ia),nexsit(ib)).gt.mxexcl) then
                ibig = max(ibig,nexsit(ia),nexsit(ib))
                safe = .false.
              else
                lexsit(ia,nexsit(ia))=ib-isite
                lexsit(ib,nexsit(ib))=ia-isite
              endif
            endif
            
          endif
          
        enddo
        
c     
c     exclude sites on basis of bond constraints
        
        do i=1,numcon(itmols)
          
          iconst=iconst+1
          ia=lstcon(iconst,1)+isite
          ib=lstcon(iconst,2)+isite
c
c     check interaction not already included

          lchk = .true.
          do jz = 1,min(nexsit(ia),mxexcl)
            if(lexsit(ia,jz).eq.ib-isite) lchk=.false.
          enddo
          if(lchk) then
            nexsit(ia)=nexsit(ia)+1
            nexsit(ib)=nexsit(ib)+1
            if(max(nexsit(ia),nexsit(ib)).gt.mxexcl) then
              ibig = max(ibig,nexsit(ia),nexsit(ib))
              safe = .false.
            else
              lexsit(ia,nexsit(ia))=ib-isite
              lexsit(ib,nexsit(ib))=ia-isite
            endif
          endif
          
        enddo
        
c     
c     exclude sites on basis of bond angles
        
        do i=1,numang(itmols)
          
          iangle=iangle+1
          if(keyang(iangle).gt.0) then
            ia=lstang(iangle,1)+isite
            ib=lstang(iangle,2)+isite
            ic=lstang(iangle,3)+isite
c
c     check if already added to lists ..
c     ia - ib interaction
            lchk = .true.
            do jz = 1,min(nexsit(ia),mxexcl)
              if(lexsit(ia,jz).eq.ib-isite) lchk=.false.
            enddo
            if(lchk) then
              nexsit(ia)=nexsit(ia)+1
              nexsit(ib)=nexsit(ib)+1
              if(max(nexsit(ia),nexsit(ib)).gt.mxexcl) then
                ibig = max(ibig,nexsit(ia),nexsit(ib))
                safe = .false.
              else
                lexsit(ia,nexsit(ia))=ib-isite
                lexsit(ib,nexsit(ib))=ia-isite
              endif
            endif
c
c     ib - ic interaction

            lchk = .true.
            do jz = 1,min(nexsit(ib),mxexcl)
              if(lexsit(ib,jz).eq.ic-isite) lchk=.false.
            enddo
            if(lchk) then
              nexsit(ib)=nexsit(ib)+1
              nexsit(ic)=nexsit(ic)+1
              if(max(nexsit(ib),nexsit(ic)).gt.mxexcl) then
                ibig = max(ibig,nexsit(ib),nexsit(ic))
                safe = .false.
              else
                lexsit(ib,nexsit(ib))=ic-isite
                lexsit(ic,nexsit(ic))=ib-isite
              endif
            endif
c
c     ia - ic interaction

            lchk = .true.
            do jz = 1,min(nexsit(ia),mxexcl)
              if(lexsit(ia,jz).eq.ic-isite) lchk=.false.
            enddo
            if(lchk) then
              nexsit(ia)=nexsit(ia)+1
              nexsit(ic)=nexsit(ic)+1
              if(max(nexsit(ia),nexsit(ic)).gt.mxexcl) then
                ibig = max(ibig,nexsit(ia),nexsit(ic))
                safe = .false.
              else
                lexsit(ia,nexsit(ia))=ic-isite
                lexsit(ic,nexsit(ic))=ia-isite
              endif
            endif

          endif
          
        enddo
        
c     
c     exclude on basis of rigid groups
        
        do i = 1,numgrp(itmols)
          
          igrp = igrp + 1
          
          id = listyp(igrp)
          
          do jj = 1,numgsit(id)-1
            
            ia = lstgst(igrp,jj)+isite
            
            do jk = jj+1,numgsit(id)
              
              ib = lstgst(igrp,jk)+isite
c
c     check interaction not already included

              lchk = .true.
              do jz = 1,min(nexsit(ia),mxexcl)
                if(lexsit(ia,jz).eq.ib-isite) lchk=.false.
              enddo
              if(lchk) then
                nexsit(ia)=nexsit(ia)+1
                nexsit(ib)=nexsit(ib)+1
                if(max(nexsit(ia),nexsit(ib)).gt.mxexcl) then
                  ibig = max(ibig,nexsit(ia),nexsit(ib))
                  safe = .false.
                else
                  lexsit(ia,nexsit(ia))=ib-isite
                  lexsit(ib,nexsit(ib))=ia-isite
                endif
              endif
              
            enddo
            
          enddo
          
        enddo
        
c     
c     exclude sites on basis of 1-4 dihedral angles
        
        do i=1,numdih(itmols)
          
          idihdr=idihdr+1
          ia=lstdih(idihdr,1)+isite
          ib=lstdih(idihdr,2)+isite
          ic=lstdih(idihdr,3)+isite
          id=lstdih(idihdr,4)+isite

c
c     check if already added to lists ..
c     ia - ib interaction

          lchk = .true.
          do jz = 1,min(nexsit(ia),mxexcl)
            if(lexsit(ia,jz).eq.ib-isite) lchk=.false.
          enddo
          if(lchk) then
            nexsit(ia)=nexsit(ia)+1
            nexsit(ib)=nexsit(ib)+1
            if(max(nexsit(ia),nexsit(ib)).gt.mxexcl) then
              ibig = max(ibig,nexsit(ia),nexsit(ib))
              safe = .false.
            else
              lexsit(ia,nexsit(ia))=ib-isite
              lexsit(ib,nexsit(ib))=ia-isite
            endif
          endif
c
c     ib - ic interaction

          lchk = .true.
          do jz = 1,min(nexsit(ib),mxexcl)
            if(lexsit(ib,jz).eq.ic-isite) lchk=.false.
          enddo
          if(lchk) then
            nexsit(ib)=nexsit(ib)+1
            nexsit(ic)=nexsit(ic)+1
            if(max(nexsit(ib),nexsit(ic)).gt.mxexcl) then
              ibig = max(ibig,nexsit(ib),nexsit(ic))
              safe = .false.
            else
              lexsit(ib,nexsit(ib))=ic-isite
              lexsit(ic,nexsit(ic))=ib-isite
            endif
          endif
c
c     ia - ic interaction

          lchk = .true.
          do jz = 1,min(nexsit(ia),mxexcl)
            if(lexsit(ia,jz).eq.ic-isite) lchk=.false.
          enddo
          if(lchk) then
            nexsit(ia)=nexsit(ia)+1
            nexsit(ic)=nexsit(ic)+1
            if(max(nexsit(ia),nexsit(ic)).gt.mxexcl) then
              ibig = max(ibig,nexsit(ia),nexsit(ic))
              safe = .false.
            else
              lexsit(ia,nexsit(ia))=ic-isite
              lexsit(ic,nexsit(ic))=ia-isite
            endif
          endif
c
c     id - ib interaction

            lchk = .true.
            do jz = 1,min(nexsit(id),mxexcl)
              if(lexsit(id,jz).eq.ib-isite) lchk=.false.
            enddo
            if(lchk) then
              nexsit(id)=nexsit(id)+1
              nexsit(ib)=nexsit(ib)+1
              if(max(nexsit(id),nexsit(ib)).gt.mxexcl) then
                ibig = max(ibig,nexsit(id),nexsit(ib))
                safe = .false.
              else
                lexsit(id,nexsit(id))=ib-isite
                lexsit(ib,nexsit(ib))=id-isite
              endif
            endif
c
c     id - ic interaction

            lchk = .true.
            do jz = 1,min(nexsit(id),mxexcl)
              if(lexsit(id,jz).eq.ic-isite) lchk=.false.
            enddo
            if(lchk) then
              nexsit(id)=nexsit(id)+1
              nexsit(ic)=nexsit(ic)+1
              if(max(nexsit(id),nexsit(ic)).gt.mxexcl) then
                ibig = max(ibig,nexsit(id),nexsit(ic))
                safe = .false.
              else
                lexsit(id,nexsit(id))=ic-isite
                lexsit(ic,nexsit(ic))=id-isite
              endif
            endif
c
c     ia - id interaction: may need to reset vdw and elec scale factors

            lchk = .true.
            do jz = 1,min(nexsit(ia),mxexcl)
              if(lexsit(ia,jz).eq.id-isite) lchk=.false.
            enddo
            if(lchk) then
              nexsit(ia)=nexsit(ia)+1
              nexsit(id)=nexsit(id)+1
              if(max(nexsit(ia),nexsit(id)).gt.mxexcl) then
                ibig = max(ibig,nexsit(ia),nexsit(id))
                safe = .false.
              else
                lexsit(ia,nexsit(ia))=id-isite
                lexsit(id,nexsit(id))=ia-isite
              endif

            else
c     
c     if already excluded reset 1..4 vdw and coulombic scale factors
            
            check = ((prmdih(idihdr,4).ne.0.d0).or.
     x        (prmdih(idihdr,5).ne.0.d0))
            
            if(check) then
              
              a1 = dble(itmols)
              a2 = dble(ia)
              a3 = dble(id)
              call warning(idnode,20,a1,a2,a3)
              
              prmdih(idihdr,4) =0.d0
              prmdih(idihdr,5) =0.d0
              
            endif
            
          endif
          
        enddo

c     
c     exclude sites on basis of inversion potentials
        
        do i=1,numinv(itmols)
          
          invers=invers+1
          ia=lstinv(invers,1)+isite
          ib=lstinv(invers,2)+isite
          ic=lstinv(invers,3)+isite
          id=lstinv(invers,4)+isite

c
c     check if already added to lists ..
c     ia - ib interaction

          lchk = .true.
          do jz = 1,min(nexsit(ia),mxexcl)
            if(lexsit(ia,jz).eq.ib-isite) lchk=.false.
          enddo
          if(lchk) then
            nexsit(ia)=nexsit(ia)+1
            nexsit(ib)=nexsit(ib)+1
            if(max(nexsit(ia),nexsit(ib)).gt.mxexcl) then
              ibig = max(ibig,nexsit(ia),nexsit(ib))
              safe = .false.
            else
              lexsit(ia,nexsit(ia))=ib-isite
              lexsit(ib,nexsit(ib))=ia-isite
            endif
          endif
c
c     ib - ic interaction

          lchk = .true.
          do jz = 1,min(nexsit(ib),mxexcl)
            if(lexsit(ib,jz).eq.ic-isite) lchk=.false.
          enddo
          if(lchk) then
            nexsit(ib)=nexsit(ib)+1
            nexsit(ic)=nexsit(ic)+1
            if(max(nexsit(ib),nexsit(ic)).gt.mxexcl) then
              ibig = max(ibig,nexsit(ib),nexsit(ic))
              safe = .false.
            else
              lexsit(ib,nexsit(ib))=ic-isite
              lexsit(ic,nexsit(ic))=ib-isite
            endif
          endif
c
c     ia - ic interaction

          lchk = .true.
          do jz = 1,min(nexsit(ia),mxexcl)
            if(lexsit(ia,jz).eq.ic-isite) lchk=.false.
          enddo
          if(lchk) then
            nexsit(ia)=nexsit(ia)+1
            nexsit(ic)=nexsit(ic)+1
            if(max(nexsit(ia),nexsit(ic)).gt.mxexcl) then
              ibig = max(ibig,nexsit(ia),nexsit(ic))
              safe = .false.
            else
              lexsit(ia,nexsit(ia))=ic-isite
              lexsit(ic,nexsit(ic))=ia-isite
            endif
          endif
c
c     id - ib interaction

            lchk = .true.
            do jz = 1,min(nexsit(id),mxexcl)
              if(lexsit(id,jz).eq.ib-isite) lchk=.false.
            enddo
            if(lchk) then
              nexsit(id)=nexsit(id)+1
              nexsit(ib)=nexsit(ib)+1
              if(max(nexsit(id),nexsit(ib)).gt.mxexcl) then
                ibig = max(ibig,nexsit(id),nexsit(ib))
                safe = .false.
              else
                lexsit(id,nexsit(id))=ib-isite
                lexsit(ib,nexsit(ib))=id-isite
              endif
            endif
c
c     id - ic interaction

            lchk = .true.
            do jz = 1,min(nexsit(id),mxexcl)
              if(lexsit(id,jz).eq.ic-isite) lchk=.false.
            enddo
            if(lchk) then
              nexsit(id)=nexsit(id)+1
              nexsit(ic)=nexsit(ic)+1
              if(max(nexsit(id),nexsit(ic)).gt.mxexcl) then
                ibig = max(ibig,nexsit(id),nexsit(ic))
                safe = .false.
              else
                lexsit(id,nexsit(id))=ic-isite
                lexsit(ic,nexsit(ic))=id-isite
              endif
            endif
c
c     ia - id interaction

            lchk = .true.
            do jz = 1,min(nexsit(ia),mxexcl)
              if(lexsit(ia,jz).eq.id-isite) lchk=.false.
            enddo
            if(lchk) then
              nexsit(ia)=nexsit(ia)+1
              nexsit(id)=nexsit(id)+1
              if(max(nexsit(ia),nexsit(id)).gt.mxexcl) then
                ibig = max(ibig,nexsit(ia),nexsit(id))
                safe = .false.
              else
                lexsit(ia,nexsit(ia))=id-isite
                lexsit(id,nexsit(id))=ia-isite
              endif

          endif
          
        enddo
c     
c     exclude sites on basis of core-shell units
        
        do i=1,numshl(itmols)
          
          ishels=ishels+1
          
          ia=lstshl(ishels,1)+isite
          ib=lstshl(ishels,2)+isite
c
c     check interaction not already included
          
          lchk = .true.
          do jz = 1,min(nexsit(ia),mxexcl)
            if(lexsit(ia,jz).eq.ib-isite) lchk=.false.
          enddo
          if(lchk) then
            nexsit(ia)=nexsit(ia)+1
            nexsit(ib)=nexsit(ib)+1
            if(max(nexsit(ia),nexsit(ib)).gt.mxexcl) then
              ibig = max(ibig,nexsit(ia),nexsit(ib))
              safe = .false.
            else
              lexsit(ia,nexsit(ia))=ib-isite
              lexsit(ib,nexsit(ib))=ia-isite
            endif
          endif

c     
c     exclude sites on basis of bonds to core-shell units

          ibonds = ibonds - numbonds(itmols)
          do kk=1,numbonds(itmols)
            
            ibonds=ibonds+1
            
            if(keybnd(ibonds).gt.0) then
              
              ia1=lstbnd(ibonds,1)+isite
              ib1=lstbnd(ibonds,2)+isite

              if(ia.eq.ia1) then
c
c     check interaction not already included
          
                lchk = .true.
                do jz = 1,min(nexsit(ib1),mxexcl)
                  if(lexsit(ib1,jz).eq.ib-isite) lchk=.false.
                enddo
                if(lchk) then
                  nexsit(ib1)=nexsit(ib1)+1
                  nexsit(ib)=nexsit(ib)+1
                  if(max(nexsit(ib1),nexsit(ib)).gt.mxexcl) then
                    ibig = max(ibig,nexsit(ib1),nexsit(ib))
                    safe = .false.
                  else
                    lexsit(ib1,nexsit(ib1))=ib-isite
                    lexsit(ib,nexsit(ib))=ib1-isite
                  endif
                endif

              endif

              if(ia.eq.ib1) then
c
c     check interaction not already included
          
                lchk = .true.
                do jz = 1,min(nexsit(ia1),mxexcl)
                  if(lexsit(ia1,jz).eq.ib-isite) lchk=.false.
                enddo
                if(lchk) then
                  nexsit(ia1)=nexsit(ia1)+1
                  nexsit(ib)=nexsit(ib)+1
                  if(max(nexsit(ia1),nexsit(ib)).gt.mxexcl) then
                    ibig = max(ibig,nexsit(ia1),nexsit(ib))
                    safe = .false.
                  else
                    lexsit(ia1,nexsit(ia1))=ib-isite
                    lexsit(ib,nexsit(ib))=ia1-isite
                  endif
                endif

              endif

              if(ib.eq.ia1) then
c
c     check interaction not already included
          
                lchk = .true.
                do jz = 1,min(nexsit(ia),mxexcl)
                  if(lexsit(ia,jz).eq.ib1-isite) lchk=.false.
                enddo
                if(lchk) then
                  nexsit(ia)=nexsit(ia)+1
                  nexsit(ib1)=nexsit(ib1)+1
                  if(max(nexsit(ia),nexsit(ib1)).gt.mxexcl) then
                    ibig = max(ibig,nexsit(ia),nexsit(ib1))
                    safe = .false.
                  else
                    lexsit(ia,nexsit(ia))=ib1-isite
                    lexsit(ib1,nexsit(ib1))=ia-isite
                  endif
                endif

              endif
              if(ib.eq.ib1) then
c
c     check interaction not already included
          
                lchk = .true.
                do jz = 1,min(nexsit(ia),mxexcl)
                  if(lexsit(ia,jz).eq.ia1-isite) lchk=.false.
                enddo
                if(lchk) then
                  nexsit(ia)=nexsit(ia)+1
                  nexsit(ia1)=nexsit(ia1)+1
                  if(max(nexsit(ia),nexsit(ia1)).gt.mxexcl) then
                    ibig = max(ibig,nexsit(ia),nexsit(ia1))
                    safe = .false.
                  else
                    lexsit(ia,nexsit(ia))=ia1-isite
                    lexsit(ia1,nexsit(ia1))=ia-isite
                  endif
                endif

              endif

            endif

          enddo
c     
c     exclude sites on basis of constraint bonds to core-shell units

          iconst = iconst - numcon(itmols)
          do kk=1,numcon(itmols)
            
            iconst=iconst+1
            
            ia1=lstcon(iconst,1)+isite
            ib1=lstcon(iconst,2)+isite

            if(ia.eq.ia1) then
c
c     check interaction not already included
          
              lchk = .true.
              do jz = 1,min(nexsit(ib1),mxexcl)
                if(lexsit(ib1,jz).eq.ib-isite) lchk=.false.
              enddo
              if(lchk) then
                nexsit(ib1)=nexsit(ib1)+1
                nexsit(ib)=nexsit(ib)+1
                if(max(nexsit(ib1),nexsit(ib)).gt.mxexcl) then
                  ibig = max(ibig,nexsit(ib1),nexsit(ib))
                  safe = .false.
                else
                  lexsit(ib1,nexsit(ib1))=ib-isite
                  lexsit(ib,nexsit(ib))=ib1-isite
                endif
              endif

            endif

            if(ia.eq.ib1) then
c
c     check interaction not already included
          
              lchk = .true.
              do jz = 1,min(nexsit(ia1),mxexcl)
                if(lexsit(ia1,jz).eq.ib-isite) lchk=.false.
              enddo
              if(lchk) then
                nexsit(ia1)=nexsit(ia1)+1
                nexsit(ib)=nexsit(ib)+1
                if(max(nexsit(ia1),nexsit(ib)).gt.mxexcl) then
                  ibig = max(ibig,nexsit(ia1),nexsit(ib))
                  safe = .false.
                else
                  lexsit(ia1,nexsit(ia1))=ib-isite
                  lexsit(ib,nexsit(ib))=ia1-isite
                endif
              endif

            endif

            if(ib.eq.ia1) then
c
c     check interaction not already included
          
              lchk = .true.
              do jz = 1,min(nexsit(ia),mxexcl)
                if(lexsit(ia,jz).eq.ib1-isite) lchk=.false.
              enddo
              if(lchk) then
                nexsit(ia)=nexsit(ia)+1
                nexsit(ib1)=nexsit(ib1)+1
                if(max(nexsit(ia),nexsit(ib1)).gt.mxexcl) then
                  ibig = max(ibig,nexsit(ia),nexsit(ib1))
                  safe = .false.
                else
                  lexsit(ia,nexsit(ia))=ib1-isite
                  lexsit(ib1,nexsit(ib1))=ia-isite
                endif
              endif

            endif
            if(ib.eq.ib1) then
c
c     check interaction not already included
          
              lchk = .true.
              do jz = 1,min(nexsit(ia),mxexcl)
                if(lexsit(ia,jz).eq.ia1-isite) lchk=.false.
              enddo
              if(lchk) then
                nexsit(ia)=nexsit(ia)+1
                nexsit(ia1)=nexsit(ia1)+1
                if(max(nexsit(ia),nexsit(ia1)).gt.mxexcl) then
                  ibig = max(ibig,nexsit(ia),nexsit(ia1))
                  safe = .false.
                else
                  lexsit(ia,nexsit(ia))=ia1-isite
                  lexsit(ia1,nexsit(ia1))=ia-isite
                endif
              endif

            endif

          enddo
c     
c     exclude sites on basis of rigid units involving  core or shell

          igrp = igrp - numgrp(itmols)
          do kk = 1,numgrp(itmols)
          
            igrp = igrp + 1
          
            id = listyp(igrp)
          
            do jj = 1,numgsit(id)
            
              ia1 = lstgst(igrp,jj)+isite
              if(ia1.eq.ia) then

                do jk = 1,numgsit(id)
            
                  if(jk.ne.jj) then
                    ib1 = lstgst(igrp,jk)+isite
c     
c     check interaction not already included
          
                    lchk = .true.
                    do jz = 1,min(nexsit(ib1),mxexcl)
                      if(lexsit(ib1,jz).eq.ib-isite) lchk=.false.
                    enddo
                    if(lchk) then
                      nexsit(ib1)=nexsit(ib1)+1
                      nexsit(ib)=nexsit(ib)+1
                      if(max(nexsit(ib1),nexsit(ib)).gt.mxexcl) then
                        ibig = max(ibig,nexsit(ib1),nexsit(ib))
                        safe = .false.
                      else
                        lexsit(ib1,nexsit(ib1))=ib-isite
                        lexsit(ib,nexsit(ib))=ib1-isite
                      endif
                    endif


                  endif

                enddo
                
              endif

              if(ia1.eq.ib) then

                do jk = 1,numgsit(id)
            
                  if(jk.ne.jj) then
                    ib1 = lstgst(igrp,jk)+isite

c
c     check interaction not already included
          
                    lchk = .true.
                    do jz = 1,min(nexsit(ia),mxexcl)
                      if(lexsit(ia,jz).eq.ib1-isite) lchk=.false.
                    enddo
                    if(lchk) then
                      nexsit(ia)=nexsit(ia)+1
                      nexsit(ib1)=nexsit(ib1)+1
                      if(max(nexsit(ia),nexsit(ib1)).gt.mxexcl) then
                        ibig = max(ibig,nexsit(ia),nexsit(ib1))
                        safe = .false.
                      else
                        lexsit(ia,nexsit(ia))=ib1-isite
                        lexsit(ib1,nexsit(ib1))=ia-isite
                      endif

                    endif

                  endif

                enddo
                
              endif

            enddo

          enddo

        enddo

        isite=isite+numsit(itmols)
        
      enddo
      
      ntpsit=isite
c
c     check for exceeded array bounds

      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) then
        call gimax(ibig,1,jj)
        if(idnode.eq.0) write(nrite,*) 'mxexcl must be at least ',ibig
        if(idnode.eq.0) write(nrite,*) 'mxexcl is currently ',mxexcl
        call error(idnode,65)
      endif
c     
c     remove redundant entries from exclusion list
c     (there shouldn't be any!)
      
      do i=1,ntpsit
        
        nlast=nexsit(i)
        do j=1,nexsit(i)-1
          
          if(j.lt.nlast)then
            
            kk=j
            do k=j+1,nexsit(i)
              
              if(lexsit(i,j).eq.lexsit(i,k))then
                
                nlast=nlast-1
                lexsit(i,k)=0
                
              else if(lexsit(i,k).gt.0)then
                
                kk=kk+1
                lexsav=lexsit(i,k)
                lexsit(i,k)=0
                lexsit(i,kk)=lexsav
                
              endif
              
            enddo
            
          endif
          
        enddo
        
        nexsit(i)=nlast
        
      enddo

#ifdef VAMPIR
      call VTEND(77, ierr)
#endif
      return
      end

