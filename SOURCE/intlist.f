      subroutine intlist
     x  (lshmov,idnode,mxnode,natms,nscons,ntangl,ntbond,ntcons,
     x  ntdihd,ntinv,ntpmls,ntteth,ntshl,ntpmf,nspmf,lashap,lishap,
     x  listang,listbnd,listcon,listdih,listinv,listshl,listin,
     x  listme,listot,lstang,lstbnd,lstcon,lstdih,lstinv,lstfrz,
     x  lsttet,listtet,lstshl,numang,numbonds,numcon,numdih,numinv,
     x  nummols,numsit,numteth,numshl,numpmf,npmf,indpmf,lstpmf,
     x  listpm,lstpmt,itest,index,kscons,msite,mconst)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for constructing the interaction lists
c     for the entire simulated system
c     
c     parallel replicated dat version : block data
      
c     copyright - daresbury laboratory 1992
c     author    - w. smith        july 1992
c     amended   - t.forester      oct 1993
c     amended   - t.forester      dec 1994 : block data
c     
c     wl
c     2001/05/30 12:40:08
c     1.7
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical safe,lshmov,safe1
      logical lchk,lfail
      
      dimension itest(mxtmls),index(mxtmls),kscons(0:mxproc-1)
      dimension msite(mxtmls),mconst(mxtmls),numinv(mxtmls)
      dimension nummols(mxtmls),numsit(mxtmls),numbonds(mxtmls)
      dimension numang(mxtmls),numdih(mxtmls),numcon(mxtmls)
      dimension lstbnd(mxtbnd,2),lstcon(mxtcon,2),lstang(mxtang,3)
      dimension lstdih(mxtdih,4),listbnd(mxbond,3),listdih(mxdihd,5)
      dimension listcon(mxcons,3),listang(mxangl,4),lstfrz(mxatms)
      dimension listme(mxatms),listin(mxatms),lashap(mxproc)
      dimension lishap(mxlshp),listot(mxatms),numshl(mxtmls)
      dimension listshl(mxshl,3),lstshl(mxtshl,2),lstinv(mxtinv,4)
      dimension numteth(mxtmls),lsttet(mxteth),listtet(msteth,2)
      dimension npmf(2),numpmf(mxtmls),indpmf(mxspmf)
      dimension lstpmf(mxspmf,mspmf),listinv(mxinv,5)
      dimension listpm(mxpmf),lstpmt(mxpmf)
#ifdef VAMPIR
      call VTBEGIN(145, ierr)
#endif
c     
c     initialise bookkeeping indices
      
      ibonds=0
      jbonds=0
      kbonds=0
      ipmf =0
      jpmf=0
      iangle=0
      jangle=0
      kangle=0
      idihed=0
      jdihed=0
      kdihed=0
      iinver=0
      jinver=0
      kinver=0
      iteths=0
      jteths=0
      kteths=0
      ishels=0
      jshels=0
      kshels=0
      safe=.true.
      safe1=.true.
c     
c     find total number of bonds,pmf constraints,bond constraints,
c     angles,dihedrals,inversions, tethers,core-shells, in system 
c     - ignoring frozen atoms
      
      ntbon0 = 0
      ntpmf0 = 0
      ntcon0 = 0
      ntang0 = 0
      ntdih0 = 0
      ntinv0 = 0
      nttet0 = 0
      ntshl0 = 0
      nscons = 0
      ntcons = 0
      
      do itmols = 1,ntpmls
        
        ntbon0 = ntbon0+nummols(itmols)*numbonds(itmols)
        ntpmf0 = ntpmf0+nummols(itmols)*numpmf(itmols)
        ntcon0 = ntcon0+nummols(itmols)*numcon(itmols)
        ntang0=  ntang0+nummols(itmols)*numang(itmols)
        ntdih0 = ntdih0+nummols(itmols)*numdih(itmols)
        ntinv0 = ntinv0+nummols(itmols)*numinv(itmols)
        nttet0 = nttet0+nummols(itmols)*numteth(itmols)
        ntshl0 = ntshl0+nummols(itmols)*numshl(itmols)
        
      enddo
      
      isite=0
      iconst=0
      jconst=0
      kconst=0
c     
c     first and last index of bonds, angles etc for this node
      
      ibnd1 = (idnode*ntbon0)/mxnode + 1
      ibnd2 = ((idnode+1)*ntbon0)/mxnode

      ipmf1 = (idnode*ntpmf0)/mxnode + 1
      ipmf2 = ((idnode+1)*ntpmf0)/mxnode
      ntpmf = ntpmf0
      nspmf = ipmf2+1-ipmf1

      iang1 = (idnode*ntang0)/mxnode + 1
      iang2 = ((idnode+1)*ntang0)/mxnode
      
      idih1 = (idnode*ntdih0)/mxnode + 1
      idih2 = ((idnode+1)*ntdih0)/mxnode
      
      iinv1 = (idnode*ntinv0)/mxnode + 1
      iinv2 = ((idnode+1)*ntinv0)/mxnode
      
      itet1 = (idnode*nttet0)/mxnode + 1
      itet2 = ((idnode+1)*nttet0)/mxnode
      
      ishl1 = (idnode*ntshl0)/mxnode + 1
      ishl2 = ((idnode+1)*ntshl0)/mxnode
      
c     
c     loop over molecule types
      
      do itmols=1,ntpmls
        
c     
c     loop over molecules in system
        
        do imols=1,nummols(itmols)
          
c     
c     construct bond constraint list later
          
c     
c     construct chemical bond interaction list
          
          do lbonds=1,numbonds(itmols)
            
            ibonds=ibonds+1
            
            if(ibonds.ge.ibnd1.and.ibonds.le.ibnd2)then
              
              jbonds=jbonds+1
              if(jbonds.le.mxbond)then
                
                listbnd(jbonds,1)=lbonds+kbonds
                listbnd(jbonds,2)=lstbnd(lbonds+kbonds,1)
     x            +isite
                listbnd(jbonds,3)=lstbnd(lbonds+kbonds,2)
     x            +isite
                
              else
                
                safe=.false.
                
              endif
              
            endif
            
          enddo
          
          if (mxnode.gt.1) call gstate(safe)
          if(.not.safe) call error(idnode,31)
c     
c     construct pmf site lists - no exclusions
          
          do lpmf=1,numpmf(itmols)
            
            ipmf=ipmf+1
            
            if(ipmf.ge.ipmf1.and.ipmf.le.ipmf2)then
              
              jpmf=jpmf+1
              if(jpmf.le.mspmf)then
                
                nnn = npmf(1)+npmf(2)
                if(nnn.le.mxspmf) then

                  do jj = 1,npmf(1)+npmf(2)
                    lstpmf(jj,jpmf)=indpmf(jj)+isite
                  enddo

                else

                  safe=.false.
                  
                endif

              else
                
                safe1=.false.
                
              endif
              
            endif
            
          enddo
          
          if (mxnode.gt.1) call gstate(safe1)
          if(.not.safe1) call error(idnode,458)

          if (mxnode.gt.1) call gstate(safe)
          if(.not.safe) call error(idnode,460)
c     
c     construct valence angle interaction list
          
          do langle=1,numang(itmols)
            
            iangle=iangle+1
            
            if(iangle.ge.iang1.and.iangle.le.iang2)then
              
              jangle=jangle+1
              if(jangle.le.mxangl)then
                
                listang(jangle,1)=langle+kangle
                listang(jangle,2)=lstang(langle+kangle,1)
     x            +isite
                listang(jangle,3)=lstang(langle+kangle,2)
     x            +isite
                listang(jangle,4)=lstang(langle+kangle,3)
     x            +isite
                
              else
                
                safe=.false.
                
              endif
              
            endif
            
          enddo
          
          if (mxnode.gt.1) call gstate(safe)
          if(.not.safe) call error(idnode,51)
          
c     
c     construct dihedral angle interaction list
          
          do ldihed=1,numdih(itmols)
            
            idihed=idihed+1
            
            if(idihed.ge.idih1.and.idihed.le.idih2)then
              
              jdihed=jdihed+1
              if(jdihed.le.mxdihd)then
                
                listdih(jdihed,1)=ldihed+kdihed
                listdih(jdihed,2)=lstdih(ldihed+kdihed,1)
     x            +isite
                listdih(jdihed,3)=lstdih(ldihed+kdihed,2)
     x            +isite
                listdih(jdihed,4)=lstdih(ldihed+kdihed,3)
     x            +isite
                listdih(jdihed,5)=lstdih(ldihed+kdihed,4)
     x            +isite
                
              else
                
                safe=.false.
                
              endif
              
            endif
            
          enddo
          
          if (mxnode.gt.1) call gstate(safe)
          if(.not.safe) call error(idnode,61)

c     
c     construct inversion potential list
          
          do linver=1,numinv(itmols)
            
            iinver=iinver+1
            
            if(iinver.ge.iinv1.and.iinver.le.iinv2)then
              
              jinver=jinver+1
              if(jinver.le.mxinv)then
                
                listinv(jinver,1)=linver+kinver
                listinv(jinver,2)=lstinv(linver+kinver,1)
     x            +isite
                listinv(jinver,3)=lstinv(linver+kinver,2)
     x            +isite
                listinv(jinver,4)=lstinv(linver+kinver,3)
     x            +isite
                listinv(jinver,5)=lstinv(linver+kinver,4)
     x            +isite
                
              else
                
                safe=.false.
                
              endif
              
            endif
            
          enddo
          
          if (mxnode.gt.1) call gstate(safe)
          if(.not.safe) call error(idnode,77)
c     
c     construct tethered atoms interaction list
          
          do lteths=1,numteth(itmols)
            
            iteths=iteths+1
            
            if(iteths.ge.itet1.and.iteths.le.itet2)then
              
              jteths=jteths+1
              if(jteths.le.msteth)then
                
                listtet(jteths,1)=lteths+kteths
                listtet(jteths,2)=lsttet(lteths+kteths)+isite
                
              else
                
                safe=.false.
                
              endif
              
            endif
            
          enddo
          
          if (mxnode.gt.1) call gstate(safe)
          if(.not.safe) call error(idnode,63)
c     
c     construct core-shell list
          
          do lshels=1,numshl(itmols)
            
            ishels=ishels+1
            
            if(ishels.ge.ishl1.and.ishels.le.ishl2)then
              
              jshels=jshels+1
              if(jshels.le.mxshl)then
                
                listshl(jshels,1)=lshels+kshels
                listshl(jshels,2)=lstshl(lshels+kshels,1)
     x            +isite
                listshl(jshels,3)=lstshl(lshels+kshels,2)
     x            +isite
                
              else
                
                safe=.false.
                
              endif
              
            endif
            
          enddo
          
          if (mxnode.gt.1) call gstate(safe)
          if(.not.safe) call error(idnode,59)
          
          isite=isite+numsit(itmols)
          
        enddo
        
        kbonds=kbonds+numbonds(itmols)
        kangle=kangle+numang(itmols)
        kdihed=kdihed+numdih(itmols)
        kinver=kinver+numinv(itmols)
        kteths=kteths+numteth(itmols)
        kshels=kshels+numshl(itmols)
        
      enddo
c     
c     store array counters for bookkeeping
      
      ntbond=ibonds
      ntangl=iangle
      ntdihd=idihed
      ntinv =iinver
      ntteth=iteths
      ntshl =ishels
c     
c     pass bond constraint information to other nodes
      
      if(ntcon0.gt.0)then
        
        ntcons = ntcon0
c     
c     find starting site no. and constraint no. for each molec. type
        
        msite(1) = 0
        mconst(1) = 0
        
        do itmols = 2,ntpmls
          msite(itmols) = msite(itmols-1) + numsit(itmols-1)*
     x      nummols(itmols-1)
          mconst(itmols) = mconst(itmols-1) + numcon(itmols-1)
        enddo
        
c     
c     sort molecules into ascending order of number of constraints
        
        do i = 1,ntpmls
          itest(i) = numcon(i)
          index(i) = 0
        enddo
        
        call shellsort(ntpmls,itest)
        
        do i = 1,ntpmls
          
          lchk =.true.
          do j = 1,ntpmls
            
            if(itest(i).eq.numcon(j)) then
              
              if(lchk) then 
                index(i) = j
                lchk=.false.
                
              endif
              
              do ii = 1,i-1
                if(index(ii).eq.j) lchk = .true.
              enddo
              
            endif
            
          enddo
          
        enddo
        
c     
c     load balance to within 10%
        
        tol = 1.0d0 + (0.10d0)/2.d0
        kcons = (ntcons)/mxnode
        ntmp =0
c     
c     find smallest constrained molecule to allocate to a node
        
        do i = 1,ntpmls
          
          if(ntmp.le.mxnode) then
            
            if(numcon(index(i)).gt.0) then
              ntmp = ntmp+nummols(index(i))
              klo = max(0,kcons - numcon(index(i))/2)
              khi = klo +  numcon(index(i)) + 1
            endif
            
          endif
          
        enddo
c     
c     reset hi/lo limits if molecules contain too many constraints
        
        if(dble(khi)/dble(max(1,klo)).gt.tol) then
          klo = nint(dble(kcons)/tol)
          khi = nint(dble(kcons)*tol)+1
        endif
c     
c     store lo value for later
        
        klo0 = klo
c     
c     begin assignment of constraints ----------------------------------
        
        ifail = -1
   50   ifail = ifail+1
        
        if(ifail.gt.ntpmls) then
          call error(idnode,432)
        endif
        
        iconst=0
        jconst=0
        kconst=0
        lconst=0
c     
c     zero running totals of constraints on each processor
 
        do id = 0,mxnode-1
          kscons(id) = 0
        enddo
        
        iloop = 0
        lfail = .false.
        iconst = 0
        jconst = 0
        nnode = 0
c     
c     assign difficult molecules in blocks
        
        if(ifail.gt.0) then
          
          nfail =0
          do i=1,ifail
            ii= ntpmls+1-i
            nfail = nfail+ nummols(index(ii))*numcon(index(ii))
          enddo
c     
c     decide on number of processors to split over
          
          nnode = int(dble(nfail)/dble(max(kcons,1))+ 1.d0/tol)
          nnode = max(2,nnode)
          nnode = min(nnode,mxnode)
c     
c     assign to processors 0..nnode-1
          
          do id=0,nnode-1
            
            nscons0 = (id*nfail)/nnode+1
            nscons1 = ((id+1)*nfail)/nnode
            
            kscons(id) = nscons1+1-nscons0
            
          enddo
c     
c     this processors block
          
          nscons0 = (idnode*nfail)/nnode+1
          nscons1 = ((idnode+1)*nfail)/nnode
          
c     
c     assign in blocks
          
          do itmols = ntpmls,ntpmls-ifail+1,-1
            
            ii = index(itmols)
            icon = numcon(ii)
            kconst = mconst(ii)
            
            do imols = 1,nummols(ii)
              
              isite = msite(ii)+ (imols-1)*numsit(ii)
c     
c     construct bond constraint list
              
              do lconst=1,numcon(ii)
                
                iconst=iconst+1
                
                if(iconst.ge.nscons0.and.iconst.le.nscons1) then
                  
                  jconst=jconst+1
                  
                  if(jconst.le.mxcons)then
                    
                    listcon(jconst,1)=lconst+kconst
                    iatom = lstcon(lconst+kconst,1)+isite
                    jatom = lstcon(lconst+kconst,2)+isite
                    
                    listcon(jconst,2)=iatom
                    listcon(jconst,3)=jatom
                    
                  else
                    
                    safe=.false.
                    
                  endif
                  
                endif
                
              enddo
              
            enddo
            
          enddo
          
        endif
c     
c     assign non-problematic molecules
        
        jdnode = mod(nnode+1,mxnode)
        
        do itmols = ntpmls-ifail,1,-1
          
          ii = index(itmols)
          icon = numcon(ii)
          kconst = mconst(ii)
          
          do imols = 1,nummols(ii)
            
            itry = 0
  100       if(kscons(jdnode)+icon.le.klo) then
              
              if(jdnode.ne.idnode) then
                kscons(jdnode) = kscons(jdnode)+icon
                jdnode = mod(jdnode+1,mxnode)
                lchk = .false.
              else
c     
c     construct bond constraint list
                
                isite = msite(ii)+ (imols-1)*numsit(ii)
                do lconst=1,numcon(ii)
                  
                  jconst=jconst+1
                  
                  if(jconst.le.mxcons)then
                    
                    listcon(jconst,1)=lconst+kconst
                    iatom = lstcon(lconst+kconst,1)+isite
                    jatom = lstcon(lconst+kconst,2)+isite
                    listcon(jconst,2)=iatom
                    listcon(jconst,3)=jatom
                    
                  else
                    
                    safe=.false.
                    
                  endif
                  
                enddo
                
                kscons(jdnode) = kscons(jdnode)+icon
                jdnode = mod(jdnode+1,mxnode)
                lchk = .false.
                
              endif
              
            else
              
              jdnode = mod(jdnode+1,mxnode)
              lchk = .true.
              itry = itry+1
              
            endif
            
            if(lchk.and.itry.gt.mxnode) then
              
              klo = kcons
              kcons = khi
              itry = 0
              iloop = iloop + 1
              
            endif
c     
c     split molecule across nodes if have to
            
            if(iloop.gt.3) then
              lfail=.true.
              kcons = ntcons/mxnode
              klo = klo0
              lchk = .false.
            endif
            
            if(lchk) goto 100
            
          enddo
          
        enddo
c     
c     check no node has less than minimum number
        
        do id = 0,mxnode-1
          if(kscons(id).lt.klo0) then 
            lfail = .true.
          endif
        enddo
        
        if(lfail) goto 50
        
        if (mxnode.gt.1) call gstate(safe)
        if(.not.safe) then
 
          if(mxnode.gt.1) call gimax(jconst,1,idum)
          if(idnode.eq.0)write(nrite,'(a,i10,a,i10)')
     x      'Number of constraints found ',jconst,'Max allowed ',mxcons

          call error(idnode,41)

        endif
        
        nscons = kscons(idnode)
        
        call passcon
     x    (lshmov,idnode,mxnode,natms,nscons,lashap,lishap,listme,
     x    listin,listot,listcon,lstfrz)
        
      endif

      if(npmf(1).gt.0) then

        call passpmf
     x    (idnode,mxnode,natms,nspmf,listpm,listin,lstpmt,lstpmf,npmf)

      endif

#ifdef VAMPIR
      call VTEND(145, ierr)
#endif
      return
      end
