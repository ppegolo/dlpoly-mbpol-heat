      subroutine quatbook
     x  (idnode,imcon,mxnode,natms,ngrp,ntpmls,ntfree,
     x  degfre,degrot,numgsit,numgrp,nummols,numsit,
     x  lstme,lstfrz,listyp,lstfre,lstgot,lstgtp,ind,
     x  lstgst,lstrgd,lstbod,xxx,yyy,zzz,weight,q0,q1,
     x  q2,q3,cell,buffer,gxx,gyy,gzz,gcmx,gcmy,gcmz,
     x  xxt,yyt,zzt,txx,tyy,tzz,gmass,rotinx,rotiny,
     x  rotinz,accum,gaxs,rotmin)
      
c**************************************************************************
c     
c     DLPOLY subroutine for setting up bookkeeping for rigid bodies
c     
c     parallel replicated data version : block data
c     
c     copyright daresbury laboratory 1993
c     author      t.forester october 1993
c     amended     t.forester dec 1994 : block data
c     
c     wl
c     2001/08/31 11:13:51
c     1.6
c     $Sate: Exp $
c     
c*************************************************************************
      
#include "dl_params.inc"
      
      logical safe,linear
      
      dimension buffer(mxbuff),accum(mxungp)
      dimension gcmx(mxgrp),gcmy(mxgrp),gcmz(mxgrp)
      dimension gmass(mxungp)
      dimension q0(mxgrp),q1(mxgrp),q2(mxgrp),q3(mxgrp)
      dimension gxx(mxungp,mxngp),gyy(mxungp,mxngp),gzz(mxungp,mxngp)
      dimension rotinx(mxungp,2),rotiny(mxungp,2),rotinz(mxungp,2)
      dimension lstfre(mxatms),lstrgd(mxgatm),numgsit(mxungp)
      dimension lstgot(mxatms),lstgtp(mxgrp),listyp(mxungp)
      dimension lstfrz(mxatms),lstbod(mxatms)
      dimension numgrp(mxtmls),nummols(mxtmls),numsit(mxtmls)
      dimension lstgst(mxungp,mxngp),lstme(mxatms)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms),cell(9)
      dimension rot(9),aa(9),rotinr(3,3),bb(9),weight(mxatms)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension ind(mxgrp,3),gaxs(mxungp,9)
      dimension rotmin(mxungp),rot1(3,3)
      
#ifdef VAMPIR
      call VTBEGIN(157, ierr)
#endif
c     
c     initialise bookkeeping indices
      
      igrp = 0
      jgrp = 0
      kgrp = 0
      isite=0
      jr =0
      jt =0
      safe=.true.
      degfre = 0.d0
      degrot = 0.d0
c     
c     rigid body identifier
      
      do i = 1,natms
        
        lstbod(i) = 0
        
      enddo
c     
c     number of rigid groups in system
      
      ngrp = 0
      do itmols=1,ntpmls
        ngrp=ngrp+nummols(itmols)*numgrp(itmols)
      enddo
c     
c     block indices for groups
      
      igrp1 = (idnode*ngrp)/mxnode + 1
      igrp2 = ((idnode+1)*ngrp)/mxnode
      
c     
c     loop over molecule types
      
      do itmols=1,ntpmls
c     
c     loop over molecules in system
        
        do imols=1,nummols(itmols)
c     
c     construct rigid body site list: each processor has a different copy
          
          do lgrp=1,numgrp(itmols)
            
            igrp=igrp+1
            
            if(igrp.le.mxgrp)then
              
              lstgtp(igrp) = listyp(lgrp+kgrp)
              id = listyp(lgrp+kgrp)
              
              if((igrp.ge.igrp1).and.(igrp.le.igrp2)) then

                jgrp = jgrp+1
                
                do jj = 1,numgsit(id)
                  
                  jr = jr +1
                  jt =jt+1
                  
                  if(jr.le.mxatms.and.jt.le.mxatms) then
                    
                    lstrgd(jr)=lstgst(id,jj)+isite
                    lstgot(jt) =lstgst(id,jj)+isite
                    lstbod(lstgst(id,jj)+isite)=igrp
                    
                  else
                    
                    safe = .false.
                    
                  endif
                  
                enddo
                
              else
                
                do jj = 1,numgsit(id)
                  
                  jt =jt+1
                  if(jt.le.mxatms) then
                    
                    lstgot(jt) =lstgst(id,jj)+isite
                    lstbod(lstgst(id,jj)+isite)=igrp
                    
                  else
                    
                    safe =.false.
                    
                  endif
                  
                enddo
                
              endif
              
            else
              
              safe =.false.
              
            endif
            
          enddo
          
          if (mxnode.gt.1) call gstate(safe)
          if(.not.safe) call error(idnode,304)
          isite=isite+numsit(itmols)
          
        enddo
        
        kgrp=kgrp+numgrp(itmols)
        
      enddo
      
      if(ngrp.eq.0) then
        
        j=0
        do i = 1,natms
          if(lstfrz(i).eq.0)then
            j=j+1
            lstfre(j) = i
          endif
        enddo
        ntfree = j
        
      else
c     
c     centre of mass of groups
c     assumes group dimensions are smaller than half box width
        
        do id = 1,mxungp
          
          gmass(id) =0.d0
          
        enddo
        
        jr = 0
        do ig = igrp1,igrp2
          
c     
c     working com is first site in group
          
          i = lstrgd(jr+1)
          txx(ig) = xxx(i)
          tyy(ig) = yyy(i)
          tzz(ig) = zzz(i)
          
          id = lstgtp(ig)
          safe = .false.
          if(gmass(id).eq.0.d0) safe=.true.
          
          do j = 1,numgsit(id)
            
            jr =jr+1
            i=lstrgd(jr)
            xxt(jr) = xxx(i) - txx(ig)
            yyt(jr) = yyy(i) - tyy(ig)
            zzt(jr) = zzz(i) - tzz(ig)
            if(safe) gmass(id) = gmass(id)+weight(i)
            
          enddo
          
        enddo
c     
c     minimum image from working com
        
        call images(imcon,0,1,jr,cell,xxt,yyt,zzt)
        
        jr = 0
        do ig = igrp1,igrp2
          
          gcmx(ig) =0.d0
          gcmy(ig) =0.d0
          gcmz(ig) =0.d0
          
          id = lstgtp(ig)
          
          do j = 1,numgsit(id)
            
            jr = jr+1
            i = lstrgd(jr)
            
            gcmx(ig)=gcmx(ig) + weight(i)*xxt(jr)
            gcmy(ig)=gcmy(ig) + weight(i)*yyt(jr)
            gcmz(ig)=gcmz(ig) + weight(i)*zzt(jr)
            
          enddo
          
          gcmx(ig) = gcmx(ig)/gmass(id) +txx(ig)
          gcmy(ig) = gcmy(ig)/gmass(id) +tyy(ig)
          gcmz(ig) = gcmz(ig)/gmass(id) +tzz(ig)
          
        enddo
        
c     
c     global communications
        
        if(mxnode.gt.1) then
          
          nbuff = mxbuff
          call merge(idnode,mxnode,ngrp,nbuff,gcmx,gcmy,gcmz,buffer)
          
        endif
c     
c     make sure all nodes have same copy of gmass
        
        if(mxnode.gt.1) then
          
          do id = 1,mxungp
            
            accum(id) = 0.d0
            if(gmass(id).gt.0.d0) accum(id) = 1.d0
            
          enddo
          
          call gdsum(gmass(1),mxungp,buffer(1))
          call gdsum(accum(1),mxungp,buffer(1))
          
          do id = 1,mxungp
            
            dnorm = max(1.d0,accum(id))
            gmass(id) = gmass(id)/dnorm
            
          enddo
          
        endif
c     
c     find a group of each type on this node to 
c     find principal axis system of the group type
        
        do id = 1,mxungp
          
          safe = .false.
          jr =0
          ij=0
          
  100     ij=ij+1
          
          if(ij.gt.ngrp) goto 110
          
          jr = jr+numgsit(lstgtp(ij))
          
          if (lstgtp(ij).eq.id) safe = .true.
          
          if(.not.safe) goto 100
          
  110     if(safe) then
c     
c     rotational inertia accumulator
            
            do k = 1,3
              do kk = 1,3
                rotinr(k,kk)=0.d0
              enddo
            enddo
            
            jr = jr - numgsit(id)
            do j = 1,numgsit(id)
              
              jr = jr +1
              i = lstgot(jr)
              
              xxt(jr) = xxx(i) -gcmx(ij)
              yyt(jr) = yyy(i) -gcmy(ij)
              zzt(jr) = zzz(i) -gcmz(ij)
              
              call images(imcon,0,1,1,cell,xxt(jr),yyt(jr),zzt(jr))
              
              rotinr(1,1)= rotinr(1,1)+weight(i)*(xxt(jr)**2)
              rotinr(1,2)= rotinr(1,2)+weight(i)*xxt(jr)*yyt(jr)
              rotinr(1,3)= rotinr(1,3)+weight(i)*xxt(jr)*zzt(jr)
              rotinr(2,2)= rotinr(2,2)+weight(i)*(yyt(jr)**2)
              rotinr(2,3)= rotinr(2,3)+weight(i)*yyt(jr)*zzt(jr)
              rotinr(3,3)= rotinr(3,3)+weight(i)*(zzt(jr)**2)
              
            enddo
            
            rotinr(2,1) = rotinr(1,2)
            rotinr(3,1) = rotinr(1,3)
            rotinr(3,2) = rotinr(2,3)
            
            call jacobi(rotinr,rot1,3)
            
            rot(1) = rot1(1,1)
            rot(4) = rot1(2,1)
            rot(7) = rot1(3,1)
            rot(2) = rot1(1,2)
            rot(5) = rot1(2,2)
            rot(8) = rot1(3,2)
            rot(3) = rot1(1,3)
            rot(6) = rot1(2,3)
            rot(9) = rot1(3,3)
            
c     
c     rotational inertia accumulators
            
            rotinx(id,1) = 0.d0
            rotiny(id,1) = 0.d0
            rotinz(id,1) = 0.d0
            
            jr =jr - numgsit(id)
            do j = 1,numgsit(id)
              
              jr = jr+1
              i = lstgot(jr)
              
c     
c     site positions in principal axis system
              
              gxx(id,j) = rot(1)*xxt(jr)+rot(4)*yyt(jr)+rot(7)*zzt(jr)
              gyy(id,j) = rot(2)*xxt(jr)+rot(5)*yyt(jr)+rot(8)*zzt(jr)
              gzz(id,j) = rot(3)*xxt(jr)+rot(6)*yyt(jr)+rot(9)*zzt(jr)
c
c              impose rounding 
              
              if(abs(gxx(id,j)).lt.1.d-8) gxx(id,j) = 0.d0
              if(abs(gyy(id,j)).lt.1.d-8) gyy(id,j) = 0.d0
              if(abs(gzz(id,j)).lt.1.d-8) gzz(id,j) = 0.d0
              
c     
c     rotational inertia tensor of group type
              
              rotinx(id,1)=rotinx(id,1)+
     x          weight(i)*(gyy(id,j)**2+gzz(id,j)**2)
              rotiny(id,1)=rotiny(id,1)+
     x          weight(i)*(gzz(id,j)**2+gxx(id,j)**2)
              rotinz(id,1)=rotinz(id,1)+
     x          weight(i)*(gxx(id,j)**2+gyy(id,j)**2)
              
            enddo
c     
c     set axis system such that: Ixx >= Iyy >= Izz
            
            rotxyz = max(rotinx(id,1),rotiny(id,1),rotinz(id,1))
            
            if(rotxyz.ge.rotinx(id,1)) then
              
              if(rotiny(id,1).ge.rotxyz) then
                
                do j = 1,numgsit(id)
                  
                  a1 = gxx(id,j)
                  gxx(id,j) = gyy(id,j)
                  gyy(id,j) = -a1
                  
                enddo
                
                rotiny(id,1) = rotinx(id,1)
                rotinx(id,1) = rotxyz
                
              elseif(rotinz(id,1).ge.rotxyz) then
                
                do j = 1,numgsit(id)
                  
                  a1 = gxx(id,j)
                  gxx(id,j) = gzz(id,j)
                  gzz(id,j) = -a1
                  
                enddo
                
                rotinz(id,1) = rotinx(id,1)
                rotinx(id,1) = rotxyz
                
              endif
              
            endif
            
            if(rotinz(id,1).gt.rotiny(id,1)) then
              
              do j = 1,numgsit(id)
                
                a1 = gyy(id,j)
                gyy(id,j) = gzz(id,j)
                gzz(id,j) = -a1
                
              enddo
              
              a1 = rotinz(id,1)
              rotinz(id,1) = rotiny(id,1)
              rotiny(id,1) = a1
              
            endif
c     
c     set up principal axis system in terms of site positions
            
c     
c     test for (near) linear unit
            
            ill = 0
            rtall = (rotinx(id,1)+rotiny(id,1)+rotinz(id,1))
            
            if(rtall.gt.1.d-5) then
              rotall = rtall
            else
              rotall = 1.d0
            endif
            
            rotmin(id) = min(rotinx(id,1),rotiny(id,1))
            rotmin(id) = min(rotmin(id),rotinz(id,1))/rotall
            
            if((rotinx(id,1)/rotall).lt.1.d-5) ill = ill+1
            if((rotiny(id,1)/rotall).lt.1.d-5) ill = ill+1
            if((rotinz(id,1)/rotall).lt.1.d-5) ill = ill+1
            
            if(ill.ge.2) then
c     
c     point particle only
              
              ind(id,1) = 1
              ind(id,2) = 1
              ind(id,3) = 1
              
              do jj = 1,9
                gaxs(id,jj) = 0.d0
              enddo
              
            elseif(ill.eq.1) then
c     
c     linear molecule
              
              ind(id,1) = 1
              ind(id,2) = 2
              ind(id,3) = 1
              
              aa(1) = gxx(id,1) - gxx(id,2)
              aa(4) = gyy(id,1) - gyy(id,2)
              aa(7) = gzz(id,1) - gzz(id,2)
              rsq = sqrt(aa(1)**2 + aa(4)**2+aa(7)**2)
              
              if(abs(aa(7)/rsq).gt.0.5d0) then
                
                rsq = sqrt(aa(4)**2+aa(7)**2)
                aa(2) =  0.d0
                aa(5) =  aa(7)/rsq
                aa(8) = -aa(4)/rsq
                
              elseif(abs(aa(4)/rsq).gt.0.5d0) then
                
                rsq = sqrt(aa(4)**2+aa(1)**2)
                aa(2) = -aa(4)/rsq
                aa(5) =  aa(1)/rsq
                aa(8) =  0.d0
                
              elseif(abs(aa(1)/rsq).gt.0.5d0) then
                
                rsq = sqrt(aa(1)**2+aa(7)**2)
                aa(2) = -aa(7)/rsq
                aa(5) =  0.d0
                aa(8) =  aa(1)/rsq
                
              endif
              
              aa(3) = aa(4)*aa(8) - aa(7)*aa(5)
              aa(6) = aa(7)*aa(2) - aa(1)*aa(8)
              aa(9) = aa(1)*aa(5) - aa(4)*aa(2)
              
              call invert(aa,bb,det)
              
              if(abs(det).lt.1.d-5) call error(idnode,306)
              
              
              do j = 1,9
                gaxs(id,j) = bb(j)
              enddo
              
            elseif(ill.eq.0) then
c     
c     non-linear molecule
              
              i1 = 1
              i2 = 1
              i3 = 1
              
  210         i2=i2+1
              i3=i2
              
  220         i3=i3+1
              
              aa(1) = gxx(id,i1) - gxx(id,i2)
              aa(4) = gyy(id,i1) - gyy(id,i2)
              aa(7) = gzz(id,i1) - gzz(id,i2)
              aa(2) = gxx(id,i1) - gxx(id,i3)
              aa(5) = gyy(id,i1) - gyy(id,i3)
              aa(8) = gzz(id,i1) - gzz(id,i3)
              aa(3) = aa(4)*aa(8) - aa(7)*aa(5)
              aa(6) = aa(7)*aa(2) - aa(1)*aa(8)
              aa(9) = aa(1)*aa(5) - aa(4)*aa(2)
c     
c     invert matrix
              
              call invert(aa,bb,det)
c
c     check on size of determinant - to see if the 3 sites are
c     too close to being linear for safety.
              
              dettest=1.d-1
              
              if(abs(det).lt.dettest.and.i3.lt.numgsit(id)) goto 220
              if(abs(det).lt.dettest.and.i2.lt.numgsit(id)-1) goto 210
              if(abs(det).lt.dettest) call error(idnode,306)
c     
c     store indices used
              
              ind(id,1) = i1
              ind(id,2) = i2
              ind(id,3) = i3
c     
c     store coefficients 
              
              do j = 1,9
                
                gaxs(id,j) = bb(j)
                
              enddo
              
            endif
            
          endif
          
        enddo

c     
c     check that rigid unit does not contain frozen atoms
        
        safe=.true.
        
        jr =0
        do ig = igrp1,igrp2
          
          id = lstgtp(ig)
          
          do j = 1,numgsit(id)
            
            jr=jr+1
            i = lstrgd(jr)
            
            if(lstfrz(i).ne.0) safe = .false.
            
          enddo
          
        enddo
c     
c     global check on error condition
        
        if(mxnode.gt.1) call gstate(safe)
        if(.not.safe) call error(idnode,360)
c     
c     quaternions for all rigid groups in system
        
        jr =0
        do ig = igrp1,igrp2
          
          id = lstgtp(ig)
          i1 = lstrgd(jr + ind(id,1))
          i2 = lstrgd(jr + ind(id,2))
          i3 = lstrgd(jr + ind(id,3))
          
          jr = jr + numgsit(id)
c     
c     group basis vectors
          
          aa(1) = xxx(i1) - xxx(i2)
          aa(4) = yyy(i1) - yyy(i2)
          aa(7) = zzz(i1) - zzz(i2)
          
          call images(imcon,0,1,1,cell,aa(1),aa(4),aa(7))
          
          if(rotmin(id).gt.1.d-5) then
            
            aa(2) = xxx(i1) - xxx(i3)
            aa(5) = yyy(i1) - yyy(i3)
            aa(8) = zzz(i1) - zzz(i3)
            
          else
            
            rsq = sqrt(aa(1)**2+aa(4)**2+aa(7)**2)
            
            if(abs(aa(7)/rsq).gt.0.5d0) then
              
              rsq = sqrt(aa(4)**2+aa(7)**2)
              aa(2) =  0.d0
              aa(5) =  aa(7)/rsq
              aa(8) = -aa(4)/rsq
              
            elseif(abs(aa(4)/rsq).gt.0.5d0) then
              
              rsq = sqrt(aa(4)**2+aa(1)**2)
              aa(2) = -aa(4)/rsq
              aa(5) =  aa(1)/rsq
              aa(8) =  0.d0
              
            elseif(abs(aa(1)/rsq).gt.0.5d0) then
              
              rsq = sqrt(aa(1)**2+aa(7)**2)
              aa(2) = -aa(7)/rsq
              aa(5) =  0.d0
              aa(8) =  aa(1)/rsq
              
            endif
            
          endif
          
          call images(imcon,0,1,1,cell,aa(2),aa(5),aa(8))
          
          aa(3) = aa(4)*aa(8) - aa(7)*aa(5)
          aa(6) = aa(7)*aa(2) - aa(1)*aa(8)
          aa(9) = aa(1)*aa(5) - aa(4)*aa(2)
          
c     
c     group rotational matrix
          
          rot(1) = gaxs(id,1)*aa(1)+gaxs(id,4)*aa(2)+gaxs(id,7)*aa(3)
          rot(2) = gaxs(id,2)*aa(1)+gaxs(id,5)*aa(2)+gaxs(id,8)*aa(3)
          rot(3) = gaxs(id,3)*aa(1)+gaxs(id,6)*aa(2)+gaxs(id,9)*aa(3)
          rot(4) = gaxs(id,1)*aa(4)+gaxs(id,4)*aa(5)+gaxs(id,7)*aa(6)
          rot(5) = gaxs(id,2)*aa(4)+gaxs(id,5)*aa(5)+gaxs(id,8)*aa(6)
          rot(6) = gaxs(id,3)*aa(4)+gaxs(id,6)*aa(5)+gaxs(id,9)*aa(6)
          rot(7) = gaxs(id,1)*aa(7)+gaxs(id,4)*aa(8)+gaxs(id,7)*aa(9)
          rot(8) = gaxs(id,2)*aa(7)+gaxs(id,5)*aa(8)+gaxs(id,8)*aa(9)
          rot(9) = gaxs(id,3)*aa(7)+gaxs(id,6)*aa(8)+gaxs(id,9)*aa(9)
          
c     
c     determine quaternions from rotational matrix
          
          aq = rot(1)+rot(5)
          bq = rot(2)-rot(4)
          cq = rot(6)-rot(8)
          dq = rot(2)+rot(4)
          eq = rot(3)+rot(7)
          fq = rot(6)+rot(8)
          gq = rot(3)-rot(7)
          hq = rot(1)-rot(5)
          
          q0(ig) = 0.5d0*sqrt(aq+sqrt(aq*aq+bq*bq))
          
          if(q0(ig).gt.1.d-4) then
            
            q1(ig) = -0.25d0*cq/q0(ig)
            q2(ig) =  0.25d0*gq/q0(ig)
            q3(ig) = -0.25d0*bq/q0(ig)
            
          else
            
            q1(ig) = 0.5d0*sqrt(hq+sqrt(hq*hq+dq*dq))
            
            if(q1(ig).gt.1.d-4) then
              
              q2(ig) = 0.25d0*dq/q1(ig)
              q3(ig) = 0.25d0*eq/q1(ig)
              
            else
              
              q2(ig) = 0.5d0*sqrt(-hq+sqrt(hq*hq+dq*dq))
              
              if(q2(ig).gt.1.d-4) then
                
                q3(ig) = 0.25d0*fq/q2(ig)
                
              else
                
                q3(ig) = 1.d0
                
              endif
              
            endif
            
          endif
c     
c     normalise quaternions
          
          rnorm = 1.d0/sqrt(q0(ig)**2+q1(ig)**2+q2(ig)**2+q3(ig)**2)
          q0(ig) =rnorm*q0(ig)
          q1(ig) =rnorm*q1(ig)
          q2(ig) =rnorm*q2(ig)
          q3(ig) =rnorm*q3(ig)
          
        enddo
c     
c     test for redundant degrees of freedom
c     and ensure rotational inertias are non-zero
        
        degrot = 0.d0
        
        do ig = 1,ngrp
          
          id = lstgtp(ig)
          rotall=1.d0/max(1.d-5,rotinx(id,1)+rotiny(id,1)+
     x      rotinz(id,1))
          
          if(rotall*rotinx(id,1).lt.1.d-5) then
            degrot = degrot -1.d0
          endif
          
          if(rotall*rotiny(id,1).lt.1.d-5) then
            degrot = degrot -1.d0
          endif
          
          if(rotall*rotinz(id,1).lt.1d-5) then
            degrot = degrot -1.d0
          endif
          
        enddo
c     
c     rotational degrees of freedom and rigid body contribution
c     to total degrees of freedom
        
        degrot = degrot + dble(ngrp)*3.d0
        degfre = degrot + dble(ngrp)*3.d0
        
c
c     summarise results

        if(idnode.eq.0) then
          
          if(gmass(1).gt.0.d0) then
            
            write(nrite,'(/,/,12x,a)') ' summary of rigid body set up'
            
            do id = 1,mxungp
              
              if(gmass(id).gt.0.d0) then
                
                write(nrite,'(/,a,i10)') ' group of type ',id
                write(nrite,'(12x,a,f20.10)') ' total mass    ',
     x            gmass(id)
                write(nrite,'(12x,a,3f20.10)')' rot. inertia  ',
     x            rotinx(id,1),rotiny(id,1),rotinz(id,1)
                write(nrite,'(/,12x,a,3(8x,a7))') ' site','a coord',
     x            'b coord','c coord'
                do j = 1,numgsit(id)
                  write(nrite,'(12x,i5,1p,3e15.5)') j,gxx(id,j),
     x              gyy(id,j),gzz(id,j)
                enddo
              endif
            enddo
          endif
        endif
c
c     find number of unique groups 

        ngp = 0
        do ig = 1,ngrp
          ngp = max(ngp,lstgtp(ig))
        enddo
c
c     calculate reciprocal of rotational inertias 

        do id = 1,ngp
          
          rotlim=max(1.d-2,rotinx(id,1)+rotiny(id,1)+
     x      rotinz(id,1))*1.d-5
          
          if(rotinx(id,1).lt.rotlim) then
            rotinx(id,2) = 0.d0
          else
            rotinx(id,2) = 1.d0/rotinx(id,1)
          endif
          
          if(rotiny(id,1).lt.rotlim) then
            rotiny(id,2) = 0.d0
          else
            rotiny(id,2) = 1.d0/rotiny(id,1)
          endif
          
          if(rotinz(id,1).lt.rotlim) then
            rotinz(id,2) = 0.d0
          else
            rotinz(id,2) = 1.d0/rotinz(id,1)
          endif
          
        enddo
c     
c     Check of quaternion set up with atomic positions
        
        jr=0
        do ig = igrp1,igrp2
c
c     group type
          
          id = lstgtp(ig)
c     
c     new rotational matrix
          
          rot(1) = q0(ig)**2+q1(ig)**2-q2(ig)**2-q3(ig)**2
          rot(2) = 2.d0*(q1(ig)*q2(ig) - q0(ig)*q3(ig))
          rot(3) = 2.d0*(q1(ig)*q3(ig) + q0(ig)*q2(ig))
          rot(4) = 2.d0*(q1(ig)*q2(ig) + q0(ig)*q3(ig))
          rot(5) = q0(ig)**2-q1(ig)**2+q2(ig)**2-q3(ig)**2
          rot(6) = 2.d0*(q2(ig)*q3(ig) - q0(ig)*q1(ig))
          rot(7) = 2.d0*(q1(ig)*q3(ig) - q0(ig)*q2(ig))
          rot(8) = 2.d0*(q2(ig)*q3(ig) + q0(ig)*q1(ig))
          rot(9) = q0(ig)**2-q1(ig)**2-q2(ig)**2+q3(ig)**2
          
          do j = 1,numgsit(id)
            
            jr = jr +1
            i = lstrgd(jr)
            
            xxt(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+
     x        rot(3)*gzz(id,j)+gcmx(ig)
            yyt(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+
     x        rot(6)*gzz(id,j)+gcmy(ig)
            zzt(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+
     x        rot(9)*gzz(id,j) +gcmz(ig)
            
            
            txx(jr) = xxx(i) - xxt(i)
            tyy(jr) = yyy(i) - yyt(i)
            tzz(jr) = zzz(i) - zzt(i)
            
          enddo
          
        enddo
        
        call images(imcon,0,1,jr,cell,txx,tyy,tzz)
c
c     set tolerance for testing quaternion setup.

        rsq = 0.d0
        tol = 1.d-2
        
        do i = 1,jr
          
          ia = lstrgd(i)
          rrr=txx(i)**2 + tyy(i)**2 + tzz(i)**2
          if(rrr.gt.tol) then 
            rsq = rrr
          endif
        enddo
c     
c     exit if error in set up
        
        safe =.true.
        if(rsq.gt.tol) safe =.false.
        if(mxnode.gt.1) call gstate(safe)
        
        if(.not.safe) call  error(idnode,310)
        
c     
c     sort lstgot into ascending order
        
        call shellsort(jt,lstgot)
c     
c     check that no site is in more than 1 rigid group
        
        safe = .true.
        i=1
  400   i=i+1
        if(i.gt.jt) goto 420
        
  410   linear = .false.
        
        if (lstgot(i).eq.lstgot(i-1)) then
          
          linear = .true.
          safe = .false.
          jt = jt-1
          
          do j=i,jt
            
            lstgot(j) = lstgot(j+1)
            
          enddo
          
        endif
        
        if(i.ge.jt) linear =.false.
        if(linear) goto 410
        
        goto 400
  420   if(.not.safe) call error(idnode,320)
        
c     
c     list of 'free' sites
        
        ii = 1
        jj = 0
        do i =1,natms
          
          if(lstgot(ii).eq.i) then
            
            ii=ii+1
            
          else
            
            if(lstfrz(i).eq.0)then
              jj=jj+1
              lstfre(jj) = i
            endif
            
          endif
          
        enddo
c     
c     number of free sites
        
        ntfree = jj
c     
c     list of atoms integrated on this node
        
        jr = 0
        do ig = igrp1,igrp2
          
          id = lstgtp(ig)
          jr = jr + numgsit(id)
          
        enddo
        
        do i = 1,jr
          lstme(i) = lstrgd(i)
        enddo
c     
c     block parameters for free atoms
        
        ifre1 = (idnode*ntfree)/mxnode + 1
        ifre2 = ((idnode+1)*ntfree)/mxnode
        
        do i = ifre1,ifre2
          
          jr = jr+1
          lstme(jr)=lstfre(i)
          
        enddo
c     
c     sort  lstme into ascending order
        
        call shellsort(jr,lstme)
c     
c     exchange quaternion data with other nodes
        
        if(mxnode.gt.1) call merge4
     x    (idnode,mxnode,ngrp,nbuff,q0,q1,q2,q3,buffer)
        
      endif

#ifdef VAMPIR
      call VTEND(157, ierr)
#endif
      return
      end

