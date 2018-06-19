      subroutine parlst_nsq
     x  (newlst,natms,idnode,mxnode,imcon,engcpe,epsq,rcut,
     x  vircpe,lexatm,nexatm,noxatm,lentry,list,lstfrz,
     x  cell,xxx,yyy,zzz,xdf,ydf,zdf,fxx,fyy,fzz,chge,stress)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for constructing the verlet neighbour
c     list based on the brode-ahlrichs decomposition .
c     Evaluates Coulombic interactions for everything outside verlet
c     shell - as part of 'secondary' forces.
c     frozen atom option included
c     
c     to be used with dlpoly_glob and multiple_nsq
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1994
c     author    - t. forester march 1994     
c     
c     stress tensor : t.forester may 1994
c     
c     wl
c     2000/01/18 14:05:51
c     1.5
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical lchk,newlst
      
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf)
      dimension lentry(msatms),list(msatms,mxlist),lstfrz(mxatms)
      dimension lexatm(msatms,mxexcl),nexatm(msatms)
      dimension noxatm(msatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension chge(mxatms),cell(9),stress(9)
#ifdef VAMPIR
      call VTBEGIN(9, ierr)
#endif
      if(newlst)then

        ibig=0
c
c     check size of work array
        if(mxxdf.lt.(natms+1)/2) then
          if(idnode.eq.0) write(nrite,*) 'mxxdf must be greater than ',
     x      (natms+1)/2
          call  error(idnode,475)
        endif
c     
c     zero force arrays
        
        do i = 1,natms
          
          fxx(i) = 0.d0
          fyy(i) = 0.d0
          fzz(i) = 0.d0
          
        enddo
c     
c     zero stress tensor

        do i = 1,9
          stress(i) = 0.d0
        enddo
c     
c     zero virial
        
        engcpe = 0.d0
        vircpe = 0.d0
        
c     
c     set control variables
        
        last=natms
        lchk=.true.
        mpm2=natms/2
        npm2=(natms-1)/2
        
c     
c     set cutoff radius - ignore border width
        
        rclim=(rcut)**2
        
c     
c     construct pair force neighbour list
        
        do i=1,msatms
          
          lentry(i)=0
          noxatm(i)=1
          
        enddo
        
c     
c     outer loop over atoms
        
        do m=1,mpm2
          
          if(m.gt.npm2)last=mpm2
          
c     
c     inner loop over atoms
          
          ii=0
          
          do i=idnode+1,last,mxnode
            
c     
c     calculate atom indices
            
            j=i+m
            if(j.gt.natms)j=j-natms
            
            if(lstfrz(i)*lstfrz(j).eq.0) then
              
              ii=ii+1
              
c     
c     calculate interatomic displacements
              
              xdf(ii)=xxx(i)-xxx(j)
              ydf(ii)=yyy(i)-yyy(j)
              zdf(ii)=zzz(i)-zzz(j)
              
            endif
            
          enddo
          
c     
c     apply minimum image convention
          
          call images(imcon,0,1,ii,cell,xdf,ydf,zdf)
c     
c     allocate atoms to neighbour list
          
          ii=0
          
          do i=idnode+1,last,mxnode
            
c     
c     calculate atom indices
            
            j=i+m
            if(j.gt.natms)j=j-natms
            
            ii=ii+1
            if((nexatm(ii).gt.0).and.(lexatm(ii,noxatm(ii)).eq.j))
     x        then
              
              noxatm(ii)=min(noxatm(ii)+1,nexatm(ii))
              
            elseif(lstfrz(i)*lstfrz(j).eq.0) then
c     
c     reject frozen atoms and calculate interatomic distance
              
              rsq=xdf(ii)**2 + ydf(ii)**2 + zdf(ii)**2
              
c     
c     running check of neighbour list array capacity
              
              if(rsq.lt.rclim)then
                
                lentry(ii)=lentry(ii)+1
                
                if(lentry(ii).gt.mxlist)then
                  
                  ibig=max(ibig,lentry(ii))
                  lchk=.false.

                endif
c     
c     compile neighbour list array
                
                if(lchk)then
                  
                  list(ii,lentry(ii))=j
                  
                endif
                
              else
                
c     
c     Coulombics for interactions outside Verlet list
                
                chgprd=chge(i)*chge(j)*r4pie0/epsq
                rrr=sqrt(rsq)
c     
c     calculate potential energy and virial
                
                coul = chgprd/rrr
                engcpe = engcpe + coul
c     
c     calculate forces
                
                fcoul = coul/rsq
                
                fxx(i) = fxx(i) + fcoul*xdf(ii)
                fyy(i) = fyy(i) + fcoul*ydf(ii)
                fzz(i) = fzz(i) + fcoul*zdf(ii)             
                
                fxx(j) = fxx(j) - fcoul*xdf(ii)
                fyy(j) = fyy(j) - fcoul*ydf(ii)
                fzz(j) = fzz(j) - fcoul*zdf(ii)  
c     
c     stress tensor

                stress(1) = stress(1) + xdf(ii)*fcoul*xdf(ii)
                stress(2) = stress(2) + xdf(ii)*fcoul*ydf(ii)
                stress(3) = stress(3) + xdf(ii)*fcoul*zdf(ii)

                stress(5) = stress(5) + ydf(ii)*fcoul*ydf(ii)
                stress(6) = stress(6) + ydf(ii)*fcoul*zdf(ii)

                stress(9) = stress(9) + zdf(ii)*fcoul*zdf(ii)

              endif
              
            endif
            
          enddo
          
        enddo
c     
c     complete stress tensor

        stress(4) = stress(2)
        stress(7) = stress(3)
        stress(8) = stress(6)
        
c     
c     terminate job if neighbour list array exceeded
        
        if(mxnode.gt.1)call gstate(lchk)
        
        if(.not.lchk) then

          call gisum(ibig,1,idum)
          if(idnode.eq.0) then
            write(nrite,*) ' mxlist must be >=  ',ibig
            write(nrite,*) ' mxlist is currenty ',mxlist
          endif
          call error(idnode,109)

        endif

c     
c     check all excluded atoms are accounted for
        
        do i = 1,ii
          
          if(nexatm(i).gt.0.and.noxatm(i).ne.nexatm(i))lchk=.false.
          
        enddo
        
        if(mxnode.gt.1)call gstate(lchk)
        
        if(.not.lchk) call error(idnode,160)
c     
c     calculate virial
        
        vircpe = -engcpe
        
      endif
#ifdef VAMPIR
      call VTEND(9, ierr)
#endif
      
      return
      end
