      subroutine parlst
     x  (newlst,natms,idnode,mxnode,imcon,rcut,delr,lexatm,nexatm,
     x  noxatm,lentry,list,lstfrz,cell,xxx,yyy,zzz,xdf,ydf,zdf)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for constructing the verlet neighbour
c     list based on the brode-ahlrichs decomposition 
c     frozen atoms taken into account
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w.smith    march 1992
c     modified  - t.forester october 1993
c     
c     
c     wl
c     2001/06/12 12:55:37
c     1.6
c     $Sate: Exp $
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical lchk,newlst,lfrzi,ldo
      
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf)
      dimension lentry(msatms),list(msatms,mxlist)
      dimension lexatm(msatms,mxexcl),nexatm(msatms),lstfrz(mxatms)
      dimension noxatm(msatms),cell(9)
#ifdef VAMPIR
      call VTBEGIN(11, ierr)
#endif
      
      if(newlst)then
        ibig = 0
c
c     check size of work array
        if(mxxdf.lt.(natms+1)/2) then
          if(idnode.eq.0) write(nrite,*) 'mxxdf must be greater than ',
     x      (natms+1)/2
          call  error(idnode,474)
        endif
c     
c     set control variables
        
        last=natms
        lchk=.true.
        mpm2=natms/2
        npm2=(natms-1)/2
        
c     
c     set cutoff radius
        
        rclim=(rcut+delr)**2
        
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
            
            ii=ii+1
            xdf(ii)=xxx(i)-xxx(j)
            ydf(ii)=yyy(i)-yyy(j)
            zdf(ii)=zzz(i)-zzz(j)
            
          enddo
          
c     
c     apply minimum image convention
          
          call images(imcon,0,1,ii,cell,xdf,ydf,zdf)
          
c     
c     allocate atoms to neighbour list
          
          ii=0
          
          do i=idnode+1,last,mxnode
            
            lfrzi=(lstfrz(i).ne.0)
c     
c     calculate atom indices
            
            j=i+m
            if(j.gt.natms)j=j-natms
            
            ii =ii+1
c     
c     reject atoms in excluded pair list
            
            if((nexatm(ii).gt.0).and.(lexatm(ii,noxatm(ii)).eq.j))
     x        then
              
              noxatm(ii)=min(noxatm(ii)+1,nexatm(ii))
c     
c     reject frozen atom pairs

            else
              
              ldo = .true.
              if(lfrzi)  ldo = (lstfrz(j).eq.0)
              
c     
              if(ldo) then
c     calculate interatomic distance
                
                if(imcon.eq.6)then

                  rsq=xdf(ii)*xdf(ii)+ydf(ii)*ydf(ii)

                else

                  rsq=xdf(ii)*xdf(ii)+ydf(ii)*ydf(ii)+zdf(ii)*zdf(ii)

                endif
c     
c     running check of neighbour list array capacity
                
                if(rsq.lt.rclim)then
                  
                  lentry(ii)=lentry(ii)+1
                  
                  if(lentry(ii).gt.mxlist)then
                    
                    lchk=.false.
                    ibig = max(lentry(ii),ibig)

                  endif
                  
c     
c     compile neighbour list array
                  
                  if(lchk)then
                    
                    list(ii,lentry(ii))=j
                    
                  endif
                  
                endif

              endif
              
            endif
            
          enddo
          
        enddo
        
c     
c     terminate job if neighbour list array exceeded
        
        if(mxnode.gt.1) call gstate(lchk)
        
        if(.not.lchk) then
          call gimax(ibig,1,idum)
          if(idnode.eq.0) then
            write(nrite,*) ' mxlist must be at least  ',ibig
            write(nrite,*) ' mxlist is currently ',mxlist
          endif
          call error(idnode,110)
        endif
c     
c     check all excluded atoms are accounted for
        
        do i = 1,ii
          
          if(nexatm(i).gt.0.and.noxatm(i).ne.nexatm(i))lchk=.false.
          
        enddo
        
        if(mxnode.gt.1) call gstate(lchk)
        
        if(.not.lchk) call error(idnode,160)
        
      endif
#ifdef VAMPIR
      call VTEND(11, ierr)
#endif
      return
      end
