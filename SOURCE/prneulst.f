      subroutine prneulst
     x  (newlst,imcon,idnode,mxnode,nneut,rprim,
     x  list,lentry,neulst,cell,xxx,yyy,zzz,xdf,ydf,zdf)
      
c***********************************************************************
c     
c     dlpoly routine to partition neutral group list into 
c     primary and secondary groups
c     loops over group ineu
c     
c     replicated data version
c     
c     copyright daresbury laboratory 1994
c     author t.forester april 1994
c     
c     wl
c     2000/01/18 14:05:53
c     1.5
c     Exp
c     
c***********************************************************************
      
#include "dl_params.inc"
      
      logical newlst,lchk,ldo
      
      dimension neulst(mxneut),lentry(msatms),list(msatms,mxlist)
      dimension cell(9),xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf)
      
#ifdef VAMPIR
      call VTBEGIN(118, ierr)
#endif
      lchk =.true.
      
      if(newlst)then
c     
c     set primary cutoff limit
        
        rclim=rprim*rprim
c     
c     set list to negative - signal for seconary shell
        
        ia =0
        do ineu = idnode+1,nneut,mxnode
          
          ia = ia +1
          
          do jj = 1,lentry(ia)
            
            list(ia,jj) = -abs(list(ia,jj))
            
          enddo
          
        enddo
c     
c     loop over neutral group ineu sites
        
        lchk = .true.
        ibig = 0
        
        ia = 0
        do ineu = idnode+1,nneut,mxnode
          
          ia = ia +1
          
          ii = 0
          do i = neulst(ineu),neulst(ineu+1)-1
            
            xi = xxx(i)
            yi = yyy(i)
            zi = zzz(i)
            
            do jj = 1,lentry(ia)
              
              jneu = -list(ia,jj)
              jj0 = neulst(jneu)
              
              if(ineu.eq.jneu) jj0 = i+1
c     
c     loop over jneu sites
              
              do j = jj0,neulst(jneu+1)-1
                
                ii=ii+1
                if(ii.le.mxxdf) then
                  xdf(ii)=xi-xxx(j)
                  ydf(ii)=yi-yyy(j)
                  zdf(ii)=zi-zzz(j)
                else
                  lchk = .false.
                  ibig = ii
                endif
              enddo
              
            enddo
            
          enddo
c     
c     apply minimum image convention
          
          ii = min(ii,mxxdf)
          call images(imcon,0,1,ii,cell,xdf,ydf,zdf)
          
c     
c     allocate groups to primary or secondary shell
c     on basis of closest atom-atom interactions
          
          ii = 0
          do i = neulst(ineu),neulst(ineu+1)-1
            
            do jj = 1,lentry(ia)
              
              rr = rclim+1.d0
              jneu = list(ia,jj)
              ldo = (jneu.lt.0)
              if(ldo) then
                jneu = -jneu
                jj0 = neulst(jneu)
                if(ineu.eq.jneu) jj0 = i+1
                
c     
c     loop over jneu sites
                
                
                do j = jj0,neulst(jneu+1)-1
                  
                  if(ldo) then 

                    ii=min(ii+1,mxxdf)

                    if(abs(xdf(ii)).lt.rprim) then
                    if(abs(ydf(ii)).lt.rprim) then
                    if(abs(zdf(ii)).lt.rprim) then
c     
c     calculate interatomic distance
                    
                    rrr=xdf(ii)**2+ydf(ii)**2+zdf(ii)**2
c     
c     put in primary list if found any interaction close enough
                    
                    if(rrr.le.rclim) then
                      ldo = .false.
                      list(ia,jj) = jneu
                    endif
                    
                    endif
                    endif
                    endif

                  endif
                  
                enddo
                
              endif
              
            enddo
            
          enddo
          
        enddo
        
        if(mxnode.gt.1) call gstate(lchk)
        if(.not.lchk) then
          call gimax(ibig,1,idum)
          if(idnode.eq.0) then
             write(nrite,*) 'mxxdf must be at least ',ibig
             write(nrite,*) 'mxxdf is currently     ',mxxdf
          endif
          call  error(idnode,477)
        endif
        
      endif
      
#ifdef VAMPIR
      call VTEND(118, ierr)
#endif
      return
      end

