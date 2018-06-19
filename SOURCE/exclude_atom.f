      subroutine exclude_atom
     x  (idnode,mxnode,natms,ntpmls,lexatm,
     x  lexsit,nexatm,nexsit,nummols,numsit)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for constructing the excluded pair
c     interaction list of the system to be simulated
c     part 2 
c
c     copyright - daresbury laboratory 1992
c     author    - w. smith        june 1992
c     
c     wl
c     2000/01/18 14:05:38
c     1.3
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension numsit(mxtmls),nummols(mxtmls)
      dimension lexatm(msatms,mxexcl),nexatm(msatms)
      dimension lexsit(mxsite,mxexcl),nexsit(mxsite)

#ifdef VAMPIR
      call VTBEGIN(78, ierr)
#endif
c     
c     construct excluded pair list for verlet neighbour correction
      
      iatom=0
      jatom=0
      lsite=0
      ksite=0
      
      do itmols=1,ntpmls
        
        do imols=1,nummols(itmols)
          
          do isite=1,numsit(itmols)
            
            iatom=iatom+1
            
            if(mod(iatom-1,mxnode).eq.idnode)then
              
              kk=0
              jatom=jatom+1
              
              do k=1,nexsit(ksite+isite)
                
                newatm=lexsit(ksite+isite,k)+lsite
c     
c     keep only brode-ahlrichs combinations of indices
                
                if(((newatm.gt.iatom).and.
     x            (newatm-iatom.le.natms/2)).or.
     x            ((newatm.lt.iatom).and.
     x            (newatm+natms-iatom.le.(natms-1)/2)))then
                  
                  kk=kk+1
                  lexatm(jatom,kk)=newatm
                  
                  if(kk.gt.1)then
                    
c     
c     sort the excluded atom list in ascending indices
                    do j=kk,2,-1
                      
                      if(lexatm(jatom,j).lt.lexatm(jatom,j-1))
     x                  then
                        latom=lexatm(jatom,j)
                        lexatm(jatom,j)=lexatm(jatom,j-1)
                        lexatm(jatom,j-1)=latom
                      endif
                      
                    enddo
                    
                  endif
                  
                endif
                
              enddo
              
              nexatm(jatom)=kk
              
            endif
            
          enddo
          
          lsite=lsite+numsit(itmols)
          
        enddo
        
        ksite=ksite+numsit(itmols)
        
      enddo
      
c     
c     final sort into brode-ahlrichs ordering
      
      ii=0
      do i=1+idnode,natms,mxnode
        
        ii=ii+1
        do j=1,nexatm(ii)
          
          if(lexatm(ii,1).lt.i)then
            
            latom=lexatm(ii,1)
            
            do k=1,nexatm(ii)-1
              
              lexatm(ii,k)=lexatm(ii,k+1)
              
            enddo
            
            lexatm(ii,nexatm(ii))=latom
            
          endif
          
        enddo
        
      enddo
      
#ifdef VAMPIR
      call VTEND(78, ierr)
#endif
      return
      end

