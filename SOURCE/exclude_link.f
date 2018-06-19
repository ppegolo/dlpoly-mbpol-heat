      subroutine exclude_link
     x  (idnode,mxnode,ntpmls,lexatm,lexsit,
     x  nexatm,nexsit,nummols,numsit,lexatm2,nexatm2)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for constructing the excluded pair
c     interaction list of the system to be simulated
c     
c     part 2 - link cell implementation
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
      dimension lexatm2(msatms,mxexcl),nexatm2(msatms)
      dimension lexsit(mxsite,mxexcl),nexsit(mxsite)

c     
c     construct excluded pair list for verlet neighbour correction
      
#ifdef VAMPIR
      call VTBEGIN(79, ierr)
#endif
      iatom=0
      jatom=0
      lsite=0
      ksite=0
      nexatm2(:)=0
      
      do itmols=1,ntpmls
        
        do imols=1,nummols(itmols)
          
          do isite=1,numsit(itmols)
            
            iatom=iatom+1
            
            if(mod(iatom-1,mxnode).eq.idnode)then
              
              kk=0
              kk2=0 
              jatom=jatom+1
              
              do k=1,nexsit(ksite+isite)
                
                newatm=lexsit(ksite+isite,k)+lsite
                
                kk=kk+1
                lexatm(jatom,kk)=newatm
                
              enddo
              
              nexatm(jatom)=kk

              do k=1,numsit(itmols)

                 if (k.ne.isite) then

                    newatm=k+lsite
                    iflag=1
              
                    do l=1,kk

                       if (newatm.eq.lexatm(jatom,l)) iflag=0

                    enddo

                    if (iflag.eq.1) then

                       kk2=kk2+1
                       lexatm2(jatom,kk2)=newatm 

                    endif

                 endif

              enddo

              nexatm2(jatom)=kk2

            endif
            
          enddo
          
          lsite=lsite+numsit(itmols)
          
        enddo
        
        ksite=ksite+numsit(itmols)
        
      enddo
      
#ifdef VAMPIR
      call VTEND(79, ierr)
#endif
      return
      end

