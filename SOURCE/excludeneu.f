      subroutine excludeneu
     x  (idnode,mxnode,nneut,lexatm,lexsit,
     x  neulst,nexatm,nexsit,nummols,numsit)

c     
c***********************************************************************
c     
c     dl_poly subroutine for constructing the excluded pair
c     interaction list of the system to be simulated
c     part 2 - neutral group implementation
c     
c     copyright - daresbury                1994
c     author    - t. forester        march 1994
c     
c     wl
c     2000/01/18 14:05:38
c     1.4
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"

      logical lchk,ldo
      
      dimension numsit(mxtmls)
      dimension nummols(mxtmls)
      dimension lexatm(msatms,mxexcl),nexatm(msatms)
      dimension lexsit(mxsite,mxexcl),nexsit(mxsite)
      dimension neulst(mxneut)
      
#ifdef VAMPIR
      call VTBEGIN(80, ierr)
#endif
c     
c     construct excluded pair list for verlet neighbour correction

      ibig = 0
      iatom=0
      jatom=0
c     
c     generate all atomic pairs and check for exclusions 
c     use Brode Ahlrichs ordering of groups
      
      last=nneut
      lchk=.true.
      mpm2=nneut/2+1
      npm2=(nneut-1)/2+1
      
c     
c     outer loop over groups
      
      do m=1,mpm2
        
        if(m.gt.npm2)last=mpm2-1
        
c     
c     inner loop over groups - include intragroup interactions
        
        ii=0
        
        do im=idnode+1,last,mxnode
          
          ii =ii+1
c     
c     first site in neutral group
          
          itmols =1
          inoff = 0
          isoff = 0
          isit = numsit(itmols)*nummols(itmols)
          iolsit = numsit(itmols)
          
c     
c     calculate j group indices
          
          jm=im+m-1
          if(jm.gt.nneut)jm=jm-nneut
          
c     
c     inner loop over neutral groups
          
          jtmols =1
          jnoff = 0
          jsoff = 0
          jsit = numsit(jtmols)*nummols(jtmols)
          jolsit = numsit(jtmols)
          
c     
c     test first sites in neutral group
          
          jatom = neulst(jm)         
c     
c     establish pointer to sets
          
  210     if(jatom.gt.jsit) then
            
            jtmols = jtmols+1
            jnoff = jsit
            jsoff = jsoff+ jolsit
            jsit = jsit + nummols(jtmols)*numsit(jtmols)
            jolsit = numsit(jtmols)
            
          endif
          
          if(jatom.gt.jsit) goto 210
          
          jn1 = jatom-jnoff
          jno1 = (jn1/jolsit)*jolsit
          jsite = jn1 -jno1
          if(jsite.eq.0) then 
            jsite = jolsit
            jno1 = jno1-jolsit
          endif
          jsite = jsite + jsoff
          jsite0 = jsite-1
          
          do iatom = neulst(im),neulst(im+1)-1
            
c     
c     establish pointer to sets
            
            if(iatom.gt.isit) then
              
  200         itmols =itmols+1
              inoff = isit
              isoff = isoff+ iolsit
              isit = isit + nummols(itmols)*numsit(itmols)
              iolsit = numsit(itmols)
              
              if(iatom.gt.isit) goto 200
              
            endif
            
            in1 = iatom-inoff
            ino1 = (in1/iolsit)*iolsit
            isite = in1 - ino1
            if(isite.eq.0) then 
              isite = iolsit
              ino1 = ino1-iolsit
            endif
            isite = isite + isoff
            
c     
c     test im and jm are neutral groups on same molecule
            
            if((jnoff.eq.inoff).and.(ino1.eq.jno1)) then
            if(abs(im-jm).lt.iolsit) then
              
              jj0 = neulst(jm)
              jsite =jsite0
              
c     
c     special case for im =jm (ie. same group)
              
              if(im.eq.jm) then 
                
                jj0 = iatom+1
                jsite=isite
                
              endif
c     
c     test for excluded interaction
              
              do jatom = jj0,neulst(jm+1)-1
                
                jsite=jsite+1
                
                do ij = 1,nexsit(isite)
                  
                  if(lexsit(isite,ij).eq.jsite-jsoff) then
                    
                    it = nexatm(ii)
                    
                    if(it+2.gt.mxexcl) then
                      
                      ibig = max(it+2,ibig)
                      nexatm(ii) = it+2
                      lchk =.false.
                      
                    else
                      
                      lexatm(ii,it+1) = iatom
                      lexatm(ii,it+2) = jatom
                      nexatm(ii) =nexatm(ii)+2
                      
                    endif

                  endif
                  
                enddo
                
              enddo
              
            endif
            endif
            
          enddo
          
        enddo
        
      enddo
c     
c     global check
      
      call gstate(lchk)
      if(.not.lchk) then
        
        if(mxnode.gt.1) call gimax(ibig,1,idum)
        if(idnode.eq.0) write(nrite,*) 'mxexcl must be at least ',ibig
        if(idnode.eq.0) write(nrite,*) 'mxexcl is currently ',mxexcl
        call error(idnode,260)

      endif

#ifdef VAMPIR
      call VTEND(80, ierr)
#endif
      return
      end
