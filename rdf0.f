      subroutine rdf0
     x  (iatm,ik,rcut,ilist,ltype,lstvdw,rdf,rsqdf)
c     
c***********************************************************************
c     
c     dl_poly subroutine for accumulating statistic for radial
c     distribution functions.
c     double precision accumulators
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1994
c     author    - t. forester    march 1994    
c     
c     wl
c     2000/01/18 14:05:54
c     1.3
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension rsqdf(mxxdf)
      dimension rdf(mxrdf,mxvdw)
      dimension ilist(mxxdf)
      dimension ltype(mxatms),lstvdw(mxvdw)

#ifdef VAMPIR
      call VTBEGIN(72, ierr)
#endif
c     
c     set cutoff condition for pair forces
      
      rcsq=rcut*rcut
c     
c     grid interval for rdf tables

      rdelr = dble(mxrdf)/rcut
c     
c     set up atom iatm type

      ai =ltype(iatm) 
c     
c     start of primary loop for rdf accumulation
      
      do m=1,ik
c     
c     atomic and potential function indices
        
        jatm=ilist(m)
        aj=ltype(jatm)
        if(ai.gt.aj) then
          ll = int(ai*(ai-1.d0)*0.5d0+aj+0.5d0)
          k = lstvdw(ll)
        else
          ll = int(aj*(aj-1.d0)*0.5d0+ai+0.5d0)
          k = lstvdw(ll)
        endif

c     
c     apply truncation of potential

        rsq=rsqdf(m)

        if(rcsq.gt.rsq)then

          rrr = sqrt(rsq)
          ll=int(rrr*rdelr+0.999999d0)
c     
c     accumulate statistic

          rdf(ll,k) = rdf(ll,k) + 1.d0

        endif
        
      enddo
      
#ifdef VAMPIR
      call VTEND(72, ierr)
#endif
      return
      end
