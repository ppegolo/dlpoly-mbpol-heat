      subroutine rdf0neu
     x  (ik,rcut,ilist,jlist,ltype,lstvdw,rdf,rsqdf)
c     
c***********************************************************************
c     
c     dl_poly subroutine for accumulating statistic for radial
c     distribution functions.
c     double precision accumulators
c     neutral group implementation
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1994
c     author    - t. forester    march 1994    
c     amended     t. forester    april 1994
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
      dimension ilist(mxxdf),jlist(mxxdf)
      dimension ltype(mxatms),lstvdw(mxvdw)
      data a0,a1,a2,a3,a4,a5/.0837557783d0,2.9399054d0,-7.8475201d0,
     x  14.1328992d0,-12.6228528d0,4.32084948d0/
      
#ifdef VAMPIR
      call VTBEGIN(75, ierr)
#endif
c     
c     set cutoff condition for pair forces
      
      rcsq=rcut**2
      rrcsq = 1.d0/rcsq
      sqlim = 0.01d0*rcsq
      
c     
c     grid interval for rdf tables
      
      rdelr = dble(mxrdf)/rcut
c     
c     start of primary loop for rdf accumulation
      
      do m=1,ik
        
c     
c     atomic and potential function indices
        
        iatm=ilist(m)
        ai=ltype(iatm)
        
        jatm=jlist(m) 
        aj=ltype(jatm)
        
        if(ai.gt.aj) then
          ll = int(ai*(ai-1.d0)*0.5d0+aj+0.5d0)
          k = lstvdw(ll)
        else
          ll = int(aj*(aj-1.d0)*0.5d0+ai+0.5d0)
          k = lstvdw(ll)
        endif
        
        rsq=rsqdf(m)
        
        if(rcsq.gt.rsq)then
          
c     determine interpolation panel for rdf table
          
          if(rsq.lt.sqlim) then
            
            rrr = sqrt(rsq)
            
          else
c     
c     interpolate square-root by polynomial plus newton-raphson
            
            sss=rsq*rrcsq
            rrr = 1.d0/
     x        (a0 +sss*(a1+sss*(a2+sss*(a3+sss*(a4+sss*a5)))))
            rrr = 0.5d0*rrr*(3.d0-sss*rrr*rrr)
            rrr = 0.5d0*rrr*(3.d0-sss*rrr*rrr)
            rrr = 0.5d0*rrr*(3.d0-sss*rrr*rrr)*sss*rcut
            
          endif
          
          ll=int(rrr*rdelr+0.999999d0)
c     
c     accumulate statistic
          
          rdf(ll,k) = rdf(ll,k) + 1.d0
          
        endif
        
      enddo
      
#ifdef VAMPIR
      call VTEND(75, ierr)
#endif
      return
      end
