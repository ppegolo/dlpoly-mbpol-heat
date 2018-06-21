      subroutine rdf1
     x  (lpgr,cfgname,unqatm,idnode,mxnode,ntpatm,ntpvdw,numrdf,
     x  rcut,volm,lstvdw,dens,rdf,buffer)

c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating radial distribution functions
c     from accumulated data.
c     double precision version
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1994
c     author    - t. forester    march 1994    
c     
c     wl
c     2003/05/08 08:45:12
c     1.5
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical lpgr,zero

      character*80 cfgname
      character*8 unqatm(mxsite)

      dimension rdf(mxrdf,mxvdw)
      dimension lstvdw(mxvdw)
      dimension dens(mxsvdw)
      dimension buffer(mxbuff)

#ifdef VAMPIR
      call VTBEGIN(169, ierr)
#endif
      if(idnode.eq.0) write(nrite,
     x  "(/,/,12X,'RADIAL DISTRIBUTION FUNCTIONS',/,/,
     x  'calculated using ',i10,' configurations')") numrdf

      if(lpgr) then

c     open RDF file and write headers

        if(idnode.eq.0)then

          open(nrdfdt,file='RDFDAT')

          write(nrdfdt,'(a80)')cfgname
          write(nrdfdt,'(2i10)')ntpvdw,mxrdf

        endif
c     
c     grid interval for rdf tables
        
        delr = rcut/dble(mxrdf)
        
        do ia = 1,ntpatm

          do ib = ia,ntpatm
            
            k = lstvdw(ib*(ib-1)/2+ia)

            if(k.le.ntpvdw) then

              if(idnode.eq.0) then

                write(nrite,
     x            "(/,'g(r)  :',2a8,/,/,8x,'r',6x,'g(r)',9x,'n(r)',/)") 
     x            unqatm(ia),unqatm(ib)
                write(nrdfdt,'(2a8)')unqatm(ia),unqatm(ib)

              endif

c     
c     global sum of data on all nodes

              if(mxnode.gt.1) call gdsum(rdf(1,k),mxrdf,buffer)
c     
c     normalisation factor 
              
              factor = volm*dens(ia)*dens(ib)*dble(numrdf)
              if((ia.eq.ib).and.(volm*dens(ia).gt.1.d0)) 
     x          factor=factor*0.5d0*dens(ia)/(dens(ia)-1.d0/volm)
c     
c     running integration of rdf
              
              sum = 0.d0
c     
c     loop over distances

              zero = .true.
              
              do j = 1,mxrdf
                
                if(zero.and.(j.lt.mxrdf-3)) 
     x            zero = (rdf(j+2,k).le.0.d0)
                
                rrr = (dble(j)-0.5d0)*delr
                dvol = 4.d0*pi*(delr*rrr*rrr+(delr**3)/12.d0)
                
                gofr = rdf(j,k)/(factor*dvol)
                sum = sum + gofr*dvol*dens(ib)
c     
c     print out information
                
                if(idnode.eq.0) then
                  
                  write(nrdfdt,"(1p,2e14.6)")rrr,gofr
                  if(.not.zero)
     x              write(nrite,"(f10.4,1p,2e14.6)")rrr,gofr,sum
                  
                endif
                
              enddo
              
            endif  

          enddo
          
        enddo

        if(idnode.eq.0)close (nrdfdt)

      endif

#ifdef VAMPIR
      call VTEND(169, ierr)
#endif
      return
      end
