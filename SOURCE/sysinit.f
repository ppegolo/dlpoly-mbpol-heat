      subroutine sysinit
     x  (lgofr,lzden,lmetal,idnode,imcon,keyfce,keyres,mxnode,
     x  natms,nstep,numacc,numrdf,ntpatm,nzden,chip,chit,chitd,
     x  conint,elrc,engunit,virlrc,rvdw,volm,lstvdw,ltpvdw,lstfrz,
     x  ltype,numtyp,numfrz,buffer,cell,dens,prmvdw,ravval,rdf,
     x  ssqval,stkval,stpval,sumval,xx0,yy0,zz0,zumval,zdens,
     x  xxs,yys,zzs,elrcm,vlrcm,eta,dipx,dipy,dipz,
     x  vdxx,vdyy,vdzz,xmsd,ymsd,zmsd,lpolar,lcp,ldpts,conintd)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for reading the REVIVE file data and 
c     defining the initial thermodynamic and structural accumulators.
c     
c     copyright - daresbury laboratory 1997
c     author    - w. smith         july 1997
c     
c     wl
c     2001/05/30 12:40:26
c     1.7
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical lgofr,lzden,lmetal,lpolar,lcp,ldpts
      
      dimension lstvdw(mxvdw),ltpvdw(mxvdw),lstfrz(mxatms)
      dimension cell(9),prmvdw(mxvdw,mxpvdw)
      dimension numtyp(mxsvdw),numfrz(mxsvdw),ltype(mxatms)
      dimension buffer(mxbuff),dens(mxsvdw),elrcm(2),vlrcm(2)
      dimension xx0(mxatms),yy0(mxatms),zz0(mxatms)
      dimension xxs(mxatms),yys(mxatms),zzs(mxatms)
      dimension stpval(mxnstk),sumval(mxnstk),ssqval(mxnstk)
      dimension zumval(mxnstk),ravval(mxnstk),stkval(mxstak,mxnstk)
      dimension eta(9),rdf(mxrdf,mxvdw),zdens(mxrdf,mxsvdw)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension vdxx(mxatms),vdyy(mxatms),vdzz(mxatms)
      dimension xmsd(mxatms),ymsd(mxatms),zmsd(mxatms)

      character*80 ccc       

#ifdef VAMPIR
      call VTBEGIN(7, ierr)
#endif
c     
c     read or initialise accumulator arrays
      
      if(keyres.eq.1.and.idnode.eq.0)then
c     
c     read accumulator data from dump file

        open(nrest,file='REVOLD',form='unformatted')
        
        read(nrest) dnstep,dnumac,dnumrd,chit,chip,conint,dnzden
        read(nrest) eta
        read(nrest) stpval
        read(nrest) sumval
        read(nrest) ssqval
        read(nrest) zumval
        read(nrest) ravval
        read(nrest) stkval
        read(nrest) xx0,yy0,zz0
        read(nrest) xxs,yys,zzs
        read(nrest) xmsd,ymsd,zmsd
        if(lpolar) read(nrest) dipx,dipy,dipz
        if(lpolar.and.lcp) read(nrest) vdxx,vdyy,vdzz
        if(lpolar.and.lcp.and.ldpts) then
           read(nrest) chitd
           read(nrest) conintd
        endif
        if(lgofr) read(nrest) rdf
        if(lzden) read(nrest) zdens

        nstep=nint(dnstep)
        numacc=nint(dnumac)
        numrdf = nint(dnumrd)
        nzden = nint(dnzden)
        close (nrest)

c        if(lpolar.and.idnode.eq.0) then
c          read(28,*)ccc
c          do i=1,mxatms
c            read(28,*)dum,dipx(i),dipy(i),dipz(i)
c          enddo
c        endif
        
      else
        
c     
c     initialise step counters
        
        nstep=0
        numacc=0
        numrdf=0
        nzden=0
c     
c     initialise temperature and pressure coupling parameters
c     and integral for conserved quantity
        
        chit = 0.d0
        chip = 0.d0
        conint = 0.d0
        do i = 1,9
          eta(i) = 0.d0
        enddo
c
c     initialize center-of-mass displacement

        do i=1,mxatms
          xmsd(i)=0.d0
          ymsd(i)=0.d0
          zmsd(i)=0.d0
        enddo
c
c     initialize dipoles

        if (lpolar) then

          do i=1,mxatms
            dipx(i)=0.d0
            dipy(i)=0.d0
            dipz(i)=0.d0
          enddo

          if (lcp) then
            do i=1,mxatms
              vdxx(i)=0.d0
              vdyy(i)=0.d0
              vdzz(i)=0.d0
            enddo
            chitd=0.d0
            conintd=0.d0
          endif

        endif

c        if(lpolar .and. keyres.ne.1) then
c          read(28,*)ccc
c          do i=1,mxatms
c            read(28,*)dum,dipx(i),dipy(i),dipz(i)
c          enddo
c        endif
c     
c     initialise accumulator arrays
        
        do i=1,mxnstk
          
          stpval(i)=0.d0
          sumval(i)=0.d0
          ssqval(i)=0.d0
          zumval(i)=0.d0
          ravval(i)=0.d0
          
        enddo
        
        do i=1,mxatms

          xx0(i)=0.d0
          yy0(i)=0.d0
          zz0(i)=0.d0
          xxs(i)=0.d0
          yys(i)=0.d0
          zzs(i)=0.d0

        enddo

        do j=1,mxnstk
          
          do i=1,mxstak
            
            stkval(i,j)=0.d0
            
          enddo
          
        enddo
        
        if(lgofr) then
          
          do i = 1,mxvdw
            
            do j = 1,mxrdf
              
              rdf(j,i) = 0.d0
              
            enddo
            
          enddo

        endif

        if(lzden) then

          do i = 1,mxsvdw
            do j = 1,mxrdf
              zdens(j,i) = 0.d0
            enddo

          enddo
          
        endif
        
      endif
c     
c     if restart then broadcast stored variables via a global sum

      if(keyres.eq.1.and.mxnode.gt.1) then

        if(mxbuff.lt.natms.or.mxbuff.lt.mxnstk*mxstak)
     x    call error(idnode,186)

        call gisum(nstep,1,idum)
        call gisum(numacc,1,idum)
        call gisum(numrdf,1,idum)
        call gdsum(chit,1,buffer)
        call gdsum(chitd,1,buffer)
        call gdsum(chip,1,buffer)
        call gdsum(conint,1,buffer)
        call gdsum(conintd,1,buffer)
        call gisum(nzden,1,idum)
        call gdsum(eta,9,buffer)
        call gdsum(stpval,mxnstk,buffer)
        call gdsum(sumval,mxnstk,buffer)
        call gdsum(ssqval,mxnstk,buffer)
        call gdsum(zumval,mxnstk,buffer)
        call gdsum(ravval,mxnstk,buffer)      
        call gdsum(stkval,mxnstk*mxstak,buffer)
        call gdsum(xx0,natms,buffer)
        call gdsum(yy0,natms,buffer)
        call gdsum(zz0,natms,buffer)
        call gdsum(xxs,natms,buffer)
        call gdsum(yys,natms,buffer)
        call gdsum(zzs,natms,buffer)
        call gdsum(dipx,natms,buffer)
        call gdsum(dipy,natms,buffer)
        call gdsum(dipz,natms,buffer)
        call gdsum(vdxx,natms,buffer)
        call gdsum(vdyy,natms,buffer)
        call gdsum(vdzz,natms,buffer)
        call gdsum(xmsd,natms,buffer)
        call gdsum(ymsd,natms,buffer)
        call gdsum(zmsd,natms,buffer)
c     
c     for rdf table - broadcast and normalise
        if(lgofr) then

          do k = 1,mxvdw
            call gdsum(rdf(1,k),mxrdf,buffer)
            
            do j = 1,mxrdf
              rdf(j,k) = rdf(j,k)/dble(mxnode)
            enddo

          enddo
          
        endif

        if(lzden) then

          do k = 1,mxsvdw
            call gdsum(zdens(1,k),mxrdf,buffer)
            do j = 1,mxrdf
              zdens(j,k) = zdens(j,k)/dble(mxnode)
            enddo
          enddo

        endif

      endif
c     
c     number densities and long-range corrections
      
      elrc = 0.d0       
      virlrc = 0.d0
      elrcm(1)=0.d0
      elrcm(2)=0.d0
      vlrcm(1)=0.d0
      vlrcm(2)=0.d0
      
      if(imcon.eq.0.or.imcon.eq.6) volm = 4.d0*pi/3.d0*rvdw**3
      
      call lrcorrect
     x  (idnode,imcon,keyfce,mxnode,natms,ntpatm,elrc,engunit,virlrc,
     x  rvdw,volm,lstvdw,ltpvdw,ltype,numtyp,numfrz,lstfrz,prmvdw,dens)
      
      if(lmetal) call lrcmetal
     x  (idnode,imcon,mxnode,natms,ntpatm,engunit,rvdw,volm,
     x  lstvdw,ltpvdw,ltype,numtyp,prmvdw,dens,elrcm,vlrcm)

      if(imcon.eq.0.or.imcon.eq.6) volm = 0.d0
#ifdef VAMPIR
      call VTEND(7, ierr)
#endif
      return
      end
