      subroutine revive
     x  (lgofr,lzden,lpolar,cfgname,atmnam,idnode,imcon,mxnode,
     x  natms,nstep,nzden,numacc,numrdf,chip,chit,chitd,conint,
     x  tstep,buffer,cell,fxx,fyy,fzz,ravval,rdf,ssqval,
     x  stkval,stpval,sumval,vxx,vyy,vzz,xxx,yyy,zumval,
     x  zzz,xx0,yy0,zz0,zdens,xxs,yys,zzs,eta,dipx,dipy,dipz,
     x  vdxx,vdyy,vdzz,xmsd,ymsd,zmsd,lcp,ldpts,conintd)
      use multibead, only: bead_suffix
c     
c***********************************************************************
c     
c     dl_poly subroutine for writing restart files at job termination
c     or at selected intervals in simulation
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith dec 1992.
c     
c     wl
c     2000/01/18 14:05:54
c     1.3
c     Exp
c     
c***********************************************************************
c     
#include "dl_params.inc"
      
      character*80 cfgname
      character*8 atmnam(mxatms)
      
      logical lgofr,lzden,lpolar,lcp,ldpts
      
      dimension rdf(mxrdf,mxvdw),eta(9)
      dimension zdens(mxrdf,mxsvdw)
      dimension cell(9),buffer(mxbuff)
      dimension stpval(mxnstk),sumval(mxnstk),ssqval(mxnstk)
      dimension zumval(mxnstk),ravval(mxnstk),stkval(mxstak,mxnstk)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension xx0(mxatms),yy0(mxatms),zz0(mxatms)
      dimension xxs(mxatms),yys(mxatms),zzs(mxatms)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension vdxx(mxatms),vdyy(mxatms),vdzz(mxatms)
      dimension xmsd(mxatms),ymsd(mxatms),zmsd(mxatms)
#ifdef VAMPIR
      call VTBEGIN(71, ierr)
#endif
c
c     define level of information in REVCON

      levcfg = 2
      
      if(mxnode.gt.1) then
c     
c     merge displacement data

        call merge(idnode,mxnode,natms,mxbuff,xx0,yy0,zz0,buffer)

c     
c     globally sum rdf information before saving
        
        if(lgofr) then
c     
c     maximum rdfs that can be summed in each step
          
          nsum = mxbuff/mxrdf
          if(nsum.eq.0) call error(idnode,200)
          
          nbuff = nsum*mxrdf
          
          do i = 1,mxvdw,nsum
            
            if((mxvdw+1-i).lt.nsum) nbuff=(mxvdw+1-i)*mxrdf
            call  gdsum(rdf(1,i),nbuff,buffer)
            
          enddo
          
        endif
c     
c     globally sum zden information before saving
        
        if(lzden) then
c     
c     maximum rdfs that can be summed in each step
          
          nsum = mxbuff/mxrdf
          if(nsum.eq.0) call error(idnode,200)
          
          nbuff = nsum*mxrdf
          
          do i =1,mxsvdw,nsum
            
            if((mxsvdw+1-i).lt.nsum) nbuff=(mxsvdw+1-i)*mxrdf
            call  gdsum(zdens(1,i),nbuff,buffer)
            
          enddo
          
        endif
        
      endif
c
c     node 0 handles i/o

      if(idnode.eq.0)then
c     
c     write configuration data to new configuration file

        open(nconf,file='REVCON'//bead_suffix,form='formatted')

        write(nconf,'(a80)') cfgname
        write(nconf,'(3i10,g20.10)') levcfg,imcon,nstep,tstep
        if(imcon.gt.0) write(nconf,'(3f20.10)') cell
        
        do i=1,natms
          
          write(nconf,'(a8,i10)') atmnam(i),i
          write(nconf,'(3g20.10)') xxx(i),yyy(i),zzz(i)
          if(levcfg.gt.0)write(nconf,'(3g20.10)')
     x      vxx(i),vyy(i),vzz(i)
          if(levcfg.gt.1)write(nconf,'(3g20.10)') 
     x      fxx(i),fyy(i),fzz(i)
          
        enddo
        
        close (nconf)
c     
c     write acccumulator data to dump file

        open(nrest,file='REVIVE'//bead_suffix,form='unformatted')
        
        write(nrest) dble(nstep),dble(numacc),dble(numrdf),chit,chip,
     x    conint,dble(nzden)
        write(nrest) eta
        write(nrest) stpval
        write(nrest) sumval
        write(nrest) ssqval
        write(nrest) zumval
        write(nrest) ravval
        write(nrest) stkval
        write(nrest) xx0,yy0,zz0
        write(nrest) xxs,yys,zzs
        write(nrest) xmsd,ymsd,zmsd
        if(lpolar) write(nrest) dipx,dipy,dipz
        if(lpolar.and.lcp) write(nrest) vdxx,vdyy,vdzz
        if(lpolar.and.lcp.and.ldpts) then
          write(nrest)chitd
          write(nrest)conintd
        endif
        if(lgofr) write(nrest) rdf
        if(lzden) write(nrest) zdens
        
        close (nrest)

      endif
c
c      if (lpolar .and. idnode.eq.0) then
c
c        open(28)
c        write(28,*)'nstep = ',nstep
c        do ii=1,natms
c           write(28,'(I6,3X,3f28.20)')
c     x      ii,dipx(ii),dipy(ii),dipz(ii)
c        enddo
c        do ii=1,natms
c           write(28,'(I6,3X,3f28.20)')
c     x      ii,vdxx(ii),vdyy(ii),vdzz(ii)
c        enddo
c        close(28)
c
c      endif
c     
c     divide rdf data between nodes
      
      rmxnode=1.d0/dble(mxnode)
      
      if(lgofr) then
        
        do i = 1,mxvdw
          
          do j = 1,mxrdf
            
            rdf(j,i) = rdf(j,i)*rmxnode
            
          enddo
          
        enddo
        
      endif
c     
c     divide zdensity data between nodes

      if(lzden) then
        
        do i = 1,mxsvdw
          
          do j = 1,mxrdf
            
            zdens(j,i) = zdens(j,i)*rmxnode
            
          enddo
          
        enddo
        
      endif
      
#ifdef VAMPIR
      call VTEND(71, ierr)
#endif
      return
      end
