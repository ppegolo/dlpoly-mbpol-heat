      subroutine static
     x  (lzeql,cfgname,idnode,intsta,imcon,keyens,natms,nstack,nstep,
     x  nsteql,ntpatm,numacc,mxnode,consv,degfre,degrot,engang,
     x  engbnd,engcpe,engdih,enginv,engke,engrot,engsrp,engunit,stpcfg,
     x  stpeng,stpeth,stpprs,stptmp,stpvir,stpvol,tstep,virbnd,engfbp,
     x  vircom,vircon,vircpe,virsrp,engfld,virfld,engtbp,virtbp,
     x  virpmf,virshl,engshl,engtet,virtet,degshl,shlke,virang,
     x  width,ltype,numtyp,buffer,cell,chge,fxx,fyy,fzz,ravval,
     x  ssqval,stkval,stpval,sumval,vxx,vyy,vzz,xxx,yyy,zzz,zumval,
     x  xx0,yy0,zz0,weight,stress,amsd,engdke,consvd,stress_pimd)
      use multibead, only: bead_suffix
c     
c***********************************************************************
c     
c     dl_poly subroutine for accumulating periodic data during the
c     molecular dynamics simulation and computing the rolling averages
c     
c     copyright daresbury laboratory 1992
c     
c     author - w. smith       august 1992
c     
c     wl
c     2001/08/31 11:13:52
c     1.11
c     $Sate: Exp $
c     
c***********************************************************************
c     

#include "dl_params.inc"
      
      logical lzeql,newjob,lfirst
      
      character*80 cfgname
      dimension ltype(mxatms),numtyp(mxsvdw)
      dimension cell(9),celprp(10),buffer(mxbuff)
      dimension stpval(mxnstk),sumval(mxnstk),ssqval(mxnstk)
      dimension zumval(mxnstk),ravval(mxnstk),stkval(mxstak,mxnstk)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms),chge(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension xx0(mxatms),yy0(mxatms),zz0(mxatms)
      dimension amsd(mxsvdw),weight(mxatms),stress(9)

!FP: added
      real(8) :: stress_pimd(9)
!FP: end

      save newjob,lfirst
      
      data newjob/.true./
      data lfirst/.true./
#ifdef VAMPIR
      call VTBEGIN(66, ierr)
#endif
c
c     write pressure tensor

c      if (lfirst .and. (idnode.eq.0)) then 
c        open(35,position='append')
c        lfirst = .false.
c      endif
c     
c     open statistics file for append
#ifndef IPI
      
      if(nstep.eq.intsta)then
        
         if(idnode.eq.0)then

           open(nstats,file='STATIS'//bead_suffix,position='append')
           
           write(nstats,'(a80)') cfgname
           if(engunit.eq.9648.530821d0) write(nstats,
     x          "(' ENERGY UNITS = electron Volts ')")
           if(engunit.eq.418.4d0)       write(nstats,
     x          "(' ENERGY UNITS = kcal/mol ')")
           if(engunit.eq.1.d2)          write(nstats,
     x          "(' ENERGY UNITS = kjoule/mol ')")
           if(engunit.eq.1.d0)          write(nstats,
     x          "(' ENERGY UNITS = DL_POLY Internal UNITS ')")
        
         endif

         newjob=.false.
        
      elseif(newjob.and.idnode.eq.0.and.mod(nstep,intsta).eq.0)then

        open(nstats,file='STATIS'//bead_suffix,position='append')
        newjob=.false.

      endif
#endif /* IPI */ 

c     
c     calculate cell volume and minimum cell half-width
      
      if(imcon.eq.0)then
        
        width=0.d0
        
        stpvol=0.d0
        do i = 1,10
           celprp(i) = 0.d0
        enddo
        
      else
        
        call dcell(cell,celprp)
        stpvol=celprp(10)
        width=min(celprp(7),celprp(8),celprp(9))/2.d0
        
        if(imcon.eq.4)then

          stpvol=0.5d0*celprp(10)
          width=sqrt(3.d0)*cell(1)/4.d0

        elseif(imcon.eq.5)then
        
          stpvol=0.5d0*celprp(10)
          width=cell(1)/2.d0

        elseif(imcon.eq.6)then

          width=min(celprp(7),celprp(8))/2.d0

        elseif(imcon.eq.7)then
        
          stpvol=0.5d0*celprp(10)

        endif
        
      endif
c     
c     energetic properties of system
      
      stpcfg=engsrp+engcpe+engbnd+engang+engdih+engfld+engtbp+engfbp
     x  +engshl+enginv
      stpvir=virsrp+vircpe+virbnd+vircon+vircom+virtbp+virang
     x  +virshl+virtet
c      stpeng=stpcfg+engke+engrot
      stpeng=stpcfg+engke+engrot+engdke
      stprot= 2.d0*(engrot)/(boltz*max(1.d0,degrot))
      stpshl= 2.d0*(shlke)/(boltz*max(1.d0,degshl))
      stptmp=2.d0*(engke+engrot)/(boltz*degfre)
      stpprs=0.d0
      if(imcon.gt.0)stpprs=(2.d0*engke-stpvir)/(3.d0*stpvol)
      stpeth=stpeng+stpprs*stpvol
c      stpcns=stpeng
c      stpcns=stpeng+consv
      stpcns=stpeng+consv+consvd
c     
c     convert pressure to units of kbar
      
      stpprs = stpprs*prsunt
c     
c     calculate mean squared displacements 
c     atomic displacements from origin of production run
      
      if((.not.lzeql).or.(nstep.gt.nsteql)) then
        call diffsn0
     x    (idnode,natms,mxnode,tstep,vxx,vyy,vzz,xx0,yy0,zz0)
      
        call  diffsn1
     x    (idnode,natms,ntpatm,mxnode,ltype,numtyp,amsd,
     x    xx0,yy0,zz0,buffer)
        
      endif

c     
c     zero statistics arrays
      
      if((nstep.le.0).or.(numacc.eq.0))then
        
        numacc=0
        
        do i=1,mxnstk
          
          stpval(i)=0.d0
          sumval(i)=0.d0
          ssqval(i)=0.d0
          
        enddo
        
        do i=1,mxatms
          
          xx0(i) = 0.d0
          yy0(i) = 0.d0
          zz0(i) = 0.d0
          
        enddo
        
      endif
      
c     
c     store current values in statistics array
      
      stpval(1) =stpcns/engunit
      stpval(2) =stptmp
      stpval(3) =stpcfg/engunit
      stpval(4) =engsrp/engunit
      stpval(5) =engcpe/engunit
      stpval(6) =engbnd/engunit
      stpval(7) =(engang+engtbp)/engunit
      stpval(8) =(engdih+enginv+engfbp)/engunit
      stpval(9) =engtet/engunit
      stpval(10)=stpeth/engunit
      stpval(11)=stprot
      stpval(12)=stpvir/engunit
      stpval(13)=virsrp/engunit
      stpval(14)=vircpe/engunit
      stpval(15)=virbnd/engunit
      stpval(16)=(virtbp+virang)/engunit
      stpval(17)=vircon/engunit
      stpval(18)=virtet/engunit
      stpval(19)=stpvol
      stpval(20)=stpshl
      stpval(21)=engshl/engunit
      stpval(22)=virshl/engunit
      stpval(23)=acos(celprp(6))*180.d0/pi
      stpval(24)=acos(celprp(5))*180.d0/pi
      stpval(25)=acos(celprp(4))*180.d0/pi
      stpval(26)=virpmf/engunit
      stpval(27)=stpprs

      iadd = 27
c     
c     mean squared displacements 
      
      if((.not.lzeql).or.(nstep.gt.nsteql)) then
        
        do k = 1,ntpatm
          
          stpval(iadd+k)=amsd(k)
          
        enddo
        
      endif

      iadd = iadd+ntpatm
      idum = iadd

#ifdef STRESS      
      if(stpvol.eq.0.d0) stpvol= 1.d0
      do i = 1,9
        stpval(iadd+i)=stress(i)*prsunt/(stpvol)
!FP: added
        stress_pimd(i) = stress(i)*prsunt/(stpvol)
!FP: end
      enddo
      iadd = iadd+9
      tdum=dble(nstep)*tstep
c      if (nstep.eq.1 .and. idnode.eq.0) 
c     x  write(35,'(3f12.6)') cell(1),cell(5),cell(9)
c      if (idnode.eq.0)   
c     x    write(35,'(f10.4,3x,6f12.6)')tdum,stpval(idum+2),
c     x        stpval(idum+3),stpval(idum+6),stpval(idum+1),
c     x        stpval(idum+5),stpval(idum+9)
#endif
      if(keyens.gt.3.and.(keyens.le.7)) then
        do i = 1,9
          stpval(iadd+i)=cell(i)
        enddo
        iadd = iadd+9
      endif
c     
c     check on number of variables for stack - 
      
      if(iadd.gt.mxnstk) call error(idnode,170)
c     
c     accumulate totals over steps
      
      numacc=numacc+1
      sclnv2=1.d0/dble(numacc)
      sclnv1=dble(numacc-1)/dble(numacc)
      
      do i=1,mxnstk
        
        ssqval(i)=sclnv1*(ssqval(i)+sclnv2*(stpval(i)-sumval(i))**2)
        sumval(i)=sclnv1*sumval(i)+sclnv2*stpval(i)
        
      enddo
      
c     
c     write statistics file
      
#ifndef IPI
      if(idnode.eq.0.and.mod(nstep,intsta).eq.0)then
        
        write(nstats,'(i10,1p,e14.6,0p,i10,/,(1p,5e14.6))')
     x    nstep,nstep*tstep,iadd,(stpval(k),k=1,iadd)
        
      endif
#endif /* IPI */
     
c     
c     zero rolling average accumulators
      
      if(nstep.le.0)then
        
        numacc=0
        
        do i=1,mxnstk
          
          zumval(i)=0.d0
          
          do j=1,mxstak
            
            stkval(j,i)=0.d0
            
          enddo
          
        enddo
        
      endif
      
c     
c     store quantities in stack
      
      kstak=mod(nstep-1,nstack)+1
      
      if(nstep.gt.nstack)then
        
        do i=1,mxnstk
          
          zumval(i)=zumval(i)-stkval(kstak,i)
          
        enddo
        
      endif
      
      do i=1,mxnstk
        
        stkval(kstak,i)=stpval(i)
        zumval(i)=zumval(i)+stpval(i)
        
      enddo
      
c     
c     calculate rolling averages
      
      zistk=min(nstack,nstep)
      
      do i=1,mxnstk
        
        ravval(i)=zumval(i)/zistk
        
      enddo
      
c     
c     zero accumulators during equilibration period
      
      if(lzeql.and.nstep.le.nsteql)then
        
        numacc=0
        do i=1,mxnstk
          
          sumval(i)=0.d0
          ssqval(i)=0.d0
          
        enddo
        
      endif
c     
c     close statistics file at regular intervals
      
      if(.not.newjob.and.mod(nstep,ndump).eq.0)then
        
        if(idnode.eq.0)close (nstats)
        newjob=.true.
        
      endif
      
#ifdef VAMPIR
      call VTEND(66, ierr)
#endif
      return
      end
