      subroutine multiple_nsq
     x  (lnsq,lgofr,lzeql,newlst,idnode,imcon,keyfce,
     x  multt,mxnode,natms,nstep,nstbgr,nsteql,numrdf,
     x  delr,dlrpot,engcpe,engsrp,engcp3,epsq,rcut,
     x  rprim,rvdw,vircpe,virsrp,vircp3,ilist,lentry,
     x  list,lstvdw,ltpvdw,ltype,buffer,cell,chge,flx,fly,
     x  flz,fxx,fyy,fzz,fpx,fpy,fpz,ggg,rdf,rsqdf,vvv,xdf,
     x  xxx,ydf,yyy,zdf,zzz,stress)
      
c***********************************************************************
c     
c     DL_POLY subroutine for multiple time step algorithm 
c     to be used with all-pairs option
c     
c     flx,fly,flz : forces from electrostatics from r > rcut
c     fpx,fpy,fpz : forces from electrostatics fron rprim < r <= rcut
c     fxx,fyy,fzz : total force
c     
c     copyright daresbury laboratory 1993
c     
c     author  t. forester,  may 1993
c     
c     keyfce = odd  ------ short range potentials calculated : srfrce
c     = 0,1  ------ no electrostatics
c     Ewald sum --- not used
c     = 4,5  ------ Distance dependent dielectric     : coul2
c     = 6,7  ------ coulombic                         : coul0
c     truncated and shifted coulombic -- not used
c     reaction field - not used
c     
c     wl
c     2000/01/18 14:05:42
c     1.4
c     Exp
c     
c****************************************************************************
      
#include "dl_params.inc"
      logical newplst,newlst,lgofr,lzeql,lgr,lnsq
      
      dimension lentry(msatms),list(msatms,mxlist),ilist(mxxdf)
      dimension ltype(mxatms),ltpvdw(mxvdw)
      dimension lstvdw(mxvdw)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms),cell(9)
      dimension chge(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension flx(mxatms),fly(mxatms),flz(mxatms)
      dimension fpx(mxatms),fpy(mxatms),fpz(mxatms)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension vvv(mxgrid,mxvdw),ggg(mxgrid,mxvdw)
      dimension buffer(mxbuff)
      dimension rdf(mxrdf,mxvdw)
      dimension stress(9),stresp(9)
      
      save engsr2,virsr2,engcp2,vircp2,nstep0,numlsts,stresp
      
      data numlsts/-1/

#ifdef VAMPIR
      call VTBEGIN(14, ierr)
#endif
      if(lnsq) then

c     
c     set extended cutoff for electrostatics - secondary shell

        rcut1 = rcut+delr
c     
c     set up primary and secondary neighbour lists if needed
        
        if(newlst) nstep0 = nstep
        
        newplst =((newlst).or.((mod(nstep-nstep0,multt).eq.0)))
        
        if (newplst) then        
          
          numlsts = numlsts + 1
c     
c     create list of primary and secondary neighbours
          
          call  primlst
     x      (idnode,mxnode,natms,imcon,rprim,lentry,list,
     x      cell,xxx,yyy,zzz,xdf,ydf,zdf)
          
c     
c     zero accumulators for secondary neighbour energies and virial
          
          engcp2 = 0.d0
          vircp2 = 0.d0
          engsr2 = 0.d0
          virsr2 = 0.d0
c     
c     zero secondary forces
          
          do i = 1,natms
            
            fpx(i) = 0.d0
            fpy(i) = 0.d0
            fpz(i) = 0.d0
            
          enddo
c     
c     zero secondary stress tensor
          
          do i = 1,9
            stresp(i) = 0.d0
          enddo
          
        endif
c     
c     zero accumulators for total energies and virials
        
        engcpe=0.d0
        engsrp=0.d0
        vircpe=0.d0
        virsrp=0.d0
        
c     
c     flag for accumulating rdfs
        
        lgr = (lgofr.and.(.not.lzeql.or.(nstep-nsteql.gt.0)))
        lgr = (lgr.and.newplst.and.(mod(numlsts,nstbgr).eq.0))
        
        if(newplst) then
c     
c     use  SECONDARY NEIGHBOURS*******************************************
          
c     
c     outer loop over atoms  **** OUTER LOOP STARTS HERE *****
          
          ii=0
          
          do i=idnode+1,natms,mxnode
            
            ii=ii+1
c     
c     calculate interatomic vectors
            
            ik = 0
            do k=1,lentry(ii)
              
              j=list(ii,k)
              
              if(j.gt.0) then
                
                ik = ik + 1
                ilist(ik) = j
                xdf(ik)=xxx(i)-xxx(j)
                ydf(ik)=yyy(i)-yyy(j)
                zdf(ik)=zzz(i)-zzz(j)
                
              endif
              
            enddo
c     
c     periodic boundary condition only for interactions greater
            
            call images(imcon,0,1,ik,cell,xdf,ydf,zdf)
c     
c     square of interatomic distances

            do k = 1,ik

              rsqdf(k) = xdf(k)**2+ydf(k)**2+zdf(k)**2
              
            enddo
c     
c     short range forces for secondary shell
            
            if((mod(keyfce,2).eq.1).and.(rvdw.gt.rprim-delr)) then
              
              call srfrce(i,ik,engacc,viracc,rvdw,dlrpot,ilist,
     x          ltype,lstvdw,ltpvdw,rsqdf,xdf,ydf,zdf,fpx,fpy,fpz,
     x          vvv,ggg,stresp)
              
              engsr2=engsr2+engacc
              virsr2=virsr2+viracc
              
            endif
c     
c     calculate coulombic force and potential terms
            
            if (keyfce/2.eq.1.or.keyfce/2.eq.6) then

              call error(idnode,424)
              
            elseif (keyfce/2.eq.2) then
c     
c     distance dependent dielectric
              
              call coul2(i,ik,engacc,viracc,rcut1,epsq,ilist,
     x          chge,rsqdf,xdf,ydf,zdf,fpx,fpy,fpz,stresp)
              
              engcp2=engcp2+engacc
              vircp2=vircp2+viracc
              
            elseif (keyfce/2.eq.3) then
c     
c     coulombic potential
              
              call coul0(i,ik,engacc,viracc,rcut1,epsq,
     x          ilist,chge,rsqdf,xdf,ydf,zdf,fpx,fpy,fpz,stresp)
              
              engcp2=engcp2+engacc
              vircp2=vircp2+viracc
              
            elseif(keyfce/2.eq.4) then
c     
c     truncated shifted coulombic potential
              
              call error(idnode,424)
              
            endif
            
c     
c     accumulate radial distribution functions : out to rcut
            
            if (lgr) call rdf0(i,ik,rcut,ilist,ltype,lstvdw,rdf,rsqdf)
            
          enddo
          
        endif
c     
c     use PRIMARY NEIGBOURS 
        
c     
c     outer loop over atoms  **** OUTER LOOP STARTS HERE *****
        
        ii=0
        
        do i=idnode+1,natms,mxnode
          
          ii=ii+1
          
c     
c     calculate interatomic distances
          
          ik = 0
          
          do k=1,lentry(ii)
            
            j= -list(ii,k)
            
            if(j.gt.0) then
              
              ik = ik + 1
              ilist(ik) = j
              xdf(ik)=xxx(i)-xxx(j)
              ydf(ik)=yyy(i)-yyy(j)
              zdf(ik)=zzz(i)-zzz(j)
              
            endif
            
          enddo
          
c     
c     periodic boundary conditions
          
          call images(imcon,0,1,ik,cell,xdf,ydf,zdf)

c     
c     square of interatomic distances

          do k = 1,ik

            rsqdf(k) = xdf(k)**2 + ydf(k)**2 + zdf(k)**2

          enddo
c     
c     accumulate radial distribution functions : out to rcut
          
          if (lgr) call rdf0(i,ik,rcut,ilist,ltype,lstvdw,rdf,rsqdf)
c     
c     calculate short range force and potential terms
          
          if(mod(keyfce,2).eq.1) then
            
            call srfrce(i,ik,engacc,viracc,rvdw,dlrpot,ilist,
     x        ltype,lstvdw,ltpvdw,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,
     x        vvv,ggg,stress)
            
            engsrp=engsrp+engacc
            virsrp=virsrp+viracc
            
          endif
c     
c     calculate coulombic force and potential terms
          
          if (keyfce/2.eq.1.or.keyfce/2.eq.6) then
            
            call error(idnode,424)

          elseif (keyfce/2.eq.2) then
c     
c     distance dependent dielectric
            
            call coul2(i,ik,engacc,viracc,rcut,epsq,
     x        ilist,chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress)
            
            engcpe=engcpe+engacc
            vircpe=vircpe+viracc
            
          elseif (keyfce/2.eq.3) then
c     
c     coulombic potential
            
            call coul0(i,ik,engacc,viracc,rcut,epsq,
     x        ilist,chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress)
            
            engcpe=engcpe+engacc
            vircpe=vircpe+viracc
            
          elseif(keyfce/2.eq.4) then

            call error(idnode,424)
            
          endif
          
        enddo
c     
c     END of PRIMARY neighbour loop
        
c     
c     counter for rdf statistics outside loop structure
        
        if(lgr) numrdf = numrdf+1
        
c     
c     add in secondary and tertiary neighbour contributions to force, energy etc
        
        do i = 1,natms
          
          fxx(i) = fxx(i) + fpx(i) + flx(i)
          fyy(i) = fyy(i) + fpy(i) + fly(i)
          fzz(i) = fzz(i) + fpz(i) + flz(i)
          
        enddo

        do i = 1,9
          stress(i) = stress(i)+stresp(i)
        enddo
        
        engsrp = engsrp + engsr2
        virsrp = virsrp + virsr2

        engcpe = engcpe + engcp2 + engcp3 
        vircpe = vircpe + vircp2 + vircp3

c     
c     sum up contributions to short range and coulombic potential
        
        if(mxnode.gt.1) then 
          
          buffer(5)=engsrp
          buffer(6)=virsrp
          buffer(7)=engcpe
          buffer(8)=vircpe
          
          call gdsum(buffer(5),4,buffer(1))
          
          engsrp=buffer(5)
          virsrp=buffer(6)
          engcpe=buffer(7)
          vircpe=buffer(8)
          
        endif

      endif

#ifdef VAMPIR
      call VTEND(14, ierr)
#endif
      return
      end
