      subroutine multipleneu
     x  (lgofr,lzeql,newlst,idnode,imcon,keyfce,multt,
     x  mxnode,natms,nneut,nstbgr,nstep,nsteql,numrdf,
     x  delr,dlrpot,engcpe,engsrp,epsq,rprim,rcut,rvdw,
     x  vircpe,virsrp,ilist,jlist,lentry,lexatm,list,
     x  lstout,lstvdw,ltype,nexatm,neulst,lstneu,lstfrz,
     x  buffer,cell,chge,fxx,fyy,fzz,flx,fly,flz,fpx,
     x  fpy,fpz,ggg,rdf,rsqdf,vvv,xdf,xxx,ydf,yyy,zdf,
     x  zzz,stress,txx,tyy,tzz,xxt,yyt,zzt,uxx,
     x  uyy,uzz)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating interatomic forces
c     using the verlet neighbour list
c     neutral group implemenation - no Ewald sum option
c     multiple timestep option
c     
c     parallel replicated data version
c     
c     fpx,fpy,fpz : forces from electrostatics fron rprim < r <= rcut
c     fxx,fyy,fzz : total force
c     
c     copyright daresbury laboratory april 1994
c     author  - t. forester april 1993
c     key:
c     
c     keyfce = odd  ------ short range potentials calculated : srfrce
c     = 0,1  ------ no electrostatics
c     = 2,3  ------ invalid
c     = 4,5  ------ distance dependent dielectric     : coul2
c     = 6,7  ------ coulombic                         : coul0
c     = 8,9  ------ invalid
c     = 10,11 ----- reaction field                    : coul3
c     
c     wl
c     2000/01/18 14:05:42
c     1.5
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical lgofr,lzeql,newlst,newplst,lgr,lchk
      
      dimension ilist(mxxdf),jlist(mxxdf),neulst(mxneut)
      dimension lentry(msatms),list(msatms,mxlist)
      dimension nexatm(msatms),lexatm(msatms,mxexcl)
      dimension lstfrz(mxatms),lstout(mxatms),lstneu(mxatms)
      dimension lstvdw(mxvdw),ltype(mxatms)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension fpx(mxatms),fpy(mxatms),fpz(mxatms)
      dimension flx(mxatms),fly(mxatms),flz(mxatms)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension vvv(mxgrid,mxvdw),ggg(mxgrid,mxvdw)
      dimension cell(9),chge(mxatms),buffer(mxbuff)
      dimension rdf(mxrdf,mxvdw)
      dimension stress(9),stresl(9),stresp(9)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension uxx(mxatms),uyy(mxatms),uzz(mxatms)
      
      save engcpl,engsrl,vircpl,virsrl,nstep0,numlsts,engcp1,vircp1,
     x  engsr1,virsr1,stresl,stresp

      data numlsts/-1/

#ifdef VAMPIR
      call VTBEGIN(18, ierr)
#endif
c     
c     error if ewald sum requested
      
      if(keyfce/2.eq.1.or.keyfce/2.eq.6) call error(idnode,250)
c     
c     set up primary and secondary neighbour lists if needed
      
      if(newlst) nstep0 = nstep
      
      newplst =(mod(nstep-nstep0,multt).eq.0)
      
      if (newplst) then        
        
        numlsts = numlsts + 1
c     
c     create list of primary and secondary neighbours
c     based on atom-atom distance
        

        call  prneulst
     x    (newlst,imcon,idnode,mxnode,nneut,rprim,
     x    list,lentry,neulst,cell,xxx,yyy,zzz,xdf,ydf,zdf)
        
      endif

      if(newplst.or.(mod(nstep-nstep0,multt).le.1)) then
c     
c     zero accumulators for secondary neighbour energies and virial
        
        engcpl = 0.d0
        vircpl = 0.d0
        engsrl = 0.d0
        virsrl = 0.d0
        
c     
c     zero secondary forces
        
        do i = 1,natms
          
          flx(i) = 0.d0
          fly(i) = 0.d0
          flz(i) = 0.d0
          
        enddo
c     
c     zero secondary stress tensor

        do i = 1,9

          stresl(i) = 0.d0

        enddo
        
      endif
c     
c     zero accumulators for total energies and virials
      
      engcpe=0.d0
      engsrp=0.d0
      vircpe=0.d0
      virsrp=0.d0
c     
c     zero force arrays
      
      do i = 1,natms
        
        fxx(i) = 0.d0
        fyy(i) = 0.d0
        fzz(i) = 0.d0
        
      enddo
c     
c     flag for accumulating rdfs
      
      lgr =.false.
      if(nstbgr.gt.0)lgr=(mod(numlsts,nstbgr).eq.0)
      lgr = (lgr.and.(newplst.and.lgofr))
      lgr = (lgr.and.((.not.lzeql).or.(nstep-nsteql.gt.0)))
c     
c     intra group vectors com
        
      do jneu = 1,nneut

         jj0 = neulst(jneu)
         jj1 = neulst(jneu+1)-1
c     
c     loop over jneu sites
            
         do j = jj0,jj1
              
            txx(j)=xxx(j)-xxx(jj0)
            tyy(j)=yyy(j)-yyy(jj0)
            tzz(j)=zzz(j)-zzz(jj0)
              
         enddo
            
      enddo
        
      call images(imcon,0,1,natms,cell,txx,tyy,tzz)
        
      do jneu = 1,nneut
         
         jj0 = neulst(jneu)
         jj1 = neulst(jneu+1)-1
c     
c     loop over jneu sites
            
         do j = jj0,jj1
              
            xxx(j)= txx(j)+xxx(jj0)
            yyy(j)= tyy(j)+yyy(jj0)
            zzz(j)= tzz(j)+zzz(jj0)
              
         enddo
c
c     centre of molecule

         uxx(jneu) = 0.d0
         uyy(jneu) = 0.d0
         uzz(jneu) = 0.d0
         anorm = 1.d0/dble(jj1-jj0+1)

         do j = jj0,jj1
              
            uxx(jneu) = uxx(jneu) + xxx(j)*anorm
            uyy(jneu) = uyy(jneu) + yyy(j)*anorm
            uzz(jneu) = uzz(jneu) + zzz(j)*anorm

         enddo
c
c     vector from site to geometric centre

         do j = jj0,jj1

            txx(j) = xxx(j) - uxx(jneu)
            tyy(j) = yyy(j) - uyy(jneu)
            tzz(j) = zzz(j) - uzz(jneu)

         enddo

      enddo
c     
c     outer loop over neutral groups
      
      lchk = .true.
      ibig = 0
      ia=0
      
      if(newplst.or.(mod(nstep-nstep0,multt).le.1)) then
c     
c     use  SECONDARY NEIGHBOURS*******************************************

c     
c     outer loop over neutral groups
        
        do ineu=idnode+1,nneut,mxnode
          
          ia=ia+1
c     
c     calculate interatomic distances
          
          isn = -1
          call neutlst
     x      (.true.,lchk,isn,imcon,idnode,mxnode,ineu,ia,ik,rcut,
     x      nexatm,lexatm,list,lentry,ilist,jlist,lstneu,lstout,
     x      lstfrz,neulst,cell,xxx,yyy,zzz,xdf,ydf,zdf,rsqdf,
     x      txx,tyy,tzz,xxt,yyt,zzt,uxx,uyy,uzz)
c
c     trap possible array bound exception 

          ibig=max(ibig,ik)
          if(ik.gt.mxxdf) ik = 0

c     
c     calculate short range force and potential terms
          
          if(mod(keyfce,2).eq.1.and.(rvdw.gt.rprim-delr)) then
            
            call srfrceneu
     x        (ik,engacc,viracc,dlrpot,rvdw,ilist,jlist,ltype,
     x        lstvdw,rsqdf,xdf,ydf,zdf,flx,fly,flz,vvv,ggg,stresl)
            
            engsrl=engsrl+engacc
            virsrl=virsrl+viracc
            
          endif
c     
c     calculate coulombic force and potential terms
          
          if (keyfce/2.eq.2) then
            
            call coul2neu
     x        (ik,engacc,viracc,epsq,ilist,jlist,
     x        chge,rsqdf,xdf,ydf,zdf,flx,fly,flz,stresl)
            
            engcpl=engcpl+engacc
            vircpl=vircpl+viracc
            
          elseif (keyfce/2.eq.3) then
            
            call coul0neu
     x        (ik,engacc,viracc,epsq,ilist,jlist,chge,
     x        rsqdf,xdf,ydf,zdf,flx,fly,flz,stresl)  
            
            engcpl=engcpl+engacc
            vircpl=vircpl+viracc

          elseif(keyfce/2.eq.4) then
            
            call error(idnode,250)

          elseif (keyfce/2.eq.5) then
            
            call coul3neu
     x        (ik,engacc,viracc,epsq,rcut,ilist,jlist,chge,
     x        rsqdf,xdf,ydf,zdf,flx,fly,flz,stresl)  
            
            engcpl=engcpl+engacc
            vircpl=vircpl+viracc

          else
            
            call error(idnode,250)
            
          endif
          
c     
c     accumulate radial distribution functions out to rcut
          
          if(lgr) then
            
            call rdf0neu(ik,rcut,ilist,jlist,ltype,lstvdw,rdf,rsqdf)
            
          endif

        enddo

      endif

c     
c     use  PRIMARY  NEIGHBOURS*******************************************
      
c     
c     outer loop over neutral groups
      
      ia=0
      
      do ineu=idnode+1,nneut,mxnode
        
        ia=ia+1
c     
c     calculate interatomic distances

        isn = 1        
        call neutlst
     x    (.true.,lchk,isn,imcon,idnode,mxnode,ineu,ia,ik,rcut,
     x    nexatm,lexatm,list,lentry,ilist,jlist,lstneu,lstout,
     x    lstfrz,neulst,cell,xxx,yyy,zzz,xdf,ydf,zdf,rsqdf,
     x    txx,tyy,tzz,xxt,yyt,zzt,uxx,uyy,uzz)
c
c     trap possible array bound exception 

        ibig=max(ibig,ik)
        if(ik.gt.mxxdf) ik = 0
c     
c     calculate short range force and potential terms
        
        if(mod(keyfce,2).eq.1) then
          
          call srfrceneu
     x      (ik,engacc,viracc,dlrpot,rvdw,ilist,jlist,ltype,
     x      lstvdw,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,vvv,ggg,stress)
          
          engsrp=engsrp+engacc
          virsrp=virsrp+viracc
          
        endif
c     
c     calculate coulombic force and potential terms
        
        if (keyfce/2.eq.2) then
          
          call coul2neu
     x      (ik,engacc,viracc,epsq,ilist,jlist,
     x      chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        elseif (keyfce/2.eq.3) then
          
          call coul0neu
     x      (ik,engacc,viracc,epsq,ilist,jlist,chge,
     x      rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress)  
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc

        elseif(keyfce/2.eq.4) then
          
          call error(idnode,250)
          
        elseif (keyfce/2.eq.5) then
          
          call coul3neu
     x      (ik,engacc,viracc,epsq,rcut,ilist,jlist,chge,
     x      rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress)  
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        endif
        
c     
c     accumulate radial distribution functions out to rcut
        
        if(lgr) then
          
          call rdf0neu(ik,rcut,ilist,jlist,ltype,lstvdw,rdf,rsqdf)
          
        endif
        
      enddo
c     
c     END of PRIMARY neighbour loop
      
c     
c     check on validity of call to neutlst

      if(mxnode.gt.1) call gstate(lchk)
      if(.not.lchk) then 
         call gimax(ibig,1,i)
         if(idnode.eq.0) write(nrite,*) 'mxxdf must be at least ',ibig
         if(idnode.eq.0) write(nrite,*) 'mxxdf is currently ',mxxdf
         call  error(idnode,479)
      endif
c     
c     counter for rdf statistics outside loop structure
      
      if(lgr) numrdf = numrdf + 1
      
c     
c     add in secondary neighbour contributions to force, energy etc
      
      if(newplst) then
        
        do i = 1,natms
          
          fpx(i) = flx(i)
          fpy(i) = fly(i)
          fpz(i) = flz(i)
          
        enddo
        
        do i = 1,9
          stresp(i) = stresl(i)
        enddo
        
        engsr1 = engsrl
        virsr1 = virsrl
        engcp1 = engcpl
        vircp1 = vircpl
        
      endif
      
      if(mod(nstep-nstep0,multt).eq.1) then
        
        do i = 1,natms
          
          flx(i) = flx(i) - fpx(i)
          fly(i) = fly(i) - fpy(i)
          flz(i) = flz(i) - fpz(i)
          
        enddo
        
        do i = 1,9
          stresl(i) = stresl(i) - stresp(i)
        enddo
        
        virsrl = virsrl - virsr1
        engsrl = engsrl - engsr1
        vircpl = vircpl - vircp1
        engcpl = engcpl - engcp1
        
      endif
      
      ann = dble(mod(nstep-nstep0,multt))
      
      do i = 1,natms
        
        fxx(i) =  fpx(i) + flx(i)*ann + fxx(i) 
        fyy(i) =  fpy(i) + fly(i)*ann + fyy(i) 
        fzz(i) =  fpz(i) + flz(i)*ann + fzz(i) 
        
      enddo

      do i = 1,9
        stress(i) = stress(i) + stresl(i)*ann + stresp(i)
      enddo
      
      engsrp = engsr1 + engsrl*ann + engsrp
      virsrp = virsr1 + virsrl*ann + virsrp
      engcpe = engcp1 + engcpl*ann + engcpe
      vircpe = vircp1 + vircpl*ann + vircpe
      
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
      
#ifdef VAMPIR
      call VTEND(18, ierr)
#endif
      return
      end

