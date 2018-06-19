      subroutine forcesneu
     x  (lgofr,lzeql,idnode,imcon,keyfce,mxnode,natms,
     x  nneut,nstbgr,nstep,nsteql,numrdf,dlrpot,engcpe,
     x  engsrp,epsq,rcut,delr,rvdw,vircpe,virsrp,ilist,jlist,
     x  lentry,lexatm,list,lstfrz,lstneu,lstout,lstvdw,
     x  ltype,nexatm,neulst,buffer,cell,chge,fxx,fyy,
     x  fzz,ggg,rdf,rsqdf,stress,vvv,xxx,yyy,zzz,xdf,
     x  ydf,zdf,txx,tyy,tzz,xxt,yyt,zzt,uxx,uyy,uzz)

c***********************************************************************
c     
c     dl_poly subroutine for calculating interatomic forces
c     using the verlet neighbour list
c     neutral group implemenation - no Ewald sum option
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     
c     modified  - t. forester april 1993
c     key:
c     
c     keyfce = odd  ------ short range potentials calculated : srfrce
c     = 0,1  ------ no electrostatics
c     = 2,3  ------ invalid
c     = 4,5  ------ distance dependent dielectric     : coul2
c     = 6,7  ------ coulombic                         : coul0
c     = 8,9  ------ invalid
c     
c     wl
c     2000/01/18 14:05:39
c     1.5
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical lgofr,lzeql,newlst,lchk
      
      dimension ilist(mxxdf),jlist(mxxdf),neulst(mxneut)
      dimension lentry(msatms),list(msatms,mxlist)
      dimension nexatm(msatms),lexatm(msatms,mxexcl)
      dimension lstfrz(mxatms),lstout(mxatms),lstneu(mxatms)
      dimension lstvdw(mxvdw),ltype(mxatms)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension vvv(mxgrid,mxvdw),ggg(mxgrid,mxvdw)
      dimension cell(9),chge(mxatms),buffer(mxbuff)
      dimension rdf(mxrdf,mxvdw),stress(9)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension uxx(mxatms),uyy(mxatms),uzz(mxatms)
      
#ifdef VAMPIR
      call VTBEGIN(17, ierr)
#endif
c     
c     initialise force and stress arrays
      
      do i=1,natms
        fxx(i)=0.d0
        fyy(i)=0.d0
        fzz(i)=0.d0
      enddo

      do i = 1,9
        stress(i) = 0.d0
      enddo
c     
c     initialise energy and virial accumulators
      
      engcpe=0.d0
      engsrp=0.d0
      vircpe=0.d0
      virsrp=0.d0

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
      
      do ineu=idnode+1,nneut,mxnode
        
        ia=ia+1
c     
c     calculate interatomic distances

        newlst =.true.

        isn = 1
        call neutlst
     x    (newlst,lchk,isn,imcon,idnode,mxnode,ineu,ia,ik,rcut,
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
c     accumulate radial distribution functions
        
        if( ((.not.lzeql).or.(nstep.gt.nsteql)).and.(lgofr).and.
     x    mod(nstep,nstbgr).eq.0) then
          
          call rdf0neu(ik,rcut,ilist,jlist,ltype,lstvdw,rdf,rsqdf)
          
        endif

      enddo
c     
c     check on validity of call to neutlst

      if(mxnode.gt.1) call gstate(lchk)
      if(.not.lchk) then 
         call gimax(ibig,1,i)
         if(idnode.eq.0) write(nrite,*) 'mxxdf must be at least ',ibig
         if(idnode.eq.0) write(nrite,*) 'mxxdf is currently ',mxxdf
         call  error(idnode,478)
      endif

      if(keyfce/2.eq.1.or.keyfce/2.eq.6) call error(idnode,250)
c     
c     counter for rdf statistics outside loop structure
      
      if( ((.not.lzeql).or.(nstep.gt.nsteql)).and.(lgofr).and.
     x  mod(nstep,nstbgr).eq.0) numrdf = numrdf + 1
      
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
      call VTEND(17, ierr)
#endif
      return
      end
