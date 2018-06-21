      subroutine neutlst
     x     (newlst,lchk,isn,imcon,idnode,mxnode,ineu,ia,ll,rcut,
     x     nexatm,lexatm,list,lentry,ilist,jlist,ifox,jfox,
     x     lstfrz,neulst,cell,xxx,yyy,zzz,xdf,ydf,zdf,rsqdf,
     x     txx,tyy,tzz,xxt,yyt,zzt,uxx,uyy,uzz)
      
c***********************************************************************
c     
c     dlpoly routine to create pair lists for neutral group
c     implementations.
c     loops over group ineu
c     
c     replicated data version
c     
c     copyright daresbury laboratory 1994
c     author t.forester march 1994
c     
c     isn = -1 => secondary neighbours
c     isn =  1 => primary neighbours - must contain excld interactions
c     
c     wl
c     2001/05/30 12:40:13
c     1.4
c     Exp
c     
c***********************************************************************
      
#include "dl_params.inc"
      
      logical newlst,lchk,lexc
      
      dimension lstfrz(mxatms),neulst(mxneut)
      dimension nexatm(msatms),lexatm(msatms,mxexcl),lentry(msatms)
      dimension list(msatms,mxlist),ilist(mxxdf),jlist(mxxdf)
      dimension cell(9),xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension jfox(mxatms),ifox(mxatms)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension uxx(mxatms),uyy(mxatms),uzz(mxatms)
      dimension fi(3),gi(3)
CDIR$ CACHE_ALIGN fi,gi
      
#ifdef VAMPIR
      call VTBEGIN(113, ierr)
#endif
      if(newlst)then

        ibig=0
c     
c     set cutoff radius
         
         ll = 0
         nexc = 1
c     
c     number of excludes found
         
         if(isn.lt.0) then
            keyexc = nexatm(ia)+2
         else
            keyexc = 1
         endif
c     
c     do centre - centre distances
         
         lenia = lentry(ia)
         
         do j = 1,lenia
            
            jneu = abs(list(ia,j))
            xxt(j) = uxx(ineu) - uxx(jneu)
            yyt(j) = uyy(ineu) - uyy(jneu)
            zzt(j) = uzz(ineu) - uzz(jneu)
            
         enddo
         
         call images(imcon,0,1,lenia,cell,xxt,yyt,zzt)
c     
c     working intragroup vectors of central group 
c     - for periodic boundaries
         
         in0 = neulst(ineu)
         in1 = neulst(ineu+1)-1
c     
c     loop over neutral groups sites of a  
         
c     
c     loop over groups in list
         
         do jj = 1,lentry(ia)
            
            jneu = list(ia,jj)*isn
            
            if(jneu.gt.0) then
               
               do i = in0,in1
                  
                  jj0 = neulst(jneu)
                  jj1 = neulst(jneu+1)-1
                  
                  if(ineu.eq.jneu) jj0 = i+1
c     
c     loop over jneu sites
                  
                  do j = jj0,jj1
                     
c     
c     reject atoms in excluded pair list
                     
                     lexc = .false.     
                     
                     if(keyexc.lt.nexatm(ia)) then
                        
                        if(lexatm(ia,keyexc).eq.i) then
                           if(lexatm(ia,keyexc+1).eq.j) then
                              lexc = .true.
                              keyexc = keyexc+2
                           endif
                        endif   
                        
                     endif
c     
c     reject frozen atom pairs
                     
                     if (lstfrz(i).ne.0) then
                        if(lstfrz(j).ne.0) lexc = .true.
                     endif
                     
                     if(.not.lexc) then
                        
                        ll = ll+1
                        if(ll.le.mxxdf) then
                           
                           xdf(ll)= txx(i) + xxt(jj) - txx(j)
                           ydf(ll)= tyy(i) + yyt(jj) - tyy(j)
                           zdf(ll)= tzz(i) + zzt(jj) - tzz(j)
                           rsqdf(ll)= xdf(ll)**2+ydf(ll)**2+zdf(ll)**2
                           ilist(ll) = i
                           jlist(ll) = j

                        else
                           
                           lchk = .false.
                           ibig = max(ibig,ll)
                           
                        endif
                        
                     endif
                     
                  enddo
                  
               enddo
               
            endif
            
         enddo
         
      endif
      
#ifdef VAMPIR
      call VTEND(113, ierr)
#endif
      return
      end
