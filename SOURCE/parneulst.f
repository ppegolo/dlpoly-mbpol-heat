      subroutine parneulst
     x     (newlst,lneut,lms,nneut,idnode,mxnode,imcon,rcut,delr,
     x     lentry,list,lstfrz,neulst,cell,xxx,yyy,zzz,xdf,ydf,zdf)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for constructing the verlet neighbour
c     list based on the brode-ahlrichs decomposition 
c     frozen atoms taken into account
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - t.forester   april   1994
c     
c     wl
c     2000/01/18 14:05:51
c     1.5
c     $Sate: Exp $
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical lchk,newlst,lneut,safe,lfrzi,ldo,lms(mxneut)
      
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf)
      dimension lentry(msatms),list(msatms,mxlist)
      dimension lstfrz(mxatms)
      dimension cell(9)
      dimension neulst(mxneut)
      dimension fi(3)
CDIR$ CACHE_ALIGN fi
#ifdef VAMPIR
      call VTBEGIN(12, ierr)
#endif
      if(newlst.and.lneut)then
c     
c     set control variables
         
         safe = .true. 
         lchk=  .true.
         mpm2=(nneut+2)/2
         npm2=(nneut+1)/2 
c     
c     set cutoff radius
         
         rcl1 =(rcut+delr)
         rclim=(rcut+delr)**2
         ibig = 0
         ill = 0
         
c     
c     construct pair force neighbour list: neutral groups
         
         do i=1,msatms
            
            lentry(i)=0
            
         enddo
         
c     
c     outer loop over groups
         
         ia=0
         
         do im=idnode+1,nneut,mxnode
            
            ia = ia+1 
            if(im.ge.mpm2) mpm2 = npm2
            
            lms(1) = .false.
            do j = 2,mpm2
               lms(j) = .true.
            enddo

            jmlast =  min(nneut,im+mpm2-1)
            jmwrap = max(0,im+mpm2-1-nneut)
c     
c     loop over atomic pairs
            
            nuei1= neulst(im)
            nuei2= neulst(im+1)-1
            
            do i = nuei1,nuei2
               
               fi(1) = xxx(i)
               fi(2) = yyy(i)
               fi(3) = zzz(i)
               lfrzi = (lstfrz(i).eq.0)
               
               ii = 0
               jm1 = 1
               do jm = im+1,jmlast
                  jm1 = jm1+1
                  if(lms(jm1)) then
                     
                     jj0 = neulst(jm)
                     jj2 = neulst(jm+1)-1
                     
                     do j = jj0,jj2
                        
                        ii=ii+1
                        if(ii.le.mxxdf) then
                           
                           xdf(ii)=fi(1)-xxx(j)
                           ydf(ii)=fi(2)-yyy(j)
                           zdf(ii)=fi(3)-zzz(j)
                           
                        else
                           
                           ibig = max(ibig,ii)
                           safe=.false.
                           
                        endif
                        
                     enddo
                  endif
                  
               enddo
               
               do jm = 1,jmwrap
                  jm1=jm1+1
                  if(lms(jm1)) then
                     
                     jj0 = neulst(jm)
                     jj2 = neulst(jm+1)-1
                     
                     do j = jj0,jj2
                        
                        ii=ii+1
                        if(ii.le.mxxdf) then
                           
                           xdf(ii)=fi(1)-xxx(j)
                           ydf(ii)=fi(2)-yyy(j)
                           zdf(ii)=fi(3)-zzz(j)
                           
                        else
                           
                           safe=.false.
                           ibig = max(ibig,ii)
                           
                        endif
                        
                     enddo
                     
                  endif
                  
               enddo
               
c     
c     apply minimum image convention
               
               ii1 = min(ii,mxxdf)
               call images(imcon,0,1,ii1,cell,xdf,ydf,zdf)
c     
c     search for those in cutoff
               
               ii=0
               jm1 = 1
               do jm = im+1,jmlast
                  jm1 = jm1+1
                  if(lms(jm1)) then
                     
                     jj0 = neulst(jm)
                     jj2 = neulst(jm+1)-1
                     
                     do j = jj0,jj2
                        
                        ii=ii+1
                        if(ii.le.mxxdf) then

                           if(lms(jm1)) then
                           
                              if(lfrzi) then
                                 if(abs(zdf(ii)).lt.rcl1) then
                                 if(abs(ydf(ii)).lt.rcl1) then
                                 if(abs(xdf(ii)).lt.rcl1) then
                                 rrr = xdf(ii)**2+ydf(ii)**2+zdf(ii)**2
                                 if(rrr.lt.rclim) lms(jm1)=.false.
                                 endif
                                 endif
                                 endif

                              elseif(lstfrz(j).eq.0) then 
                                 if(abs(zdf(ii)).lt.rcl1) then
                                 if(abs(ydf(ii)).lt.rcl1) then
                                 if(abs(xdf(ii)).lt.rcl1) then
                                 rrr = xdf(ii)**2+ydf(ii)**2+zdf(ii)**2
                                 if(rrr.lt.rclim) lms(jm1)=.false.
                                 endif
                                 endif
                                 endif
                              
                              endif
                           
                           endif

                        endif
                        
                     enddo


                  endif
                  
               enddo
               
               do jm = 1,jmwrap
                  jm1=jm1+1
                  if(lms(jm1)) then
                     
                     jj0 = neulst(jm)
                     jj2 = neulst(jm+1)-1
                     
                     do j = jj0,jj2
                        
                        ii=ii+1
                        if(ii.le.mxxdf) then

                           if(lms(jm1)) then

                              if(lfrzi) then
                                 if(abs(zdf(ii)).lt.rcl1) then
                                 if(abs(ydf(ii)).lt.rcl1) then
                                 if(abs(xdf(ii)).lt.rcl1) then
                                 rrr = xdf(ii)**2+ydf(ii)**2+zdf(ii)**2
                                 if(rrr.lt.rclim) lms(jm1)=.false.
                                 endif
                                 endif
                                 endif
                              
                              elseif(lstfrz(j).eq.0) then
                                 if(abs(zdf(ii)).lt.rcl1) then
                                 if(abs(ydf(ii)).lt.rcl1) then
                                 if(abs(xdf(ii)).lt.rcl1) then
                                 rrr = xdf(ii)**2+ydf(ii)**2+zdf(ii)**2
                                 if(rrr.lt.rclim) lms(jm1)=.false.
                                 endif
                                 endif
                                 endif
                              
                              endif
                           
                           endif

                        endif
                        
                     enddo
                     
                  endif
                  
               enddo
               
            enddo
c     
c     compile neighbour list for ia
c     with running check of neighbour list array capacity
            
            jm1 = 0
            do jm = im,jmlast
               jm1 = jm1+1
               if(.not.lms(jm1)) then
                  
                  lentry(ia)=lentry(ia)+1
                  if(lentry(ia).le.mxlist)then
                     
                     list(ia,lentry(ia))=jm
                     
                  else
                     
                     ill =max(ill,lentry(ia))
                     lchk=.false.
                     
                  endif
                  
               endif
               
            enddo
            
            do jm = 1,jmwrap
               jm1 = jm1+1
               if(.not.lms(jm1)) then
                  
                  lentry(ia)=lentry(ia)+1
                  if(lentry(ia).le.mxlist)then
                     
                     list(ia,lentry(ia))=jm
                     
                  else
                     
                     ill =max(ill,lentry(ia))
                     lchk=.false.
                     
                  endif
                  
               endif
               
            enddo
            
         enddo
         
c     
c     terminate job if neighbour list array exceeded

         if(mxnode.gt.1) call gstate(lchk)
         if(.not.lchk) then

           call gimax(ill,1,idum)
           if(idnode.eq.0) then
             write(nrite,*) ' mxlist must be at least  ',ill
             write(nrite,*) ' mxlist is currently ',mxlist
           endif
           call error(idnode,108)

         endif   
c
c     terminate job if work arrays exceeded

         if(mxnode.gt.1) call gstate(safe)
         if(.not.safe) then
            call gimax(ibig,1,idum)
            if(idnode.eq.0) then
              write(nrite,*)'mxxdf must be at least ',ibig
              write(nrite,*)'mxxdf is currently ',mxxdf
            endif
            call  error(idnode,476)
         endif
         
      endif
#ifdef VAMPIR
      call VTEND(12, ierr)
#endif
      return
      end
