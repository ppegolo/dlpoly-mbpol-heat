      subroutine shellsort(n,list)

c***********************************************************************
c     
c     dlpoly shell sort routine. 
c     Sorts an array of integers into ascending order
c     
c     copyright daresbury laboratory 1993
c     author - t.forester   november 1993
c     
c     wl
c     1996/02/15 14:33:18
c     1.1.1.1
c     Exp
c     
c***********************************************************************

      implicit real*8(a-h,o-z)

      dimension list(*)

c     
c     set up sort

      if(n.gt.1) then

c     number of lists
        nl = n/2

c     iterate shell sort

   10   do nn = 1,nl
c     
c     begin insertion sort on nnth list
          
          do i = nn+nl,n,nl

            imax = list(i)
            ix = i
c     
c     find location for insertion
            
            j = i
  100       j = j-nl

            if(j.lt.1) goto 110

            if (list(j).gt.imax) then
              
              ix = j
              
            else
              
              j =1

            endif
            
            goto 100
  110       continue
            
c     
c     insert in index array

            do j = i,ix+nl,-nl
              
              list(j) = list(j-nl)
              
            enddo
            
            list(ix) = imax

          enddo

        enddo
        
        nl = nl/2
        if(nl.gt.0) goto 10
        
      endif

      return
      end
