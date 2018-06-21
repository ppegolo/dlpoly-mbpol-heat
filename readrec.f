      subroutine readrec(record,codeline,carr,iarr,rarr)
c
c     This subroutine reads a record and extracts the entries
c     from it according to the code of the record
c
      implicit real*8 (a-h,o-z)
      character*256 record
      character*(*) codeline
      character*8 carr(20)
      character*1 space
      integer iarr(20)
      real*8 rarr(20)
      data space /' '/
      
      do i=1,20
         carr(i) = ''
         iarr(i) = 0
         rarr(i) = 0.d0
      enddo

      len = 256
      i_carr = 0
      i_iarr = 0
      i_rarr = 0

c      do i=1,len
c         if (codeline(i:i).eq.space) then
c            ispace = i
c            exit
c         endif
c      enddo

      ispace = 1
      do while(codeline(ispace:ispace).ne.space)
         ispace = ispace + 1
      enddo
      
      n = ispace - 1

      do i=1,n

         call strip(record,len)

         if (codeline(i:i).eq.'c') then

c           character string
c           ----------------
            do k=1,len
               if (record(k:k).eq.space) then
                  ispace = k
                  exit
               endif
            enddo
            i_carr = i_carr + 1
            carr(i_carr) = record(1:ispace-1)
            record = record(ispace:)
            len = len - ispace + 1

         elseif (codeline(i:i).eq.'i') then

c           integer
c           -------
            i_iarr = i_iarr + 1
            iarr(i_iarr) = intstr(record,len,idum)
            record = record(idum+1:)
            len = len - idum

         elseif (codeline(i:i).eq.'r') then

c           float
c           -----
            i_rarr = i_rarr + 1
            rarr(i_rarr) = dblstr(record,len,idum)
            record = record(idum+1:)
            len = len - idum

         else

            stop 'Error in readrec: unknown entry code...'

         endif

      enddo

      return
      end
