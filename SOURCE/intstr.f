      function intstr(word,len,lst)

c     
c***********************************************************************
c     
c     dl_poly function for extracting integers from a 
c     character string
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith may 1994.
c     
c     parameters:
c     word   - input character string
c     len    - working length of character string
c     lst    - location of space character at end of
c     integer string
c     
c     wl
c     1996/02/15 14:32:58
c     1.1.1.1
c     Exp
c     
c***********************************************************************
c     
      
      character*1 n(0:9),word(len),ksn
      logical flag,final
      
      data n/'0','1','2','3','4','5','6','7','8','9'/

      isn=1
      ksn='+'
      intstr=0
      final=.false.
      
      do lst=1,len
        
        flag=.false.

        do j=0,9
          
          if(n(j).eq.word(lst))then
            
            intstr=10*intstr+j
            flag=.true.
            
          endif
          
        enddo

        if(flag.and.ksn.eq.'-')isn=-1
        ksn=word(lst)

        if(flag)then

          final=.true.

        else

          if(final)then

            intstr=isn*intstr
            return

          endif

        endif

      enddo

      intstr=isn*intstr
      
      return
      
      end




