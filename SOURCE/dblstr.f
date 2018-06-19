      function dblstr(word,len,lst)
c     
c***********************************************************************
c     
c     dl_poly function for extracting double precisions from a 
c     character string. 
c     modified from dl_poly function intstr
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith may 1994.
c     modified  - t. forester april 1994
c     
c     parameters:
c     word   - input character string
c     len    - working length of character string
c     lst    - location of space character at end of
c     double precision string
c     
c     wl
c     1996/02/20 08:49:46
c     1.2
c     Exp
c     
c***********************************************************************
c     
      
      implicit real*8 (a-h,o-z)
      
      character*1 n(0:9),word(len),ksn,dot,d,e,work(256)
      logical flag,ldot,start
      
      data n/'0','1','2','3','4','5','6','7','8','9'/
      data dot/'.'/
      data d/'d'/
      data e/'e'/
      
      sn=1.d0
      ksn='+'
      ten = 10.d0
      one = 1.d0
      
      dblstr=0.d0
      iexp = 0
      idum =0
      start=.false.
      ldot = .false.
      
      do lst=1,len
        
        flag=.false.
        
        do j=0,9
          
          if(n(j).eq.word(lst))then
            
            dblstr= ten*dblstr+one*dble(j)
            flag=.true.
            start=.true.
            
          endif
          
        enddo
        
        if(dot.eq.word(lst)) then
          
          flag=.true.
          ten = 1.d0
          ldot =.true.
          start=.true.
          
        endif

        if(flag.and.ksn.eq.'-') sn=-1.d0
        
        ksn=word(lst)
        if(ksn.eq."D")ksn="d"
        if(ksn.eq."E")ksn="e"
        
        if(ldot) then
          
          one = one/10.d0
          
        endif
        
        if(start) then
          
          if(d.eq.ksn.or.e.eq.ksn) then
            
            do i = 1,len-lst
              work(i) = word(i+lst)
            enddo
            iexp = intstr(work,len-lst,idum)
            go to 100

          endif

          if(.not.flag)go to 100
          
        endif
        
      enddo
      
  100 dblstr = sn*dblstr*(10.d0**iexp)
      lst = lst +idum
      
      return
      end
