      function duni()
c     
c*********************************************************************
c     
c     dl_poly random number generator based on the universal
c     random number generator of marsaglia, zaman and tsang
c     (stats and prob. lett. 8 (1990) 35-39.) it must be
c     called once to initialise parameters u,c,cd,cm
c     
c     copyright daresbury laboratory 1992
c     author -  w.smith         july 1992
c     
c     wl
c     2001/05/30 12:40:03
c     1.3
c     Exp
c     
c*********************************************************************
c     
      real*8 duni
      logical new
      real*4 u(97)
      save u,c,cd,cm,uni,ir,jr,new
      data new/.true./
#ifdef VAMPIR
      call VTBEGIN(134, ierr)
#endif
      if(new)then
c     
c     initial values of i,j,k must be in range 1 to 178 (not all 1)
c     initial value of l must be in range 0 to 168.
        i=12
        j=34
        k=56
        l=78
c     
        ir=97
        jr=33
        new=.false.
        do 200 ii=1,97
          s=0.0
          t=0.5
          do 100 jj=1,24
            m=mod(mod(i*j,179)*k,179)
            i=j
            j=k
            k=m
            l=mod(53*l+1,169)
            if(mod(l*m,64).ge.32)s=s+t
            t=0.5*t
  100     continue
          u(ii)=s
  200   continue
        c =  362436.0/16777216.0
        cd= 7654321.0/16777216.0
        cm=16777213.0/16777216.0
        duni=0.d0
      else
c     
c     calculate random number
        uni=u(ir)-u(jr)
        if(uni.lt.0.0)uni=uni+1.0
        u(ir)=uni
        ir=ir-1
        if(ir.eq.0)ir=97
        jr=jr-1
        if(jr.eq.0)jr=97
        c=c-cd
        if(c.lt.0.0)c=c+cm
        uni=uni-c
        if(uni.lt.0.0)uni=uni+1.0
        duni=dble(uni)
      endif
#ifdef VAMPIR
      call VTEND(134, ierr)
#endif
      return
      end
