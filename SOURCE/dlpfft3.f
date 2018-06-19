      subroutine dlpfft3
     x  (ind,isw,ndiv1,ndiv2,ndiv3,key1,key2,key3,
     x  ww1,ww2,ww3,aaa)

c***********************************************************************
c     
c     dl-poly 3D fast fourier transform routine (in place)
c     
c     copyright daresbury laboratory 1998
c     
c     author w smith july 1998
c     
c     wl
c     2002/05/31 13:58:07
c     1.2
c     Exp
c     
c***********************************************************************
      
      implicit real*8(a-h,o-z)
      
      logical lkx,lky,lkz
      dimension key1(ndiv1),key2(ndiv2),key3(ndiv3)
      complex*16 ww1(ndiv1),ww2(ndiv2),ww3(ndiv3)
      complex*16 ttt,aaa(ndiv1,ndiv2,ndiv3)
      save nu1,nu2,nu3

      data tpi/6.283185307179586d0/

      if(ind.gt.0)then

c     check FFT array dimensions

        idm=1
        lkx=.true.
        lky=.true.
        lkz=.true.

        do i=1,30
          
          idm=2*idm

          if(idm.eq.ndiv1)then

            lkx=.false.
            nu1=i

          endif
          if(idm.eq.ndiv2)then

            lky=.false.
            nu2=i

          endif
          if(idm.eq.ndiv3)then

            lkz=.false.
            nu3=i

          endif
          
        enddo
        
        if(lkx.or.lky.or.lkz)then
          
          write(*,*)'error - FFT array not 2**N'
          stop
          
        endif
        
c     set reverse bit address arrays
        
        do kkk=1,ndiv1

          iii=0
          jjj=kkk-1

          do j=1,nu1

            jj2=jjj/2
            iii=2*(iii-jj2)+jjj
            jjj=jj2

          enddo

          key1(kkk)=iii+1

        enddo

        do kkk=1,ndiv2

          iii=0
          jjj=kkk-1

          do j=1,nu2

            jj2=jjj/2
            iii=2*(iii-jj2)+jjj
            jjj=jj2

          enddo

          key2(kkk)=iii+1

        enddo

        do kkk=1,ndiv3

          iii=0
          jjj=kkk-1

          do j=1,nu3

            jj2=jjj/2
            iii=2*(iii-jj2)+jjj
            jjj=jj2

          enddo

          key3(kkk)=iii+1

        enddo

c     initialise complex exponential factors
        
        ww1(1)=(1.d0,0.d0)

        do i=1,ndiv1/2

          arg=(tpi/dble(ndiv1))*dble(i)
          ww1(i+1)=cmplx(cos(arg),sin(arg),kind=8)
          ww1(ndiv1+1-i)=conjg(ww1(i+1))

        enddo
        
        ww2(1)=(1.d0,0.d0)

        do i=1,ndiv2/2

          arg=(tpi/dble(ndiv2))*dble(i)
          ww2(i+1)=cmplx(cos(arg),sin(arg),kind=8)
          ww2(ndiv2+1-i)=conjg(ww2(i+1))

        enddo
        
        ww3(1)=(1.d0,0.d0)

        do i=1,ndiv3/2

          arg=(tpi/dble(ndiv3))*dble(i)
          ww3(i+1)=cmplx(cos(arg),sin(arg),kind=8)
          ww3(ndiv3+1-i)=conjg(ww3(i+1))

        enddo

        return

      endif

c     take conjugate of exponentials if required
      
      if(isw.lt.0)then
        
        do i=1,ndiv1

          ww1(i)=conjg(ww1(i))

        enddo

        do i=1,ndiv2

          ww2(i)=conjg(ww2(i))

        enddo

        do i=1,ndiv3

          ww3(i)=conjg(ww3(i))

        enddo
        
      endif

c     perform fourier transform in X direction
      
      kkk=0
      num=ndiv1/2

      do l=1,nu1

        do while(kkk.lt.ndiv1)

          do i=1,num
            
            iii=key1(kkk/num+1)
            kk1=kkk+1
            k12=kk1+num
            
            do j=1,ndiv2
              
              do k=1,ndiv3
                
                ttt=aaa(k12,j,k)*ww1(iii)
                aaa(k12,j,k)=aaa(kk1,j,k)-ttt
                aaa(kk1,j,k)=aaa(kk1,j,k)+ttt
                
              enddo
              
            enddo
            
            kkk=kkk+1
            
          enddo
          
          kkk=kkk+num
          
        enddo
        
        kkk=0
        num=num/2

      enddo

c     unscramble the fft using bit address array
      
      do kkk=1,ndiv1

        iii=key1(kkk)

        if(iii.gt.kkk)then

          do j=1,ndiv2

            do k=1,ndiv3

              ttt=aaa(kkk,j,k)
              aaa(kkk,j,k)=aaa(iii,j,k)
              aaa(iii,j,k)=ttt

            enddo

          enddo

        endif

      enddo

c     perform fourier transform in Y direction
      
      kkk=0
      num=ndiv2/2

      do l=1,nu2

        do while(kkk.lt.ndiv2)

          do i=1,num
            
            iii=key2(kkk/num+1)
            kk1=kkk+1
            k12=kk1+num
            
            do j=1,ndiv1
              
              do k=1,ndiv3
                
                ttt=aaa(j,k12,k)*ww2(iii)
                aaa(j,k12,k)=aaa(j,kk1,k)-ttt
                aaa(j,kk1,k)=aaa(j,kk1,k)+ttt
                
              enddo
              
            enddo
            
            kkk=kkk+1
            
          enddo
          
          kkk=kkk+num
          
        enddo

        kkk=0
        num=num/2

      enddo

c     unscramble the fft using bit address array
      
      do kkk=1,ndiv2

        iii=key2(kkk)

        if(iii.gt.kkk)then

          do j=1,ndiv1

            do k=1,ndiv3

              ttt=aaa(j,kkk,k)
              aaa(j,kkk,k)=aaa(j,iii,k)
              aaa(j,iii,k)=ttt

            enddo

          enddo

        endif

      enddo

c     perform fourier transform in Z direction
      
      kkk=0
      num=ndiv3/2

      do l=1,nu3

        do while(kkk.lt.ndiv3)

          do i=1,num

            iii=key3(kkk/num+1)
            kk1=kkk+1
            k12=kk1+num
            
            do j=1,ndiv1
              
              do k=1,ndiv2
                
                ttt=aaa(j,k,k12)*ww3(iii)
                aaa(j,k,k12)=aaa(j,k,kk1)-ttt
                aaa(j,k,kk1)=aaa(j,k,kk1)+ttt
                
              enddo
              
            enddo
            
            kkk=kkk+1
            
          enddo
          
          kkk=kkk+num
          
        enddo

        kkk=0
        num=num/2

      enddo

c     unscramble the fft using bit address array
      
      do kkk=1,ndiv3

        iii=key3(kkk)

        if(iii.gt.kkk)then

          do j=1,ndiv1

            do k=1,ndiv2

              ttt=aaa(j,k,kkk)
              aaa(j,k,kkk)=aaa(j,k,iii)
              aaa(j,k,iii)=ttt

            enddo

          enddo

        endif

      enddo

c     restore exponentials to unconjugated values if necessary
      
      if(isw.lt.0)then
        
        do i=1,ndiv1

          ww1(i)=conjg(ww1(i))

        enddo
        
        do i=1,ndiv2

          ww2(i)=conjg(ww2(i))

        enddo
        
        do i=1,ndiv3

          ww3(i)=conjg(ww3(i))

        enddo
        
      endif
      
      return
      end
