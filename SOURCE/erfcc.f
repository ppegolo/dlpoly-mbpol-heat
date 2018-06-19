      function erfcc(x)
      real*8 erfcc,x
      real*8 t,z
      z=abs(x)
      t=1.d0/(1.d0+0.5d0*z)
      erfcc=t*dexp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+t*
     $     (.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t*
     $     (1.48851587+t*(-0.82215223+t*.17087277)))))))))
      if(x.lt.0.d0) erfcc=2.d0-erfcc
      return
      end

