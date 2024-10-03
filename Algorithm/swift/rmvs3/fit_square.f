      SUBROUTINE fit_square(x,y,a,b,siga,sigb)
c numerical recipes in f77 p559 
c b**2 + a*x +c
       include '../swift.inc'

      real*8 a,b,siga,sigb,y(NB_DT),x(NB_DT)
      integer i
      real*8 ss,sx,sx2,sx3,sx4,sy,sxy,sx2y
      real*8 delta2,delta3,delta4,R1,R2



      sx=0.d0
      sx2=0.d0
      sx3=0.d0
      sx4=0.d0
      sy=0.d0
      sxy=0.d0
      sx2y=0.d0    
      b=0.
      do i=1,NB_DT
         sx= sx + x(i)
         sx2= sx2 + x(i)**2
         sx3= sx3 + x(i)**3
         sx4= sx4 + x(i)**4
         sy = sy + y(i)
         sxy = sxy + x(i)*y(i)
         sx2y = sx2y + x(i)**2*y(i) 
      enddo 
      ss=float(NB_DT)
      sx=sx/ss
      sx2=sx2/ss
      sx3=sx3/ss
      sx4=sx4/ss
      sy=sy/ss
      sxy=sxy/ss
      sx2y=sx2y/ss

      delta2=sx2-sx**2
      delta3=sx3-sx*sx2
      delta4=sx4-sx2**2
 
      R1=sxy-sy*sx
      R2=sx2y-sy*sx2
 
      b=(R1*delta3-R2*delta2)/(delta3**2-delta4*delta2)
      a=(R1-b*delta3)/delta2

      sigb=sqrt(NB_DT*delta2/(delta2*delta4-delta3**2))
      siga=sigb*sqrt(delta4/delta2)
     
      return
      END     
