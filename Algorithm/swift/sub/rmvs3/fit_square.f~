      SUBROUTINE fit(x,y,a,b,siga,sigb,chi2)
c numerical recipes in f77 p559 

       include '../swift.inc'

      real*8 a,b,siga,sigb,chi2,y(NB_DT),x(NB_DT)
      integer i
      real*8 sigdat,ss,st2,sx,sxoss,sy,t

      sx=0.d0
      sy=0.d0
      st2=0.d0
      b=0.
      do i=1,NB_DT
         sx= sx + x(i)
         sy = sy + y(i)
      enddo 
      ss=float(NB_DT)
      sxoss=sx/ss
      do i=1,NB_DT
        t=x(i)-sxoss
        st2=st2+t*t
        b=b+t*y(i)
      enddo
      b=b/st2
      a=(sy-sx*b)/ss
      siga=sqrt((1.d0+sx*sx/(ss*st2))/ss)
      sigb=sqrt(1.d0/st2)
      chi2=0.d0
      do i=1,NB_DT
        chi2=chi2+(y(i)-a-b*x(i))**2
      enddo
      sigdat=sqrt(chi2/(NB_DT-2))
      siga=siga*sigdat
      sigb=sigb*sigdat
      
      return
      END     
