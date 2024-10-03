      program essai

      implicit none

      real*8 a,b,c,m,s
      integer*4 i

      m = 10.0d0
      s = 2.0d0
      do i=1,30
        call rgauss(m,s,a)
        print*,a
      end do
      end
