C**************************************************************************
C	    		        HELIOQ_GETF
C**************************************************************************
c This is the subroutine that calculates the transition F function
c
c             Input:
c              nbod      ==> number of planets  (int scalar)
c             xh,yh,zh   ==> Current heliocentric position (real arrays)
c             rtrans     ==>  solar encounter radii (real array)
c
c             Output:
c 	         ffunc   ==>  the transition function (real scalar)
c                             (real scalar)
c 	         dffunc   ==> the derivative of the transition function 
c                             (real aray)
c
c Remarks:  
c Authors:  Hal Levison
c Date:    2/2/2K
c Last revision: 

      subroutine helioq_getf(nbod,xh,yh,zh,rtrans,ffunc,dffunc)

      include '../swift.inc'
      include 'helioq.inc'

c...  Inputs Only: 
      integer nbod
      real*8 xh(nbod),yh(nbod),zh(nbod),rtrans(2)

c...  Outputs:
      real*8 ffunc,dffunc(nbod)

c.... Internals
      integer i
      real*8 rr,ffuncl(NTPMAX),r2

cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c      real*8 rr_tmp(NTPMAX)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c----
c...  Executable code 

      ffunc = 1.0d0

      do i=2,nbod

         r2 = xh(i)**2 + yh(i)**2 + zh(i)**2

         if (r2.lt.rtrans(1)) then
            ffuncl(i) = 1.0d0
            dffunc(i) = 0.0d0
         else if (r2.lt.rtrans(2)) then
            rr = (rtrans(2)-r2)/(rtrans(2)-rtrans(1))
            ffuncl(i) = 10.0d0*(rr**3) - 15.0d0*(rr**4) + 6.0d0*(rr**5)
            dffunc(i) = (-30.0d0*(rr**2) + 60.0d0*(rr**3) - 
     &           30.0d0*(rr**4))
     &           /(rtrans(2)-rtrans(1))
         else
            ffuncl(i) = 0.0d0
            dffunc(i) = 0.0d0
         endif

         ffunc = ffunc * (1.0d0 - ffuncl(i))
      enddo
      ffunc = 1.0d0 - ffunc

      if( (ffunc.eq.1.0d0) .or. (ffunc.eq.0.0d0) ) then
         do i=2,nbod
            dffunc(i) = 0.0d0
         enddo
      else
         do i=2,nbod
            if(dffunc(i).ne.0.0) then
               dffunc(i) = dffunc(i) * (1.0d0-ffunc)/(1.0d0-ffuncl(i))
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c               if(abs(dffunc(i)).gt.1.0d6) then
c                  write(*,*) 'Here #F1',i,dffunc(i),ffunc,ffuncl(i),
c     &                 rr_tmp(i),rtrans
c               endif
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            endif 
         enddo
      endif

      return
      end   ! helio_getf
c--------------------------------------------------------------------------


