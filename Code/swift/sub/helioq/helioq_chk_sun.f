c*************************************************************************
c                            HELIOQ_CHK_SUN.F
c*************************************************************************
c This subroutine checks to see if anything is getting close to the Sun
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                                    (real arrays)
c                 xh,yh,zh      ==>  position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  velocity in bary coord 
c                                    (real arrays)
c                 rtrans        ==>  solar encounter radii (real array)
c                 dt            ==>  time step (real scalar)
c             Output:
c                 isun          ==>  Flag for solar encounter (int scalar)
c                                    = 1  if there are ecnounters
c                                    = 0 if not
c                                       
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/2/2K
c Last revision: 2/5/2K

      subroutine helioq_chk_sun(nbod,xh,yh,zh,vxb,vyb,vzb,rtrans,
     &     dt,isun)

      include '../swift.inc'
      include 'helioq.inc'

c...  Inputs Only: 
      integer nbod
      real*8 xh(nbod),yh(nbod),zh(nbod),rtrans(2)
      real*8 vxb(NTPMAX),vyb(NTPMAX),vzb(NTPMAX),dt

c...  Outputs:
      integer isun

c.... Internals
      integer i
      real*8 xtmp(NTPMAX),ytmp(NTPMAX),ztmp(NTPMAX),rho2
      real*8 rtransl(2)
      real*8 ffunc,dffunc(NTPMAX)

c----
c...  Executable code 

      rtransl(1) = rtrans(1)
      rtransl(2) = OVERSIZE*rtrans(2)

      call helioq_getf(nbod,xh,yh,zh,rtransl,ffunc,dffunc)

      if(ffunc.ne.0.0) then
         isun = 1
         return                  ! <========= NOTE !!!!
      endif

      do i=2,nbod
         xtmp(i) = xh(i) + dt*vxb(i)
         ytmp(i) = yh(i) + dt*vyb(i)
         ztmp(i) = zh(i) + dt*vzb(i)
      enddo

      call helioq_getf(nbod,xtmp,ytmp,ztmp,rtransl,ffunc,dffunc)

      if(ffunc.ne.0.0) then
         isun = 1
      else
         isun = 0
      endif

      return
      end   ! helioq_chk_sun
c------------------------------------------------------------------------
