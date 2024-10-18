c*************************************************************************
c                            STEP_KDK_MIG.F
c*************************************************************************
c This subroutine takes a step in helio coord.  
c both massive and test particles
c

      subroutine step_kdk_mig_fix(i1st,time,nbod,ntp,mass,j2rp2,j4rp4,
     &     xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,
     &     istat,rstat,dt,avar)	

      include '../../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,time,j2rp2,j4rp4

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 rstat(NTPMAX,NSTATR)
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)
      real*8 avar
  
c...  Internals
      integer i1sttp,i
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)

c----
c...  Executable code 

      i1sttp = i1st

c...  remember the current position of the planets
      do i=1,nbod
         xbeg(i) = xh(i)
         ybeg(i) = yh(i)
         zbeg(i) = zh(i)
      enddo


c...  first do the planets
      call step_kdk_mig_fix_pl(i1st,nbod,mass,j2rp2,j4rp4,
     &     xh,yh,zh,vxh,vyh,vzh,dt,avar)

      if(ntp.ne.0) then

c...     now remember these positions
         do i=1,nbod
            xend(i) = xh(i)
            yend(i) = yh(i)
            zend(i) = zh(i)
         enddo

c...     next the test particles
         call step_kdk_mig_par_tp(i1sttp,nbod,ntp,mass,j2rp2,j4rp4,
     &        xbeg,ybeg,zbeg,xend,yend,zend,
     &        xht,yht,zht,vxht,vyht,vzht,istat,dt)	

      endif

      return
      end   ! step_kdk
c------------------------------------------------------------------------

