c*************************************************************************
c                            STEP_KDK_MIG_TP.F
c*************************************************************************
c This subroutine takes a step in helio coord.  
c Does a KICK than a DRIFT than a KICK.
c ONLY DOES TEST PARTICLES
c

      subroutine step_kdk_mig_tp(i1st,nbod,ntp,mass,j2rp2,j4rp4,
     &              xbeg,ybeg,zbeg,xend,yend,zend,
     &              xht,yht,zht,vxht,vyht,vzht,istat,dt)	

      include '../../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4  
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Internals:
c      integer i
      real*8 dth 
      real*8 axht(NTPMAX),ayht(NTPMAX),azht(NTPMAX)

      save axht,ayht,azht     ! Note this !!

c----
c...  Executable code 

      dth = 0.5d0*dt

      if(i1st.eq.0) then 
c...      Get the accelerations in helio frame.
          call getacch_tp(nbod,ntp,mass,j2rp2,j4rp4,xbeg,ybeg,zbeg,
     &                  xht,yht,zht,istat,axht,ayht,azht)
          i1st = 1    ! turn this off
      endif

c...  Apply a heliocentric kick for a half dt 
      call kickvh_tp(ntp,vxht,vyht,vzht,axht,ayht,azht,istat,dth) 

c...  Take a drift forward full step
      call drift_mig_tp(ntp,mass(1),xht,yht,zht,vxht,vyht,vzht,dt,istat)
	

c...  Get the accelerations in helio frame.
      call getacch_tp(nbod,ntp,mass,j2rp2,j4rp4,xend,yend,zend,
     &                  xht,yht,zht,istat,axht,ayht,azht)

c...  Apply a heliocentric kick for a half dt 
      call kickvh_tp(ntp,vxht,vyht,vzht,axht,ayht,azht,istat,dth) 

      return
      end   ! step_kdk_tp
c---------------------------------------------------------------------

