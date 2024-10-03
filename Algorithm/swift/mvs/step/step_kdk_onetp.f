c*************************************************************************
c                            STEP_KDK_ONETP.F
c*************************************************************************
c This subroutine takes a step in helio coord.  
c Does a KICK than a DRIFT than a KICK.
c ONLY DOES 1 TEST PARTICLE, AND DOES NOT SAVE ACCEL'S (for parallelization)
c
c             Input:
c                 i1st           ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod           ==>  number of massive bodies (int scalar)
c                 mass           ==>  mass of bodies (real array)
c                 j2rp2,j4rp4    ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xbeg,ybeg,zbeg ==>  massive part position at beginning of dt
c                                       (real arrays)
c                 xend,yend,zend ==>  massive part position at end of dt
c                                       (real arrays)
c                 xht,yht,zht    ==>  initial part position in helio coord 
c                                      (real scalars)
c                 vxht,vyht,vzht ==>  initial velocity in helio coord 
c                                        (real scalars)
c                 axht,ayht,azht ==>  initial accel's in helio coord 
c                                        (real scalars)
c                 istat           ==>  status of the test particle
c                                      (1d integer array)
c                                      istat(1) = 0 ==> active:  = 1 not
c                                      istat(2) = -1 ==> Danby did not work
c                 dt             ==>  time step
c             Output:
c                 xht,yht,zht    ==>  final position in helio coord 
c                                       (real scalars)
c                 vxht,vyht,vzht ==>  final position in helio coord 
c                                       (real scalars)
c                 axht,ayht,azht ==>  initial accel's in helio coord 
c                                        (real scalars)      
c
c Remarks: Adapted from step_kdk_tp.f
c Authors:  Hervé Beust
c Date:    Apr 20, 2023

      subroutine step_kdk_onetp(i1st,nbod,mass,j2rp2,j4rp4,
     &            xbeg,ybeg,zbeg,xend,yend,zend,
     &            xht,yht,zht,vxht,vyht,vzht,axht,ayht,azht,istat,dt)	

      include '../../swift.inc'

c...  Inputs Only: 
      integer nbod,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4  
      real*8 xbeg(nbod),ybeg(nbod),zbeg(nbod)
      real*8 xend(nbod),yend(nbod),zend(nbod)

c...  Inputs and Outputs:
      integer istat(NSTAT)
      real*8 xht,yht,zht,vxht,vyht,vzht,axht,ayht,azht

c...  Internals:
      integer i,iflg
      real*8 dth 

c----
c...  Executable code 

      dth = 0.5d0*dt

      if(i1st.eq.0) then 
c...      Get the accelerations in helio frame.
          call getacch_onetp(nbod,mass,j2rp2,j4rp4,xbeg,ybeg,zbeg,
     &                  xht,yht,zht,istat(1),axht,ayht,azht)
          i1st = 1    ! turn this off
      endif

c...  Apply a heliocentric kick for a half dt 
      call kickvh_tp(1,vxht,vyht,vzht,axht,ayht,azht,istat(1),dth) 

c...  Take a drift forward full step
      call drift_one(mass(1),xht,yht,zht,vxht,vyht,vzht,dt,iflg)	
      if (iflg.ne.0) then
            istat(1) = 1
            istat(2) = -1
      end if
      
c...  Get the accelerations in helio frame.
      call getacch_onetp(nbod,mass,j2rp2,j4rp4,xend,yend,zend,
     &                  xht,yht,zht,istat(1),axht,ayht,azht)

c...  Apply a heliocentric kick for a half dt 
      call kickvh_tp(1,vxht,vyht,vzht,axht,ayht,azht,istat(1),dth) 

      return
      end   ! step_kdk_onetp
c---------------------------------------------------------------------

