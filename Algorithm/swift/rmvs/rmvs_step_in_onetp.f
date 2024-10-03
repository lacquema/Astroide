c*************************************************************************
c                            RMVS_STEP_IN_ONETP.F
c*************************************************************************
c This subroutine takes a step for 1 tp in inner region.
c  and don't save accels arrays (parallelization)
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
c                 istat           ==>  status of the test paricles
c                                      (1d integer array)
c                                      istat(1) = 0 ==> active:  = 1 not
c                                      istat(2) = -1 ==> Danby did not work
c                 dt             ==>  time step
c                 ipl            ==>  the planet that is currently in the 
c                                      center (integer scalar)
c           aoblxb,aoblyb,aoblzb ==> acceleration of the Sun on the central pl
c                                    at beginning of dt due to J2 and J4
c                                      (real scalars)
c           aoblxe,aoblye,aoblze ==> acceleration of the Sun on the central pl
c                                    at end of dt  due to J2 and J4
c                                         (real scalars)
c             Input & Output:
c                 xht,yht,zht    ==>  Position in helio coord 
c                                       (real scalars)
c                 vxht,vyht,vzht ==>  Velocity in helio coord 
c                                      (real scalars)
c                 axht,ayht,azht ==>  Acceleration in helio coord 
c                                       (real scalars)      
c
c Remarks: Taken from rmvs_step_in_tp.f
c Authors:  Hervé Beust
c Date:    Apr 18, 2023
c Last revision: 

      subroutine rmvs_step_in_onetp(i1st,nbod,mass,j2rp2,j4rp4,
     &      xbeg,ybeg,zbeg,xend,yend,zend,
     &      xht,yht,zht,vxht,vyht,vzht,axht,ayht,azht,
     &      istat,dt,ipl,aoblxb,aoblyb,aoblzb,aoblxe,aoblye,aoblze)

      include '../swift.inc'
      include 'rmvs.inc'

c...  Inputs Only: 
      integer nbod,i1st,ipl
      real*8 mass(nbod),dt,j2rp2,j4rp4  
      real*8 xbeg(nbod),ybeg(nbod),zbeg(nbod)
      real*8 xend(nbod),yend(nbod),zend(nbod)
      real*8 aoblxb,aoblyb,aoblzb,aoblxe,aoblye,aoblze

c...  Inputs and Outputs:
      integer istat(NSTAT)
      real*8 xht,yht,zht,vxht,vyht,vzht,axht,ayht,azht

c...  Internals:
      real*8 dth,j2rp2i,j4rp4i
      integer iflg
      
c----
c...  Executable code 

      dth = 0.5d0*dt

c...  disable J2 and J4 terms in the getacch_tp calls
      j2rp2i = 0.0d0
      j4rp4i = 0.0d0

      if(i1st.eq.0) then 
c...      Get the accelerations in helio frame.
          call getacch_onetp(nbod,mass,j2rp2i,j4rp4i,xbeg,ybeg,zbeg,
     &                  xht,yht,zht,istat(1),axht,ayht,azht)
          call rmvs_obl_acc_one(nbod,ipl,mass,j2rp2,j4rp4,xbeg,
     &         ybeg,zbeg,aoblxb,aoblyb,aoblzb,xht,yht,zht,istat,
     &         axht,ayht,azht)
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
      call getacch_onetp(nbod,mass,j2rp2i,j4rp4i,xend,yend,zend,
     &                  xht,yht,zht,istat(1),axht,ayht,azht)
      call rmvs_obl_acc_one(nbod,ipl,mass,j2rp2,j4rp4,xend,yend,
     &     zend,aoblxe,aoblye,aoblze,xht,yht,zht,istat,
     &     axht,ayht,azht)

c...  Apply a heliocentric kick for a half dt 
      call kickvh_tp(1,vxht,vyht,vzht,axht,ayht,azht,istat(1),dth) 

      return
      end   ! rmvs_step_in_onetp
c---------------------------------------------------------------------

