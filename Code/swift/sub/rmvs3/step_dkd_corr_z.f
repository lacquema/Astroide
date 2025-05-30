c*************************************************************************
c                            STEP_DKD_CORR_Z.F
c*************************************************************************
c This subroutine takes a CORRECTOR step in helio coordinates 
c both massive and test particles
c Does Z(a*dt,b*dt)=Chi(a*dt,b*dt).Chi(-a*dt,-b*dt) Wisdom 2006, AJ 131, 2294
c
c             Input:
c                 a,b          ==> The corrector coefficients
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                                    not used here !!!
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 xht,yht,zht    ==>  initial part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 rstat           ==>  status of the test paricles
c                                      (2d real array)
c                 dt            ==>  time step
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  final velocity in helio coord 
c                                       (real arrays)
c                 xht,yht,zht    ==>  final position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final position in helio coord 
c                                       (real arrays)
c
c
c Authors:  Hervé Beust
c Date:    Feb 13, 2023

      subroutine step_dkd_corr_z(a,b,i1st,time,nbod,ntp,mass,
     &     j2rp2,j4rp4,xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,
     &     vxht,vyht,vzht,istat,rstat,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 a,b
      real*8 mass(nbod),dt,time,j2rp2,j4rp4

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 rstat(NTPMAX,NSTATR)
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Internals
      integer i1sttp
      real*8 xhm(NPLMAX),yhm(NPLMAX),zhm(NPLMAX)

c----
c...  Executable code 

      i1sttp = i1st

c...  First step, Chi(-a*dt,-b*dt)

c...  first do the planets
      call step_dkd_pl_corr(-a,-b,i1st,nbod,mass,j2rp2,j4rp4,
     &                 xh,yh,zh,vxh,vyh,vzh,dt,xhm,yhm,zhm)

      if(ntp.ne.0) then
c...     next the test particles
         call step_dkd_tp_corr(-a,-b,i1sttp,nbod,ntp,mass,j2rp2,j4rp4,
     &     xhm,yhm,zhm,xht,yht,zht,vxht,vyht,vzht,istat,dt)	

      endif

c...  Second step, Chi(a*dt,b*dt)      
            
c...  first do the planets
      call step_dkd_pl_corr(a,b,i1st,nbod,mass,j2rp2,j4rp4,
     &                 xh,yh,zh,vxh,vyh,vzh,dt,xhm,yhm,zhm)	

c...  next the test particles
      if(ntp.ne.0) then
         call step_dkd_tp_corr(a,b,i1sttp,nbod,ntp,mass,j2rp2,j4rp4,
     &     xhm,yhm,zhm,xht,yht,zht,vxht,vyht,vzht,istat,dt)	
      endif

      return
      end   ! step_dkd_corr_z
c------------------------------------------------------------------------
