c*************************************************************************
c                            STEP_DKD_TP_CORR.F
c*************************************************************************
c This subroutine takes a CORRECTOR step in helio coords  
c Does a DRIFT for a*dt then a KICK for b*dt then a DRIFT for -adt.
c Does Chi(a*dt,b*dt) cf Wisdom 2006, AJ 131, 2294
c ONLY DOES TEST PARTICLES
c
c             Input:
c                 a,b            ==>  The corrector coefficients
c                 i1st           ==>  = 0 if first step; = 1 not (int scalar)
c                                     not used here !!!
c                 nbod           ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 mass           ==>  mass of bodies (real array)
c                 j2rp2,j4rp4    ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xhm,yhm,zhm   ==>  massive part position in middle of step
c                                       (real arrays)
c                 xht,yht,zht    ==>  initial part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 dt             ==>  time step
c             Output:
c                 xht,yht,zht    ==>  final position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final position in helio coord 
c                                       (real arrays)
c
c Remarks: Adapted from step_dkd_tp.f
c Authors:  Hervé Beust
c Date:    Feb 13, 2023

      subroutine step_dkd_tp_corr(a,b,i1st,nbod,ntp,mass,j2rp2,j4rp4,
     &     xhm,yhm,zhm,xht,yht,zht,vxht,vyht,vzht,istat,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 a,b
      real*8 mass(nbod),dt,j2rp2,j4rp4
      real*8 xhm(nbod),yhm(nbod),zhm(nbod)

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Internals:
      integer istati(NTPMAX,NSTAT)
      real*8 dta,dtb
      real*8 axht(NTPMAX),ayht(NTPMAX),azht(NTPMAX)

c----
c...  Executable code 

      dta = a*dt
      dtb = b*dt
      istati = istat
      
c...  Take a drift forward -dta to begin the step (the first kick = 0).
      call drift_tp(ntp,mass(1),xht,yht,zht,vxht,vyht,vzht,-dta,istati)	

c...  Get the accelerations in helio frame.
      call getacch_tp(nbod,ntp,mass,j2rp2,j4rp4,xhm,yhm,zhm,
     &                xht,yht,zht,istati,axht,ayht,azht)

c...  Apply a heliocentric kick for dtb 
      call kickvh_tp(ntp,vxht,vyht,vzht,axht,ayht,azht,istati,dtb) 

      istati = istat
c..   Drift again in Jacobi coords for the final +dta step	  
      call drift_tp(ntp,mass(1),xht,yht,zht,vxht,vyht,vzht,dta,istati)	

      i1st = 1
      return
      end   ! step_dkd_tp_corr
c---------------------------------------------------------------------

