c*************************************************************************
c                            STEP_DKD_PL_HJS_TID.F
c*************************************************************************
c This subroutine takes a step in generalized Jacobi coord (HJS)  
c Does a DRIFT, a KICK then a DRIFT.
c ONLY DOES MASSIVE PARTICLES
c                  ==> takes tides & relativity into account
c                  ==> Also computes the evolution of the spin vectors         
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 tid         ==> 0 (no tide), 1 (tide) or 2 (averaged tide)
c                 mass          ==>  mass of bodies (real array)
c                 umat,mat      ==>  Conversion matrixes (Jacobi - Barycentric)
c                                     (2D real arrays)
c                 inert       ==> Moments of inertia (real array)
c                 krpl5       ==> k(apsidal consants) * r(radius)^5
c                 Qlov        ==> Q constant's (array)
c                 oloc          ==>  Link Bodies <--> Orbits (2D int array)
c                 eta,mu        ==>  Masses for centers & sats (real arrays)
c                 xj,yj,zj      ==>  initial position in Gen. Jacobi coord 
c                                    (real arrays)
c                 vxj,vyj,vzj   ==>  initial velocity in Gen. Jacobi coord 
c                                    (real arrays)
c                 dt            ==>  time step
c            Input & Output
c                 xbbeg,ybbeg,zbbeg   ==>  initial position in bary coord 
c                                    (real arrays)
c                 axbbeg,aybbeg,azbbeg   ==>  initial accel. in bary coord 
c                                    (real arrays)
c                 xj,yj,zj      ==>  final position in Gen. Jacobi coord 
c                                       (real arrays)
c                 vxj,vyj,vzj   ==>  final velocity in Gen. Jacobi coord 
c                                       (real arrays)
c                 rotx,roty,rotz ==> Spin vectors of the bodies (real arrays)
c                 ir3jbeg,ir3jend ==> Inverse radii^3 beg & end (real arrays)
c                 xbend,ybend,zbend ==> Final position in bary coord 
c                                    (real arrays)
c                 axbend,aybend,azbend ==> Final accel in bary coord 
c                                    (real arrays)
c                 vxjh,vyjh,vzjh==>  Middle velocity in Gen. Jacobi coord 
c                                       (real arrays)
c
c Remarks: Adapted from step_kdk_pl_hjs_tid.f
c Authors:  Herve Beust
c Date:    Mar 10, 2010

      subroutine step_dkd_pl_hjs_tid(i1st,nbod,tid,oloc,umat,mat,mass,
     &      inert,krpl5,Qlov,eta,mu,xj,yj,zj,vxj,vyj,vzj,rotx,roty,rotz,
     &      xjmid,yjmid,zjmid,vxjmid,vyjmid,vzjmid,
     &      ir3jmid,xbmid,ybmid,zbmid,axbmid,aybmid,azbmid,dt)	


      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,i1st,oloc(NPLMAX,NPLMAX)
      integer tid(nbod)
      real*8 inert(nbod),krpl5(nbod),Qlov(nbod)
      real*8 mass(nbod),dt,eta(nbod),mu(nbod)
      real*8 umat(NPLMAX,NPLMAX),mat(NPLMAX,NPLMAX)

c...  Inputs and Outputs:
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)
      real*8 rotx(nbod),roty(nbod),rotz(nbod)
      real*8 vxjh(nbod),vyjh(nbod),vzjh(nbod)
      real*8 xbmid(nbod),ybmid(nbod),zbmid(nbod)
      real*8 axbmid(nbod),aybmid(nbod),azbmid(nbod)
      real*8 xjmid(nbod), yjmid(nbod), zjmid(nbod)
      real*8 vxjmid(nbod), vyjmid(nbod), vzjmid(nbod)    
      real*8 ir3jmid(nbod)


c...  Internals:
      integer i,ialpha
      real*8 dth,a,e,q 
      real*8 axj(NPLMAX),ayj(NPLMAX),azj(NPLMAX)
      real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)
c      real*8 drotx(NPLMAX),droty(NPLMAX),drotz(NPLMAX)

c      save axj,ayj,azj   ! Note this !!

c----
c...  Executable code 

       dth = 0.5d0*dt

      if (i1st.eq.0) i1st = 1
c...   Drift in Jacobi coords for a half step
      call drift_hjs(nbod,mass,eta,mu,xj,yj,zj,vxj,vyj,vzj,dth)

c...   Do the tides for a half step
      call kickvj_tid(nbod,tid,oloc,mass,inert,krpl5,Qlov,
     &     eta,mu,xj,yj,zj,vxj,vyj,vzj,dth,rotx,roty,rotz)

c...     Compute barycentric coords
      call  coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
     &      xbmid,ybmid,zbmid,vxb,vyb,vzb)
      do i=1,nbod
        xjmid(i) = xj(i)
        yjmid(i) = yj(i)
        zjmid(i) = zj(i)
        vxjmid(i) = vxj(i)
        vyjmid(i) = vyj(i)
        vzjmid(i) = vzj(i)
      end do
c...     Get the Jacobi accels
      call getaccj_hjs(nbod,oloc,mass,eta,mu,xj,yj,zj,
     &     xbmid,ybmid,zbmid,ir3jmid,axj,ayj,azj)
c...     Compute barycentric accels.
      call coord_g2b(nbod,umat,mass,vxj,vyj,vzj,axj,ayj,azj,
     &     vxb,vyb,vzb,axbmid,aybmid,azbmid)
c...       Apply a Jacobi kick for a full dt 
      call kickvh(nbod,vxj,vyj,vzj,axj,ayj,azj,dt)
c...   Do the tides for a half step
      call kickvj_tid(nbod,tid,oloc,mass,inert,krpl5,Qlov,
     &     eta,mu,xj,yj,zj,vxj,vyj,vzj,dth,rotx,roty,rotz)
c...   Drift in Jacobi coords for a half step
      call drift_hjs(nbod,mass,eta,mu,xj,yj,zj,vxj,vyj,vzj,dth)

      return
      end   ! step_dkd_pl_hjs_tid
c---------------------------------------------------------------------

