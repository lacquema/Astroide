c*************************************************************************
c                            STEP_DKD_PL_CORR_HJS.F
c*************************************************************************
c This subroutine takes a CORRECTOR step in Generalized Jacobi coords (HJS)  
c Does a DRIFT for a*dt then a KICK for b*dt then a DRIFT for -adt.
c Does Chi(a*dt,b*dt) cf Wisdom 2006, AJ 131, 2294
c ONLY DOES MASSIVE PARTICLES
c
c             Input:
c                 a,b           ==>  the corrector coefficients
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 umat,mat      ==>  Conversion matrixes (Jacobi - Barycentric)
c                                     (2D real arrays)
c                 oloc          ==>  Link Bodies <--> Orbits (2D int array)
c                 eta,mu        ==>  Masses for centers & sats (real arrays)
c                 xj,yj,zj      ==>  initial position in Gen. Jacobi coord 
c                                    (real arrays)
c                 vxj,vyj,vzj   ==>  initial velocity in Gen. Jacobi coord 
c                                    (real arrays)
c                 dt            ==>  time step
c            Input & Output
c                 xbmid,ybmid,zbmid   ==> middle position in bary coord 
c                                    (real arrays)
c                 axbmid,aybmid,azbmid   ==>  middle accel. in bary coord 
c                                    (real arrays)
c                 xj,yj,zj      ==>  final position in Gen. Jacobi coord 
c                                       (real arrays)
c                 vxj,vyj,vzj   ==>  final velocity in Gen. Jacobi coord 
c                                       (real arrays)
c                 ir3jmid ==> Inverse radii^3 at middle point (real arrays)
c                 xjmid,yjmid,zjmid ==> Jac. positions at middle point
c                                    (real arrays)
c                 vxjmid,vyjmid,vzjmid==> Middle velocity in Gen. Jacobi coord 
c                                       (real arrays)
c
c Remarks: Adapted from step_dkd_pl_hjs.f
c Authors:  Herve Beust
c Date:   Aug 31, 2006


      subroutine step_dkd_pl_corr_hjs(a,b,i1st,nbod,oloc,umat,mat,
     &     mass,eta,mu,xj,yj,zj,vxj,vyj,vzj,xjmid,yjmid,zjmid,
     &     vxjmid,vyjmid,vzjmid,ir3jmid,xbmid,ybmid,zbmid,
     &     axbmid,aybmid,azbmid,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,i1st,oloc(NPLMAX,NPLMAX)
      real*8 a,b
      real*8 mass(nbod),dt,eta(nbod),mu(nbod)
      real*8 umat(NPLMAX,NPLMAX),mat(NPLMAX,NPLMAX)

c...  Inputs and Outputs:
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)
      real*8 xbbeg(nbod),ybbeg(nbod),zbbeg(nbod)
      real*8 axbbeg(nbod),aybbeg(nbod),azbbeg(nbod)
      real*8 ir3jbeg(nbod),ir3jend(nbod)
      real*8 xbend(nbod),ybend(nbod),zbend(nbod)
      real*8 axbend(nbod),aybend(nbod),azbend(nbod)

c...  Internals:
      integer i
      real*8 dta,dtb 
      real*8 axj(NPLMAX),ayj(NPLMAX),azj(NPLMAX)
      real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)
      real*8 xbmid(NPLMAX),ybmid(NPLMAX),zbmid(NPLMAX)
      real*8 axbmid(NPLMAX),aybmid(NPLMAX),azbmid(NPLMAX)
      real*8 xjmid(NPLMAX),yjmid(NPLMAX),zjmid(NPLMAX)
      real*8 vxjmid(NPLMAX),vyjmid(NPLMAX),vzjmid(NPLMAX)
      real*8 ir3jmid(NPLMAX)


c      save axj,ayj,azj     ! Note this !!

c----
c...  Executable code 

      dta = a*dt
      dtb = b*dt

c...   Drift in Jacobi coordinates for a half step
      call drift_hjs(nbod,mass,eta,mu,xj,yj,zj,vxj,vyj,vzj,-dta)

c...   Compute barycentric coords
      call coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
     &      xbmid,ybmid,zbmid,vxb,vyb,vzb)
      do i=1,nbod
        xjmid(i) = xj(i)
        yjmid(i) = yj(i)
        zjmid(i) = zj(i)
      end do

c...    Get the accelerations in bary frame. if frist time step
      call getacch_hjs(nbod,oloc,mass,eta,mu,xj,yj,zj,
     &     xbmid,ybmid,zbmid,ir3jmid,axbmid,aybmid,azbmid)

c...     Convert bary accels to Jacobi accels (use vel. procedure)
      call coord_vb2vg(nbod,mat,mass,axbmid,aybmid,azbmid,
     &                                               axj,ayj,azj)

c...       Apply a Jacobi kick for a full dt 
      call kickvh(nbod,vxj,vyj,vzj,axj,ayj,azj,dtb)

c...   Drift in Jacobi coordinates for a half step
      call drift_hjs(nbod,mass,eta,mu,xj,yj,zj,vxj,vyj,vzj,dta)

      return
      end   ! step_dkd_pl_corr_hjs
c---------------------------------------------------------------------

