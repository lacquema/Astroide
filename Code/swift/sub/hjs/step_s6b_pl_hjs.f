c*************************************************************************
c                            STEP_S6B_PL_HJS.F
c*************************************************************************
c This subroutine takes a step in generalized Jacobi coord (HJS)  
c    Uses S6B 6th order symplectic algorithm (Chambers & Murison 2000)
c ONLY DOES MASSIVE PARTICLES
c
c             Input:
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
c                 xjpl,yjpl,zjpl ==>  Bodies Jacobi pos. at substeps of dt
c                                       (real 2D arrays) 1 = begin; 4 = end 
c                 vxjpl,vyjpl,vzjpl ==> Bodies Jacobi vels. at substeps of dt
c                                       (real 2D arrays) 1 = begin; 4 = end 
c                 ir3j   ==>  Bodies inverse Jac. radii^3 at substeps of dt
c                                       (real 2D arrays) 1 = begin; 4 = end 
c                 xbpl,ybpl,zbpl ==> Bodies Bary. pos. at substeps of dt
c                                       (real 2D arrays) 1 = begin; 4 = end 
c                 axbpl,aybpl,azbpl ==> Bodies Bary accels. at substeps of dt
c                                       (real 2D arrays) 1 = begin; 4 = end 
c                 xj,yj,zj      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxj,vyj,vzj   ==>  final velocity in helio coord 
c                                       (real arrays)
c
c     xbbeg,ybbeg,zbbeg   ==>  initial position in bary coord 
c                                    (real arrays)
c                 axbbeg,aybbeg,azbbeg   ==>  initial accel. in bary coord 
c                                    (real arrays)
c                 xj,yj,zj      ==>  final position in Gen. Jacobi coord 
c                                       (real arrays)
c                 vxj,vyj,vzj   ==>  final velocity in Gen. Jacobi coord 
c                                       (real arrays)
c
c     Remarks: Adapted from step_kdk_pl_hjs.f
c Authors:  Herve Beust
c Last modif:    Feb 25, 2025

      subroutine step_s6b_pl_hjs(i1st,nbod,oloc,umat,mat,mass,eta,mu,
     &     xjpl,yjpl,zjpl,vxjpl,vyjpl,vzjpl,ir3j,xbpl,ybpl,zbpl,
     &     axbpl,aybpl,azbpl,xj,yj,zj,vxj,vyj,vzj,dt)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,i1st,oloc(NPLMAX,NPLMAX)
      real*8 mass(nbod),dt,eta(nbod),mu(nbod)
      real*8 umat(NPLMAX,NPLMAX),mat(NPLMAX,NPLMAX)

c...  Inputs and Outputs:
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)
      real*8 xjpl(nbod,4),yjpl(nbod,4),zjpl(nbod,4),ir3j(nbod,4)
      real*8 vxjpl(nbod,4),vyjpl(nbod,4),vzjpl(nbod,4)
      real*8 xbpl(nbod,4),ybpl(nbod,4),zbpl(nbod,4)
      real*8 axbpl(nbod,4),aybpl(nbod,4),azbpl(nbod,4)

c...  Internals:
      integer i
      real*8 dt1,dt2,dt3,dt4
      real*8 axj(NPLMAX),ayj(NPLMAX),azj(NPLMAX)
      real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)

c      real*8, parameter :: kk = 1.2599210498948731647672106d0  ! 2^(1/3)
c      real*8, parameter :: cc = 2.0d0-kk
      real*8, parameter :: sq5 = 2.236067977499789696409174d0
      real*8, parameter :: t1 = 1.d0/12.d0
      real*8, parameter :: t2 = 0.5d0*(1.d0-1.d0/sq5)
      real*8, parameter :: t3 = 5.d0/12.d0
      real*8, parameter :: t4 = 1.d0/sq5


      save axj,ayj,azj     ! Note this !!

c----
c...  Executable code 

c      dt1 = 0.5d0*dt/cc
c      dt2 = dt/cc
c      dt3 = 0.5d0*dt*(1.0d0-kk)/cc
c      dt4 = -dt*kk/cc
      dt1 = t1*dt
      dt2 = t2*dt
      dt3 = t3*dt
      dt4 = t4*dt

c...   Store beginning positions (#1)  (for tp's after)
      xjpl(1:nbod,1) = xj(1:nbod)
      yjpl(1:nbod,1) = yj(1:nbod)
      zjpl(1:nbod,1) = zj(1:nbod)      
      vxjpl(1:nbod,1) = vxj(1:nbod)
      vyjpl(1:nbod,1) = vyj(1:nbod)
      vzjpl(1:nbod,1) = vzj(1:nbod)      

      
      if (i1st.eq.0) then
c...     Compute barycentric coords
        call  coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
     &      xbpl(:,1),ybpl(:,1),zbpl(:,1),vxb,vyb,vzb)
c...     Get the Jacobi accels if first time step
         call getaccj_hjs(nbod,oloc,mass,eta,mu,xj,yj,zj,
     &     xbpl(:,1),ybpl(:,1),zbpl(:,1),ir3j(:,1),axj,ayj,azj)
         call coord_g2b(nbod,umat,mass,vxj,vyj,vzj,axj,ayj,azj,
     &     vxb,vyb,vzb,axbpl(:,1),aybpl(:,1),azbpl(:,1))
        i1st = 1    ! turn this off
      endif
      
c...  Apply a Jacobi kick for dt1 
      call kickvh(nbod,vxj,vyj,vzj,axj,ayj,azj,dt1)

c...   First Drift in Jacobi coords for dt2
      call drift_hjs(nbod,mass,eta,mu,xj,yj,zj,vxj,vyj,vzj,dt2)

c...  After drift, compute bary. xb and vb for acceleration calculations
      call coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
     &         xbpl(:,2),ybpl(:,2),zbpl(:,2),vxb,vyb,vzb)
      call getaccj_hjs(nbod,oloc,mass,eta,mu,xj,yj,zj,
     &         xbpl(:,2),ybpl(:,2),zbpl(:,2),ir3j(:,2),axj,ayj,azj)
      call coord_g2b(nbod,umat,mass,vxj,vyj,vzj,axj,ayj,azj,
     &         vxb,vyb,vzb,axbpl(:,2),aybpl(:,2),azbpl(:,2))

c...  Store intermediate values #2 (for tp's after)
      xjpl(1:nbod,2) = xj(1:nbod)
      yjpl(1:nbod,2) = yj(1:nbod)
      zjpl(1:nbod,2) = zj(1:nbod)     
      vxjpl(1:nbod,2) = vxj(1:nbod)
      vyjpl(1:nbod,2) = vyj(1:nbod)
      vzjpl(1:nbod,2) = vzj(1:nbod)     

      
c...  Apply a second Jacobi kick for dt3 
      call kickvh(nbod,vxj,vyj,vzj,axj,ayj,azj,dt3)

c...  Second Drift in Jacobi coords for dt4
      call drift_hjs(nbod,mass,eta,mu,xj,yj,zj,vxj,vyj,vzj,dt4)

c...  After drift, compute bary. xb and vb for acceleration calculations
      call coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
     &                   xbpl(:,3),ybpl(:,3),zbpl(:,3),vxb,vyb,vzb)
      call getaccj_hjs(nbod,oloc,mass,eta,mu,xj,yj,zj,
     &     xbpl(:,3),ybpl(:,3),zbpl(:,3),ir3j(:,3),axj,ayj,azj)
      call coord_g2b(nbod,umat,mass,vxj,vyj,vzj,axj,ayj,azj,
     &     vxb,vyb,vzb,axbpl(:,3),aybpl(:,3),azbpl(:,3))

c...  Store intermediate values #3 (for tp's after)
      xjpl(1:nbod,3) = xj(1:nbod)
      yjpl(1:nbod,3) = yj(1:nbod)
      zjpl(1:nbod,3) = zj(1:nbod)     
      vxjpl(1:nbod,3) = vxj(1:nbod)
      vyjpl(1:nbod,3) = vyj(1:nbod)
      vzjpl(1:nbod,3) = vzj(1:nbod)     

c...  Apply a third Jacobi kick for dt3 
      call kickvh(nbod,vxj,vyj,vzj,axj,ayj,azj,dt3)

c...  Third Drift in Jacobi coords for dt2
      call drift_hjs(nbod,mass,eta,mu,xj,yj,zj,vxj,vyj,vzj,dt2)

c...  After drift, compute bary. xb and vb for acceleration calculations
      call coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
     &                   xbpl(:,4),ybpl(:,4),zbpl(:,4),vxb,vyb,vzb)
      call getaccj_hjs(nbod,oloc,mass,eta,mu,xj,yj,zj,
     &     xbpl(:,4),ybpl(:,4),zbpl(:,4),ir3j(:,4),axj,ayj,azj)
      call coord_g2b(nbod,umat,mass,vxj,vyj,vzj,axj,ayj,azj,
     &     vxb,vyb,vzb,axbpl(:,4),aybpl(:,4),azbpl(:,4))

c...  Store final values #4 (for tp's after)
      xjpl(1:nbod,4) = xj(1:nbod)
      yjpl(1:nbod,4) = yj(1:nbod)
      zjpl(1:nbod,4) = zj(1:nbod)     
      vxjpl(1:nbod,4) = vxj(1:nbod)
      vyjpl(1:nbod,4) = vyj(1:nbod)
      vzjpl(1:nbod,4) = vzj(1:nbod)     

c...  Apply a fourth Jacobi kick for dt1
      call kickvh(nbod,vxj,vyj,vzj,axj,ayj,azj,dt1)

      return
      end   ! step_s6b_pl_hjs
c---------------------------------------------------------------------

