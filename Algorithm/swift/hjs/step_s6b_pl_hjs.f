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
c                 xbbeg,ybbeg,zbbeg   ==>  initial position in bary coord 
c                                    (real arrays)
c                 axbbeg,aybbeg,azbbeg   ==>  initial accel. in bary coord 
c                                    (real arrays)
c                 xj,yj,zj      ==>  final position in Gen. Jacobi coord 
c                                       (real arrays)
c                 vxj,vyj,vzj   ==>  final velocity in Gen. Jacobi coord 
c                                       (real arrays)
c                 ir3jbeg      ==> Inverse radii^3 beg & end (real arrays)
c                 xj1,yj1,zj1 ==> Bodies Jac. position at substep 1          
c                 ir3j1                (real arrays)                         
c                 xb1,yb1,zb1 ==> Bodies bary position at subsetep 1         
c                                       (real arrays)                        
c                 vxj1,vyj1,vzj1 ==> Bodies Jac. veloc. at substep 1         
c                                       (real arrays)                        
c                 axb1,ayb1,azb1 ==> Bodies bary accs. at substep 1          
c                                       (real arrays)                        
c                 xj2,yj2,zj2 ==> Bodies Jac. position at substep 2          
c                 ir3j2                (real arrays)                         
c                 xb2,yb2,zb2 ==> Bodies bary position at subsetep 2         
c                                       (real arrays)                        
c                 vxj2,vyj2,vzj2 ==> Bodies Jac. veloc. at substep 2         
c                                       (real arrays)                        
c                 axb2,ayb2,azb2 ==> Bodies bary accs. at substep 2          
c                                       (real arrays)                        
c                 xj3,yj3,zj3 ==> Bodies Jac. position at substep 3          
c                 ir3j3                (real arrays)                         
c                 xb3,yb3,zb3 ==> Bodies bary position at subsetep 3         
c                                       (real arrays)                        
c                 vxj3,vyj3,vzj3 ==> Bodies Jac. veloc. at substep 3         
c                                       (real arrays)                        
c                 axb3,ayb3,azb3 ==> Bodies bary accs. at substep 3          
c                                       (real arrays)         
c
c Remarks: Adapted from step_kdk_pl_hjs.f
c Authors:  Herve Beust
c Date:    Mar 24, 2010

      subroutine step_s6b_pl_hjs(i1st,nbod,oloc,umat,mat,mass,eta,mu,
     &              ir3jbeg,xbbeg,ybbeg,zbbeg,axbbeg,aybbeg,azbbeg,
     &              xj1,yj1,zj1,vxj1,vyj1,vzj1,
     &              ir3j1,xb1,yb1,zb1,axb1,ayb1,azb1,
     &              xj2,yj2,zj2,vxj2,vyj2,vzj2,
     &              ir3j2,xb2,yb2,zb2,axb2,ayb2,azb2,
     &              xj3,yj3,zj3,vxj3,vyj3,vzj3,
     &              ir3j3,xb3,yb3,zb3,axb3,ayb3,azb3,
     &              xj,yj,zj,vxj,vyj,vzj,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,i1st,oloc(NPLMAX,NPLMAX)
      real*8 mass(nbod),dt,eta(nbod),mu(nbod)
      real*8 umat(NPLMAX,NPLMAX),mat(NPLMAX,NPLMAX)

c...  Inputs and Outputs:
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)
      real*8 vxjh(nbod),vyjh(nbod),vzjh(nbod)
      real*8 xbbeg(nbod),ybbeg(nbod),zbbeg(nbod)
      real*8 axbbeg(nbod),aybbeg(nbod),azbbeg(nbod)
      real*8 ir3jbeg(nbod)

      real*8 xj1(nbod),yj1(nbod),zj1(nbod),ir3j1(nbod)
      real*8 vxj1(nbod),vyj1(nbod),vzj1(nbod)         
      real*8 xb1(nbod),yb1(nbod),zb1(nbod)            
      real*8 axb1(nbod),ayb1(nbod),azb1(nbod)         

      real*8 xj2(nbod),yj2(nbod),zj2(nbod),ir3j2(nbod)
      real*8 vxj2(nbod),vyj2(nbod),vzj2(nbod)         
      real*8 xb2(nbod),yb2(nbod),zb2(nbod)            
      real*8 axb2(nbod),ayb2(nbod),azb2(nbod)         

      real*8 xj3(nbod),yj3(nbod),zj3(nbod),ir3j3(nbod)
      real*8 vxj3(nbod),vyj3(nbod),vzj3(nbod)         
      real*8 xb3(nbod),yb3(nbod),zb3(nbod)            
      real*8 axb3(nbod),ayb3(nbod),azb3(nbod)         

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

      if (i1st.eq.0) then
c...     Compute barycentric coords
        call  coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
     &      xbbeg,ybbeg,zbbeg,vxb,vyb,vzb)
c...     Get the Jacobi accels if first time step
         call getaccj_hjs(nbod,oloc,mass,eta,mu,xj,yj,zj,
     &     xbbeg,ybbeg,zbbeg,ir3jbeg,axj,ayj,azj)
         call coord_g2b(nbod,umat,mass,vxj,vyj,vzj,axj,ayj,azj,
     &     vxb,vyb,vzb,axbbeg,aybbeg,azbbeg)
        i1st = 1    ! turn this off
      endif
c...       Apply a Jacobi kick for dt1 
      call kickvh(nbod,vxj,vyj,vzj,axj,ayj,azj,dt1)

c...   First Drift in Jacobi coords for dt2
      call drift_hjs(nbod,mass,eta,mu,xj,yj,zj,vxj,vyj,vzj,dt2)

c...  After drift, compute bary. xb and vb for acceleration calculations
      call coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
     &                   xb1,yb1,zb1,vxb,vyb,vzb)
      call getaccj_hjs(nbod,oloc,mass,eta,mu,xj,yj,zj,
     &     xb1,yb1,zb1,ir3j1,axj,ayj,azj)
      call coord_g2b(nbod,umat,mass,vxj,vyj,vzj,axj,ayj,azj,
     &     vxb,vyb,vzb,axb1,ayb1,azb1)

c...  Store intermediate values #1 (for tp's after)
      do i=1,nbod
        xj1(i) = xj(i)
        yj1(i) = yj(i)
        zj1(i) = zj(i)
        vxj1(i) = vxj(i)
        vyj1(i) = vyj(i)
        vzj1(i) = vzj(i)
      end do

c...  Apply a second Jacobi kick for dt3 
      call kickvh(nbod,vxj,vyj,vzj,axj,ayj,azj,dt3)

c...  Second Drift in Jacobi coords for dt4
      call drift_hjs(nbod,mass,eta,mu,xj,yj,zj,vxj,vyj,vzj,dt4)

c...  After drift, compute bary. xb and vb for acceleration calculations
      call coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
     &                   xb2,yb2,zb2,vxb,vyb,vzb)
      call getaccj_hjs(nbod,oloc,mass,eta,mu,xj,yj,zj,
     &     xb2,yb2,zb2,ir3j2,axj,ayj,azj)
      call coord_g2b(nbod,umat,mass,vxj,vyj,vzj,axj,ayj,azj,
     &     vxb,vyb,vzb,axb2,ayb2,azb2)

c...  Store intermediate values #2 (for tp's after)
      do i=1,nbod
        xj2(i) = xj(i)
        yj2(i) = yj(i)
        zj2(i) = zj(i)
        vxj2(i) = vxj(i)
        vyj2(i) = vyj(i)
        vzj2(i) = vzj(i)
      end do

c...  Apply a third Jacobi kick for dt3 
      call kickvh(nbod,vxj,vyj,vzj,axj,ayj,azj,dt3)

c...  Third Drift in Jacobi coords for dt2
      call drift_hjs(nbod,mass,eta,mu,xj,yj,zj,vxj,vyj,vzj,dt2)

c...  After drift, compute bary. xb and vb for acceleration calculations
      call coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
     &                   xb3,yb3,zb3,vxb,vyb,vzb)
      call getaccj_hjs(nbod,oloc,mass,eta,mu,xj,yj,zj,
     &     xb3,yb3,zb3,ir3j3,axj,ayj,azj)
      call coord_g2b(nbod,umat,mass,vxj,vyj,vzj,axj,ayj,azj,
     &     vxb,vyb,vzb,axb3,ayb3,azb3)

c...  Store intermediate values #3 (for tp's after)
      do i=1,nbod
        xj3(i) = xj(i)
        yj3(i) = yj(i)
        zj3(i) = zj(i)
        vxj3(i) = vxj(i)
        vyj3(i) = vyj(i)
        vzj3(i) = vzj(i)
      end do

c...  Apply a fourth Jacobi kick for dt1
      call kickvh(nbod,vxj,vyj,vzj,axj,ayj,azj,dt1)

      return
      end   ! step_s6b_pl_hjs
c---------------------------------------------------------------------

