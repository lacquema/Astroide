c*************************************************************************
c                            STEP_KDK_PL_HJS_TMRP.F
c*************************************************************************
c This subroutine takes a step in generalized Jacobi coord (HJS)  
c Does a KICK than a DRIFT than a KICK.
c ONLY DOES MASSIVE PARTICLES, with Poincare time regularization
c         Source : Mikkola, 1997, CeMDA 67, 147
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 umat,mat      ==>  Conversion matrixes (Jacobi - Barycentric)
c                                     (2D real arrays)
c                 oloc          ==>  Link Bodies <--> Orbits (2D int array)
c                 coef          ==> The coefficients of time transformation
c                 eta,mu        ==>  Masses for centers & sats (real arrays)
c                 p0          ==> The conjugate momentum of time ( = -H)
c                 xj,yj,zj      ==>  initial position in Gen. Jacobi coord 
c                                    (real arrays)
c                 vxj,vyj,vzj   ==>  initial velocity in Gen. Jacobi coord 
c                                    (real arrays)
c                 h            ==>  dependant variable step
c            Input & Output
c                 xbbeg,ybbeg,zbbeg   ==>  initial position in bary coord 
c                                    (real arrays)
c                 axbbeg,aybbeg,azbbeg   ==>  initial accel. in bary coord 
c                                    (real arrays)
c                 xj,yj,zj      ==>  final position in Gen. Jacobi coord 
c                                       (real arrays)
c                 vxj,vyj,vzj   ==>  final velocity in Gen. Jacobi coord 
c                                       (real arrays)
c                 ir3jbeg,ir3jend ==> Inverse radii^3 beg & end (real arrays)
c                 xbend,ybend,zbend ==> Final position in bary coord 
c                                    (real arrays)
c                 axbend,aybend,azbend ==> Final accel in bary coord 
c                                    (real arrays)
c                 vxjh,vyjh,vzjh==>  Middle velocity in Gen. Jacobi coord 
c                                       (real arrays)
c                 dt            ==>  time step       
c
c Remarks: Adapted from step_kdk_pl_hjs.f
c Authors:  Herve Beust
c Date:    Sep. 28, 2006

      subroutine step_kdk_pl_hjs_tmrp(i1st,nbod,oloc,coef,umat,mat,
     &     mass,eta,mu,p0,
     &     xj,yj,zj,vxj,vyj,vzj,gbeg,xbbeg,ybbeg,zbbeg,axbbeg,aybbeg,
     &     azbbeg,ir3jbeg,vxjh,vyjh,vzjh,gend,xbend,ybend,zbend,
     &     axbend,aybend,azbend,ir3jend,h,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,i1st,oloc(NPLMAX,NPLMAX)
      real*8 mass(nbod),eta(nbod),mu(nbod),coef(nbod)
      real*8 umat(NPLMAX,NPLMAX),mat(NPLMAX,NPLMAX)
      real*8 p0,h

c...  Inputs and Outputs:
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)
      real*8 vxjh(nbod),vyjh(nbod),vzjh(nbod)
      real*8 xbbeg(nbod),ybbeg(nbod),zbbeg(nbod)
      real*8 axbbeg(nbod),aybbeg(nbod),azbbeg(nbod)
      real*8 ir3jbeg(nbod),ir3jend(nbod)
      real*8 xbend(nbod),ybend(nbod),zbend(nbod)
      real*8 axbend(nbod),aybend(nbod),azbend(nbod)
      real*8 dt,gbeg,gend

c...  Internals:
      integer i,ialpha
      real*8 a,e,inc,capom,omega,capm,q
      real*8 hh,irj(NPLMAX),eps,vj2,g,hr,pp0
      real*8 axj(NPLMAX),ayj(NPLMAX),azj(NPLMAX)
      real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)

      save axj,ayj,azj,irj,g,hr,pp0     ! Note this !!

c----
c...  Executable code 

      hh = 0.5d0*h
c      do i=2,nbod
c          call orbel_xv2el(xj(i),yj(i),zj(i),vxj(i),vyj(i),vzj(i),
c     &     eta(i)+mu(i),ialpha,a,e,inc,capom,omega,capm)
c           print*,'h=',h
c           print*,'step_kdk_pl_hjs_tmrp orb ',i,sngl(a),sngl(e),
c     &      sngl(inc*180d0/PI),sngl(capom*180d0/PI),
c     &      sngl(omega*180d0/PI),sngl(capm*180d0/PI)
c        end do



      if (i1st.eq.0) then
c        print*,'on est bien la'
c...     Compute barycentric coords
        call  coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
     &      xbbeg,ybbeg,zbbeg,vxb,vyb,vzb)
c...     Get the accelerations in both frames if frist time step
        call getaccj_hjs_tmrp(nbod,mat,coef,oloc,mass,eta,mu,p0,
     &     xj,yj,zj,vxj,vyj,vzj,xbbeg,ybbeg,zbbeg,ir3jbeg,irj,
     &     axbbeg,aybbeg,azbbeg,axj,ayj,azj,gbeg,hr)
        i1st = 1    ! turn this off
      else
        gbeg = g
      endif

c...       Apply a Jacobi kick for a half h
      call kickvh(nbod,vxj,vyj,vzj,axj,ayj,azj,hh)

c...  Re-Compute epsilon (velocities changed)
      eps = 0.d0
      do i=2,nbod
        vj2 = vxj(i)*vxj(i)+vyj(i)*vyj(i)+vzj(i)*vzj(i)
        eps = eps + eta(i)*mu(i)*(0.5d0*vj2/(eta(i)+mu(i))-irj(i))
c         call orbel_xv2aeq(xj(i),yj(i),zj(i),vxj(i),vyj(i),vzj(i),
c     &                     eta(i)+mu(i),ialpha,a,e,q)
c         eps = eps-0.5d0*mu(i)*eta(i)/a
      end do
      eps = gbeg*(eps + p0)

c...   Drift in Jacobi coords for the full step 
      call drift_multi(nbod,coef,eta,mu,eps,
     &                          xj,yj,zj,vxj,vyj,vzj,h,dt)

c...  After drift, compute bary. xb and vb for acceleration calculations
      call coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
     &                   xbend,ybend,zbend,vxb,vyb,vzb)
c...  Save Jacobi velocities at middle point (for tp's after)
      do i=1,nbod
        vxjh(i) = vxj(i)
        vyjh(i) = vyj(i)
        vzjh(i) = vzj(i)
      end do

c... Compute now accelerations in both frames
        call getaccj_hjs_tmrp(nbod,mat,coef,oloc,mass,eta,mu,p0,
     &     xj,yj,zj,vxj,vyj,vzj,xbend,ybend,zbend,ir3jend,irj,
     &     axbend,aybend,azbend,axj,ayj,azj,gend,hr)

      g = gend

c...  Re-Compute epsilon (velocities changed)
c      eps = 0.d0
c      do i=2,nbod
c        vj2 = vxj(i)*vxj(i)+vyj(i)*vyj(i)+vzj(i)*vzj(i)
c        eps = eps + eta(i)*mu(i)*(0.5d0*vj2/(eta(i)+mu(i))-irj(i))
c      end do
c      pp0 = -eps-hr

c...  Apply a Jacobi kick for a half h 
      call kickvh(nbod,vxj,vyj,vzj,axj,ayj,azj,hh)


c... Compute barycentric coordinates
c      call coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
c     &             xb,yb,zb,vxb,vyb,vzb)

      return
      end   ! step_kdk_pl_hjs_tmrp.f
c---------------------------------------------------------------------

