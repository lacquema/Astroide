c*************************************************************************
c                            STEP_S6B_TP_HJS.F
c*************************************************************************
c This subroutine takes a step in Generalized Jacobi coords (HJS)  
c    Uses S6B 6th order symplectic algorithm (Chambers & Murison 2000)
c ONLY DOES TEST PARTICLES
c
c             Input:
c                 i1st           ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod           ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of test bodies (int scalar)
c                 mass           ==>  mass of bodies (real array)
c                 matp,umatp     ==>  Conversion vectors for tp's
c                                      (2D real arrays)
c                 oloct          ==>  Link between tp's and orbits
c                                      (2D integer array)
c                 eta,mu,etatp   ==>  Masses of centers & sats (real arrays)
c                 xjbeg,yjbeg,zjbeg ==> Bodies Jac. position at beginning
c                 ir3jbeg                (real arrays)
c                 xbbeg,ybbeg,zbbeg ==> Bodies bary position at beginning
c                                       (real arrays)
c                 vxjbeg,vyjbeg,vzjbeg ==> Bodies Jac. veloc. at beginning
c                                       (real arrays)
c                 axbbeg,aybbeg,azbbeg ==> Bodies bary accs. at beginning
c                                       (real arrays)
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
c                 xjt,yjt,zjt    ==>  initial part position in Jacobi coord 
c                                      (real arrays)
c                 vxjt,vyjt,vzjt ==>  initial velocity in Jacobi coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 dt             ==>  time step
c             Output:
c                 xjt,yjt,zjt    ==>  final position in Jacobi coord 
c                                       (real arrays)
c                 vxjt,vyjt,vzjt ==>  final position in Jacobi coord 
c                                       (real arrays)
c
c Remarks: Adapted from step_kdk_tp.f
c Authors:  Herve Beust
c Date:    Mar 24, 2010

      subroutine step_s6b_tp_hjs(i1st,nbod,ntp,matp,umatp,oloct,mass,
     &              eta,mu,etatp,
     &              xjbeg,yjbeg,zjbeg,vxjbeg,vyjbeg,vzjbeg,
     &              ir3jbeg,xbbeg,ybbeg,zbbeg,axbbeg,aybbeg,azbbeg,
     &              xj1,yj1,zj1,vxj1,vyj1,vzj1,
     &              ir3j1,xb1,yb1,zb1,axb1,ayb1,azb1,
     &              xj2,yj2,zj2,vxj2,vyj2,vzj2,
     &              ir3j2,xb2,yb2,zb2,axb2,ayb2,azb2,
     &              xj3,yj3,zj3,vxj3,vyj3,vzj3,
     &              ir3j3,xb3,yb3,zb3,axb3,ayb3,azb3,
     &              xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,dt)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st,oloct(NPLMAX,NTPMAX)
      real*8 umatp(NPLMAX,NTPMAX),matp(NPLMAX,NTPMAX)
      real*8 mass(nbod),eta(nbod),mu(nbod),dt,etatp(ntp)

      real*8 xjbeg(nbod),yjbeg(nbod),zjbeg(nbod),ir3jbeg(nbod)
      real*8 vxjbeg(nbod),vyjbeg(nbod),vzjbeg(nbod)
      real*8 xbbeg(nbod),ybbeg(nbod),zbbeg(nbod)
      real*8 axbbeg(nbod),aybbeg(nbod),azbbeg(nbod)

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

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xjt(ntp),yjt(ntp),zjt(ntp)
      real*8 vxjt(ntp),vyjt(ntp),vzjt(ntp)

c...  Internals:
      integer j,k,istati(NSTAT)
      real*8 dt1,dt2,dt3,dt4
      real*8 axbttp,aybttp,azbttp,axjttp,ayjttp,azjttp 
      real*8 axjt(NTPMAX),ayjt(NTPMAX),azjt(NTPMAX)
      real*8 xbttp,zbttp,ybttp,vxbttp,vybttp,vzbttp 

c      real*8, parameter :: kk = 1.2599210498948731647672106d0  ! 2^(1/3)
c      real*8, parameter :: cc = 2.0d0-kk
      real*8, parameter :: sq5 = 2.236067977499789696409174d0
      real*8, parameter :: t1 = 1.d0/12.d0
      real*8, parameter :: t2 = 0.5d0*(1.d0-1.d0/sq5)
      real*8, parameter :: t3 = 5.d0/12.d0
      real*8, parameter :: t4 = 1.d0/sq5
      save axjt,ayjt,azjt     ! Note this !!

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

c...  loop over all tp's
      do j = 1,ntp
        if (istat(j,1).eq.0) then
          if (i1st.eq.0) then 
c...      Convert to barycentric coordinates
            call coord_g2b_tp(nbod,umatp(1,j),xjbeg,yjbeg,zjbeg,
     &          vxjbeg,vyjbeg,vzjbeg,xjt(j),yjt(j),zjt(j),vxjt(j),
     &          vyjt(j),vzjt(j),xbttp,ybttp,zbttp,vxbttp,vybttp,vzbttp)

c...      Get the acceleration in bary frame.
            call getacch_tp_hjs(nbod,mass,eta,mu,xjbeg,yjbeg,
     &           zjbeg,xbbeg,ybbeg,zbbeg,ir3jbeg,oloct(1,j),
     &           etatp(j),xbttp,ybttp,zbttp,xjt(j),yjt(j),zjt(j),
     &           axbttp,aybttp,azbttp)
c...        Convert bary accels to Jacobi accels (use vel. procedure)
            call coord_vb2vg_tp(nbod,matp(1,j),axbbeg,aybbeg,azbbeg,
     &                 axbttp,aybttp,azbttp,axjttp,ayjttp,azjttp)
          else
            axjttp = axjt(j)
            ayjttp = ayjt(j)
            azjttp = azjt(j)
          end if
c...  Apply a Jacobi kick for a dt1
          call kickvh_tp(1,vxjt(j),vyjt(j),vzjt(j),
     &                            axjttp,ayjttp,azjttp,istat(j,1),dt1) 

c...  Take a first Jacobi drift for dt2
          do k=1,NSTAT
            istati(k) = istat(j,k)
          end do
          call drift_tp(1,etatp(j),xjt(j),yjt(j),zjt(j),vxjt(j),
     &                    vyjt(j),vzjt(j),dt2,istati)	

c...  After drift, compute bary pos. and vels. 
          call coord_g2b_tp(nbod,umatp(1,j),xj1,yj1,zj1,
     &          vxj1,vyj1,vzj1,xjt(j),yjt(j),zjt(j),vxjt(j),vyjt(j),
     &          vzjt(j),xbttp,ybttp,zbttp,vxbttp,vybttp,vzbttp)

c...  Get the acceleration in bary frame.
          call getacch_tp_hjs(nbod,mass,eta,mu,xj1,yj1,
     &           zj1,xb1,yb1,zb1,ir3j1,oloct(1,j),
     &           etatp(j),xbttp,ybttp,zbttp,xjt(j),yjt(j),zjt(j),
     &           axbttp,aybttp,azbttp)
c...        Convert bary accels to Jacobi accels (use vel. procedure)
          call coord_vb2vg_tp(nbod,matp(1,j),axb1,ayb1,azb1,
     &                 axbttp,aybttp,azbttp,axjttp,ayjttp,azjttp)
c...  Apply a second Jacobi kick for dt3 
          call kickvh_tp(1,vxjt(j),vyjt(j),vzjt(j),
     &                           axjttp,ayjttp,azjttp,istat(j,1),dt3) 
c...  Second drift in Jacobi coords for dt4
          do k=1,NSTAT
            istati(k) = istat(j,k)
          end do
          call drift_tp(1,etatp(j),xjt(j),yjt(j),zjt(j),vxjt(j),
     &                    vyjt(j),vzjt(j),dt4,istati)	

c...  After drift, compute bary pos. and vels. 
          call coord_g2b_tp(nbod,umatp(1,j),xj2,yj2,zj2,
     &          vxj2,vyj2,vzj2,xjt(j),yjt(j),zjt(j),vxjt(j),vyjt(j),
     &          vzjt(j),xbttp,ybttp,zbttp,vxbttp,vybttp,vzbttp)

c...  Get the acceleration in bary frame.
          call getacch_tp_hjs(nbod,mass,eta,mu,xj2,yj2,
     &           zj2,xb2,yb2,zb2,ir3j2,oloct(1,j),
     &           etatp(j),xbttp,ybttp,zbttp,xjt(j),yjt(j),zjt(j),
     &           axbttp,aybttp,azbttp)
c...        Convert bary accels to Jacobi accels (use vel. procedure)
          call coord_vb2vg_tp(nbod,matp(1,j),axb2,ayb2,azb2,
     &                 axbttp,aybttp,azbttp,axjttp,ayjttp,azjttp)

c...  Apply a third Jacobi kick for dt3 
          call kickvh_tp(1,vxjt(j),vyjt(j),vzjt(j),
     &                           axjttp,ayjttp,azjttp,istat(j,1),dt3) 

c...  Third drift in Jacobi coords for dt2
          do k=1,NSTAT
            istati(k) = istat(j,k)
          end do
          call drift_tp(1,etatp(j),xjt(j),yjt(j),zjt(j),vxjt(j),
     &                    vyjt(j),vzjt(j),dt2,istati)	

c...  After drift, compute bary pos. and vels. 
          call coord_g2b_tp(nbod,umatp(1,j),xj3,yj3,zj3,
     &          vxj3,vyj3,vzj3,xjt(j),yjt(j),zjt(j),vxjt(j),vyjt(j),
     &          vzjt(j),xbttp,ybttp,zbttp,vxbttp,vybttp,vzbttp)

c...  Get the acceleration in bary frame.
          call getacch_tp_hjs(nbod,mass,eta,mu,xj3,yj3,
     &           zj3,xb3,yb3,zb3,ir3j3,oloct(1,j),
     &           etatp(j),xbttp,ybttp,zbttp,xjt(j),yjt(j),zjt(j),
     &           axbttp,aybttp,azbttp)
c...        Convert bary accels to Jacobi accels (use vel. procedure)
          call coord_vb2vg_tp(nbod,matp(1,j),axb3,ayb3,azb3,
     &                 axbttp,aybttp,azbttp,axjttp,ayjttp,azjttp)

c...  Apply a fourth Jacobi kick for dt1 
          call kickvh_tp(1,vxjt(j),vyjt(j),vzjt(j),
     &                           axjttp,ayjttp,azjttp,istat(j,1),dt1) 

          axjt(j) = axjttp
          ayjt(j) = ayjttp
          azjt(j) = azjttp
        end if
      end do
      i1st = 1

      return
      end   ! step_s6b_tp_hjs
c---------------------------------------------------------------------

