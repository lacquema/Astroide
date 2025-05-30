c*************************************************************************
c                            STEP_KDK_HJS_TID.F
c*************************************************************************
c This subroutine takes a step in generalized Jacobi coordinates (HJS)  
c both massive and test particles
c
c                  ==> takes tides & relativity into account
c                  ==> Also computes the evolution of the spin vectors         
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 tid         ==> 0 (no tide), 1 (tide) or 2 (averaged tide)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 oloc          ==>  Link bodies - orbits (2D int array)
c                 oloct         ==>  Link orbits - tp's (2D int array)
c                 mass          ==>  mass of bodies (real array)
c                 inert       ==> Moments of inertia (real array)
c                 krpl5       ==> k(apsidal consants) * r(radius)^5
c                 Qlov        ==> Q constant's (array)
c                 mat,umat      ==>  Conversion matrixes for bodies
c                                     (2D real arrays)
c                 matp,umatp    ==>  Conversion vectors for tp's
c                                     (2D real arrays)
c                 eta,mu        ==>  Masses for center & satellites for bods.
c                                     (real arrays)
c                 etatp         ==>  Masses for centers for tp's (real array)
c                 xj,yj,zj      ==>  initial position in jacobi coord 
c                                    (real arrays)
c                 vxj,vyj,vzj   ==>  initial velocity in jacobi coord 
c                                    (real arrays)
c                 xjt,yjt,zjt    ==>  initial part position in jacobi coord 
c                                      (real arrays)
c                 vxjt,vyjt,vzjt ==>  initial velocity in jacobi coord 
c                                        (real arrays)
c                 rotx,roty,rotz ==> Spin vectors of the bodies (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 rstat           ==>  status of the test paricles
c                                      (2d real array)
c                 dt            ==>  time step
c             Output:
c                 xj,yj,zj      ==>  final position in jacobi coord 
c                                       (real arrays)
c                 vxj,vyj,vzj   ==>  final velocity in jacobi coord 
c                                       (real arrays)
c                 xjt,yjt,zjt    ==>  final position in jacobi coord 
c                                       (real arrays)
c                 vxjt,vyjt,vzjt ==>  final position in jacobi coord 
c                                       (real arrays)
c
c
c Remarks: Adopted from step_kdk_hjs.f
c Authors:  Herve Beust 
c Date:    Jun 6, 2008, revised Jun 19

      subroutine step_kdk_hjs_tid(i1st,time,nbod,tid,ntp,oloc,oloct,
     &     mat,umat,matp,umatp,mass,inert,krpl5,Qlov,eta,mu,etatp,
     &     xj,yj,zj,vxj,vyj,vzj,rotx,roty,rotz,xjt,yjt,zjt,vxjt,vyjt,
     &     vzjt,istat,rstat,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,time
      real*8 mu(nbod),eta(nbod),etatp(ntp)
      real*8 umat(NPLMAX,NPLMAX),mat(NPLMAX,NPLMAX)
      real*8 umatp(NPLMAX,NTPMAX),matp(NPLMAX,NTPMAX)
      real*8 krpl5(nbod),Qlov(nbod),inert(nbod)
      integer tid(nbod)
      integer oloc(NPLMAX,NPLMAX),oloct(NPLMAX,NTPMAX)

c...  Inputs and Outputs
      integer istat(NTPMAX,NSTAT)
      real*8 rstat(NTPMAX,NSTATR)
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)
      real*8 xjt(ntp),yjt(ntp),zjt(ntp)
      real*8 vxjt(ntp),vyjt(ntp),vzjt(ntp)
      real*8 rotx(nbod),roty(nbod),rotz(nbod)

c...  Internals
      integer i1sttp,i
      real*8 xbbeg(NPLMAX),ybbeg(NPLMAX),zbbeg(NPLMAX)
      real*8 axbbeg(NPLMAX),aybbeg(NPLMAX),azbbeg(NPLMAX)
      real*8 xjbeg(NPLMAX),yjbeg(NPLMAX),zjbeg(NPLMAX)
      real*8 xjend(NPLMAX),yjend(NPLMAX),zjend(NPLMAX)
      real*8 vxjbeg(NPLMAX),vyjbeg(NPLMAX),vzjbeg(NPLMAX)
      real*8 ir3jbeg(NPLMAX),ir3jend(NPLMAX)
      real*8 vxjh(NPLMAX),vyjh(NPLMAX),vzjh(NPLMAX)
      real*8 xbend(NPLMAX),ybend(NPLMAX),zbend(NPLMAX)
      real*8 axbend(NPLMAX),aybend(NPLMAX),azbend(NPLMAX)

      save xbbeg,ybbeg,zbbeg,axbbeg,aybbeg,azbbeg,ir3jbeg ! Note this !
c----
c...  Executable code 

      i1sttp = i1st

c...  remember the current position & velocities of the massive bodies
      do i=1,nbod
         xjbeg(i) = xj(i)
         yjbeg(i) = yj(i)
         zjbeg(i) = zj(i)
         vxjbeg(i) = vxj(i)
         vyjbeg(i) = vyj(i)
         vzjbeg(i) = vzj(i)
      enddo

c...  first do the planets

      call step_kdk_pl_hjs_tid(i1st,nbod,tid,oloc,umat,mat,mass,
     &      inert,krpl5,Qlov,eta,mu,xj,yj,zj,vxj,vyj,vzj,rotx,roty,rotz,
     &      xbbeg,ybbeg,zbbeg,axbbeg,aybbeg,azbbeg,ir3jbeg,vxjh,vyjh,
     &      vzjh,xbend,ybend,zbend,axbend,aybend,azbend,ir3jend,dt)	

      if(ntp.ne.0) then

c...     now remember these positions
         do i=1,nbod
            xjend(i) = xj(i)
            yjend(i) = yj(i)
            zjend(i) = zj(i)
         enddo

c...     next the test particles
         call step_kdk_tp_hjs(i1sttp,nbod,ntp,matp,umatp,oloct,mass,
     &              eta,mu,etatp,xjbeg,yjbeg,zjbeg,vxjbeg,vyjbeg,vzjbeg,
     &              ir3jbeg,xbbeg,ybbeg,zbbeg,axbbeg,aybbeg,azbbeg,
     &              vxjh,vyjh,vzjh,xjend,yjend,zjend,ir3jend,
     &              xbend,ybend,zbend,axbend,aybend,azbend,
     &              xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,dt)

c...     store things for next step
         do i=1,nbod
           xbbeg(i) = xbend(i) 
           ybbeg(i) = ybend(i)
           zbbeg(i) = zbend(i)
           axbbeg(i) = axbend(i) 
           aybbeg(i) = aybend(i)
           azbbeg(i) = azbend(i)
           ir3jbeg(i) = ir3jend(i)
         end do

      endif

      return
      end   ! step_kdk_hjs_tid
c------------------------------------------------------------------------

