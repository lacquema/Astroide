c*************************************************************************
c                            STEP_DKD_HJS_TID.F
c*************************************************************************
c This subroutine takes a step in generalized Jacobi coordinates (HJS)  
c both massive and test particles
c
c                  ==> takes tides & relativity into account
c                  ==> Also computes the evolution of the spin vectors         
c                  ==> Does correction steps
c
c             Input:
c                 nsta          ==> # of correction stages
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
c Remarks: Adopted from step_kdk_hjs_tid.f & step_dkd_hjs.f
c Authors:  Herve Beust 
c Date:   Mar 10, 2010

      subroutine step_dkd_hjs_tid(nsta,i1st,time,nbod,tid,ntp,
     &     oloc,oloct,mat,umat,matp,umatp,mass,inert,krpl5,Qlov,
     &     eta,mu,etatp,xj,yj,zj,vxj,vyj,vzj,rotx,roty,rotz,
     &     xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,rstat,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nsta,nbod,ntp,i1st
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
      integer i1sttp,i,s
      real*8 xbmid(NPLMAX),ybmid(NPLMAX),zbmid(NPLMAX)
      real*8 axbmid(NPLMAX),aybmid(NPLMAX),azbmid(NPLMAX)
      real*8 xjmid(NPLMAX),yjmid(NPLMAX),zjmid(NPLMAX)
      real*8 vxjmid(NPLMAX),vyjmid(NPLMAX),vzjmid(NPLMAX)
      real*8 ir3jmid(NPLMAX)
      real*8 vxjh(NPLMAX),vyjh(NPLMAX),vzjh(NPLMAX)

      save xbmid,ybmid,zbmid,axbmid,aybmid,azbmid,ir3jmid ! Note this !
c----
c...  Executable code 

c... If nsta>1 first do a positive (s=1) correction step

      if (nsta.gt.1) then
        s = 1
        call  step_dkd_corr_hjs_tid(nsta,s,i1st,time,nbod,tid,ntp,
     &     oloc,oloct,mat,umat,matp,umatp,mass,inert,krpl5,Qlov,
     &     eta,mu,etatp,xj,yj,zj,vxj,vyj,vzj,rotx,roty,rotz,
     &     xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,rstat,dt)	
      end if

c...  Do nsta integration steps. If nsta=1, this is just a normal single
c...  step
      do i = 1,nsta
        i1sttp = i1st

c...  first do the planets

        call step_dkd_pl_hjs_tid(i1st,nbod,tid,oloc,umat,mat,mass,
     &      inert,krpl5,Qlov,eta,mu,xj,yj,zj,vxj,vyj,vzj,rotx,roty,rotz,
     &      xjmid,yjmid,zjmid,vxjmid,vyjmid,vzjmid,
     &      ir3jmid,xbmid,ybmid,zbmid,axbmid,aybmid,azbmid,dt)	

        if(ntp.ne.0) then

c...     next the test particles

          call step_dkd_tp_hjs(i1sttp,nbod,ntp,matp,umatp,oloct,mass,
     &              eta,mu,etatp,xjmid,yjmid,zjmid,vxjmid,vyjmid,vzjmid,
     &              ir3jmid,xbmid,ybmid,zbmid,axbmid,aybmid,azbmid,
     &              xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,dt)

        end if
      end do

c...  If nsta>1, do finally a negative (s=-1) corrector step

      if (nsta.gt.1) then
        s = -1
        call  step_dkd_corr_hjs_tid(nsta,s,i1st,time,nbod,tid,ntp,
     &     oloc,oloct,mat,umat,matp,umatp,mass,inert,krpl5,Qlov,
     &     eta,mu,etatp,xj,yj,zj,vxj,vyj,vzj,rotx,roty,rotz,
     &     xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,rstat,dt)	
      end if

      return
      end   ! step_dkd_hjs_tid
c------------------------------------------------------------------------

