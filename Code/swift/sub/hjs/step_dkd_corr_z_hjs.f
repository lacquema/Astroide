c*************************************************************************
c                            STEP_DKD_CORR_Z_HJS.F
c*************************************************************************
c This subroutine takes a CORRECTOR step in generalized Jacobi 
c coordinates (HJS). both massive and test particles
c Does Z(a*dt,b*dt)=Chi(a*dt,b*dt).Chi(-a*dt,-b*dt) Wisdom 2006, AJ 131, 2294
c
c             Input:
c                 a,b           ==>  The corrector coefficients
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 oloc          ==>  Link bodies - orbits (2D int array)
c                 oloct         ==>  Link orbits - tp's (2D int array)
c                 mass          ==>  mass of bodies (real array)
c                 mat,umat      ==>  Conversion matrixes for bodies
c                                     (2D real arrays)
c                 matp,umatp    ==>  Conversion vectors for tp's
c                                     (2D real arrays)
c                 eta,mu        ==>  Masses for center & satellites for bods.
c                                     (real arrays)
c                 etatp         ==>  Masses for centers for tp's (real array)
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
c Remarks: Adopted from step_dkd_hjs.f
c Authors:  Herve Beust 
c Date:    Aug 31, 2006


      subroutine step_dkd_corr_z_hjs(a,b,i1st,time,nbod,ntp,oloc,oloct,
     &     mat,umat,matp,umatp,mass,eta,mu,etatp,xj,yj,zj,vxj,vyj,vzj,
     &     xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,rstat,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 a,b
      real*8 mass(nbod),dt,time
      real*8 mu(nbod),eta(nbod),etatp(ntp)
      real*8 umat(NPLMAX,NPLMAX),mat(NPLMAX,NPLMAX)
      real*8 umatp(NPLMAX,NTPMAX),matp(NPLMAX,NTPMAX)
      integer oloc(NPLMAX,NPLMAX),oloct(NPLMAX,NTPMAX)

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 rstat(NTPMAX,NSTATR)
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)
      real*8 xjt(ntp),yjt(ntp),zjt(ntp)
      real*8 vxjt(ntp),vyjt(ntp),vzjt(ntp)

c...  Internals
      integer i1sttp,i
      real*8 xbmid(NPLMAX),ybmid(NPLMAX),zbmid(NPLMAX)
      real*8 axbmid(NPLMAX),aybmid(NPLMAX),azbmid(NPLMAX)
      real*8 xjmid(NPLMAX),yjmid(NPLMAX),zjmid(NPLMAX)
      real*8 vxjmid(NPLMAX),vyjmid(NPLMAX),vzjmid(NPLMAX)
      real*8 ir3jmid(NPLMAX)

c      save xbbeg,ybbeg,zbbeg,axbbeg,aybbeg,azbbeg,ir3jbeg ! Note this !
c----
c...  Executable code 

      i1sttp = i1st

c...  First step, Chi(-a*dt,-b*dt)

c...  first do the planets
      call step_dkd_pl_corr_hjs(-a,-b,i1st,nbod,oloc,umat,mat,mass,
     &     eta,mu,xj,yj,zj,vxj,vyj,vzj,xjmid,yjmid,zjmid,vxjmid,
     &     vyjmid,vzjmid,ir3jmid,xbmid,ybmid,zbmid,axbmid,aybmid,
     &                                                   azbmid,dt)

      if(ntp.ne.0) then
c...     next the test particles
         call step_dkd_tp_corr_hjs(-a,-b,i1sttp,nbod,ntp,matp,umatp,
     &              oloct,mass,eta,mu,etatp,xjmid,yjmid,zjmid,
     &              vxjmid,vyjmid,vzjmid,ir3jmid,xbmid,ybmid,zbmid,
     &              axbmid,aybmid,azbmid,xjt,yjt,zjt,vxjt,vyjt,vzjt,
     &              istat,dt)
      endif

c...  Second step, Chi(a*dt,b*dt)

c...  first do the planets
      call step_dkd_pl_corr_hjs(a,b,i1st,nbod,oloc,umat,mat,mass,
     &     eta,mu,xj,yj,zj,vxj,vyj,vzj,xjmid,yjmid,zjmid,vxjmid,vyjmid,
     &     vzjmid,ir3jmid,xbmid,ybmid,zbmid,axbmid,aybmid,azbmid,dt)

      if(ntp.ne.0) then
c...     next the test particles
         call step_dkd_tp_corr_hjs(a,b,i1sttp,nbod,ntp,matp,umatp,
     &              oloct,mass,eta,mu,etatp,xjmid,yjmid,zjmid,
     &              vxjmid,vyjmid,vzjmid,ir3jmid,xbmid,ybmid,zbmid,
     &              axbmid,aybmid,azbmid,xjt,yjt,zjt,vxjt,vyjt,vzjt,
     &              istat,dt)
      endif

      return
      end   ! step_dkd_corr_z_hjs
c------------------------------------------------------------------------

