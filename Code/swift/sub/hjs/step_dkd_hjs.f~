c*************************************************************************
c                            STEP_DKD_HJS.F
c*************************************************************************
c This subroutine takes nsta dt*step in generalized Jacobi coordinates (HJS)  
c both massive and test particles
c It performs eventually a corrector step.
c
c             Input:
c                 nsta          ==> The number of stages for the correction
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
c Remarks: Adopted from step_kdk_hjs.f
c Authors:  Herve Beust 
c Date:    Sep 1, 2006


      subroutine step_dkd_hjs(nsta,i1st,time,nbod,ntp,oloc,oloct,
     &     mat,umat,matp,umatp,mass,eta,mu,etatp,xj,yj,zj,vxj,vyj,vzj,
     &     xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,rstat,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st,nsta
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
      integer i1sttp,i,s
      real*8 xbmid(NPLMAX),ybmid(NPLMAX),zbmid(NPLMAX)
      real*8 axbmid(NPLMAX),aybmid(NPLMAX),azbmid(NPLMAX)
      real*8 xjmid(NPLMAX),yjmid(NPLMAX),zjmid(NPLMAX)
      real*8 vxjmid(NPLMAX),vyjmid(NPLMAX),vzjmid(NPLMAX)
      real*8 ir3jmid(NPLMAX)
      real*8 vxjh(NPLMAX),vyjh(NPLMAX),vzjh(NPLMAX)

c      save xbbeg,ybbeg,zbbeg,axbbeg,aybbeg,azbbeg,ir3jbeg ! Note this !
c----
c...  Executable code 

c... If nsta>1 first do a positive (s=1) correction step

      if (nsta.gt.1) then
        s = 1
        call  step_dkd_corr_hjs(nsta,s,i1st,time,nbod,ntp,oloc,
     &     oloct,mat,umat,matp,umatp,mass,eta,mu,etatp,xj,yj,zj,
     &     vxj,vyj,vzj,xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,rstat,dt)
      end if

c...  Do nsta integration steps. If nsta=1, this is just a normal single
c...  step
      do i = 1,nsta
        i1sttp = i1st

c...  first do the planets

        call step_dkd_pl_hjs(i1st,nbod,oloc,umat,mat,mass,eta,mu,
     &     xj,yj,zj,vxj,vyj,vzj,xjmid,yjmid,zjmid,vxjmid,vyjmid,vzjmid,
     &     ir3jmid,xbmid,ybmid,zbmid,axbmid,aybmid,azbmid,dt)

        if(ntp.ne.0) then

c...     next the test particles

           call step_dkd_tp_hjs(i1sttp,nbod,ntp,matp,umatp,oloct,mass,
     &              eta,mu,etatp,xjmid,yjmid,zjmid,vxjmid,vyjmid,vzjmid,
     &              ir3jmid,xbmid,ybmid,zbmid,axbmid,aybmid,azbmid,
     &              xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,dt)

        endif
      end do

c...  If nsta>1, do finally a negative (s=-1) corrector step

      if (nsta.gt.1) then
        s = -1
        call  step_dkd_corr_hjs(nsta,s,i1st,time,nbod,ntp,oloc,
     &     oloct,mat,umat,matp,umatp,mass,eta,mu,etatp,xj,yj,zj,
     &     vxj,vyj,vzj,xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,rstat,dt)
      end if


      return
      end   ! step_kdk_hjs
c------------------------------------------------------------------------

