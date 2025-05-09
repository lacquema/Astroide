c*************************************************************************
c                            STEP_S6B_HJS.F
c*************************************************************************
c This subroutine takes a step in generalized Jacobi coordinates (HJS)  
c both massive and test particles
c    Uses S6B 6th order symplectic algorithm (Chambers & Murison 2000)
c
c             Input:
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
c Date:    Mar 24, 2010

      subroutine step_s6b_hjs(i1st,time,nbod,ntp,oloc,oloct,mat,umat,
     &     matp,umatp,mass,eta,mu,etatp,xj,yj,zj,vxj,vyj,vzj,
     &     xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,rstat,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
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
      real*8 xbbeg(NPLMAX),ybbeg(NPLMAX),zbbeg(NPLMAX)
      real*8 axbbeg(NPLMAX),aybbeg(NPLMAX),azbbeg(NPLMAX)
      real*8 xjbeg(NPLMAX),yjbeg(NPLMAX),zjbeg(NPLMAX)
      real*8 xjend(NPLMAX),yjend(NPLMAX),zjend(NPLMAX)
      real*8 vxjbeg(NPLMAX),vyjbeg(NPLMAX),vzjbeg(NPLMAX)
      real*8 ir3jbeg(NPLMAX)

      real*8 xj1(NPLMAX),yj1(NPLMAX),zj1(NPLMAX),ir3j1(NPLMAX)
      real*8 vxj1(NPLMAX),vyj1(NPLMAX),vzj1(NPLMAX)         
      real*8 xb1(NPLMAX),yb1(NPLMAX),zb1(NPLMAX)            
      real*8 axb1(NPLMAX),ayb1(NPLMAX),azb1(NPLMAX)         

      real*8 xj2(NPLMAX),yj2(NPLMAX),zj2(NPLMAX),ir3j2(NPLMAX)
      real*8 vxj2(NPLMAX),vyj2(NPLMAX),vzj2(NPLMAX)         
      real*8 xb2(NPLMAX),yb2(NPLMAX),zb2(NPLMAX)            
      real*8 axb2(NPLMAX),ayb2(NPLMAX),azb2(NPLMAX)         

      real*8 xj3(NPLMAX),yj3(NPLMAX),zj3(NPLMAX),ir3j3(NPLMAX)
      real*8 vxj3(NPLMAX),vyj3(NPLMAX),vzj3(NPLMAX)         
      real*8 xb3(NPLMAX),yb3(NPLMAX),zb3(NPLMAX)            
      real*8 axb3(NPLMAX),ayb3(NPLMAX),azb3(NPLMAX)         

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

      call step_s6b_pl_hjs(i1st,nbod,oloc,umat,mat,mass,eta,mu,
     &              ir3jbeg,xbbeg,ybbeg,zbbeg,axbbeg,aybbeg,azbbeg,  
     &              xj1,yj1,zj1,vxj1,vyj1,vzj1,                      
     &              ir3j1,xb1,yb1,zb1,axb1,ayb1,azb1,                
     &              xj2,yj2,zj2,vxj2,vyj2,vzj2,                      
     &              ir3j2,xb2,yb2,zb2,axb2,ayb2,azb2,                
     &              xj3,yj3,zj3,vxj3,vyj3,vzj3,                      
     &              ir3j3,xb3,yb3,zb3,axb3,ayb3,azb3,                
     &              xj,yj,zj,vxj,vyj,vzj,dt)    

      if(ntp.ne.0) then

c...     next the test particles
         call step_s6b_tp_hjs(i1sttp,nbod,ntp,matp,umatp,oloct,mass,
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

c...     store things for next step
         do i=1,nbod
           xbbeg(i) = xb3(i) 
           ybbeg(i) = yb3(i)
           zbbeg(i) = zb3(i)
           axbbeg(i) = axb3(i) 
           aybbeg(i) = ayb3(i)
           azbbeg(i) = azb3(i)
           ir3jbeg(i) = ir3j3(i)
         end do

      endif

      return
      end   ! step_s6b_hjs
c------------------------------------------------------------------------

