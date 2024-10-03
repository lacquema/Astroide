c*************************************************************************
c                            STEP_DKD_TP_CORR_HJS.F
c*************************************************************************
c This subroutine takes a CORRECTOR step in Generalized Jacobi coords (HJS)  
c Does a DRIFT for a*dt then a KICK for b*dt then a DRIFT for -adt.
c Does Chi(a*dt,b*dt) cf Wisdom 2006, AJ 131, 2294
c ONLY DOES TEST PARTICLES
c
c             Input:
c                 a,b            ==>  The corrector coefficients
c                 i1st           ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod           ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of test bodies (int scalar)
c                 mass           ==>  mass of bodies (real array)
c                 matp,umatp     ==>  Conversion vectors for tp's
c                                      (2D real arrays)
c                 oloct          ==>  Link between tp's and orbits
c                                      (2D integer array)
c                 eta,mu,etatp   ==>  Masses of centers (real arrays)
c                 xjmid,yjmid,zjmid ==> Bodies Jac. position at middle point
c                                       (real arrays)
c                 xbmid,ybmid,zbmid ==> Bodies bary position at middle point
c                                       (real arrays)
c                 vxjmid,vyjmid,vzjmid ==> Bodies Jac. veloc. at middle point
c                                       (real arrays)
c                 ir3jmid ==> 1/rj^3 at middle point (real arrays)
c                 axbmid,aybmid,azbmid ==> Bodies bary accs. at middle point
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
c                 xht,yht,zht    ==>  final position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final position in helio coord 
c                                       (real arrays)
c
c Remarks: Adapted from step_dkd_tp_hjs.f
c Authors:  Herve Beust
c Date:    Sep 1, 2006 


      subroutine step_dkd_tp_corr_hjs(a,b,i1st,nbod,ntp,matp,umatp,
     &              oloct,mass,eta,mu,etatp,xjmid,yjmid,zjmid,
     &              vxjmid,vyjmid,vzjmid,ir3jmid,xbmid,ybmid,zbmid,
     &              axbmid,aybmid,azbmid,xjt,yjt,zjt,vxjt,vyjt,vzjt,
     &              istat,dt)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st,oloct(NPLMAX,NTPMAX)
      real*8 a,b
      real*8 umatp(NPLMAX,NTPMAX),matp(NPLMAX,NTPMAX)
      real*8 mass(nbod),eta(nbod),mu(nbod),dt
      real*8 xjmid(nbod),yjmid(nbod),zjmid(nbod)
      real*8 vxjmid(nbod),vyjmid(nbod),vzjmid(nbod)
      real*8 xbmid(nbod),ybmid(nbod),zbmid(nbod)
      real*8 axbmid(nbod),aybmid(nbod),azbmid(nbod)
      real*8 ir3jmid(nbod),etatp(ntp)

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xjt(ntp),yjt(ntp),zjt(ntp)
      real*8 vxjt(ntp),vyjt(ntp),vzjt(ntp)

c...  Internals:
      integer j,k,istati(NSTAT)
      real*8 dta,dtb,axbttp,aybttp,azbttp,axjttp,ayjttp,azjttp 
      real*8 axjt(NTPMAX),ayjt(NTPMAX),azjt(NTPMAX)
      real*8 xbttp,zbttp,ybttp,vxbttp,vybttp,vzbttp 

c      save axjt,ayjt,azjt     ! Note this !!

c----
c...  Executable code 

      dta = a*dt
      dtb = b*dt

c...  loop over all tp's
      do j = 1,ntp
        if (istat(j,1).eq.0) then

c...  Take a drift forward half step
          do k=1,NSTAT
            istati(k) = istat(j,k)
          end do
          call drift_tp(1,etatp(j),xjt(j),yjt(j),zjt(j),vxjt(j),
     &                    vyjt(j),vzjt(j),-dta,istati)	

c...      Convert to barycentric coordinates
          call coord_g2b_tp(nbod,umatp(1,j),xjmid,yjmid,zjmid,
     &          vxjmid,vyjmid,vzjmid,xjt(j),yjt(j),zjt(j),vxjt(j),
     &          vyjt(j),vzjt(j),xbttp,ybttp,zbttp,vxbttp,vybttp,vzbttp)

c...      Get the acceleration in bary frame.
          call getacch_tp_hjs(nbod,mass,eta,mu,xjmid,yjmid,
     &           zjmid,xbmid,ybmid,zbmid,ir3jmid,oloct(1,j),
     &           etatp(j),xbttp,ybttp,zbttp,xjt(j),yjt(j),zjt(j),
     &           axbttp,aybttp,azbttp)
c...      Convert bary accels to Jacobi accels (use vel. procedure)
          call coord_vb2vg_tp(nbod,matp(1,j),axbmid,aybmid,azbmid,
     &                 axbttp,aybttp,azbttp,axjttp,ayjttp,azjttp)

c...  Apply a Jacobi kick for a full dt 
          call kickvh_tp(1,vxjt(j),vyjt(j),vzjt(j),
     &                            axjttp,ayjttp,azjttp,istat(j,1),dtb) 

c...  Take a drift forward half step
          do k=1,NSTAT
            istati(k) = istat(j,k)
          end do
          call drift_tp(1,etatp(j),xjt(j),yjt(j),zjt(j),vxjt(j),
     &                    vyjt(j),vzjt(j),dta,istati)	

        end if
      end do
      i1st = 1

      return
      end   ! step_dkd_tp_corr_hjs
c---------------------------------------------------------------------

