c*************************************************************************
c                            STEP_S6B_ONETP_HJS.F
c*************************************************************************
c This subroutine takes a step in Generalized Jacobi coords (HJS) for 1 tp 
c    Uses S6B 6th order symplectic algorithm (Chambers & Murison 2000)
c ONLY DOES ONE TEST PARTICLE
c
c       Input:
c          i1st        ==>  = 0 if first step; = 1 not (int scalar)
c          nbod        ==>  number of massive bodies (int scalar)
c          mass        ==>  mass of bodies (real array)
c          matp,umatp  ==>  Conversion vectors for tp's (real arrays)
c          oloct       ==>  Link between tp's and orbits (integer array)
c          eta,mu,etatp   ==>  Masses of centers & sats (real arrays+scalar)
c          xjpl,yjpl,zjpl,ir3j ==> Bodies Jacobi positions at substeps of dt
c                                       (real 2D arrays) 1 = begin; 4 = end 
c          vxjpl,vyjpl,vzjpl ==> Bodies Jacobi velocities at substeps of dt
c                                       (real 2D arrays) 1 = begin; 4 = end 
c          xbpl,ybpl,zbpl ==> Bodies Barycentric positions at substeps of dt
c                                       (real 2D arrays) 1 = begin; 4 = end 
c          axbpl,aybpl,azbpl ==> Bodies Barycentric accels. at substeps of dt
c                                       (real 2D arrays) 1 = begin; 4 = end 
c                 istat        ==>  status of the test paricle (integer array)
c                                   istat(1) = 0 ==> active:  = 1 not
c                                   istat(2) = -1 ==> Danby did not work
c                 dt             ==>  time step
c             Output:
c                 xjt,yjt,zjt    ==>  final position in Jacobi coord (reals)
c                 vxjt,vyjt,vzjt ==>  final position in Jacobi coord (reals)
c
c Remarks: Adapted from step_s6b_onetp.f & step_s6b_tp_hjs.f
c Authors:  Herve Beust
c Date:    Feb. 13, 2025

      subroutine step_s6b_onetp_hjs(i1st,nbod,matp,umatp,oloct,mass,
     &              eta,mu,etatp,xjpl,yjpl,zjpl,ir3j,
     &              vxjpl,vyjpl,vzjpl,xbpl,ybpl,zbpl,axbpl,aybpl,azbpl,
     &              xjt,yjt,zjt,vxjt,vyjt,vzjt,axjt,ayjt,azjt,istat,dt)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st,oloct(nbod)
      real*8 umatp(nbod),matp(nbod)
      real*8 mass(nbod),eta(nbod),mu(nbod),dt,etatp

      real*8 xjpl(nbod,4),yjpl(nbod,4),zjpl(nbod,4),ir3j(nbod,4)
      real*8 vxjpl(nbod,4),vyjpl(nbod,4),vzjpl(nbod,4)      
      real*8 xbpl(nbod,4),ybpl(nbod,4),zbpl(nbod,4)
      real*8 axbpl(nbod,4),aybpl(nbod,4),azbpl(nbod,4)
      
c...  Inputs and Outputs:
      integer istat(NSTAT)
      real*8 xjt,yjt,zjt,vxjt,vyjt,vzjt,axjt,ayjt,azjt

c...  Internals:
      integer iflg
      real*8 dt1,dt2,dt3,dt4
      real*8 axbt,aybt,azbt
      real*8 xbt,zbt,ybt,vxbt,vybt,vzbt 

c      real*8, parameter :: kk = 1.2599210498948731647672106d0  ! 2^(1/3)
c      real*8, parameter :: cc = 2.0d0-kk
      real*8, parameter :: sq5 = 2.236067977499789696409174d0 ! sqrt(5)
      real*8, parameter :: t1 = 1.d0/12.d0
      real*8, parameter :: t2 = 0.5d0*(1.d0-1.d0/sq5)
      real*8, parameter :: t3 = 5.d0/12.d0
      real*8, parameter :: t4 = 1.d0/sq5

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
c...      Convert to barycentric coordinates
         call coord_g2b_tp(nbod,umatp,xjpl(:,1),yjpl(:,1),zjpl(:,1),
     &          vxjpl(:,1),vyjpl(:,1),vzjpl(:,1),xjt,yjt,zjt,vxjt,
     &          vyjt,vzjt,xbt,ybt,zbt,vxbt,vybt,vzbt)

c...      Get the acceleration in bary frame.
         call getacch_tp_hjs(nbod,mass,eta,mu,xjpl(:,1),yjpl(:,1),
     &        zjpl(:,1),xbpl(:,1),ybpl(:,1),zbpl(:,1),ir3j(:,1),
     &        oloct,etatp,xbt,ybt,zbt,xjt,yjt,zjt,axbt,aybt,azbt)
c...        Convert bary accels to Jacobi accels (use vel. procedure)
         call coord_vb2vg_tp(nbod,matp,axbpl(:,1),aybpl(:,1),azbpl(:,1),
     &                 axbt,aybt,azbt,axjt,ayjt,azjt)
         i1st = 1    ! turn this off        
      end if

c...  Apply a Jacobi kick for a dt1
      call kickvh_tp(1,vxjt,vyjt,vzjt,axjt,ayjt,azjt,istat(1),dt1) 

c...  Take a first Jacobi drift for dt2
      call drift_one(etatp,xjt,yjt,zjt,vxjt,vyjt,vzjt,dt2,iflg)
      if (iflg.ne.0) then
            istat(1) = 1
            istat(2) = -1
      end if
      
c...  After drift, compute bary pos. and vels. 
      call coord_g2b_tp(nbod,umatp,xjpl(:,2),yjpl(:,2),zjpl(:,2),
     &       vxjpl(:,2),vyjpl(:,2),vzjpl(:,2),xjt,yjt,zjt,vxjt,
     &       vyjt,vzjt,xbt,ybt,zbt,vxbt,vybt,vzbt)

c...  Get the acceleration in bary frame.
      call getacch_tp_hjs(nbod,mass,eta,mu,xjpl(:,2),yjpl(:,2),
     &       zjpl(:,2),xbpl(:,2),ybpl(:,2),zbpl(:,2),ir3j(:,2),
     &       oloct,etatp,xbt,ybt,zbt,xjt,yjt,zjt,axbt,aybt,azbt)
c...        Convert bary accels to Jacobi accels (use vel. procedure)
      call coord_vb2vg_tp(nbod,matp,axbpl(:,2),aybpl(:,2),azbpl(:,2),
     &                 axbt,aybt,azbt,axjt,ayjt,azjt)
c...  Apply a second Jacobi kick for dt3 
      call kickvh_tp(1,vxjt,vyjt,vzjt,axjt,ayjt,azjt,istat(1),dt3) 

c...  Second drift in Jacobi coords for dt4
      call drift_one(etatp,xjt,yjt,zjt,vxjt,vyjt,vzjt,dt4,iflg)	
      if (iflg.ne.0) then
            istat(1) = 1
            istat(2) = -1
      end if
      
c...  After drift, compute bary pos. and vels. 
      call coord_g2b_tp(nbod,umatp,xjpl(:,3),yjpl(:,3),zjpl(:,3),
     &       vxjpl(:,3),vyjpl(:,3),vzjpl(:,3),xjt,yjt,zjt,vxjt,
     &       vyjt,vzjt,xbt,ybt,zbt,vxbt,vybt,vzbt)

c...  Get the acceleration in bary frame.
      call getacch_tp_hjs(nbod,mass,eta,mu,xjpl(:,3),yjpl(:,3),
     &       zjpl(:,3),xbpl(:,3),ybpl(:,3),zbpl(:,3),ir3j(:,3),
     &       oloct,etatp,xbt,ybt,zbt,xjt,yjt,zjt,axbt,aybt,azbt)
c...  Convert bary accels to Jacobi accels (use vel. procedure)
      call coord_vb2vg_tp(nbod,matp,axbpl(:,3),aybpl(:,3),azbpl(:,3),
     &                 axbt,aybt,azbt,axjt,ayjt,azjt)

c...  Apply a third Jacobi kick for dt3 
      call kickvh_tp(1,vxjt,vyjt,vzjt,axjt,ayjt,azjt,istat(1),dt3) 

c...  Third drift in Jacobi coords for dt2
      call drift_one(etatp,xjt,yjt,zjt,vxjt,vyjt,vzjt,dt2,iflg)	
      if (iflg.ne.0) then
            istat(1) = 1
            istat(2) = -1
      end if

c...  After drift, compute bary pos. and vels. 
      call coord_g2b_tp(nbod,umatp,xjpl(:,4),yjpl(:,4),zjpl(:,4),
     &       vxjpl(:,4),vyjpl(:,4),vzjpl(:,4),xjt,yjt,zjt,vxjt,
     &       vyjt,vzjt,xbt,ybt,zbt,vxbt,vybt,vzbt)

c...  Get the acceleration in bary frame.
      call getacch_tp_hjs(nbod,mass,eta,mu,xjpl(:,4),yjpl(:,4),
     &       zjpl(:,4),xbpl(:,4),ybpl(:,4),zbpl(:,4),ir3j(:,4),
     &       oloct,etatp,xbt,ybt,zbt,xjt,yjt,zjt,axbt,aybt,azbt)
c...  Convert bary accels to Jacobi accels (use vel. procedure)
      call coord_vb2vg_tp(nbod,matp,axbpl(:,4),aybpl(:,4),azbpl(:,4),
     &                 axbt,aybt,azbt,axjt,ayjt,azjt)

c...  Apply a fourth Jacobi kick for dt1 
      call kickvh_tp(1,vxjt,vyjt,vzjt,axjt,ayjt,azjt,istat(1),dt1) 

      return
      end   ! step_s6b_onetp_hjs
c---------------------------------------------------------------------

