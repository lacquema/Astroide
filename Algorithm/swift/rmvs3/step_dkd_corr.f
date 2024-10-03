c*************************************************************************
c                            STEP_DKD_CORR.F
c*************************************************************************
c This subroutine takes a nth-stage CORRECTOR step in heliocentric 
c coordinates. both massive and test particles
c Does C(dt) cf  Wisdom 2006, AJ 131, 2294
c     
c             Input:
c                 nsta,s        ==>  Number of stages and sign
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                                    not used here !!!
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
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
c Authors:  Hervé Beust
c Date:    Feb 13, 2023

      subroutine step_dkd_corr(nsta,s,i1st,time,nbod,ntp,mass,
     &     j2rp2,j4rp4,xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,
     &     vxht,vyht,vzht,istat,rstat,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      integer s,nsta
      real*8 mass(nbod),dt,time,j2rp2,j4rp4

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 rstat(NTPMAX,NSTATR)
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Internals
      integer i,level
      real*8, parameter :: alpha = 0.41833001326703777399d0  ! sqrt(7/40)
      real*8, parameter :: beta = 1.d0/(48.d0*alpha)
      real*8, parameter :: z=0.d0
      real*8, dimension(10,4), parameter :: a = reshape ( (/
     &       -alpha,alpha,z,z,z,z,z,z,z,z,
     &       2.d0*alpha,alpha,-alpha,-2.d0*alpha,z,z,z,z,z,z,
     &       3.d0*alpha,2.d0*alpha,alpha,-alpha,-2.d0*alpha,-3.d0*alpha,
     &                                                         z,z,z,z,
     &       5.d0*alpha,4.d0*alpha,3.d0*alpha,2.d0*alpha,alpha,
     &         -alpha,-2.d0*alpha,-3.d0*alpha,4.d0*alpha,5.d0*alpha /),
     &                                                    (/ 10,4 /) ) 
      real*8, dimension(10,4), parameter :: b = reshape ( (/
     &       -0.5d0*beta,0.5d0*beta,z,z,z,z,z,z,z,z,
     &       -beta/6.d0,5.d0*beta/6.d0,-5.d0*beta/6.d0,beta/6.d0,
     &                                                   z,z,z,z,z,z,
     &       12361.d0*beta/246960.d0,-22651.d0*beta/61740.d0,
     &       53521.d0*beta/49392.d0,-53521.d0*beta/49392.d0,
     &       22651.d0*beta/61740.d0,-12361.d0*beta/246960.d0,z,z,z,z,
     &       2798927.d0*beta/684573120.d0,-329447.d0*beta/6985440.d0,
     &       895249.d0*beta/3622080.d0,-14556229.d0*beta/19015920.d0,
     &       3394141.d0*beta/2328480.d0,-3394141.d0*beta/2328480.d0,
     &       14556229.d0*beta/19015920.d0,-895249.d0*beta/3622080.d0,
     &       329447.d0*beta/6985440.d0,-2798927.d0*beta/684573120.d0 /),
     &                                                  (/ 10,4 /) )
      real*8 aa(10),bb(10)

c----
c...  Executable code 

c...  Select level of correction
      select case ( nsta )
        case (2)
          level = 1
        case (4)
          level = 2
        case (6)
          level = 3 
        case(10)
          level = 4
      end select

c...  Fix coefficients depending whether we are doing C or C^-1

      if (s.eq.1) then
        do i=1,nsta
          aa(i) = a(nsta+1-i,level)
          bb(i) = b(nsta+1-i,level)
        end do
      else if (s.eq.-1) then
        do i=1,nsta
          aa(i) = a(i,level)
          bb(i) = -b(i,level)
        end do
      end if


c...  Do the nsta stages of correction

      do i = 1,nsta
         call step_dkd_corr_z(aa(i),bb(i),i1st,time,nbod,ntp,
     &        mass,j2rp2,j4rp4,xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,
     &        vxht,vyht,vzht,istat,rstat,dt)	
      end do

      return
      end                       ! step_dkd_corr
c------------------------------------------------------------------------
