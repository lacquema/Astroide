c*************************************************************************
c                            SYMBA5Q_STEP_HELIO.F
c*************************************************************************
c This subroutine takes a step in helio coord.  
c Does a KICK than a DRIFT than a KICK.
c ONLY DOES MASSIVE PARTICLES
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 dt            ==>  time step
c                 rtrans        ==>  solar encounter radii (real array)
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  final velocity in helio coord 
c                                       (real arrays)
c Remarks: Based on symba5_step_helio.f but uses the BS routines.
c Authors:  Hal Levison 
c Date:    2/2/2K
c Last revision: 12/13/2K 

      subroutine symba5q_step_helio(nbod,nbodm,mass,j2rp2,
     &     j4rp4,xh,yh,zh,vxh,vyh,vzh,dt,rtrans)

      include '../swift.inc'
      include '../symba5/symba5.inc'
      include '../helioq/helioq.inc'

c...  Inputs Only: 
      integer nbod,nbodm
      real*8 mass(nbod),dt,j2rp2,j4rp4,rtrans(2)

c...  Inputs and Outputs:
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)

c...  Internals:
      integer iz
      real*8 dth 
      real*8 axh(NTPMAX),ayh(NTPMAX),azh(NTPMAX)
      real*8 vxb(NTPMAX),vyb(NTPMAX),vzb(NTPMAX),msys
      real*8 ptxb,ptyb,ptzb            ! Not used here
      real*8 ptxe,ptye,ptze

      save axh,ayh,azh     ! Note this !!
      save vxb,vyb,vzb     ! Note this !!

c----
c...  Executable code 

      dth = 0.5d0*dt

c...  Convert vel to bery to jacobi coords
      call coord_vh2b(nbod,mass,vxh,vyh,vzh,vxb,vyb,vzb,msys)

c...  Do the linear drift due to momentum of the Sun
      call helioq_setbs(1,nbod,mass,xh,yh,zh,vxb,vyb,vzb,dth,
     &     rtrans)

c...  Get the accelerations in helio frame. if frist time step
      iz = 0
      call symba5_helio_getacch(iz,nbod,nbodm,mass,j2rp2,j4rp4,
     &     xh,yh,zh,axh,ayh,azh)

c...  Apply a heliocentric kick for a half dt 
      call kickvh(nbod,vxb,vyb,vzb,axh,ayh,azh,dth)

c..   Drift in helio coords for the full step 
      call helioq_setbs(2,nbod,mass,xh,yh,zh,vxb,vyb,vzb,dt,
     &     rtrans)

c...  Get the accelerations in helio frame. if frist time step
      iz = 0
      call symba5_helio_getacch(iz,nbod,nbodm,mass,j2rp2,j4rp4,
     &     xh,yh,zh,axh,ayh,azh)

c...  Apply a heliocentric kick for a half dt 
      call kickvh(nbod,vxb,vyb,vzb,axh,ayh,azh,dth)

c...  Do the linear drift due to momentum of the Sun
      call helioq_setbs(1,nbod,mass,xh,yh,zh,vxb,vyb,vzb,dth,
     &     rtrans)

c...  convert back to helio velocities
      call coord_vb2h(nbod,mass,vxb,vyb,vzb,vxh,vyh,vzh)

      return
      end   ! symba5q_step_helio
c---------------------------------------------------------------------

