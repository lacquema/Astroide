c*************************************************************************
c                            HELIOQ_STEP_PL.F
c*************************************************************************
c This subroutine takes a step in helio coord.  
c Does a KICK than a DRIFT than a KICK.
c ONLY DOES MASSIVE PARTICLES
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 dt            ==>  time step
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  final velocity in helio coord 
c                                       (real arrays)
c Remarks: Based on helio_step_pl.f but does not pass the intermediate
c          positions and velocities back for the TP to use.
c Authors:  Hal Levison 
c Date:    2/2/2K
c Last revision: 

      subroutine helioq_step_pl(i1st,nbod,mass,j2rp2,
     &     j4rp4,xh,yh,zh,vxh,vyh,vzh,dt)

      include '../swift.inc'
      include 'helioq.inc'

c...  Inputs Only: 
      integer nbod,i1st,nbodm
      real*8 mass(nbod),dt,j2rp2,j4rp4

c...  Inputs and Outputs:
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)

c...  Internals:
      real*8 dth 
      real*8 axh(NTPMAX),ayh(NTPMAX),azh(NTPMAX)
      real*8 vxb(NTPMAX),vyb(NTPMAX),vzb(NTPMAX),msys
      real*8 ptxb,ptyb,ptzb            ! Not used here
      integer isun,i1stloc
      real*8 rtrans(2)

      data i1stloc/0/

      save axh,ayh,azh     ! Note this !!
      save vxb,vyb,vzb     ! Note this !!
      save rtrans,i1stloc

c----
c...  Executable code 

      dth = 0.5d0*dt

      if(i1stloc.eq.0) then
         write(*,*) 'Input the inner and outer transition distance'
         read(*,*) rtrans
         rtrans(1) = rtrans(1)**2
         rtrans(2) = rtrans(2)**2
         i1stloc = 1    ! turn this off
      endif


      if(i1st.eq.0) then
c...      Convert vel to bery to jacobi coords
          call coord_vh2b(nbod,mass,vxh,vyh,vzh,vxb,vyb,vzb,msys)
         i1st = 1    ! turn this off
      endif

c..   check to see if we are getting close to the Sun
      call helioq_chk_sun(nbod,xh,yh,zh,vxb,vyb,vzb,rtrans,dt,isun)
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c      write(77,*) isun,' = isun' 
c      isun = 1
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c...  Do the linear drift due to momentum of the Sun
      if(isun.eq.0) then
         call helio_lindrift(nbod,mass,vxb,vyb,vzb,dth,
     &        xh,yh,zh,ptxb,ptyb,ptzb)
      else
         call helioq_setbs(1,nbod,mass,xh,yh,zh,vxb,vyb,vzb,dth,
     &        rtrans)
      endif

c...  Get the accelerations in helio frame. 
      call helio_getacch(nbod,mass,j2rp2,j4rp4,xh,yh,zh,axh,ayh,azh)

c...  Apply a heliocentric kick for a half dt 
      call kickvh(nbod,vxb,vyb,vzb,axh,ayh,azh,dth)

c..   Drift in helio coords for the full step 
      if(isun.eq.0) then
         call helio_drift(nbod,mass,xh,yh,zh,vxb,vyb,vzb,dt)
      else
         call helioq_setbs(2,nbod,mass,xh,yh,zh,vxb,vyb,vzb,dt,
     &        rtrans)
      endif

c...  Get the accelerations in helio frame. if frist time step
      call helio_getacch(nbod,mass,j2rp2,j4rp4,xh,yh,zh,axh,ayh,azh)

c...  Apply a heliocentric kick for a half dt 
      call kickvh(nbod,vxb,vyb,vzb,axh,ayh,azh,dth)

c...  Do the linear drift due to momentum of the Sun
      if(isun.eq.0) then
         call helio_lindrift(nbod,mass,vxb,vyb,vzb,dth,
     &        xh,yh,zh,ptxb,ptyb,ptzb)
      else
         call helioq_setbs(1,nbod,mass,xh,yh,zh,vxb,vyb,vzb,dth,
     &        rtrans)
      endif

c...  convert back to helio velocities
      call coord_vb2h(nbod,mass,vxb,vyb,vzb,vxh,vyh,vzh)

      return
      end   ! symba5_step_helio
