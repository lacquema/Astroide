c*************************************************************************
c                            HELIOQ_SETBS.F
c*************************************************************************
c This sets up the BS integrator for helioq
c
c             Input:
c                 iwhich        ==> determines which H to integrate
c                                       (int scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxb,vyb,vzb   ==>  initial velocity in bary coord 
c                                    (real arrays)
c                 dt            ==>  time step (real scalar)
c                 rtrans        ==>  solar encounter radii (real array)
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxb,vyb,vzb   ==>  final velocity in bary coord 
c                                       (real arrays)
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/2/2K
c Last revision: 2/18/2K

      subroutine helioq_setbs(iwhich,nbod,mass,xh,yh,zh,vxb,
     &        vyb,vzb,dt,rtrans)

      include '../swift.inc'
      include 'helioq.inc'
      
c...  Inputs Only: 
      integer nbod,iwhich
      real*8 mass(nbod),dt
      real*8 vxb(nbod),vyb(nbod),vzb(nbod)
      real*8 rtrans(2)

c...  Inputs and Outputs:
      real*8 xh(nbod),yh(nbod),zh(nbod)

c...  Internals:
      integer i,ism
      real*8 ybs(6,NTPMAX),tfake,dttmp,dtleft,dtnext

c----
c...  Executable code 

c...  copy to the big array
      do i=1,nbod-1
         ybs(1,i) = xh(i+1)
         ybs(2,i) = yh(i+1)
         ybs(3,i) = zh(i+1)
         ybs(4,i) = vxb(i+1)
         ybs(5,i) = vyb(i+1)
         ybs(6,i) = vzb(i+1)
      enddo

      tfake = 0.0d0
      dttmp = dt

      ism = 0
      do while( (abs(tfake-dt)/dt) .gt. EPS )    ! just to be real safe
         dtleft = dt-tfake

         call helioq_bs(iwhich,nbod,mass,tfake,dttmp,ybs,rtrans,dtnext)

cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c         write(*,*) iwhich,dt,tfake,dttmp
c         write(77,*) 'Here #setbs1 ',iwhich,dttmp,dt,tfake
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

         if( ((dttmp/dtleft) .lt. 1.0d-6 ) .and. (ism.eq.0)) then 
            write(*,*) 'WARNING:'
            write(*,*) 'BS integrator is taking very small timesteps'
            write(*,*) 'local dt = ',dttmp
            ism = 1
         endif

c         dttmp = dt - tfake
c         dttmp = min(dttmp,(dt - tfake))
         dttmp = min(dtnext,(dt - tfake))

      enddo

cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c      write(77,*) 'Here #setbs100 done'
c      write(*,*) 'done'
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c..... put things back
      do i=1,nbod-1
         xh(i+1) = ybs(1,i)
         yh(i+1) = ybs(2,i)
         zh(i+1) = ybs(3,i)
         vxb(i+1) = ybs(4,i)
         vyb(i+1) = ybs(5,i)
         vzb(i+1) = ybs(6,i)
      enddo

      return
      end   ! helioq_setbs
c-------------------------------------------------------------------------
