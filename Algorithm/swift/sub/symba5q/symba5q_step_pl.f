c*************************************************************************
c                            SYMBA5Q_STEP_PL.F
c*************************************************************************
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 nbodm         ==>  location of the last massie body 
c                                    (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 dt            ==>  time step
c                 lclose        ==> .true. --> marge particles if they
c                                    get too close. Read in that 
c                                    distance in io_init_pl
c                                      (logical*2 scalar)
c                 rpl           ==>  physical size of a planet.
c                                    (real array)
c                 eoff          ==>  Energy offset (real scalar)
c                 rhill         ==>  size of planet's hills sphere
c                                    (real array)
c                 mtiny         ==>  Small mass  (real array)
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  final velocity in helio coord 
c                                       (real arrays)
c                 rpl           ==>  Recalculated physical size of a planet.
c                                    if merger happened (real array)
c                 nbod          ==>  Recalculated number of massive bodies 
c                                    if merger happened (int scalar)
c                 mass          ==>  Recalculated mass of bodies 
c                                    if merger happened (real array)
c                 isenc         ==>  0 --> No encounter during last dt
c                                    1 --> There was encounters
c                                     (integer scalar)
c                 mergelst      ==>  list of mergers (int array)
c                 mergecnt      ==>  count of mergers (int array)
c                 iecnt         ==>  Number of encounters (int*2 array)
c                 eoff          ==>  Energy offset (real scalar)
c                 rhill         ==>  size of planet's hills sphere
c                                    (real array)
c
c Remarks: Based on symba5_step_pl.f 
c Authors:  Hal Levison
c Date:    02/3/2K
c Last revision: 

      subroutine symba5q_step_pl(i1st,time,nbod,nbodm,mass,j2rp2,
     &     j4rp4,xh,yh,zh,vxh,vyh,vzh,dt,lclose,rpl,isenc,
     &     mergelst,mergecnt,iecnt,eoff,rhill,mtiny)

      include '../swift.inc'
      include '../symba5/symba5.inc'

c...  Inputs Only: 
      integer nbod,i1st,nbodm
      real*8 mass(nbod),dt,time,j2rp2,j4rp4,mtiny
      logical*2 lclose 

c...  Inputs and Outputs:
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 rpl(nbod),eoff,rhill(NTPMAX)

c...  Outputs only
      integer isenc
      integer*2 iecnt(NTPMAX),ielev(NTPMAX)
      integer mergelst(2,NTPMAX),mergecnt

c...  Internals
      integer i,j,ieflg,irec,i1stloc,isun
      logical*1 svdotr            ! Not used in the routine
      integer*2 ielst(2,NENMAX),ielc
      real*8 rtrans(2)

      data i1stloc/0/

      save rtrans,i1stloc

c----
c...  Executable code 

      if(i1stloc.eq.0) then
         write(*,*) 'Input the inner and outer transition distance'
         read(*,*) rtrans
         rtrans(1) = rtrans(1)**2
         rtrans(2) = rtrans(2)**2
         i1stloc = 1    ! turn this off
      endif

      do i=2,nbod
         iecnt(i) = 0
         ielev(i) = -1
      enddo
      isenc = 0

c...  check for encounters
      irec = 0
      ielc = 0
      do j=2,nbodm
         do i=j+1,nbod
            ieflg = 0
            call symba5_chk(rhill,nbod,i,j,mass,xh,yh,zh,vxh,
     &           vyh,vzh,dt,irec,ieflg,svdotr)
            if(ieflg.ne.0) then
               isenc = 1
               iecnt(i) = iecnt(i) + 1
               iecnt(j) = iecnt(j) + 1
               ielev(i) = 0
               ielev(j) = 0
               ielc = ielc + 1
               if(ielc.gt.NENMAX) then
                  write(*,*) 'ERROR: encounter matrix is filled.'
                  write(*,*) 'STOPPING'
                  call util_exit(1)
               endif
               ielst(1,ielc) = i
               ielst(2,ielc) = j
            endif
         enddo
      enddo

c      write(*,*) 'isenc = ',isenc

c..   check to see if we are getting clonse to the Sun
      call helioq_chk_sun(nbod,xh,yh,zh,vxh,vyh,vzh,rtrans,dt,isun)

cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c      isun = 1
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c...  do a step
      if(isenc.eq.0) then
         if(isun.eq.0) then
            call symba5_step_helio(i1st,nbod,nbodm,mass,j2rp2,
     &           j4rp4,xh,yh,zh,vxh,vyh,vzh,dt)
         else
            call symba5q_step_helio(nbod,nbodm,mass,j2rp2,
     &           j4rp4,xh,yh,zh,vxh,vyh,vzh,dt,rtrans)
            i1st = 0  
         endif
         mergecnt=0
      else
         if(isun.eq.0) then
            call symba5_step_interp(time,iecnt,ielev,nbod,
     &           nbodm,mass,rhill,j2rp2,j4rp4,lclose,rpl,xh,yh,zh,
     &           vxh,vyh,vzh,dt,mergelst,mergecnt,eoff,ielc,ielst,
     &           mtiny)
         else
            call symba5q_step_interp(time,iecnt,ielev,nbod,
     &           nbodm,mass,rhill,j2rp2,j4rp4,lclose,rpl,xh,yh,zh,
     &           vxh,vyh,vzh,dt,mergelst,mergecnt,eoff,ielc,ielst,
     &           mtiny,rtrans)
         endif
         i1st = 0  
      endif

      return
      end                       ! symba5q_step_pl.f
c-----------------------------------------------------------

