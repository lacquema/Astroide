c*************************************************************************
c                            SYMBA6_STEP_PL.F
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
c
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
c Authors:  Craig Agnor & Hal Levison
c Date:    7/20/99
c Last revision: 

      subroutine symba6_step_pl(i1st,time,nbod,nbodm,mass,j2rp2,
     &     j4rp4,xh,yh,zh,vxh,vyh,vzh,dt,lclose,rpl,isenc,
     &     mergelst,mergecnt,iecnt,eoff,rhill,mtiny)

      include '../swift.inc'
      include 'symba6.inc'

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
      integer i,j,ieflg,irec
      logical*1 svdotr            ! Not used in the routine
      integer*2 ielst(2,NENMAX),ielc
      integer ip,jp,ids,ie
      logical notstored
      integer ichck(NSECTMAX), sect_chklst(NSECTMAX,NTPMAX) 
      integer*2 ielst1(NENMAX1,NTPMAX)

      integer nsect,nrz,naz,i1stin
      real*8 rbnd1,dr,dtheta
      character*(chlen) inznfile

      data i1stin/0/

      save nsect,nrz,naz,rbnd1,dr,dtheta,i1stin

c----
c...  Executable code 

      if(i1stin.eq.0) then
         write(*,*) ' '
         write(*,*) 'Enter name of search zones data file : '
         read(*,fmt='(a)') inznfile
         call io_symba6_init_zones(inznfile,nrz,naz,nsect)
         call symba6_sect_bndrs(nrz,naz,nbod,xh,yh,rbnd1,dr,dtheta)
         i1stin = 1
      endif

      do i=2,nbod
         iecnt(i) = 0
         ielev(i) = -1
      enddo
      isenc = 0

c...  check for encounters  
      call symba6_sector_chklist(ichck,sect_chklst,nbod,dt,nsect,
     &     naz,nrz,xh,yh,vxh,vyh,rhill,rbnd1,dr,dtheta)
      irec = 0
      ielc = 0

      do ids = 1,nsect
         i = 1
         do while( (i.le.ichck(ids)) .and. 
     &        (sect_chklst(ids,i).le.nbodm) )
            do j = i+1,ichck(ids)

               ip = sect_chklst(ids,j) ! note ip > jp always
               jp = sect_chklst(ids,i)

               ieflg = 0
               call symba6_chk(rhill,nbod,ip,jp,mass,xh,yh,zh,vxh,
     &              vyh,vzh,dt,irec,ieflg,svdotr)

c.      search encounter list for prior detections of this encounter
               if(ieflg .ne. 0) then
                  notstored = .true.
                  if( (iecnt(ip).ne.0) .and. (iecnt(jp).ne.0)) then
                     do ie=1,iecnt(ip)
                        if(jp.eq.ielst1(ie,ip)) then
                           notstored = .false.
                        endif
                     enddo
                  endif
                  
                  if (notstored) then
                     isenc = 1
                     iecnt(ip) = iecnt(ip) + 1
                     iecnt(jp) = iecnt(jp) + 1
                     ielev(ip) = 0
                     ielev(jp) = 0
                     ielc = ielc + 1
                     if( (ielc.gt.NENMAX) 
     &                    .or. (iecnt(ip).gt.NENMAX1)
     &                    .or. (iecnt(jp).gt.NENMAX1) ) then
                        write(*,*) 'ERROR: encounter new search failed.'
                        write(*,*) 'STOPPING'
                        call util_exit(1)
                     endif
                     ielst(1,ielc) = ip
                     ielst(2,ielc) = jp
                     ielst1(iecnt(jp),jp) = ip
                     ielst1(iecnt(ip),ip) = jp
                  endif
               endif
            enddo
            i = i + 1
         enddo
      enddo

cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c      write(77,*) time,ielc
c      do i=1,ielc
c         write(77,*) '   ',i,ielst(1,i),ielst(2,i)
c      enddo
c      write(77,*) '------------------------------'
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c...  do a step
      if(isenc.eq.0) then
         call symba6_step_helio(i1st,nbod,nbodm,mass,j2rp2,
     &     j4rp4,xh,yh,zh,vxh,vyh,vzh,dt)
         mergecnt=0
      else
         call symba6_step_interp(time,iecnt,ielev,nbod,
     &     nbodm,mass,rhill,j2rp2,j4rp4,lclose,rpl,xh,yh,zh,
     &     vxh,vyh,vzh,dt,mergelst,mergecnt,eoff,ielc,ielst,
     &     mtiny)
         i1st = 0  
      endif

      return
      end ! symba6_step_pl.f
c-----------------------------------------------------------
