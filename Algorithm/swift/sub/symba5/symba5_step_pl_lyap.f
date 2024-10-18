c*************************************************************************
c                            SYMBA5_STEP_PL_LYAP.F
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
c Remarks: Based on symba5_step_pl.f & lyap_step.f 
c Authors:  H. Beust
c Date:    06/08/2007
c Last revision: 

      subroutine symba5_step_pl_lyap(diro,fopenstat,
     &     i1st,time,nbod,nbodm,mass,j2rp2,
     &     j4rp4,xh,yh,zh,vxh,vyh,vzh,dt,lclose,rpl,isenc,
     &     mergelst,mergecnt,iecnt,eoff,rhill,mtiny)

      include '../swift.inc'
      include 'symba5.inc'

c...  Inputs Only: 
      integer nbod,i1st,nbodm
      real*8 mass(nbod),dt,time,j2rp2,j4rp4,mtiny
      logical*2 lclose 
      character*(chlen) fopenstat,diro

c...  Inputs and Outputs:
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 rpl(nbod),eoff,rhill(NTPMAX)

c...  Outputs only
      integer isenc
      integer*2 iecnt(NTPMAX),ielev(NTPMAX)
      integer mergelst(2,NTPMAX),mergecnt

c...  Internals
      integer isencsh,iflgchk
      integer, save :: i1stsh
      integer*2 iecntsh(NTPMAX),ielevsh(NTPMAX)
      integer mergelstsh(2,NTPMAX),mergecntsh
      integer :: i1stin = 0
      integer :: iul = 50
      real*8, save :: xsh(NTPMAX),ysh(NTPMAX),zsh(NTPMAX)
      real*8, save :: vxsh(NTPMAX),vysh(NTPMAX),vzsh(NTPMAX)
      real*8, save :: dist0(NTPMAX),lrsum(NTPMAX),logpr(NTPMAX)
      integer i,j,ieflg,irec
      logical*1 svdotr            ! Not used in the routine
      integer*2 ielst(2,NENMAX),ielc,ialpha
      real*8, save :: tin,tnorm,dtnorm
      character*(chlen) inshfile,dataname
      real*8, save :: rhills(NTPMAX)
      real*8 r2hills(NTPMAX)
      real*8 a,e,q,as,es,qs

c----
c...  Executable code 

c...  set things up if this is the initial call
      if(i1stin.eq.0) then
	write(*,*) 'Enter name of shadow data file : '
	read(*,999) inshfile
999 	format(a)
	iul = 50
	call io_lyap_init(diro,inshfile,nbod,xh,yh,zh,vxh,vyh,vzh,
     &            xsh,ysh,zsh,vxsh,vysh,vzsh,dist0,lrsum,dtnorm,iul)
        tin = 0.0
        i1stin = 1
        i1stsh = 0
        tnorm = dtnorm

        call util_hills(nbod,mass,xsh,ysh,zsh,vxsh,vysh,vzsh,r2hills) 
        do i=2,nbod
          rhills(i) = sqrt(r2hills(i))
        end do

        if ((fopenstat(1:3).eq.'new')
     &       .or.(fopenstat(1:3).eq.'NEW')) then
          dataname = 'cp '//trim(inshfile)//' '//trim(diro)
          call system(dataname)

c        iflgchk = 0
c        j2rp2 = 0.0d0
c        j4rp4 = 0.0d0
c        call io_dump_pl_symba('copy.in',nbod,mass,xsh,ysh,zsh,
c     &     vxsh,vysh,vzsh,lclose,iflgchk,rpl,rhill,j2rp2,j4rp4)
      end if

c Initialize initial



        write(*,*) ' CONTINUE: '
      endif

      do i=2,nbod 
         call orbel_xv2aeq(xh(i),yh(i),zh(i),vxh(i),vyh(i),vzh(i),
     &     mass(1)+mass(i),ialpha,a,e,q)
        call orbel_xv2aeq(xsh(i),ysh(i),zsh(i),vxsh(i),vysh(i),vzsh(i),
     &     mass(1)+mass(i),ialpha,as,es,qs)
c        print*,i,'a=',sngl(a),sngl(as),'e=',sngl(e),sngl(es)
c        if (es.gt.1.) then
c          print*,'x',xh(i),xsh(i)
c          print*,'y',yh(i),ysh(i)
c          print*,'z',zh(i),zsh(i)
c          print*,'vx',vxh(i),vxsh(i)
c          print*,'vy',vyh(i),vysh(i)
c          print*,'vz',vzh(i),vzsh(i)
c          stop
c        end if
      end do
c...  now the normal step

      do i=2,nbod
         iecnt(i) = 0
         ielev(i) = -1
      enddo
      isenc = 0

c...  check for encounters
      irec = 0
      ielc = 0
c      print*,'check encounters'
      do j=2,nbodm
         do i=j+1,nbod
            ieflg = 0
            call symba5_chk(rhill,nbod,i,j,mass,xh,yh,zh,vxh,
     &           vyh,vzh,dt,irec,ieflg,svdotr)
            if(ieflg.ne.0) then
c               print*,'part enc',i,j,irec
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

c...  do a step
c      print*,'integrate solution'
      if(isenc.eq.0) then
         call symba5_step_helio(i1st,nbod,nbodm,mass,j2rp2,
     &     j4rp4,xh,yh,zh,vxh,vyh,vzh,dt)
         mergecnt=0
      else
         call symba5_step_interp(time,iecnt,ielev,nbod,
     &     nbodm,mass,rhill,j2rp2,j4rp4,lclose,rpl,xh,yh,zh,
     &     vxh,vyh,vzh,dt,mergelst,mergecnt,eoff,ielc,ielst,
     &     mtiny)
         i1st = 0  
      endif

c...  and the same for the shadow bodies (use the same temporary variables)
c...                                       when possible
      do i=2,nbod
         iecntsh(i) = 0
         ielevsh(i) = -1
      enddo
      isencsh = 0

c      print*,'check encounters for shadws'
c...  check for encounters
      irec = 0
      ielc = 0
      do j=2,nbodm
         do i=j+1,nbod
            ieflg = 0
            call symba5_chk(rhills,nbod,i,j,mass,xsh,ysh,zsh,vxsh,
     &           vysh,vzsh,dt,irec,ieflg,svdotr)
            if(ieflg.ne.0) then
               isencsh = 1
c               print*,'enc sha',i,j,irec
c               print*,'x',sngl(xh(i)),sngl(xh(j)),
c     &                    sngl(xsh(i)),sngl(xsh(j))
c               print*,'y',sngl(yh(i)),sngl(yh(j)),
c     &                    sngl(ysh(i)),sngl(ysh(j))
               iecntsh(i) = iecntsh(i) + 1
               iecntsh(j) = iecntsh(j) + 1
               ielevsh(i) = 0
               ielevsh(j) = 0
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
      
c...  do a step
c      print*,'inetgrate shadows'
      if(isencsh.eq.0) then
         call symba5_step_helio_lyap(i1stsh,nbod,nbodm,mass,j2rp2,
     &     j4rp4,xsh,ysh,zsh,vxsh,vysh,vzsh,dt)
         mergecntsh=0
      else
         call symba5_step_interp_lyap(time,iecntsh,ielevsh,nbod,
     &     nbodm,mass,rhills,j2rp2,j4rp4,lclose,rpl,xsh,ysh,zsh,
     &     vxsh,vysh,vzsh,dt,mergelstsh,mergecntsh,eoff,ielc,ielst,
     &     mtiny)
         i1stsh = 0  
      endif

      tin = tin + dt

c...  now do the lyap anal stuff

      if(tin .ge. tnorm) then 
c      print*,'lyap anal.......................................'
         do i=2,nbod
           call lyap_renorm(xh(i),yh(i),zh(i),
     &              vxh(i),vyh(i),vzh(i),xsh(i),ysh(i),zsh(i),
     &              vxsh(i),vysh(i),vzsh(i),dist0(i),lrsum(i))
           logpr(i) = dlog10(lrsum(i)/tin) 
         enddo

         i1stsh = 0
         call io_lyap_write(diro,iul,tin,logpr,nbod)

         tnorm = tnorm + dtnorm
      endif

      return
      end ! symba5_step_pl_lyap.f
c-----------------------------------------------------------

