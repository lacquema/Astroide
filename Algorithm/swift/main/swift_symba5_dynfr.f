c**********************************************************************
c		      SWIFT_SYMBA5.F
c**********************************************************************
c
c                 To run, need 2 input files. The code prompts for
c                 the file names, but examples are :
c
c                   parameter file like       param.in
c		    planet file like          pl.in
c
c  NOTE:  No test particles in this code and the massive bodies 
c         are dimensioned at NTPMAX
c
c Authors:  Hal Levison \& Martin Duncan
c Date:    11/21/96
c Last revision: 12/27/96

     
      include 'swift.inc'

      real*8 mass(NTPMAX),j2rp2,j4rp4
      real*8 xh(NTPMAX),yh(NTPMAX),zh(NTPMAX)
      real*8 vxh(NTPMAX),vyh(NTPMAX),vzh(NTPMAX)

      real*8 xht(1),yht(1),zht(1)       ! Dummy for the io
      real*8 vxht(1),vyht(1),vzht(1)
      integer ntp,istat(1)

      integer nbod,i1st,i,nbodm,nbodo
      integer iflgchk,iub,iuj,iud,iue,ium
      
      real*8 t0,tstop,dt,dtout,dtdump
      real*8 t,tout,tdump,tfrac,eoff
      real*8 rpl(NTPMAX),rhill(NTPMAX)

      real*8 rmin,rmax,rmaxu,qmin,mtiny
      logical*2 lclose 
      integer isenc,ihills
      integer mergelst(2,NTPMAX),mergecnt
      integer*2 iecnt(NTPMAX)

      character*(chlen) outfile,inparfile,inplfile,fopenstat

c-----
c...  Executable code 

      ntp = 0

c...  print version number
      call util_version

c Get data for the run and the test particles
      write(*,*) 'Enter name of parameter data file : '
      read(*,999) inparfile
      call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,
     &     iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)

c Prompt and read name of planet data file
      write(*,*) ' '
      write(*,*) 'Enter name of planet data file : '
      read(*,999) inplfile
 999  format(a)
      call io_init_pl_symba(inplfile,lclose,iflgchk,nbod,mass,
     &     xh,yh,zh,vxh,vyh,vzh,rpl,rhill,j2rp2,j4rp4)

      write(*,*) 'Enter the smallest mass to self gravitate :'
      read(*,*) mtiny
      write(*,*) ' mtiny = ',mtiny

c Initialize initial time and times for first output and first dump
      t = t0
      tout = t0 + dtout
      tdump = t0 + dtdump

      iub = 20
      iuj = 30
      iud = 40
      iue = 60
      ium = 21

c...    Do the initial io write
      if(btest(iflgchk,0))  then ! bit 0 is set
         call io_write_frame(t0,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
     &        xht,yht,zht,vxht,vyht,vzht,istat,outfile,iub,fopenstat)
         call io_write_mass(t0,nbod,mass,outfile,ium,fopenstat)
      endif
      if(btest(iflgchk,1))  then ! bit 1 is set
         call io_write_frame_r(t0,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
     &        xht,yht,zht,vxht,vyht,vzht,istat,outfile,iub,fopenstat)
         call io_write_mass_r(t0,nbod,mass,outfile,ium,fopenstat)
      endif
      if(btest(iflgchk,2))  then ! bit 2 is set
         eoff = 0.0d0
         call anal_energy_write(t0,nbod,mass,j2rp2,j4rp4,xh,yh,zh,vxh,
     &        vyh,vzh,iue,fopenstat,eoff)
      endif

c...  must initize discard io routine
      if(btest(iflgchk,4))  then ! bit 4 is set
         call io_discard_mass(0,t,0,mass(1),rpl(1),xh(1),yh(1),zh(1),
     &        vxh(1),vyh(1),vzh(1),iud,-1,fopenstat)
      endif

c...  Calculate the location of the last massive particle
      call symba5_nbodm(nbod,mass,mtiny,nbodm)

      ihills = 0
      i1st = 0

c***************here's the big loop *************************************
      write(*,*) ' ************** MAIN LOOP ****************** '

      do while ( (t .le. tstop) .and. (nbod.gt.1) )

         call symba5_step_pl(i1st,t,nbod,nbodm,mass,j2rp2,j4rp4,
     &        xh,yh,zh,vxh,vyh,vzh,dt,lclose,rpl,isenc,
     &        mergelst,mergecnt,iecnt,eoff,rhill,mtiny)

         t = t + dt

         if(btest(iflgchk,4))  then ! bit 4 is set
            nbodo = nbod
            call discard_massive5(t,dt,nbod,mass,xh,yh,zh,
     &           vxh,vyh,vzh,rmin,rmax,rmaxu,qmin,lclose,
     &           rpl,rhill,isenc,mergelst,mergecnt,
     &           iecnt,eoff,i1st)
            if(nbodo.ne.nbod) then
               call symba5_nbodm(nbod,mass,mtiny,nbodm)
            endif
         endif


c if it is time, output orb. elements, 
         if(t .ge. tout) then 

            if(btest(iflgchk,0))  then ! bit 0 is set
               call  io_write_frame(t,nbod,ntp,mass,xh,yh,zh,vxh,
     &              vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,outfile,
     &              iub,fopenstat)
               call io_write_mass(t,nbod,mass,outfile,ium,fopenstat)
            endif
            if(btest(iflgchk,1))  then ! bit 1 is set
               call  io_write_frame_r(t,nbod,ntp,mass,xh,yh,zh,vxh,
     &              vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,outfile,
     &              iub,fopenstat)
               call io_write_mass_r(t,nbod,mass,outfile,ium,fopenstat)
            endif

	    tout = tout + dtout
         endif

c If it is time, do a dump
         if(t.ge.tdump) then

            tfrac = (t-t0)/(tstop-t0)
            write(*,998) t,tfrac,nbod
 998        format(' Time = ',1p1e12.5,': fraction done = ',0pf5.3,
     &            ': Number of bodies =',i4)
            call io_dump_pl_symba('dump_pl.dat',nbod,mass,xh,yh,zh,
     &           vxh,vyh,vzh,lclose,iflgchk,rpl,rhill,j2rp2,j4rp4)
            call io_dump_param('dump_param.dat',t,tstop,dt,dtout,
     &           dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)
            tdump = tdump + dtdump

            if(btest(iflgchk,2))  then ! bit 2 is set
               call anal_energy_write(t,nbod,mass,j2rp2,j4rp4,
     &              xh,yh,zh,vxh,vyh,vzh,iue,fopenstat,eoff)
            endif
            
	  endif

	enddo
c********** end of the big loop from time 't0' to time 'tstop'

c Do a final dump for possible resumption later 

	call io_dump_pl_symba('dump_pl.dat',nbod,mass,xh,yh,zh,
     &            vxh,vyh,vzh,lclose,iflgchk,rpl,rhill,j2rp2,j4rp4)
	call io_dump_param('dump_param.dat',t,tstop,dt,dtout,
     &         dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)

        call util_exit(0)
        end    ! swift_symba5.f
c---------------------------------------------------------------------


c*************************************************************************
c                            SYMBA5_STEP_INTERP.F
c*************************************************************************
c
c             Input:
c                 time          ==> Current time (real scalar)
c                 iecnt         ==>  The number of objects that each planet 
c                                    is encountering (int*2 array)
c                 ielev         ==>  The level that this particle should go
c                                             (int*2 array)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 nbodm         ==>  Location of last massive body(int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 rhill         ==>  Radius of hill sphere (real array)
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
c                ielc           ==>  number of encounters (integer*2 scalar)
c                ielst          ==>  list of ecnounters (2D integer*2 array)
c                mtiny          ==>  Small mass  (real array)
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  final velocity in helio coord 
c                                       (real arrays)
c                 rpl           ==>  Recalculated physical size of a planet.
c                                    if merger happened (real array)
c                 nbod          ==>  Recalculated number of massive bodies 
c                                    if merger happened (int scalar)
c                 nbodm         ==>  Location of last massive body(int scalar)
c                 mass          ==>  Recalculated mass of bodies 
c                                    if merger happened (real array)
c                 mergelst      ==>  list of mergers (int array)
c                 mergecnt      ==>  count of mergers (int array)
c                 eoff          ==>  Energy offset (real scalar)
c Remarks: 
c Authors:  Hal Levison 
c Date:    11/21/96
c Last revision: 5/13/99
c temptative modification by morby to output close encounters between massive 
c bodies (nbod < nbodm)

      subroutine symba5_step_interp(time,iecnt,ielev,nbod,
     &     nbodm,mass,rhill,j2rp2,j4rp4,lclose,rpl,xh,yh,zh,vxh,
     &     vyh,vzh,dt,mergelst,mergecnt,eoff,ielc,ielst,mtiny)

      include 'swift.inc'
      include '../symba5/symba5.inc'

c...  Inputs Only: 
      real*8 mass(NTPMAX),dt,j2rp2,j4rp4,time,mtiny
      integer*2 iecnt(NTPMAX),ielev(NTPMAX)
      logical*2 lclose 
      integer*2 ielst(2,NENMAX),ielc

c...  Inputs and Outputs:
      integer nbod,nbodm
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 rpl(nbod),eoff
      real*8 rhill(nbod)

c...  Outputs
      integer mergelst(2,NTPMAX),mergecnt

c...  Internals:
      integer irec,ilevl(NTPMAX),i 
      real*8 dth
      real*8 axh(NTPMAX),ayh(NTPMAX),azh(NTPMAX)
      real*8 vxb(NTPMAX),vyb(NTPMAX),vzb(NTPMAX),msys
      real*8 ptxb,ptyb,ptzb            ! Not used here
      real*8 ptxe,ptye,ptze
      real*8 vri(NENMAX),vro(NENMAX)
      real*8 xr,yr,zr,vxr,vyr,vzr
      real*8 gm
      real*8 mass_orig(NTPMAX)
      integer ialpha,ie,j
      real*8 a,e,inc,capom,omega,capm,peri
      real*8 s,b,bmax,dU,U,bhill,Vesc2,fr
      PARAMETER (bmax=3.d0) ! maximal multiple of a Hill sphere used for drag 
      logical*1 svdotr(NENMAX)  ! Used by symba_step_recur


      save axh,ayh,azh     ! Note this !!
      save vxb,vyb,vzb     ! Note this !!

c----
c...  Executable code 

      dth = 0.5d0*dt

c--   compute r cdot v for each couple of massive - non massive encountering bodies
c     NOTE: j is massive by construction (see swift_symba5_pl.f)
      do ie=1,ielc
         i = ielst(1,ie)
         j = ielst(2,ie) 
         if(mass(i).lt.mtiny)then
            vri(ie)=(xh(i)-xh(j))*(vxh(i)-vxh(j))+
     $           (yh(i)-yh(j))*(vyh(i)-vyh(j))+
     $           (zh(i)-zh(j))*(vzh(i)-vzh(j))
         endif
      enddo
            

c...  Convert vel to bery to jacobi coords
      call coord_vh2b(nbod,mass,vxh,vyh,vzh,vxb,vyb,vzb,msys)

c...  Do the linear drift due to momentum of the Sun
      call helio_lindrift(nbod,mass,vxb,vyb,vzb,dth,
     &     xh,yh,zh,ptxb,ptyb,ptzb)

c...  Get the accelerations in helio frame. For each object
c...     only include those guys that it is not encountering with. 
      call symba5_getacch(nbod,nbodm,mass,j2rp2,j4rp4,
     &     xh,yh,zh,axh,ayh,azh,mtiny,ielc,ielst)

c...  Apply a heliocentric kick for a half dt 
      call kickvh(nbod,vxb,vyb,vzb,axh,ayh,azh,dth)

c..   Do a recursion step for full dt
      irec = -1
      call symba5_helio_drift(nbod,ielev,irec,mass,xh,
     &     yh,zh,vxb,vyb,vzb,dt)
c     If object i is not massive and 
c     involved in an encounter, put its mass to zero for the recursive procedure.
      do ie=1,ielc
         i = ielst(1,ie)
         j = ielst(2,ie) 
         if( (mass(i).lt.mtiny) .and. (mass(i).gt.(1.1d0*TINY) ) )then
            mass_orig(i)=mass(i)
            mass(i)=TINY
         endif
         if( (mass(j).lt.mtiny) .and. (mass(j).gt.(1.1d0*TINY) ) )then
            mass_orig(j)=mass(j)
            mass(j)=TINY
         endif
      enddo
      irec = 0
      do i=2,nbod
         ilevl(i) = 0
      enddo
      mergecnt = 0
      call symba5_step_recur(time,nbod,nbodm,mass,irec,ilevl,
     &     iecnt,ielev,rhill,xh,yh,zh,vxb,vyb,vzb,lclose,
     &     rpl,mergelst,mergecnt,dt,eoff,svdotr,ielc,ielst)
c     reset the correct values of the masses of the small bodies
      do ie=1,ielc
         i = ielst(1,ie)
         j = ielst(2,ie) 
         if(mass(i).lt.mtiny) then
            mass(i)=mass_orig(i)
         endif
      enddo
c...  Get the accelerations in helio frame. For each object
c...     only include those guys that it is not encountering with. 
      call symba5_getacch(nbod,nbodm,mass,j2rp2,j4rp4,
     &     xh,yh,zh,axh,ayh,azh,mtiny,ielc,ielst)

c...  Apply a heliocentric kick for a half dt 
      call kickvh(nbod,vxb,vyb,vzb,axh,ayh,azh,dth)

c...  Do the linear drift due to momentum of the Sun
      call helio_lindrift(nbod,mass,vxb,vyb,vzb,dth,
     &     xh,yh,zh,ptxe,ptye,ptze)

c...  convert back to helio velocities
      call coord_vb2h(nbod,mass,vxb,vyb,vzb,vxh,vyh,vzh)

c--   compute r cdot v for each couple of massive -- non massive encountering bodies
      do ie=1,ielc
         i = ielst(1,ie)
         j = ielst(2,ie)
         if(mass(i).lt.mtiny)then
            xr=(xh(i)-xh(j))
            vxr=(vxh(i)-vxh(j))
            yr=(yh(i)-yh(j))
            vyr=(vyh(i)-vyh(j))
            zr=(zh(i)-zh(j))
            vzr=(vzh(i)-vzh(j))
            vro(ie)=xr*vxr+yr*vyr+zr*vzr
            if(vri(ie).lt.0..and.vro(ie).ge.0.)then
c     a closest approach distance occurred. Compute planetocentric 
c     orbital elements.
               gm=mass(j)    ! +TINY
               call orbel_xv2el(xr,yr,zr,vxr,vyr,vzr,gm,
     &              ialpha,a,e,inc,capom,omega,capm)
               if(ialpha.gt.0)then
c     if hyperbolic, get the relative velocity at -infinity and impact parameter b
                  capm=-1.d6
                  call orbel_el2xv(gm,ialpha,a,e,inc,capom,omega,capm,
     &                xr,yr,zr,vxr,vyr,vzr)
                  U=sqrt(vxr**2+vyr**2+vzr**2)
                  s=-(xr*vxr+yr*vyr+zr*vzr)/U**2
                  xr=xr+vxr*s
                  yr=yr+vyr*s
                  zr=zr+vzr*s
                  b=sqrt(xr**2+yr**2+zr**2)
                  bhill=bmax*rhill(j)
                  if(b.lt.bhill)then
c     compute the fraction of the beam that hits the planet
                    Vesc2=2.d0*mass(j)/rpl(j)
                    fr=((rpl(j)/bhill)**2)*(1.d0+Vesc2/U**2)
                    if(fr.gt.1.d0)fr=1.d0
c     compute the velocity kick to the planet due to the particles that 
c     flyby but not collide
                    dU=2.d0*mass(i)*mass(j)/bhill**2/U**3
     %                    *(log(1.d0+((bhill**2)*(U**4)/(mass(j)**2)))-
     %                    log(1.d0+((fr*bhill**2)*(U**4)/(mass(j)**2))))
c     add the velocity kick due au bodies that do collide
                    dU=dU+U*fr*mass(i)/mass(j)
c     this formula comes from Binney and Tremaine 7-13b by substituting 
c     m*f(v_m)*d^3v_m*(v_m-v_M) with  mass(i)/(PI*bhill**2), G with 1, 
c     M+m with mass(j) and b_max with bhill 
c                     write(99,*)time,j,i,b/bhill,dU ! for later checks 
c     attribute the fraction of mass that hits the planet 
c     to the body j and take it away from tracer i
                     mass(j)=mass(j)+fr*mass(i)
                     mass(i)=mass(i)*(1.d0-fr)
                  else
                     dU=0.d0
                     U=1.d0
                  endif
               else
                  dU=0.0 
                  U=1.d0
               endif
c     now apply a velocity change dU to the massive guy j, directed as (vxr,vyr,vzr)
               vxh(j)=vxh(j)+dU*vxr/U
               vyh(j)=vyh(j)+dU*vyr/U
               vzh(j)=vzh(j)+dU*vzr/U
            endif
         endif
      enddo

      return
      end   ! symba5_step_interp
c---------------------------------------------------------------------

