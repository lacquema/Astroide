c converts binary file to ascii file

	include 'swift.inc'

	real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

	real*8 mass(NPLMAX),dr,peri
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

	integer istat(NTPMAX,NSTAT)
        real*8 rstat(NTPMAX,NSTATR)
	integer nbod,ntp,ierr,ifol,istep
	integer iflgchk,iu,nleft,i,id
        integer io_read_hdr,io_read_line
        integer io_read_hdr_r,io_read_line_r

	real*8 t0,tstop,dt,dtout,dtdump
	real*8 t,tmax,DumpTime,DumpTimeStep

	real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
        logical*2 lclose
        real*8 a,e,inc,capom,omega,capm,j2rp2,j4rp4
        real*8 elh,elk,elp,elq,apo

	character*(chlen) outfile,inparfile,inplfile,intpfile,fopenstat

c Get data for the run and the test particles
	read(*,999) inparfile
        call io_open(7,inparfile,'old','formatted',ierr)
        rewind(7)
	read(7,*) t0,tstop,dt
	read(7,*) dtout,dtdump
        read(7,*) (lflg(i),i=IO_NBITS-1,0,-1)

        iflgchk=0
        do i=0,IO_NBITS-1
           if(lflg(i)) then
              iflgchk = ibset(iflgchk,i)
           endif
        enddo

        if(btest(iflgchk,0) .and. btest(iflgchk,1))  then 
           write(*,*) ' SWIFT ERROR: in io_init_param:'
           write(*,*) '    Invalid logical flags '
           write(*,*) '    You cannot request that both a real and ',
     &                '       an integer binary file be written '
           call util_exit(1)
        endif

        if(btest(iflgchk,4))  then ! bit 4 is set
           read(7,*) rmin,rmax,rmaxu,qmin,lclose
     &          rmin,rmax,rmaxu,qmin,lclose
        else
           rmin = -1.0
           rmax = -1.0
           rmaxu = -1.0
           qmin = -1.0
           lclose = .false.
        endif

	close(unit = 7)


c Prompt and read name of planet data file
c	write(*,*) 'Enter name of planet data file : '
	read(*,999) inplfile
999 	format(a)
        call io_open(7,inplfile,'old','formatted',ierr)

c Read number of planets
	read(7,*) nbod

        if(nbod.gt.NPLMAX) then
           write(*,*) ' SWIFT ERROR: in io_init_pl: '
           write(*,*) '   The number of massive bodies,',nbod,','
           write(*,*) '   is too large, it must be less than',NPLMAX
           call util_exit(1)
        endif

c For each planet read mass, 
c and helioc. position and vel .
        if(btest(iflgchk,5))  then ! bit 5 is set
           read(7,*) mass(1),j2rp2,j4rp4
        else
           read(7,*) mass(1)
           j2rp2 = 0.0d0
           j4rp4 = 0.0d0
        endif
        read(7,*) xh(1),yh(1),zh(1)
        read(7,*) vxh(1),vyh(1),vzh(1)
        rplsq(1) = 0.0d0

        if(  (xh(1).ne.0.0d0) .or.
     &       (yh(1).ne.0.0d0) .or.
     &       (zh(1).ne.0.0d0) .or.
     &       (vxh(1).ne.0.0d0) .or.
     &       (vyh(1).ne.0.0d0) .or.
     &       (vzh(1).ne.0.0d0) ) then
           write(*,*) ' SWIFT ERROR: in io_init_pl: '
           write(*,*) '   Input MUST be in heliocentric coordinates '
           write(*,*) '   Position and Vel. of Massive body 1 .ne. 0'
           call util_exit(1)
        endif

	do j=2,nbod
           if(lclose) then
              read(7,*) mass(j),rpl
              rplsq(j) = rpl*rpl
           else
              read(7,*) mass(j)
           endif
	  read(7,*) xh(j),yh(j),zh(j)
	  read(7,*) vxh(j),vyh(j),vzh(j)
	enddo
	close(unit = 7)
c Get data for the run and the test particles
c	write(*,*) 'Enter name of test particle data file : '
	read(*,999) intpfile
        
        ntp = 0     
        call io_open(7,infile,'old','formatted',ierr)
        rewind(7)
	read(7,*) ntp

        if(ntp.gt.NTPMAX) then
           write(*,*) ' SWIFT ERROR: in io_init_tp: '
           write(*,*) '   The number of test bodies,',ntp,','
           write(*,*) '   is too large, it must be less than',NTPMAX
           call util_exit(1)
        endif

        if(ntp.eq.0) then
           close(unit = 7)
           end
        endif               ! <===== NOTE


	  do  i=1,ntp
            read(7,*) beta 
	    read(7,*) xht(i),yht(i),zht(i)
	    read(7,*) vxht(i),vyht(i),vzht(i)
	    read(7,*) (istat(i,j),j=1,NSTAT)
	    read(7,*) (rstat(i,j),j=1,NSTATR)
	  enddoiu = 20

        dr = 180.0/PI

        if(btest(iflgchk,0)) then
c           write(*,*) ' Reading an integer*2 binary file '
        else if(btest(iflgchk,1)) then
c           write(*,*) ' Reading an real*4 binary file '
        else
c           write(*,*) ' ERROR: no binary file format specified '
c           write(*,*) '        in param file '
           stop
        endif

        read(*,*) DumpTimeStep
c        write(*,*) ' Input the particle number to follow '
c        read(*,*) ifol
c        write(*,*) ' Following particle ',ifol

        open(unit=iu, file=outfile, status='old',form='unformatted')
c        open(unit=7,file='follow.out')
c        open(unit=8,file='follow_hkpq.out')
c
c        write(*,*) '1 2 3  4    5     6    7    8    9 '
c        write(*,*) 't,a,e,inc,capom,omega,capm,peri,apo'

        tmax = t0
        dumpTime = t0
 1      continue
             if(btest(iflgchk,0))  then ! bit 0 is set
                ierr = io_read_hdr(iu,t,nbod,nleft) 
             else
                ierr = io_read_hdr_r(iu,t,nbod,nleft) 
             endif


             write(*,*) 'time: ',t

             if(ierr.ne.0) goto 2

             if (t.ne.dumpTime) then
               do i=2,nbod
                if(btest(iflgchk,0))  then ! bit 0 is set
                   ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm) 
                else
                   ierr = io_read_line_r(iu,id,a,e,inc,capom,omega,capm) 
                endif
               enddo
               do i=1,nleft
                if(btest(iflgchk,0))  then ! bit 0 is set
                   ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm) 
                else
                   ierr = io_read_line_r(iu,id,a,e,inc,capom,omega,capm) 
                endif
               enddo
               goto 1   
             else
               dumpTime = t + DumpTimeStep*dtout
             endif 

             istep = 0
             do i=2,nbod
                if(btest(iflgchk,0))  then ! bit 0 is set
                   ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm) 
                else
                   ierr = io_read_line_r(iu,id,a,e,inc,capom,omega,capm) 
                endif
                if(ierr.ne.0) goto 2
c                if(id.eq.ifol) then
                   istep = 1
                   elh = e*cos(omega+capom)
                   elk = e*sin(omega+capom)
                   elp = sin(inc/2.0)*cos(capom)
                   elq = sin(inc/2.0)*sin(capom)
                   inc = inc*dr
                   capom = capom*dr
                   omega = omega*dr
                   capm = capm*dr
                   peri = a*(1.0d0-e)
                   apo = a*(1.0d0+e)
                   write(*,*) id,a,e,inc,capom,omega,capm,peri,apo
                    tmax = t
c               endif
             enddo

             do i=1,nleft
                if(btest(iflgchk,0))  then ! bit 0 is set
                   ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm) 
                else
                   ierr = io_read_line_r(iu,id,a,e,inc,capom,omega,capm) 
                endif
                if(ierr.ne.0) goto 2
c                if(id.eq.ifol) then
                   istep = 1
                   elh = e*cos(omega+capom)
                   elk = e*sin(omega+capom)
                   elp = sin(inc/2.0)*cos(capom)
                   elq = sin(inc/2.0)*sin(capom)
                   inc = inc*dr
                   capom = capom*dr
                   omega = omega*dr
                   capm = capm*dr
                   peri = a*(1.0d0-e)
                   apo = a*(1.0d0+e)
                   write(*,*) id,a,e,inc,capom,omega,capm,peri,apo
                   tmax = t
c                endif
             enddo
c             if(istep.eq.0) goto 2     ! did not find particle this times step

        goto 1

 2      continue

        stop
        end
