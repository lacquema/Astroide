c************************************************************************
c                          IO_INIT_PARAM.F
c************************************************************************
c INIT_PARAM reads in the parameters for the integration. 

        subroutine io_init_param_old(infile,t0,tstop,dt,dtout,dtdump,
     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,
     &         outfile,fopenstat)

	include '../swift.inc'
	include 'io.inc'

c...    Input
	character*(chlen) infile

c...  Outputs: 
	integer iflgchk
	real*8 t0,tstop,dt
	real*8 dtout,dtdump
	real*8 rmin,rmax,rmaxu,qmin
        logical*2 lclose
	character*(chlen) outfile,fopenstat,dataname

c...  Internals
        logical*1 lflg(0:IO_NBITS-1)
        integer i,ierr

c-----
c...  Executable code 

c	write(*,*) 'Parameter data file is ',infile
        call io_open(7,infile,'old','formatted',ierr)
        rewind(7)
	read(7,*) t0,tstop,dt
c	write(*,*) 't0,tstop,dt : ',t0,tstop,dt
	read(7,*) dtout,dtdump
c	write(*,*) 'dtout,dtdump : ',dtout,dtdump
        read(7,*) (lflg(i),i=IO_NBITS-1,0,-1)

        iflgchk=0
        do i=0,IO_NBITS-1
           if(lflg(i)) then
              iflgchk = ibset(iflgchk,i)
           endif
        enddo

c        write(*,*) (lflg(i),i=IO_NBITS-1,0,-1),' = ',iflgchk

        if(btest(iflgchk,0) .and. btest(iflgchk,1))  then 
           write(*,*) ' SWIFT ERROR: in io_init_param_old:'
           write(*,*) '    Invalid logical flags '
           write(*,*) '    You cannot request that both a real and ',
     &                '       an integer binary file be written '
           call util_exit(1)
        endif

        if(btest(iflgchk,4))  then ! bit 4 is set
           read(7,*) rmin,rmax,rmaxu,qmin,lclose
c           write(*,*) 'rmin,rmax,rmaxu,qmin,lclose :',
c     &          rmin,rmax,rmaxu,qmin,lclose
        else
           rmin = -1.0
           rmax = -1.0
           rmaxu = -1.0
           qmin = -1.0
           lclose = .false.
        endif

 999       format(a)
       
        if(btest(iflgchk,0) .or. btest(iflgchk,1))  then 
           read(7,999) outfile
c           write(*,*) 'outfile : ', outfile
c           write(*,*) ' '
        endif

        read(7,999) fopenstat
        if(  (fopenstat(1:3).ne.'new') .and. 
     &       (fopenstat(1:3).ne.'NEW') .and.
     &       (fopenstat(1:7).ne.'unknown') .and. 
     &       (fopenstat(1:7).ne.'UNKNOWN') .and.
     &       (fopenstat(1:6).ne.'append') .and. 
     &       (fopenstat(1:6).ne.'APPEND') ) then
           write(*,*) ' SWIFT ERROR: in io_init_param_old:'
           write(*,*) '    Invalid status flag:',fopenstat(1:7),':'
           call util_exit(1)
        endif
        
	close(unit = 7)

	return
	end     ! io_init_param
c____________________________________________________________________________
c
c

