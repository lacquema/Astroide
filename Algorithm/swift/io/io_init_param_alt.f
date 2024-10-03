c************************************************************************
c                          IO_INIT_PARAM.F
c************************************************************************
c INIT_PARAM reads in the parameters for the integration. 

        subroutine io_init_param_alt(infile,t0,tstop,dt,dtout,dtdump,
     &         diskMass,
     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,diro,
     &         outfile)

	include '../swift.inc'
	include 'io.inc'

c...    Input
	character*120 infile

c...  Outputs: 
	integer iflgchk
	real*8 t0,tstop,dt
	real*8 dtout,dtdump,diskMass
	real*8 rmin,rmax,rmaxu,qmin
        logical*2 lclose
	character*120 diro,outfile

c...  Internals
        logical*1 lflg(0:IO_NBITS-1)
        integer i,ierr

c-----
c...  Executable code 

	write(*,*) 'Parameter data file is ',infile
        call io_open(7,infile,'old','formatted',ierr)
        rewind(7)
	read(7,*) t0,tstop,dt
	write(*,*) 't0,tstop,dt : ',t0,tstop,dt
	read(7,*) dtout,dtdump
	write(*,*) 'dtout,dtdump : ',dtout,dtdump
        read(7,*) DiskMass
        write(*,*) 'disk mass: ',DiskMass 
        read(7,*) (lflg(i),i=IO_NBITS-1,0,-1)

        
        iflgchk=0
        do i=0,IO_NBITS-1
           if(lflg(i)) then
              iflgchk = ibset(iflgchk,i)
           endif
        enddo

        write(*,*) (lflg(i),i=IO_NBITS-1,0,-1),' = ',iflgchk

        if(btest(iflgchk,0) .and. btest(iflgchk,1))  then 
           write(*,*) ' SWIFT ERROR: in io_init_param:'
           write(*,*) '    Invalid logical flags '
           write(*,*) '    You cannot request that both a real and ',
     &                '       an integer binary file be written '
           call util_exit(1)
        endif

        if(btest(iflgchk,4))  then ! bit 4 is set
           read(7,*) rmin,rmax,rmaxu,qmin,lclose
           write(*,*) 'rmin,rmax,rmaxu,qmin,lclose :',
     &          rmin,rmax,rmaxu,qmin,lclose
        else
           rmin = -1.0
           rmax = -1.0
           rmaxu = -1.0
           qmin = -1.0
           lclose = .false.
        endif

 999       format(a)
       

           read(7,999) diro
c           call system(dataname)
           read(7,999) outfile
           outfile=trim(diro)//outfile 
           write(*,*) 'outfile : ',outfile
           write(*,*) ' '

         
	close(unit = 7)

	return
	end     ! io_init_param
c____________________________________________________________________________
c
c

