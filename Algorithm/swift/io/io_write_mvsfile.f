c************************************************************************
c                          IO_WRITE_MVSFILE
c************************************************************************
c IO_WRITE_MVSFILE dumps out the mvs files for further continuation 
c       
c      Input:
c            diro,genfile  ==>  Directory & generic file (character*(chlen))
c
c Authors:  Hervé Beust
c Date:    Feb 14, 2023
c 

	subroutine io_write_mvsfile(diro,genfile)	

	include '../swift.inc'
	include 'io.inc'

c...   Inputs
	character*(*) diro,genfile

c...  Internals
        character*(chlen) mvsfile
        integer ierr

c-----
c...  Executable code 

        mvsfile = trim(diro)//'/mvs.in'

        call io_open(7,mvsfile,'unknown','formatted',ierr)

        write(7,'(a)')trim(genfile)
        write(7,'(a)')trim(diro)//'/dump_param.dat'
        write(7,'(a)')trim(diro)//'/dump_pl.dat'
        write(7,'(a)')trim(diro)//'/dump_tp.dat'
	write(7,'(a)')'1d-8'

	close(unit = 7)

	return
	end     ! io_write_mvsfile
c____________________________________________________________________________
c
c

