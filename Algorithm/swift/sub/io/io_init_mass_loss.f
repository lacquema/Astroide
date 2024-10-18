c**********************************************************************
c			IO_INIT_MASS_LOSS.F
c**********************************************************************

	subroutine io_init_mass_loss(time,indice,TabulatedTime,
     &                                              TabulatedMass)

	include '../swift.inc'
	include 'io.inc'


c...    Output
	real*8 time

c...   Internal
	integer i,imax,indice
	real*8 TabulatedTime(INDICEMAX),TabulatedMass(INDICEMAX)	

c-----
c...  Executable code 
      open(unit=1000, file='mass.dat', status='old',
     &        form='formatted')
      rewind(1000)
	read(1000,*) imax
	i=1
        do while(i.le.imax)
         read(1000,*) TabulatedTime(i),TabulatedMass(i)
	 i=i+1
	enddo
        indice=1 
	do while(time.ge.TabulatedTime(indice+1))
         indice = indice + 1
        enddo

	end    ! io_init_tp.f
c-----------------------------------------------------------------

