c************************************************************************
c                         IO_DUMP_PL_HJS_TID.F
c************************************************************************
c Dumps the data for the massive bodies, HJS case 
c
c
c             Input:
c                 dplfile       ==>  Name of file to write to (character*(chlen))
c                 nbod          ==>  number of massive bodies (int scalar)
c                 tid         ==> 0 (no tide), 1 (tide) or 2 (averaged tide)
c                 oloc          ==>  Link matrix between bodies & orbits
c                         oloc(j,i)=1  : body #i is a satellite in orbit #j
c                         oloc(j,i)=-1 : body #i is a center in orbit #j
c                                         (2D integer array)
c                 mass          ==>  mass of bodies (real array)
c                 inert         ==>  moments of inertia (real array)
c                 Qlov          ==> tidal constants (real array)
c                 krpl5         ==> love number x radius^5 (real array)
c                 umat          ==>  Conversion matrix Gen. Jacobi => Bary
c                                         (2D real array)
c                 xj,yj,zj      ==>  initial position in Jacobi coord 
c                                    (real arrays)
c                 vxj,vyj,vzj   ==>  initial position in Jacobi coord 
c                                    (real arrays)
c                 rotx,roty,rotz ==> spin vectors (real arrays)
c                 lclose        ==> .true. --> discard particle if it gets 
c                                    too close to a planet. Read in that 
c                                    distance in io_init_pl
c                                      (logical*2 scalar)
c                 iflgchk       ==>  bit 5 set ==>  include J2 and J4 terms
c                 rplsq         ==>  min distance^2 that a tp can get from pl
c                                    (real array)
c
c Remarks: based on io_dump_pl_hjs.f
c Authors:  Hervé Beust
c Date:    Jun 8, 2008; revised Jun 19

	subroutine io_dump_pl_hjs_tid(dplfile,nbod,tid,oloc,
     &     mass,umat,xj,yj,zj,vxj,vyj,vzj,rotx,roty,rotz,inert,
     &     krpl5,Qlov,lclose,iflgchk,rplsq)

	include '../swift.inc'
	include '../io/io.inc'

c...    Input
        integer nbod, iflgchk, oloc(NPLMAX,NPLMAX)
	real*8 mass(nbod),rplsq(nbod)
	real*8 xj(nbod),yj(nbod),zj(nbod)
	real*8 vxj(nbod),vyj(nbod),vzj(nbod)
        real*8 rotx(nbod),roty(nbod),rotz(nbod)
        real*8 umat(NPLMAX,NPLMAX)
        real*8 inert(nbod),krpl5(nbod),Qlov(nbod)
	character*(*) dplfile
        logical*2 lclose
        integer tid(nbod)

c...   Internal
	integer j,ierr,i
	real*8 xb(NPLMAX),yb(NPLMAX),zb(NPLMAX)
	real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)
        real*8 rpl,rgyr,klov

c-----
c...  Executable code      

        call io_open(7,dplfile,'unknown','formatted',ierr)

	write(7,*) nbod

        call coord_g2b(nbod,umat,mass,xj,yj,zj,vxj,vyj,vzj,
     &                                   xb,yb,zb,vxb,vyb,vzb) 

	do j=1,nbod
           if(lclose) then
              rpl = sqrt(rplsq(j))
              if (rplsq(j).ne.0.d0) then
                rgyr = sqrt(inert(j)/mass(j)/rplsq(j))
                klov = krpl5(j)/rplsq(j)/rplsq(j)/rpl
              else
                rgyr = 0.d0
                klov = 0.d0
              end if
              write(7,123) mass(j),rpl,klov,Qlov(j),rgyr
           else
              write(7,123) mass(j)
           endif
           write(7,123) xb(j),yb(j),zb(j)
           write(7,123) vxb(j),vyb(j),vzb(j)
           write(7,123) rotx(j),roty(j),rotz(j)
	enddo
 123    format(3(1p1e23.16,1x))

        do j = 2,nbod
          write(7,*)(oloc(j,i),i=1,nbod)   
        end do

        write(7,*)tid(2:nbod)
	close(unit = 7)
	return
	end    ! io_dump_pl_hjs.f
c--------------------------------------------------------------------------

