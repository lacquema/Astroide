c**********************************************************************
c		            symba6_radzone.f
c**********************************************************************
c
c  Purpose:  This function takes the cylindrical radial coordinate 
c            and returns the integer id number for the corresponding
c            radial zone that it resides in.
c
c  Input:
c            rho   -->   cylindrical radial coordinate (real scalar)
c            rbnd1 -->   outer boundary of inner most radial zone  
c                            (real scalar)
c            dr    -->   width of radial zones  (real scalar)
c            nrz   -->   number of radial zones  (real scalar)
c
c  Output:
c         symba6_radzone ==> radial zone (int scalar)
c
c Remarks: 
c Author:  cba
c Date:    5/4/99
c Last revision: 

      integer function symba6_radzone(rho,rbnd1,dr,nrz)

      include '../swift.inc'
      include 'symba6.inc'

c...  Inputs Only: 
      integer nrz
      real*8  rho,rbnd1,dr

c...  Internals
      integer ir

c----
c...  Executable code 

      ir = (rho - rbnd1)/dr

      if (ir .lt. 1) then
          symba6_radzone = 1
      elseif (ir .gt. nrz) then
          symba6_radzone = nrz
      else
          symba6_radzone = ir
      endif

      return
      end   ! symba6_radzone
c----------------------------------------------------------------------


