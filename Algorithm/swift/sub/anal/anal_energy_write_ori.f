c*************************************************************************
c                          ANAL_ENERGY_WRITE_ORI.F
c*************************************************************************
c Writes the energy of the total system (massive bodies) wrt time.
c
c      Input:
c            t             ==>  current time
c            nbod          ==>  number of massive bodies (int scalar)
c            mass          ==>  mass of bodies (real array)
c            j2rp2         ==>  scaled value of j2 moment (real*8 scalar)
c            j4rp4         ==>  scaled value of j4 moment (real*8 scalar)
c            xh,yh,zh      ==>  current position in helio coord 
c                               (real arrays)
c            vxh,vyh,vzh   ==>  current velocity in helio coord 
c                               (real arrays)
c            iu            ==>  unit to write to (int scalar)
c            fopenstat     ==>  The status flag for the open 
c                                statements of the output files.  
c                                      (character*(chlen))
c            eoff          ==> An energy offset that is added to the energy
c                                      (real*8 scalar)
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    3/4/93
c Last revision: Oct 10, 2006

      subroutine anal_energy_write_ori(t,nbod,mass,j2rp2,j4rp4,xh,yh,zh,
     &     vxh,vyh,vzh,iu,diro,fopenstat,eoff)

      include '../swift.inc'

c...  Inputs: 
      integer nbod,iu
      real*8 mass(nbod),t,j2rp2,j4rp4,eoff
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      character*(chlen) fopenstat,diro

c...  Internals
      integer i1st
      real*8 energy,eltot(3),ke,pot

      data i1st/0/
      save i1st

c----
c...  Executable code 

c Compute and print initial ke,pot,energy and ang. mom.
      call anal_energy(nbod,mass,j2rp2,j4rp4,xh,yh,zh,
     &           vxh,vyh,vzh,ke,pot,energy,eltot)

      energy = energy + eoff

      call io_energy_write(i1st,t,energy,eltot,iu,diro,fopenstat)

      if(i1st.eq.0) then
         i1st=1
      endif

      return
      end     !anal_energy_write
c--------------------------------------------------------------------------

