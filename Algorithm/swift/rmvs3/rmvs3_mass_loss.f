c*************************************************************************
c                            RMVS3_MASS_LOSS.F
c*************************************************************************

      subroutine rmvs3_mass_loss(time,indice,nbod,mass,
     &                           TabulatedTime,TabulatedMass)

      include '../swift.inc'
      include '../rmvs/rmvs.inc'

c...  Inputs Only: 
      integer nbod,indice
      real*8 time
      real*8 TabulatedTime(INDICEMAX),TabulatedMass(INDICEMAX)

c...  Inputs and Outputs:
      real*8 mass(NPLMAX)

c...  Internals


c----
c...  Executable code 

      do while(time.ge.TabulatedTime(indice+1))
         indice = indice + 1
      enddo

      mass(1)=TabulatedMass(indice)-TabulatedMass(indice+1)
      mass(1)=mass(1)/(TabulatedTime(indice)-TabulatedTime(indice+1))
      mass(1)=mass(1)*(time-TabulatedTime(indice))
      mass(1)=mass(1)+TabulatedMass(indice)
      mass(1)=mass(1)*3.9478417604357354d1
      end   ! step_enc
c------------------------------------------------------------------------






