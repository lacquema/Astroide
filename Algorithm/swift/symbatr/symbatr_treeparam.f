c**************************
c  SYMBATR_TREEPARAM.F
c**************************
c  Gets the parameters that are used in the
c  treecode part of swift_symbatr


      subroutine symbatr_treeparam(intreefile)

      include './inc/treedefs.f'

      character*80 intreefile
      
      open(unit=10, file=intreefile)
      read(10,*) EPS, USQUAD, THETA, treefreq

      print *,'treefreq=',treefreq


      return
      end




