c*************************
cSYMBATR_INTERFACE.F
c*************************
c This subroutine acts as an interface between 
c Symba and the treecode stuff



      subroutine symbatr_interface(xh,yh,zh,axh,ayh,azh,
     &     mass,nbod,nbodm)

      include '../swift.inc'
      include 'symbatr.inc'
      include './inc/treedefs.f'

c...  Inputs only:
      integer nbod,nbodm
      real*8 mass(nbod)
      real*8 xh(nbod),yh(nbod),zh(nbod)

c...  Inputs and outputs:
      real*8 axh(nbod),ayh(nbod),azh(nbod)

c...  Internals:
      integer i,j
c      print *,'in symbatr_interface'
c     Determine number of bodies to be processed
c     Note that nbodm is actually one less than
c     the number of non-massive bodies
      NBODY = nbod-nbodm

c     Determine whether treecode needs to be called,
c     or whether previous accelerations can be re-used:

c      print *,'accnum, callcount, treefreq ='
c      print *,accnum,callcount,treefreq

      if((accnum.eq.2).or.(tree_i1st.eq.1)) then
         if(callcount.eq.treefreq) then


c      if(((accnum.eq.2).and.(callcount.eq.treefreq))
c     &     .or.(tree_i1st.eq.1)) then
c         print *,'calculating tree accell'
         

c     Initialize box coordinate system:
         RSIZE = 4.0



c     Load positions and masses into treecode's arrays:


         do i = 1, NBODY
            POS(i,1) = xh(nbodm + i)
            POS(i,2) = yh(nbodm + i)
            POS(i,3) = zh(nbodm + i)
            MASS_T(i) = mass(nbodm + i)
         enddo



c     TEST
c      print *,'calling treecode'

         call symbatr_tree_MKTREE
c     TEST: short-circuit that mother
c      goto 100
         call symbatr_tree_ACCEL

c 100  continue


      endif
      callcount = callcount + 1
      if(callcount.gt.treefreq) callcount = 1
      endif

c     Add the treecode-calculated accelerations  to those
c     calculated in Symba for the M<mtiny bodies:
c     TEST
c
c      j=1
      write(50,100) time_global, ACC(1,1),ACC(1,2),ACC(1,3)
      write(51,100) time_global, ACC(10,1),ACC(10,2),ACC(10,3)
      write(52,100) time_global, ACC(50,1),ACC(50,2),ACC(50,3)
      write(53,100) time_global, ACC(100,1),ACC(100,2),ACC(100,3)
      write(54,100) time_global, ACC(500,1),ACC(500,2),ACC(500,3)

 100  format(1x,e15.7,1x,e15.7,1x,e15.7,1x,e15.7)




      do i = 1, NBODY
c     TEST
c         print *,'b: axh(',i+nbodm,')=',axh(i+nbodm)
c         print *,'ACC(',i,'1)=',ACC(i,1)
         axh(i + nbodm) = axh(i + nbodm) + ACC(i,1)
c         print *,'a: axh(',i+nbodm,')=',axh(i+nbodm)
         ayh(i + nbodm) = ayh(i + nbodm) + ACC(i,2)
         azh(i + nbodm) = azh(i + nbodm) + ACC(i,3)
      enddo
      if(tree_i1st.eq.1) tree_i1st = 0
c      callcount = callcount + 1
c      if(callcount.gt.treefreq) callcount = 1

c     TEST
c      print *,'RETURNING'

      return

      end

