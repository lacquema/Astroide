c*************************************************************************
c                            HELIOQ_BS.F
c*************************************************************************
c This is the BS integrator for helioq
c
c             Input:
c                 iwhich        ==> determines which H to integrate
c                                       (int scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                  x            ==> initial value independent variable 
c                                        (real scalar)
c 	           h0           ==> stepsize  (real scalar)
c 	           y            ==> initial value dependent variables  
c                                     (real array)
c                 rtrans        ==>  solar encounter radii (real array)
c             Output:
c 	          y  ==> final value dependent variables  (real array)
c                 x  ==> final value independent variable (real scalar)
c 	          h0 ==> recommended stepsize for the next call (real scalar)
c              hnext ==> suggested timestep for the next time (real scalar)
c
c Remarks: Based on numerical recipe routine
c Authors:  Hal Levison 
c Date:    2/2/2K
c Last revision:  2/18/2K 

      subroutine helioq_bs(iwhich,nbod,mass,x,h0,y,rtrans,hnext)

      include '../swift.inc'
      include 'helioq.inc'

c...  Inputs Only: 
      integer nbod,iwhich
      real*8 mass(nbod)
      real*8 rtrans(2)

c...  Inputs and Outputs:
      real*8 y(N6DBS),x,h0,hnext

c...  Internals:
      real*8 yseq(N6DBS),dy(N6DBS),yscal(N6DBS)
      real*8 yerr(N6DBS),ysav(N6DBS),dysav(N6DBS)
      integer nseq
      real*8 xsav,xest,errmax,hdid,h,sum
      integer i,n,j,ic
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c      integer ib
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c----
c...  Executable code 

c...  set things up for the BS step

      n = 6*(nbod-1)

      if(iwhich.eq.1) then
         call helioq_bs_der1(rtrans,nbod,mass,y,dy)
      else
         call helioq_bs_der2(rtrans,nbod,mass,y,dy)
      endif

      do i=1,2*(nbod-1)
         sum = 0.0d0
         do j = 1,3
            ic = 3*(i-1) + j
            sum = sum + abs(y(ic))+abs(h0*dy(ic)) + TINY
         enddo
         sum = sum/3.0d0
         do j = 1,3
            ic = 3*(i-1) + j
            yscal(ic) = sum
         enddo
      enddo

      h=h0
      xsav=x
      do i=1,n
        ysav(i)=y(i)
        dysav(i)=dy(i)
      enddo

      do while(.true.)
         do i=1,NTRYS
            nseq = 2*i
            call helioq_mmid(ysav,dysav,nbod,xsav,h,nseq,
     &           yseq,iwhich,rtrans,mass)
            xest=(h/nseq)**2

            call helioq_extr(i,xest,yseq,y,yerr,n)

            errmax=0.
            do j=1,n
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c               if(errmax.lt.abs(yerr(j)/yscal(j))) then
c                  ib = j
c               endif
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
               errmax=max(errmax,abs(yerr(j)/yscal(j)))
            enddo
            errmax=errmax/EPS
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c            write(55,*) i,ib,errmax,y(ib),dy(ib),yerr(ib),
c     &              yscal(ib)
c            ib = 4693
c            write(56,*) i,y(ib),dy(ib),yerr(ib),yscal(ib),
c     &           yseq(ib),ysav(ib)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

            if(errmax.lt.1.0d0) then
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c               write(*,*) 'Here BS1',ib,y(ib),dy(ib),yerr(ib),
c     &              yscal(ib),errmax
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
               x=x+h
               hdid=h
               if(i.eq.NUSE) then
                  hnext=h*SHRINK
               else if(i.eq.NUSE-1) then
                  hnext=h*GROW
               else
                  hnext=(h*2*(NUSE-1))/nseq
               endif
               h0 = hdid
               return           ! <========= NOTE RETURN
            endif
         enddo

         h=0.25*h/2**((NTRYS-NUSE)/2)
         if(x+h.eq.x) then
            write(*,*) 'ERROR: Step size underflow in helioq_bs.'
            call util_exit(1)
         endif
      enddo
c
      end      !  helioq_bs
c------------------------------------------------------------------------
