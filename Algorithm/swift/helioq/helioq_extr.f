c*************************************************************************
c                            HELIOQ_EXTR.F
c*************************************************************************
c This is extrapolates the values to a dt=h=0.
c
c             Input:
c                 iest          ==> the level of interation on dt
c                                       (int scalar)
c                 xest          ==> estimate of the independant var
c                                       (real scalar)
c                 yest          ==> estimate of the dependant var
c                                       (real array)
c                 nv             ==>  Number of independant var (int scalar)
c
c             Output:
c
c                 yz            ==> extrapolated values of the dependant var
c                                       (real array)
c                 dy            ==> extimated error on dy
c                                       (real array)
c Remarks: Based on numerical recipe routine pzextr.f
c Authors:  Hal Levison 
c Date:    2/8/2K
c Last revision: 

      subroutine helioq_extr(iest,xest,yest,yz,dy,nv)

      include '../swift.inc'
      include 'helioq.inc'

c...  Inputs Only: 
      integer iest,nv
      real*8 xest,yest(N6DBS)

c...  Outputs:
      real*8 yz(N6DBS),dy(N6DBS)

c...  Internals:
      real*8 x(ntrys),qcol(N6DBS,NUSE),d(N6DBS)
      real*8 delta,f1,f2,q
      integer j,k,m1,k1

      save x,d,qcol

c----
c...  Executable code 

      x(iest)=xest

      do j=1,nv
         dy(j)=yest(j)
         yz(j)=yest(j)
      enddo

      if(iest.eq.1) then

        do j=1,nv
          qcol(j,1)=yest(j)
       enddo

      else

         m1=min(iest,NUSE)

         do j=1,nv
            d(j)=yest(j)
         enddo

         do k1=1,m1-1
            delta=1.0d0/(x(iest-k1)-xest)
            f1=xest*delta
            f2=x(iest-k1)*delta
            do j=1,nv
               q=qcol(j,k1)
               qcol(j,k1)=dy(j)
               delta=d(j)-q
               dy(j)=f1*delta
               d(j)=f2*delta
               yz(j)=yz(j)+dy(j)
            enddo
         enddo

         do j=1,nv
            qcol(j,m1)=dy(j)
         enddo

      endif

      return
      end    ! helioq_extr.f
c---------------------------------------------------------
