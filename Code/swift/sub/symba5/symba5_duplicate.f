c*************************************************************************
c                            SYMBA5_DUPLICATE_.F
c*************************************************************************
c This subroutine checks to see if a massive body enter at less than 1.
c times the periastre of the planet. And then it duplicate the particle to
c have a significative effect of the particles if they comes from an
c asteroid belt
c
c             Input:
c                 time          ==>  current time (real scalar)
c                 dt            ==>  time step  (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>   position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>   pl vel in helio coord 
c                                    (real arrays)
c                 Mmin,Mmax     ==>  maximum and min mass to create
c                                     random particles
c             Output:
c                 nbod          ==>  recalculated number of massive bodies 
c                                       (int scalar)
c                 mass          ==>  recalculated mass of bodies (real array)
c                 xh,yh,zh      ==>  recalculated position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  recalculated pl vel in helio coord 
c                                    (real arrays)
c
c
c Remarks: 
c
c Authors:  Herve Beust & Julien Vandeportal
c Date:    10 oct. 2009
c Last revision: 21 oct. 2009

C-----------------------------------------------------------------
C         Generateur aleatoire Gaussien normalise
C------------------------------------------------------------------
C
        SUBROUTINE GASDEV(S)

        IMPLICIT NONE

        REAL*8          S,R,V(2),S1,S2   ! Parametres internes
        INTEGER*4       J                ! aleatoire 1 ou 2


        CALL RANDOM_NUMBER(S1)
        J=1+FLOOR(2.d0*S1)
        R=2.d0
        DO WHILE(R.GE.1.d0)
          CALL RANDOM_NUMBER(S1)
          CALL RANDOM_NUMBER(S2)
          V(1)=2.d0*S1-1.d0
          V(2)=2.d0*S2-1.d0
          R=V(1)*V(1)+V(2)*V(2)
        END DO
        S = V(J)*SQRT(-2.d0*LOG(R)/R)
        END
C-----------------------------------------------------------------
C         Generateur de direction aleatoire
C------------------------------------------------------------------
C
        SUBROUTINE RANDOM_DIRECTION(theta,phi)

        IMPLICIT NONE

        REAL*8          theta,phi   ! colatitude,longitude
        REAL*8          s                   ! aleatoire 
        REAL*8, PARAMETER :: PI=3.14159265358979323846d-8

        CALL RANDOM_NUMBER(s)
        phi = 2.d0*Pi*s
        call random_number(s)
        theta = ACOS(2.0d0*s-1.d0)
        END

c---------------------------------------------------------------------
C
c...  give a random mass between Mmin and Mmax with the classical
c...  distribution in -3.5

      REAL*8 FUNCTION rand_mass(Mmin,Mmax)

      implicit none 

      REAL*8 Mmin,Mmax,s
      real*8, parameter :: z = -5.d0/6.d0, uz = -6.d0/5.d0

      call random_number(s)
      rand_mass = ((Mmax**z-Mmin**z)*s+Mmin**z)**uz

      RETURN

      END

c -----------------------------------------------------------------


      SUBROUTINE symba5_duplicate(nbod,nbodm,iouter,msmall,mtiny,
     &                 mass,xh,yh,zh,vxh,vyh,vzh,rpl,rhill)

      include '../swift.inc'

c...  Inputs: 
      real*8 time,dt,msmall,mtiny
      integer iouter,nbodm

c...  Input and Output

      integer nbod
      real*8 mass(NTPMAX),xh(NTPMAX),yh(NTPMAX),zh(NTPMAX)
      real*8 vxh(NTPMAX),vyh(NTPMAX),vzh(NTPMAX)
      real*8 rpl(NTPMAX),rhill(NTPMAX)

c...  internal

      integer ialpha,ic
      real*8 a,e,q,apoiouter
      logical*1 needdup
      real*8 xtemp,ytemp,ztemp,vxtemp,vytemp,vztemp,masstemp,mbottom,s
      real*8 rpltemp,rhilltemp,vtemp,dv,theta,phi,sth,cth,sphi,cphi
      real*8 xadd(NTPMAX),yadd(NTPMAX),zadd(NTPMAX),massadd(NTPMAX)
      real*8 vxadd(NTPMAX),vyadd(NTPMAX),vzadd(NTPMAX),mtot
      real*8 rhilladd(NTPMAX),rpladd(NTPMAX)
      real*8 rand_mass
      real*8, parameter :: third = 1.d0/3.d0
      integer isperih(NTPMAX)
      integer nbodtemp,nbcopy,i,j,k,l

      mbottom = msmall**2/mtiny

      call orbel_xv2aeq(xh(iouter),yh(iouter),zh(iouter),
     &             vxh(iouter),vyh(iouter),vzh(iouter),
     &             mass(iouter)+mass(1),ialpha,a,e,q)
      apoiouter = a*(1.d0+e) 

      nbodtemp = nbod
      ic = nbodm+1
      do i = nbodm+1,nbod

c.... we check if the little body needs to be duplicated
        needdup = .false.
        if (mass(ic).gt.msmall) then
          call orbel_xv2aeq(xh(ic),yh(ic),zh(ic),
     &       vxh(ic),vyh(ic),vzh(ic),mass(ic)+mass(1),ialpha,a,e,q)
          needdup =  (q.lt.(1.1d0*apoiouter))
c          print*,i,ic,sngl(apoiouter),sngl(q/(1.1d0*apoiouter))
        end if

c...  Duplicate it if necessary
        if (needdup) then
          xtemp = xh(ic)
          ytemp = yh(ic)
          ztemp = zh(ic)
          vxtemp = vxh(ic)
          vytemp = vyh(ic)
          vztemp = vzh(ic)
          masstemp = mass(ic)
          rpltemp = rpl(ic)
          rhilltemp = rhill(ic)
          vtemp = sqrt(vxtemp*vxtemp + vytemp*vytemp + vztemp*vztemp)

c.... supressing the particle duplicated

          call discard_mass_reorder5(ic,nbodtemp,mass,xh,yh,zh,
     &           vxh,vyh,vzh,rpl,rhill,isperih)

c.... Creation of the list of particles to add

          mtot = 0
          nbcopy = 0
          do while ((mtot.LT.masstemp).and.
     &           (nbodtemp+nbcopy.LT.NTPMAX))

            nbcopy = nbcopy+1
            call random_direction(theta,phi)
            call gasdev(s)
            dv = 0.01d0*vtemp*s
            cth = cos(theta)
            sth = sin(theta)
            cphi = cos(phi)
            sphi = sin(phi)
            vxadd(nbcopy) = vxtemp + dv*sth*cphi
            vyadd(nbcopy) = vytemp + dv*sth*sphi
            vzadd(nbcopy) = vztemp + dv*cth
            xadd(nbcopy)= xtemp
            yadd(nbcopy)= ytemp
            zadd(nbcopy)= ztemp
            massadd(nbcopy) = rand_mass(mbottom,Msmall)
            mtot = mtot + massadd(nbcopy) 
          end do

c             print*,i,ic,'duplique en',nbcopy,'total ',nbodtemp
c     &    sngl(masstemp/mtiny),sngl(masstemp/msmall),
c     &    sngl(msmall/mbottom)
c             print*,'--->',sngl(massadd(1:nbcopy)/mbottom)
          if (nbodtemp+nbcopy.GE.NTPMAX)then
            write(*,*) "TOO MUCH PARTICLES: impossible to duplicate"
            return
          end if
c.... rescaling masses

          if (mtot.NE.masstemp) then
            do j = 1,nbcopy
              massadd(j) = massadd(j)*(masstemp/mtot)
            end do
          end if

c.... Calculated radii and Hill spheres             
          do j = 1, nbcopy
            rpladd(j) = rpltemp*(massadd(j)/masstemp)**third
            call util_hills1(mass(1),massadd(j),
     &             xadd(j),yadd(j),zadd(j),
     &             vxadd(j),vyadd(j),vzadd(j),rhilladd(j))
          end do

c.... incorporating new particles 

          do j = 1,nbcopy
            k = nbodtemp

c.... checking where to add particle j
            do while (massadd(j).GT.mass(k))
              k = k-1
            end do
c.... we have to add in k+1
            nbodtemp = nbodtemp+1

c.... we move particles
            do l = nbodtemp,k+2,-1
              mass(l)= mass(l-1)
              xh(l)  = xh(l-1)
              yh(l)  = yh(l-1)
              zh(l)  = zh(l-1)
              vxh(l) = vxh(l-1)
              vyh(l) = vyh(l-1)
              vzh(l) = vzh(l-1)
              rpl(l) = rpl(l-1)
              rhill(l) = rhill(l-1)
            end do
         
c.... we add the particle j

            mass(k+1) = massadd(j)
            xh(k+1) = xadd(j)
            yh(k+1) = yadd(j)
            zh(k+1) = zadd(j)
            vxh(k+1) = vxadd(j)
            vyh(k+1) = vyadd(j)
            vzh(k+1) = vzadd(j)
            rpl(k+1) = rpladd(j)
            rhill(k+1) = rhilladd(j)
          end do      
        else
c... We don't need to duplicate the particle. Move to the next one
          ic = ic+1
        end if
      end do

      nbod = nbodtemp

      end SUBROUTINE

