c**********************************************************************
c                        SYMBA6_SECT_BNDRS.F
c**********************************************************************
c
c  Purpose:  This subroutine determines the boundaries of the search
c            sectors.  Azimuthal zones are divided evenly in angles.
c            Radial zones are divided into annuli with equal radial
c            spacing (dr)
c
c    Input:
c            nrz        ==>  Number of radial search zones (int scalar)  
c            naz        ==>  Number of azimuthal search zones  (int scalar)  
c            nbod       ==>  number of massive bodies (int scalar)
c            xh,yh      ==>  position in helio coord 
c                                    (real arrays)
c
c   Output:
c            rbnd1      ==> outer boudary of innermost zone (real*8 scalar)
c            dr         ==> width of radial zones (real*8 scalar)
c            dtheta     ==> extent of azimuthal zones (real*8 scalar)
c
c Remarks:
c Author:  cba
c Date:    5/4/99
c Last revision: 
c

      subroutine symba6_sect_bndrs(nrz,naz,nbod,xh,yh,rbnd1,dr,dtheta)

      include '../swift.inc'
      include 'symba6.inc'

c.. Inputs
      integer nrz,naz,nbod
      real*8 xh(nbod),yh(nbod)

c.. Outputs
      real*8 rbnd1,dr,dtheta

c.. Internals
      integer iaz, irz, i,j
      real*8 rho(NTPMAX), rhotmp

c----
c.. Executable code

c  calculate the radial cylindrical coordinate of each body

      do i = 2,nbod 
         rho(i)=sqrt( xh(i)**2 + yh(i)**2 )
      enddo

c  sort the radial coordinate in ascending order
    
      rho(1) = 0.d0

      do i = 2,nbod   
         do j= i,nbod
             if (rho(i) .gt. rho(j)) then
                 rhotmp = rho(i)
                 rho(i) = rho(j)
                 rho(j) = rhotmp
             endif
         enddo
      enddo


c  set the radial boundaries for each zone


      dr = (rho(nbod) - rho(2))/float(nrz)  ! width of radial zones

      rbnd1 = rho(2) + dr              ! outer boudary of innermost zone      

      write(*,*) ' radial sector boundaries '

      do irz = 1,nrz
        write(*,*) irz, rbnd1 + float(irz-1)*dr   
      enddo

c  set the azimuthal boundaries for each zone

      dtheta = TWOPI/float(naz)      

      write(*,*) 'azimuthal sector boundaries'

      do iaz = 1,naz                      ! azimuthal zone index

         write(*,*) iaz, float(iaz)*dtheta   

      enddo

      return
      end

