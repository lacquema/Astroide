c**********************************************************************
c                      SYMBA6_SECTOR_CHKLIST.F
c**********************************************************************
c
c  Purpose:  This subroutine gets the encounter search list for each
c            sector in the grid.
c
c  Inputs:
c            nbod       ==>  number of massive bodies (int scalar)
c            dt         ==>  time step
c            nsect      ==>  Total number of sectors(int scalar)
c            nrz        ==>  Number of radial search zones (int scalar)  
c            naz        ==>  Number of azimuthal search zones  (int scalar)  
c            xh,yh      ==>  position in helio coord 
c                                    (real arrays)
c            vxh,vyh    ==>  velocity in helio coord 
c                                    (real arrays)
c            rhill      ==>  size of planet's hills sphere
c                                    (real array)
c            rbnd1      ==> outer boudary of innermost zone (real*8 scalar)
c            dr         ==> width of radial zones (real*8 scalar)
c            dtheta     ==> extent of azimuthal zones (real*8 scalar)
c             
c   Output:
c         sect_chklst   ==> what particles in a sectors check list 
c                           (2D int array)
c         ichck         ==> number of particles to be checked in sector 
c                           (int array)
c
c  Author:  cba
c  Date:  3/5/99
c  Last Modified:
c
c
      subroutine symba6_sector_chklist(ichck,sect_chklst,nbod,dt,
     &          nsect,naz,nrz,xh,yh,vxh,vyh,rhill,rbnd1,dr,dtheta)

       include '../swift.inc'
       include 'symba6.inc'

c..Inputs
       integer nbod,naz,nrz,nsect
       real*8 xh(nbod),yh(nbod),vxh(nbod),vyh(nbod)
       real*8 rhill(nbod),dt,rbnd1,dr,dtheta

c..Outputs     
       integer ichck(nsect), sect_chklst(NSECTMAX,NTPMAX)
     
c..Internals
       real*8 rho,theta,xhnext,yhnext
       real*8 rcrit,rhonext,rhomax,rhomin
       real*8 thetacrit,thetanext,thetamax,thetamin
       integer irmax,irmin,iazmax,iazmin
       integer ir,ia,ids,i,ja

c..Functions
       integer symba6_radzone

c----
c.. begin executable code

       do ids=1,nsect
          ichck(ids) = 0
       enddo


       do i = 2, nbod

          xhnext = xh(i) + vxh(i)*dt
          yhnext = yh(i) + vyh(i)*dt

          rho = sqrt( xh(i)**2 + yh(i)**2 )
          rhonext = sqrt(xhnext**2 + yhnext**2)  


c  getting the id numbers of the outer most and inner most radial
c  zones encountered

          rcrit = rhill(i)*RHSCALE

          if (rhonext .gt. rho) then    ! moving radially out
             rhomax = rhonext + rcrit
             rhomin = rho - rcrit
          else                          ! moving radially in
             rhomax = rho + rcrit
             rhomin = rhonext - rcrit
          endif

          irmax = symba6_radzone(rhomax,rbnd1,dr,nrz) 
          irmin = symba6_radzone(rhomin,rbnd1,dr,nrz)

c  getting the id numbers of the left most and right most azimuthal
c  zones encountered

          if (rcrit .ge. min(rho,rhonext)) then  ! can encounter origin       
             iazmin = 1                          ! and all azimuthal zones
             iazmax = naz
          else

             theta = atan2(yh(i),xh(i))
             thetanext = atan2(yhnext,xhnext)
             thetacrit = asin(rcrit/min(rho,rhonext))

             if (thetanext .gt. theta) then   ! moving counter clockwise
                thetamax = thetanext + thetacrit
                thetamin = theta - thetacrit
             else                             ! moving clockwise
                thetamax = theta + thetacrit
                thetamin = thetanext - thetacrit
             endif

             thetamax = mod(thetamax,TWOPI)
             if(thetamax.lt.0.0d0) then
                thetamax = thetamax + TWOPI
             endif
             thetamin = mod(thetamin,TWOPI)
             if(thetamin.lt.0.0d0) then
                thetamin = thetamin + TWOPI
             endif

             iazmax = thetamax/dtheta + 1
             iazmin = thetamin/dtheta + 1

c.. if the azimuthal zones encountered span the 0:2PI boundary

             if (iazmin .gt. iazmax) iazmax = iazmax + naz 

          endif 


c.. add body 'i' to the checklists for all encountered sectors

          do ir = irmin,irmax

             do ja = iazmin, iazmax

                ia = mod(ja,naz)
                if(ia .eq. 0) ia = naz

                ids = ia + (ir-1)*naz

                if (ids .le. 0) then
                   write(*,*) 'ERROR in symba6_sector_chklst'
                   call util_exit(1)
                endif

                ichck(ids) = ichck(ids) + 1
                sect_chklst(ids,ichck(ids)) = i

             enddo
          enddo


      enddo

      return
      end   ! symba6_sector_chklist
c------------------------------------------------------------

