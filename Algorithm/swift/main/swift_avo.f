c**********************************************************************
c		      SWIFT_AV.F
c**********************************************************************
c  This code is a realization of Hierarchical Jacobi Symplectic
c  n-body mapping method for hierarchical stellar systems
c  (Beust 2003)
c
c                 TAKES TIDES AND RELATIVITY INTO ACCOUNT
c
c                 NO CLOSE ENCOUNTERS
c                 To run, need 3 input files. The code prompts for
c                 the file names, but examples are :
c
c                   parameter file like       param.in
c		    planet file like          pl.in
c                   test particle file like   tp.in
c
c Authors:  Herve Beust
c Date:    Jun 6, 2008
c Remark : based on swift_hjs.f

     
        program swift_av 

	include 'swift.inc'

	real*8 mass(3),inert(3)
        real*8 rotx(3),roty(3),rotz(3)
        real*8 vecte(3),vectep(3),vecth(3),vecthp(3)

	integer nbod,ntp,nleft
	integer iub,iuj,iud,iue,i,j

	real*8 t0,tstop,dt,dtout,dtdump
	real*8 t,tout,tdump,tfrac,eoff

	real*8 rmin,rmax,rmaxu,qmin,rplsq(3),krpl5(3)
        real*8 Qlov(3)
        integer tid(3)
        logical*2 lclose 

	character*(chlen) outfile,inparfile,inplfile,intpfile,
     &             fopenstat,diro,dirs,gname,dataname,genfile

        integer ialpha,fverb
        real*8 a,e,inc,capom,omega,capm
        logical ok

c-----
c...    Executable code 

        call io_init_av(mass,vecte,vecth,vectep,vecthp,
     &         rotx,roty,rotz,inert,krpl5,Qlov,rplsq,
     &         t0,tstop,dt,dtout,dtdump,fverb,diro,dirs,gname,
     &         outfile,fopenstat)

c Copy initial files input into work directory
        if ((fopenstat(1:3).eq.'new')
     &       .or.(fopenstat(1:3).eq.'NEW')) then
          dataname = 'cp swift_av.sh '//trim(diro)
          call system(dataname)
        end if


c Initialize initial time and times for first output and first dump
	t = t0
	tout = t0 + dtout
	tdump = t0 + dtdump

        iub = 20
        iuj = 30
        iud = 40
        iue = 60

c...    Do the initial io write

        call io_write_frame_av(t0,fverb,mass,vecte,vecth,
     &                  vectep,vecthp,rotx,roty,rotz,
     &                  trim(diro)//'/'//outfile,iub,fopenstat)

c        if(btest(iflgchk,2))  then    ! bit 2 is set
c           eoff = 0.0d0
c           call anal_energy_write_hjs(t,nbod,umat,mass,xj,yj,zj,
c     &                         vxj,vyj,vzj,iue,fopenstat,diro,eoff)
c        endif
c        if(btest(iflgchk,3))  then    ! bit 3 is set
c           call anal_jacobi_write(t0,nbod,ntp,mass,xh,yh,zh,vxh,
c     &        vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,2,iuj,fopenstat)
c        endif            No Jacobi intagral in HJS...

c...    must initize discard io routine


c***************here's the big loop *************************************
        write(*,*) ' ************** MAIN LOOP ****************** '

          ok = .true.
c	  do while ( (t .le. tstop) .and. 
c     &       ((ntp.eq.0).or.(nleft.gt.0)) )
          do while (ok)

            call kickav_tid(mass,inert,krpl5,Qlov,
     &     vecte,vecth,vectep,vecthp,rotx,roty,rotz,dt)

	    t = t + dt
            ok = (ok.and.(t.le.tstop))

c if it is time, output orb. elements, 
            if (t.ge.tout) then 

              call io_write_frame_av(t,fverb,mass,vecte,vecth,
     &                  vectep,vecthp,rotx,roty,rotz,
     &                  trim(diro)//'/'//outfile,iub,fopenstat)
	      tout = tout + dtout
	    endif

c If it is time, do a dump
          if(t.ge.tdump) then

             tfrac = (t-t0)/(tstop-t0)
c             write(*,998) t,tfrac
 998         format(' Time = ',1p1e12.5,': fraction done = ',0pf5.3)

	     call io_dump_pl_av(trim(diro)//'/'//'swift_cont.sh',
     &              t,tstop,dt,dtout,dtdump,
     &              mass,vecte,vecth,vectep,vecthp,rotx,roty,rotz,
     &              inert,krpl5,Qlov,rplsq,fverb,dirs,gname,outfile)

	     tdump = tdump + dtdump

c             if (btest(iflgchk,2))  then    ! bit 2 is set
c               call anal_energy_write_hjs(t,nbod,umat,mass,xj,yj,zj,
c     &                           vxj,vyj,vzj,iue,fopenstat,diro,eoff)
c             endif
c        if(btest(iflgchk,3))  then    ! bit 3 is set
c           call anal_jacobi_write(t0,nbod,ntp,mass,xh,yh,zh,vxh,
c     &        vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,2,iuj,fopenstat)
c        endif            No Jacobi intagral in HJS...

	  endif

	enddo
c********** end of the big loop from time 't0' to time 'tstop'

c Do a final dump for possible resumption later 

	call io_dump_pl_av(trim(diro)//'/'//'swift_cont.sh',
     &              t,tstop,dt,dtout,dtdump,
     &              mass,vecte,vecth,vectep,vecthp,rotx,roty,rotz,
     &              inert,krpl5,Qlov,rplsq,fverb,dirs,gname,outfile)

        call util_exit(0)
        end    ! swift_hjs_tid.f
c---------------------------------------------------------------------








c*************************************************************************
c                        KICKAV_TID.F
c*************************************************************************
c This subroutine tidally moves the orbit and to the rotation
c axes of the solid bodies (averaged)


      subroutine kickav_tid(mass,inert,krpl5,Qlov,
     &     vecte,vecth,vectep,vecthp,rotx,roty,rotz,dt)

      include 'swift.inc'

c...  Inputs: 
      real*8 inert(3),krpl5(3),Qlov(3)
      real*8 mass(3)
      real*8 dt

c...  Inputs & outputs
      real*8 vecte(3),vectep(3),vecth(3),vecthp(3)
      real*8 rotx(3),roty(3),rotz(3)

c...  Internals:
      integer i,j,k,kk,iflg,ialpha,ic
      real*8 s,adt
      real*8 eps,tfake,dttmp,ybs(18)

c----
c...  Executable code 

c...  get thr r^-3's

      eps = 1d-8

      ybs(1:3) = vecte(1:3)
      ybs(4:6) = vecth(1:3)
      ybs(7:9) = vectep(1:3)
      ybs(10:12) = vecthp(1:3)  
      ybs(13) = rotx(1)
      ybs(14) = roty(1)
      ybs(15) = rotz(1)
      ybs(16) = rotx(2)
      ybs(17) = roty(2)
      ybs(18) = rotz(2)

      tfake = 0.0d0
      adt = abs(dt)
      s = 1.d0
      if (dt.lt.0.d0) s = -s
      dttmp = adt
c          print*,'bs'

c      do while(tfake.lt.dt)
c          ic = 0
      do while( (abs(tfake-adt)/adt) .gt. 1.0d-7 )  ! just to be real safe
            
c            print*,'ybs avnat',sngl(ybs)
        call bs_av_tid(s,mass,inert,krpl5,Qlov,
     &                                    tfake,dttmp,ybs,eps)
c        call rkuttaqs(mass,inert,krpl5,Qlov,tfake,dttmp,ybs,eps)
c            print*,'ybs apres',sngl(ybs)
        dttmp = adt - tfake
c            ic = ic+1
      end do
c          print*,'ic',ic
c          print*,'bs fini'

      vecte(1:3) = ybs(1:3)
      vecth(1:3) = ybs(4:6)
      vectep(1:3) = ybs(7:9)
      vecthp(1:3) = ybs(10:12)
      rotx(1) = ybs(13)
      roty(1) = ybs(14)
      rotz(1) = ybs(15)
      rotx(2) = ybs(16)
      roty(2) = ybs(17)
      rotz(2) = ybs(18)

      return
      end      ! kickav_tid

c---------------------------------------------------------------------

      subroutine bs_av_tid(s,mass,inert,krpl5,Qlov,x,h0,y,eps)

      include 'swift.inc'

c...  Inputs Only: 
      integer nbod,ntp
      real*8 mass(3),inert(3),krpl5(3),Qlov(3),h0,eps,j2rp2,j4rp4,s

c...  Input & Output
      real*8 x,y(18)

c...  Internals
      real*8 tp(12*18),dy(18),d(6),alt(10),lt(10)
      integer idiv,i,ii,m,l,m1,k,mmk,i1,i1max,ik,n,ic
      real*8 xa,xb,varm,fl,h,hd,flt,eta2,dta,yb,c,b1
      real*8 den,dtn,b,var,varma
      real*8 gm,mu,ff(10)
      integer, parameter :: NTRYS = 30
      logical lbreak

      data lt/1,2,3,4,6,8,12,16,24,32/
      data alt/1.d0,2.d0,3.d0,4.d0,6.d0,8.d0,12.d0,16.d0,24.d0,32.d0/

      save lt,alt

c----
c...  Executable code 


c     ic = 0
      gm = mass(1)+mass(2)
      ff(1) = (1.d0+mass(2)/mass(1))*krpl5(1)
      ff(2) = (1.d0+mass(1)/mass(2))*krpl5(2)
      ff(3) = krpl5(1)*mass(2)/mass(1)/Qlov(1)
      ff(4) = krpl5(2)*mass(1)/mass(2)/Qlov(2)
      mu =  mass(2)*mass(1)/gm
      ff(5) = mu/gm  
      ff(6) = -mu/inert(1)
      ff(8) = -mu/inert(2)
      ff(7) = gm
      ff(9) = ff(1)*mass(2)
      ff(10) = ff(2)*mass(1)

      n = 18

      xa=x
      call derivs_tid(mass,ff,y,dy)
c      ic=ic+1
      if (s.lt.0.d0) dy = s*dy
      do i=1,n
         ii=12*i
         tp(ii-1)=dabs(y(i))
c
         if(tp(ii-1).lt.eps) then
            tp(ii-1)=eps
         endif
c
         tp(ii-4)=dy(i)
         tp(ii)=y(i)
      enddo

      do idiv=0,NTRYS

         xb=h0+xa
c
c        successive extrapolations
c
c         do m=1,10
         m = 1
         lbreak = .true.
         do while( (m.le.10) .and. lbreak )

            l=lt(m)
            fl=alt(m)
            varm=0.d0
            m1=min0(m-1,6)
c
c           calculation of d(k)=(h(m-k)/h(m))**2 for equation (6)
c
            if(m1.ne.0) then
               do k=1,m1
                  mmk=m-k
                  flt=alt(mmk)
                  d(k)=(fl/flt)**2
               enddo
            endif
            h=h0/fl
            hd=0.5d0*h

c           integration
c
            do i=1,n
               ii=12*i
               tp(ii-3)=tp(ii) 
               y(i)=tp(ii)+hd*tp(ii-4) !    equation (3b)
            enddo
            i1max=2*l-1
            x=xa
         
            do i1=1,i1max
               x=x+hd
               call derivs_tid(mass,ff,y,dy)
c               ic=ic+1
               if (s.lt.0.d0) dy = s*dy
               do i=1,n
                  ii=12*i
                  tp(ii-1)=dmax1(tp(ii-1),dabs(y(i)))
                  eta2=tp(ii-3)+h*dy(i)
                  tp(ii-3)=y(i)
                  y(i)=eta2
               enddo 
            enddo
         
            call derivs_tid(mass,ff,y,dy)
c            ic=ic+1
            if (s.lt.0.d0) dy = s*dy
            do i=1,n
               ii=12*i
               dta=tp(ii-11)
               yb=(tp(ii-3)+y(i)+hd*dy(i))/2.d0 !    equation (3d)
c     
c              extrapolated values
c
               c=yb             ! equation (6b)
               tp(ii-11)=yb     ! equation (6a)

               if(m1.ne.0) then
                  do k=1,m1
                     b1=d(k)*dta
                     den=b1-c
                     dtn=dta
                     if(den.ne.0.d0) then
                        b=(c-dta)/den
                        dtn=c*b       !   equation (6c)
                        c=b1*b        !    equation (6d)
                     endif
                     ik=ii-11+k
                     dta=tp(ik)
                     tp(ik)=dtn       !   equation (6c)
                     yb=yb+dtn        !     equation (6f)
                  enddo
                  var=dabs(tp(ii-2)-yb)/tp(ii-1)
                  varm=dmax1(varm,var)
               endif
               tp(ii-2)=yb
            enddo
            
            if(m.gt.3) then
               if(varm.le.eps) then       !   return results to calling program
                  x=xb
                  do i=1,n
                     ii=12*i
                     y(i)=tp(ii-2)
                  enddo
                  h0=h0*1.5d0*0.6d0**(m-1-m1)    !   recommend a new step size
c                  if (ic.ne.21) print*,'ic',ic
                  return
               endif

               if(varm.ge.varma) then !  calculation did not converge
c                                        start again with half the step size
                  lbreak = .false.
               endif
            endif
            varma=varm
            m = m + 1
         enddo                  ! m

         h0=h0/2.d0
      enddo      ! idiv

      write(*,*) ' ERROR (b_int): lack of convergence !!!! '

c
      end      !  bs_av_tid
c-----------------------------------------------------------------------------

C**************************************************************************
C	    		        BS_DER_TID
C**************************************************************************
c This is the subroutine that calculates the tidal derivatives 
c
c             Input:
c              ff(10)  ==> various constants
c              ybs(18)   ==> Variables 
c
c             Output:
c 	         dy  ==> derivatives of the independant var (real array)
c
c Authors:  Herve Beust
c Date:    2/17/2009 (Based on bs_der)
c

      subroutine derivs_tid(mass,ff,ybs,dy)
       
      include 'swift.inc'

c...  Inputs Only: 
      real*8 ff(10)
      real*8 ybs(18)
      real*8 mass(3)

c...  Output
      real*8 dy(18)

c...  Internals:
      integer i,j,k,kk,ialpha,iflg
      real*8 a,ap,e,q,gm,gmp,mu,f2m2,f2m1
      real*8 f1,f2,f3,ux,uy,uz,nq,nn,nnp
      real*8 f2f3,f4,f5,f6,f7,f8,f9,f10,f11
      real*8 exq,exq2,exq3,ex4,ex6,ex2,exq4
      real*8 excp,excp2,exqp,exqp2,exqp3,exqp4,exqp5,exqp6,exqp7,exqp8
      real*8 alp,alpha,alpha2,alpha3,alpha4,alpha5
      real*8 ep,ep2,ep3,ep5,ep6,ep12,ep52,ep72,ep92,ep132
      real*8 q2,q3,q5,q8,z332,z232
      real*8 gqdx(2),gqdy(2),gqdz(2),gtfx(2),gtfy(2),gtfz(2)
      real*8 dqdx(2),dqdy(2),dqdz(2),dtfx(2),dtfy(2),dtfz(2)
      real*8 grelx,grely,grelz
      real*8 e1(3),e2(3),e3(3),rx(2),ry(2),rz(2)
      real*8 e1p(3),e2p(3),e3p(3),z(3,3)
      real*8 vecth(3),vecte(3),h
      real*8 vecthp(3),vectep(3),hp
      real*8 hh1,hh2,hh3,hh4,hh5,hh6,hh7
      real*8 vdotr,posx,posy,vx,vy,vz,cu,su,capm
      real*8 gtotx,gtoty,gtotz,dtotx,dtoty,dtotz,mrel,dcapm
      real*8 rotx(2),roty(2),rotz(2),mqd(2) 
      real*8 dvecte(3),dvecth(3),drotx(2),droty(2),drotz(2)
      real*8 dvectpe(3),dvectph(3)
      real*8 dhdtx(5),dhdty(5),dhdtz(5),dhpdtx(5),dhpdty(5),dhpdtz(5)
      real*8 drundtx(5),drundty(5),drundtz(5)
      real*8 drunpdtx(5),drunpdty(5),drunpdtz(5)
      real*8 zz0,zz1,zz2,zz3,zz4,zz5,zz6,zz7,zz8,zz9,zz10,zz11,dd
      real*8, parameter :: 
     &           cl = 2.9979248d8*86400.d0*365.2422d0/1.4959787061d11
      real*8, parameter :: cl2 = cl*cl
      integer, parameter :: imax = 4

c----
c...  Executable code 

c...  move things so that I can deal with it
      vecte(1:3) = ybs(1:3)
      vecth(1:3) = ybs(4:6)
      vectep(1:3) = ybs(7:9)
      vecthp(1:3) = ybs(10:12)
      rotx(1) = ybs(13)
      roty(1) = ybs(14)
      rotz(1) = ybs(15)
      rotx(2) = ybs(16)
      roty(2) = ybs(17)
      rotz(2) = ybs(18)

      e = sqrt(vecte(1)*vecte(1) + vecte(2)*vecte(2)
     &                                  + vecte(3)*vecte(3))
      h = sqrt(vecth(1)*vecth(1) + vecth(2)*vecth(2)
     &                                  + vecth(3)*vecth(3))
      excp = sqrt(vectep(1)*vectep(1) + vectep(2)*vectep(2)
     &                                  + vectep(3)*vectep(3))
      hp = sqrt(vecthp(1)*vecthp(1) + vecthp(2)*vecthp(2)
     &                                  + vecthp(3)*vecthp(3))

c..............................

      gm = ff(7)
c      gm = mass(1)+mass(2)
      gmp = gm + mass(3)
c      mu = mass(1)*mass(2)/(mass(1)+mass(2))

      e3(1) = vecth(1)/h
      e3(2) = vecth(2)/h
      e3(3) = vecth(3)/h

      if (e*e.gt.TINY) then
        e1(1) = vecte(1)/e
        e1(2) = vecte(2)/e
        e1(3) = vecte(3)/e
      else
        i = 1
        if (abs(vecth(2)).gt.abs(vecth(i))) i = 2
        if (abs(vecth(3)).gt.abs(vecth(i))) i = 3
        j = mod(i,3)+1
        k = mod(j,3)+1
        f2 = 1.d0/sqrt(1.d0-e3(k)*e3(k))
        e1(i) = -e3(j)*f2
        e1(j) = e3(i)*f2
        e1(k) = 0.d0
      end if 

      e2(1) = e3(2)*e1(3) - e3(3)*e1(2)
      e2(2) = e3(3)*e1(1) - e3(1)*e1(3)
      e2(3) = e3(1)*e1(2) - e3(2)*e1(1)

c..............................

      e3p(1) = vecthp(1)/hp
      e3p(2) = vecthp(2)/hp
      e3p(3) = vecthp(3)/hp

      if (excp*excp.gt.TINY) then
        e1p(1) = vectep(1)/excp
        e1p(2) = vectep(2)/excp
        e1p(3) = vectep(3)/excp
      else
        i = 1
        if (abs(vecthp(2)).gt.abs(vecthp(i))) i = 2
        if (abs(vecthp(3)).gt.abs(vecthp(i))) i = 3
        j = mod(i,3)+1
        k = mod(j,3)+1
        f2 = 1.d0/sqrt(1.d0-e3p(k)*e3p(k))
        e1p(i) = -e3p(j)*f2
        e1p(j) = e3p(i)*f2
        e1p(k) = 0.d0
      end if 

      e2p(1) = e3p(2)*e1p(3) - e3p(3)*e1p(2)
      e2p(2) = e3p(3)*e1p(1) - e3p(1)*e1p(3)
      e2p(3) = e3p(1)*e1p(2) - e3p(2)*e1p(1)

c.......................................
c      if (e.gt.0.6) then
c      print*,e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3)
c      print*,e1(1)*e2(1)+e1(2)*e2(2)+e1(3)*e2(3)
c      print*,e1(1)*e3(1)+e1(2)*e3(2)+e1(3)*e3(3)
c      print*,e2(1)*e1(1)+e2(2)*e1(2)+e2(3)*e1(3)
c      print*,e2(1)*e2(1)+e2(2)*e2(2)+e2(3)*e2(3)
c      print*,e2(1)*e3(1)+e2(2)*e3(2)+e2(3)*e3(3)
c      print*,e3(1)*e1(1)+e3(2)*e1(2)+e3(3)*e1(3)
c      print*,e3(1)*e2(1)+e3(2)*e2(2)+e3(3)*e2(3)
c      print*,e3(1)*e3(1)+e3(2)*e3(2)+e3(3)*e3(3)
c       print*,'---------'
c      print*,e1p(1)*e1p(1)+e1p(2)*e1p(2)+e1p(3)*e1p(3)
c      print*,e1p(1)*e2p(1)+e1p(2)*e2p(2)+e1p(3)*e2p(3)
c      print*,e1p(1)*e3p(1)+e1p(2)*e3p(2)+e1p(3)*e3p(3)
c      print*,e2p(1)*e1p(1)+e2p(2)*e1p(2)+e2p(3)*e1p(3)
c      print*,e2p(1)*e2p(1)+e2p(2)*e2p(2)+e2p(3)*e2p(3)
c      print*,e2p(1)*e3p(1)+e2p(2)*e3p(2)+e2p(3)*e3p(3)
c      print*,e3p(1)*e1p(1)+e3p(2)*e1p(2)+e3p(3)*e1p(3)
c      print*,e3p(1)*e2p(1)+e3p(2)*e2p(2)+e3p(3)*e2p(3)
c      print*,e3p(1)*e3p(1)+e3p(2)*e3p(2)+e3p(3)*e3p(3)
c      stop
c      end if

      z(1,1) = e1(1)*e1p(1)+e1(2)*e1p(2)+e1(3)*e1p(3)
      z(1,2) = e1(1)*e2p(1)+e1(2)*e2p(2)+e1(3)*e2p(3)
      z(1,3) = e1(1)*e3p(1)+e1(2)*e3p(2)+e1(3)*e3p(3)      
      z(2,1) = e2(1)*e1p(1)+e2(2)*e1p(2)+e2(3)*e1p(3)
      z(2,2) = e2(1)*e2p(1)+e2(2)*e2p(2)+e2(3)*e2p(3)
      z(2,3) = e2(1)*e3p(1)+e2(2)*e3p(2)+e2(3)*e3p(3)      
      z(3,1) = e3(1)*e1p(1)+e3(2)*e1p(2)+e3(3)*e1p(3)
      z(3,2) = e3(1)*e2p(1)+e3(2)*e2p(2)+e3(3)*e2p(3)
      z(3,3) = e3(1)*e3p(1)+e3(2)*e3p(2)+e3(3)*e3p(3)      
     
c      do i=1,3
c        do j=1,3
c          z2(i,j) = z(i,j)*z(i,j)
c          z3(i,j) = z2(i,j)*z(i,j)
c        end do
c      end do
c      z2 = z*z
c      z3 = z2*z
c      z4 = z2*z2

c.......................................

      q = h*h/gm/(1.d0+e)
      q2 = q*q
      q3 = q2*q
      q5 = q3*q2
      q8 = q5*q3
c      print*,'yb',ybs

      a = q/(1.d0-e)
      ap = hp*hp/gmp/(1.d0-excp*excp)
      alpha = a/ap
      alpha2 = alpha*alpha
      alpha3 = alpha2*alpha
      alpha4 = alpha3*alpha
      alpha5 = alpha4*alpha
      nn = sqrt(gm/(a*a*a))
      nnp = sqrt(gmp/(ap*ap*ap))
      nq = sqrt(gm/q3)

      ex2 = e*e      
      excp2 = excp*excp
      exq2 = 1.d0-ex2
      exqp2 = 1.d0-excp2
      exq = sqrt(exq2)
      exqp = sqrt(exqp2)
      exq3 = exq2*exq
      exqp3 = exqp2*exqp
      exqp4 = exqp2*exqp2
      exqp5 = exqp4*exqp
      exqp6 = exqp4*exqp2 
      exqp7 = exqp6*exqp
      exqp8 = exqp6*exqp2
      exq4 = exq2*exq2
      ex4 = ex2*ex2
      ex6 = ex4*ex2
      ep = 1.d0+e
      ep12 = sqrt(ep)
      ep2 = ep*ep
      ep3 = ep2*ep
      ep5 = ep2*ep3
      ep6 = ep3*ep3 
      ep52 = ep2*ep12
      ep72 = ep3*ep12
      ep92 = ep2*ep52
      ep132 = ep6*ep12

      hh1 = 30.d0*(1.d0+1.5d0*ex2+0.125d0*ex4)
      hh2 = 54.d0*(1.d0+3.75d0*ex2+1.875d0*ex4+0.078125d0*ex6)
      hh3 = 3.d0*(1.d0+1.5d0*ex2+0.125d0*ex4)
      hh4 = 6.d0*(1.d0+7.5d0*ex2+5.625d0*ex4+0.3125d0*ex6)
      hh5 = 3.d0*(1.d0+4.5d0*ex2+0.625d0*ex4)
      hh6 = 6.d0*(1.d0+3.d0*ex2+0.375d0*ex4)
      hh7 = 6.d0*(1.d0-4.5d0*ex2-0.875d0*ex4)

c........
      dd = mass(3)*alpha2/ap/exqp3
      z332 = z(3,3)*z(3,3)
      z232 = z(2,3)*z(2,3)
      dhdtx(2) = -0.75d0*dd*z(3,3)*z(2,3)*exq2
      dhdty(2) = +0.75d0*dd*z(1,3)*z(3,3)*(4.d0*ex2+1.d0)
      dhdtz(2) = -3.75d0*dd*z(1,3)*z(2,3)*ex2

      dd = mass(3)*alpha3*nn*e*exq/exqp3/gm
      drundtx(2) = 3.75d0*dd*z(2,3)*z(1,3)
      drundty(2) = 0.75d0*dd*(-3.d0+5.d0*z232+4.d0*z332)
      drundtz(2) = +0.75d0*dd*z(3,3)*z(2,3)

      dd = gmp*ff(5)*alpha2/ap/exqp3
      dhpdtx(2) = -0.75d0*dd*(z(3,3)*z(3,2)+
     &          (5.d0*z(2,3)*z(2,2)+4.d0*z(3,3)*z(3,2))*ex2)
      dhpdty(2) = 0.75d0*dd*(z(3,3)*z(3,1)+
     &          (5.d0*z(2,3)*z(2,1)+4.d0*z(3,3)*z(3,1))*ex2)
      dhpdtz(2) = 0.d0

      dd = ff(5)*alpha2*nnp*excp/exqp4
      drunpdtx(2) = 0.d0
      drunpdty(2) = 0.375d0*dd*(-1.d0+3.d0*z332+
     &          (-9.d0+15.d0*z232+12.d0*z332)*ex2)
      drunpdtz(2) = 0.75d0*dd*(z(3,3)*z(3,2)+
     &          (5.d0*z(2,3)*z(2,2)+4.d0*z(3,3)*z(3,2))*ex2)

c...............

      if (imax.gt.2)
     &  call ordre3(mass,ap,gm,gmp,nn,nnp,alpha3,alpha4,e,ex2,
     &  exq,exq2,excp,excp2,exqp4,exqp5,exqp6,z,dhdtx(3),
     &  dhdty(3),dhdtz(3),drundtx(3),drundty(3),drundtz(3),dhpdtx(3),
     &  dhpdty(3),dhpdtz(3),drunpdtx(3),drunpdty(3),drunpdtz(3))
       
      if (imax.gt.3)
     &  call ordre4(mass,ap,gm,gmp,nn,nnp,alpha4,alpha5,e,ex2,exq,exq2,
     &  exq4,excp,excp2,exqp2,exqp6,exqp7,exqp8,z,dhdtx(4),
     &  dhdty(4),dhdtz(4),drundtx(4),drundty(4),drundtz(4),dhpdtx(4),
     &  dhpdty(4),dhpdtz(4),drunpdtx(4),drunpdty(4),drunpdtz(4))
       
      if (imax.gt.4)
     &  call ordre5(mass,ap,gm,gmp,nn,nnp,alpha4,alpha5,e,ex2,
     &  exq,exq2,exq4,excp,excp2,exqp2,exqp6,exqp7,exqp8,z,dhdtx(5),
     &  dhdty(5),dhdtz(5),drundtx(5),drundty(5),drundtz(5),dhpdtx(5),
     &  dhpdty(5),dhpdtz(5),drunpdtx(5),drunpdty(5),drunpdtz(5))

c      dd = 2.d0*e*drundtx(imax)/exq2+2.d0*dhdtz(imax)/h
c       dd = dhdtx(3)*e+drundtz(3)*h

c      if (e.gt.0.6) print*,'daa4 = ',sngl(dd),
c     &             sngl(dhdtx(imax)*e),sngl(drundtz(imax)*h)     
c     &   sngl(2.d0*e*drundtx(imax)/exq2),sngl(2.d0*dhdtz(imax)/h)

c      dd = 2.d0*excp*drunpdtx(imax)/exqp2+2.d0*dhpdtz(imax)/hp

c      if (e.gt.0.6) print*,'daap4 = ',sngl(dd),
c     &             sngl(dhdtx(3)*e),sngl(drundtz(3)*h)     
c     & sngl(2.d0*excp*drunpdtx(imax)/exqp2),sngl(2.d0*dhpdtz(imax)/hp)

c........

      rx(1) = rotx(1)*e1(1) + roty(1)*e1(2) + rotz(1)*e1(3)
      ry(1) = rotx(1)*e2(1) + roty(1)*e2(2) + rotz(1)*e2(3)
      rz(1) = rotx(1)*e3(1) + roty(1)*e3(2) + rotz(1)*e3(3)
      rx(2) = rotx(2)*e1(1) + roty(2)*e1(2) + rotz(2)*e1(3)
      ry(2) = rotx(2)*e2(1) + roty(2)*e2(2) + rotz(2)*e2(3)
      rz(2) = rotx(2)*e3(1) + roty(2)*e3(2) + rotz(2)*e3(3)

      f4 = 1.d0/(nq*ep72*q5)
      f5 = 1.d0/(nq*ep132*q8)
      f7 = nq/(ep92*q3)
      f8 = e*hh3/(ep5*q5)
      f9 = nq*e*hh2/(ep132*q5)
      f10 = nq*nq*hh4/(ep6*q3)
      f11 = e*exq3*f4

      f2 = ff(1)
      f2m2 = ff(9)
      f3 = (rx(1)*rx(1)+ry(1)*ry(1)-2.d0*rz(1)*rz(1))
      f2f3 = f2*f3
      f6 = f2*exq3/(q3*ep3)

      gqdx(1) = 0.d0
      gqdy(1) = -0.5d0*f2f3*f11
      gqdz(1) = f2*ry(1)*rz(1)*f11
      gqdy(1) = gqdy(1) + f2m2*hh1*exq3*e*f5

      dqdx(1) = -f6*rz(1)*ry(1)
      dqdy(1) = f6*rz(1)*rx(1)
      dqdz(1) = 0.d0

c      mqd(1) = -(0.5d0*f2f3*f4 + f2m2*hh7*f5)*exq4

      f2 = ff(2)
      f2m1 = ff(10)
      f3 = (rx(2)*rx(2)+ry(2)*ry(2)-2.d0*rz(2)*rz(2))
      f2f3 = f2*f3
      f6 = f2*exq3/(q3*ep3)

      gqdx(2) = 0.d0
      gqdy(2) = -0.5d0*f2f3*f11
      gqdz(2) = f2*ry(2)*rz(2)*f11
      gqdy(2) = gqdy(2) + f2m1*hh1*exq3*e*f5
          
      dqdx(2) = -f6*rz(2)*ry(2)
      dqdy(2) = f6*rz(2)*rx(2)
      dqdz(2) = 0.d0
       
c      mqd(2) = -(0.5d0*f2f3*f4 + f2m1*hh7*f5)*exq4

      f2 = ff(3)

      gtfx(1) = -f2*f9
      gtfx(1) = gtfx(1) + 11.d0*rz(1)*f2*f8
      gtfz(1) = -rx(1)*f2*f8
      gtfy(1) = 0.d0

      dtfz(1) = -f10*f2
      dtfx(1) = f2*rx(1)*hh3*f7
      dtfy(1) = f2*ry(1)*hh5*f7
      dtfz(1) = dtfz(1) + f2*rz(1)*hh6*f7

      f2 = ff(4)

      gtfx(2) = -f2*f9    
      gtfx(2) = gtfx(2) + 11.d0*rz(2)*f2*f8
      gtfz(2) = -rx(2)*f2*f8
      gtfy(2) = 0.d0

      dtfz(2) = -f10*f2  
      dtfx(2) = f2*rx(2)*hh3*f7
      dtfy(2) = f2*ry(2)*hh5*f7
      dtfz(2) = dtfz(2) + f2*rz(2)*hh6*f7

      f4 = nq*nq*nq*q2/(cl2*ep52)

      grelx = 0.d0
      grely = 3.d0*f4*exq3*e
      grelz = 0.d0
        
      f2 = ff(5)

c      mrel = -(15.d0-6.d0*exq-9.d0*f2+7.d0*exq*f2)*f4*exq4

c          print*,sngl(gqdx),sngl(dqdx)
c          print*,sngl(gqdy),sngl(dqdy)
c          print*,sngl(gqdz),sngl(dqdz)
c          print*,sngl(gtfx),sngl(dtfx)
c          print*,sngl(gtfy),sngl(dtfy)
c          print*,sngl(gtfz),sngl(dtfy)
c          print*,sngl(grelx),sngl(grely),sngl(grelz)
c
      gtotx = gqdx(1) + gqdx(2) + gtfx(1) + gtfx(2) + grelx
      gtoty = gqdy(1) + gqdy(2) + gtfy(1) + gtfy(2) + grely
      gtotz = gqdz(1) + gqdz(2) + gtfz(1) + gtfz(2) + grelz

      dtotx = dqdx(1) + dqdx(2) + dtfx(1) + dtfx(2)
      dtoty = dqdy(1) + dqdy(2) + dtfy(1) + dtfy(2)
      dtotz = dqdz(1) + dqdz(2) + dtfz(1) + dtfz(2)

c      dcapm = mqd(1) + mqd(2) + mrel
c          print*,gtotx,e1(1),gtoty,e2(1),gtotz,e3(1)
c          print*,gtotx,e1(2),gtoty,e2(2),gtotz,e3(2)
c          print*,gtotx,e1(3),gtoty,e2(3),gtotz,e3(3)
          
      dvecte(1) = gtotx*e1(1) + gtoty*e2(1) + gtotz*e3(1)
      dvecte(2) = gtotx*e1(2) + gtoty*e2(2) + gtotz*e3(2)
      dvecte(3) = gtotx*e1(3) + gtoty*e2(3) + gtotz*e3(3)

      dvecth(1) = dtotx*e1(1) + dtoty*e2(1) + dtotz*e3(1)
      dvecth(2) = dtotx*e1(2) + dtoty*e2(2) + dtotz*e3(2)
      dvecth(3) = dtotx*e1(3) + dtoty*e2(3) + dtotz*e3(3)

c      dvecte(1:3) = 0.d0    !!!!!!!
c      dvecth(1:3) = 0.d0    !!!!!!!
      dvectpe(1:3) = 0.d0
      dvectph(1:3) = 0.d0

      do i = 2,imax
        dvecth(1) = dvecth(1) + dhdtx(i)*e1(1) + dhdty(i)*e2(1)
     &                                         + dhdtz(i)*e3(1) 
        dvecth(2) = dvecth(2) + dhdtx(i)*e1(2) + dhdty(i)*e2(2)
     &                                         + dhdtz(i)*e3(2) 
        dvecth(3) = dvecth(3) + dhdtx(i)*e1(3) + dhdty(i)*e2(3)
     &                                         + dhdtz(i)*e3(3) 
        dvecte(1) = dvecte(1) + drundtx(i)*e1(1) + drundty(i)*e2(1)
     &                                         + drundtz(i)*e3(1) 
        dvecte(2) = dvecte(2) + drundtx(i)*e1(2) + drundty(i)*e2(2)
     &                                         + drundtz(i)*e3(2) 
        dvecte(3) = dvecte(3) + drundtx(i)*e1(3) + drundty(i)*e2(3)
     &                                         + drundtz(i)*e3(3) 
        dvectph(1) = dvectph(1) + dhpdtx(i)*e1p(1) + dhpdty(i)*e2p(1)
     &                                         + dhpdtz(i)*e3p(1) 
        dvectph(2) = dvectph(2) + dhpdtx(i)*e1p(2) + dhpdty(i)*e2p(2)
     &                                         + dhpdtz(i)*e3p(2) 
        dvectph(3) = dvectph(3) + dhpdtx(i)*e1p(3) + dhpdty(i)*e2p(3)
     &                                         + dhpdtz(i)*e3p(3) 
        dvectpe(1) = dvectpe(1) + drunpdtx(i)*e1p(1) 
     &                      + drunpdty(i)*e2p(1) + drunpdtz(i)*e3p(1) 
        dvectpe(2) = dvectpe(2) + drunpdtx(i)*e1p(2)
     &                      + drunpdty(i)*e2p(2) + drunpdtz(i)*e3p(2) 
        dvectpe(3) = dvectpe(3) + drunpdtx(i)*e1p(3)
     &                      + drunpdty(i)*e2p(3) + drunpdtz(i)*e3p(3) 
      end do

c       if (e.gt.0.6) then
c         print*,'-------------'
c         print*, 'h',sngl(mass(1)*mass(2)/(mass(1)+mass(2))*dvecth)
c         print*, 'hp',sngl(
c     &  (mass(1)+mass(2))*mass(3)/(mass(1)+mass(2)+mass(3))*dvectph)
c         print*, 'sum',sngl(mass(1)*mass(2)/(mass(1)+mass(2))*dvecth
c     &    +(mass(1)+mass(2))*mass(3)/(mass(1)+mass(2)+mass(3))*dvectph)
c       end if
       
c       dd = 2.d0/h/h*(vecth(1)*dvecth(1)+vecth(2)*dvecth(2)
c     &                                +vecth(3)*dvecth(3))
c     &     +2.d0/exq2*(vecte(1)*dvecte(1)+vecte(2)*dvecte(2)
c     &                                +vecte(3)*dvecte(3))
c         if (e.gt.0.6) print*,'-------------',sngl(dd),
c     & sngl(2.d0/h/h*(vecth(1)*dvecth(1)+vecth(2)*dvecth(2)
c     &                                +vecth(3)*dvecth(3))),
c     & sngl(+2.d0/exq2*(vecte(1)*dvecte(1)+vecte(2)*dvecte(2)
c     &                                +vecte(3)*dvecte(3)))
c       if (e.gt.0.6) then
c         print*,'-----------------'
c         print*,'dh',sngl(vecth(1)*dvecth(1)+vecth(2)*dvecth(2)
c     &    +vecth(3)*dvecth(3)),sngl(h*(dhdtz(2)+0*dhdtz(3)))
c         print*,'de',sngl(vecte(1)*dvecte(1)+vecte(2)*dvecte(2)
c     &    +vecte(3)*dvecte(3)),sngl(e*(drundtx(2)+0*drundtx(3)))
c       end if

c      dvectpe(1:3)=0.d0
c      dvectph(1:3)=0.d0
c      dvecte(1:3)=0.d0
c      dvecth(1:3)=0.d0
      dtotx = dqdx(1) + dtfx(1)
      dtoty = dqdy(1) + dtfy(1)
      dtotz = dqdz(1) + dtfz(1)

      alp = ff(6)
      drotx(1) = (dtotx*e1(1) + dtoty*e2(1) + dtotz*e3(1))*alp
      droty(1) = (dtotx*e1(2) + dtoty*e2(2) + dtotz*e3(2))*alp
      drotz(1) = (dtotx*e1(3) + dtoty*e2(3) + dtotz*e3(3))*alp

      dtotx = dqdx(2) + dtfx(2)
      dtoty = dqdy(2) + dtfy(2)
      dtotz = dqdz(2) + dtfz(2)

      alp = ff(8)
      drotx(2) = (dtotx*e1(1) + dtoty*e2(1) + dtotz*e3(1))*alp
      droty(2) = (dtotx*e1(2) + dtoty*e2(2) + dtotz*e3(2))*alp
      drotz(2) = (dtotx*e1(3) + dtoty*e2(3) + dtotz*e3(3))*alp

      dy(1:3) = dvecte(1:3)
      dy(4:6) = dvecth(1:3)
      dy(7:9) = dvectpe(1:3)
      dy(10:12) = dvectph(1:3)
      dy(13) = drotx(1)
      dy(14) = droty(1)
      dy(15) = drotz(1)
      dy(16) = drotx(2)
      dy(17) = droty(2)
      dy(18) = drotz(2)
c      dy = 0.d0
c      print*,'dy',dy

      return
      end     ! derivs_tid

c--------------------------------------------------------------------


C
C-----------------------------------------------------------------------------
C       Determination des elements d'orbite     
C-----------------------------------------------------------------------------
C                                                                             
        SUBROUTINE ORBITE(GM,VECTE,VECTH,A,E,I,O,OM)                    

        IMPLICIT NONE

        REAL*8          GM,             ! Constante de couplage
     &                  E3(3),           !           
     &                  H,H2,vecth(3),  !         
     &                  E,Vecte(3),     ! Vecteur excentricite (et exc)
     &                  I,SI,           ! Inclinaison et sinus         
     &                  O,              ! Longidude du noeud ascendant 
     &                  OM,             ! Argument du pericentre       
     &                  A               ! Demi-grand axe                       

        e = sqrt(vecte(1)*vecte(1)+vecte(2)*vecte(2)
     &                                     +vecte(3)*vecte(3))
        h2 = vecth(1)*vecth(1)+vecth(2)*vecth(2)+vecth(3)*vecth(3)
        a = h2/gm/(1.d0-e*e)
        h = sqrt(h2)
        e3(1) = vecth(1)/h
        e3(2) = vecth(2)/h
        e3(3) = vecth(3)/h
        si=sqrt(e3(1)*e3(1)+e3(2)*e3(2))
        i=atan2(si,e3(3))        
        if (si.eq.0.0d0) then           
          o=0.0d0                       
          if (e.eq.0.0d0) then          
            om=0.0d0       
          else                          
            om=atan2(vecte(2),vecte(1))
          end if
        else
          o=atan2(e3(1),-e3(2))
          if (e.eq.0.0d0) then
            om=0.0d0
          else
            om=atan2(vecte(3),e3(1)*vecte(2)-e3(2)*vecte(1))
          end if
        end if
        end

c.................................................................

       subroutine io_write_frame_av(time,fverb,mass,
     &   vecte,vecth,vectep,vecthp,rotx,roty,rotz,oname,iu,fopenstat)

       include 'swift.inc'

c...  Inputs: 
       integer iu,fverb
       real*8 mass(3),time
       real*8 vecte(3),vecth(3),vectep(3),vecthp(3)
       real*8 rotx(3),roty(3),rotz(3)                  
       character*(*) oname,fopenstat                            

c...  Internals
       integer id,ierr,iverb
       logical ok,verbose

       real*8 a,e,inc,capom,omega,capm,n,rrot
       real*8 ap,ep,incp,capomp,omegap,capmp
       real*8 gm,gmp,ci,co,si,so,rot(2),irot(2),orot(2)
       real*8 istat(1),irel
       integer i1st    ! =0 first time through; =1 after
       data i1st/0/
       save i1st,iverb

       if(i1st.eq.0) then             
         call io_open(iu,oname,fopenstat,'UNFORMATTED',ierr)
         if(ierr.ne.0) then                                 
           write(*,*) ' SWIFT ERROR: in io_write_frame: '   
           write(*,*) '     Could not open binary output file:'
           call util_exit(1)                                   
         end if                                                 
         i1st = 1
         iverb = 0
       else
         call io_open(iu,oname,'append','UNFORMATTED',ierr)
       endif

       verbose = (mod(iverb,fverb).eq.0)
       call io_write_hdr_r(iu,time,3,0,istat)

       if (verbose) print*,'t=',time

       gm = mass(1)+mass(2)
       gmp = gm + mass(3)
       call orbite(gm,vecte,vecth,a,e,inc,capom,omega)
       call orbite(gmp,vectep,vecthp,ap,ep,incp,capomp,omegap)
       capm = 0.d0
       capmp = 0.d0

       irel = acos(cos(inc)*cos(incp)
     &           +cos(capom-capomp)*sin(inc)*sin(incp))*180./Pi
       if (verbose) then
         print*,'a = ',sngl(a),'e = ',sngl(e)
         print*,'ap = ',sngl(ap),'ep = ',sngl(ep),'irel = ',sngl(irel)
       end if

c       print*, mass(1)*mass(2)/(mass(1)+mass(2))*vecth
c     &    +(mass(1)+mass(2))*mass(3)/(mass(1)+mass(2)+mass(3))*vecthp

c       print*,'------>',sqrt(1.-e*e)*cos(irel*Pi/180.)
c       print*,inc*180/PI,capom*180/PI,omega*180/PI
c       print*,incp*180/PI,capomp*180/PI,omegap*180/PI
       id = -2
       call io_write_line_r(iu,id,a,e,inc,capom,omega,capm)   
       n = sqrt((mass(1)+mass(2))/a/a/a)    
       ci = cos(inc)                                            
       si = sin(inc)                                            
       co = cos(capom)                                          
       so = sin(capom)     
       rrot = sqrt(rotx(1)*rotx(1)+roty(1)*roty(1)+rotz(1)*rotz(1))
       rot(1) = rrot/n
       irot(1) = acos((si*(rotx(1)*so-roty(1)*co)+ci*rotz(1))/rrot)
       orot(1) = atan2(co*rotx(1)+so*roty(1),
     &              ci*(so*rotx(1)-co*roty(1))+si*rotz(1))
       rrot = sqrt(rotx(2)*rotx(2)+roty(2)*roty(2)+rotz(2)*rotz(2))
       rot(2) = rrot/n
       irot(2) = acos((si*(rotx(2)*so-roty(2)*co)+ci*rotz(2))/rrot)
       orot(2) = atan2(co*rotx(2)+so*roty(2),
     &              ci*(so*rotx(2)-co*roty(2))+si*rotz(2))
c           print*,'irot',irot(1)*180./PI,irot(2)*180./PI
       if (verbose) print*,'rot/n',sngl(rot(1:2)),'irot',
     &                        sngl(irot(1:2)*180.d0/Pi)
       call io_write_line_r(iu,id,rot(1),irot(1),orot(1),rot(2),
     &                                            irot(2),orot(2))
       id = -3
       call io_write_line_r(iu,id,ap,ep,incp,capomp,omegap,capmp)
       close(iu)

       iverb = iverb+1
       end

c************************************************************************

	subroutine io_dump_pl_av(dplfile,time,tstop,dt,dtout,dtdump,
     &              mass,vecte,vecth,vectep,vecthp,rotx,roty,rotz,
     &              inert,krpl5,Qlov,rplsq,fverb,dirs,gname,outfile)

	include 'swift.inc'

c...    Input
	real*8 mass(3),rplsq(3)
	real*8 vecte(3),vecth(3),vectep(3),vecthp(3)
        real*8 rotx(3),roty(3),rotz(3)
        real*8 inert(3),krpl5(3),Qlov(3)
        real*8 time,tstop,dt,dtout,dtdump
	character*(*) dplfile,dirs,gname,outfile
        integer fverb

c...   Internal
	integer j,ierr,i
        real*8 rpl,rgyr,klov
        real*8 a,e,inc,capom,omega,capm,n,rrot
        real*8 ap,ep,incp,capomp,omegap,capmp
        real*8 gm,gmp,ci,co,si,so,rot(2),irot(2),orot(2)

        REAL*8, PARAMETER :: SMASSYR = TWOPI*TWOPI,
     &          DR = 1.7453292519943294d-2,      ! 180/pi 
     &          EPOCH = 2449101.0d0,                      
     &          YEAR = 365.2422d0,                        
     &          AU = 1.495978707d8               ! km     


c-----
c...  Executable code      

        call io_open(7,dplfile,'unknown','formatted',ierr)

	write(7,'(a)')'cd /home/beust/Documents/gl436'
        write(7,'(a)')'./swift_av <<!'
        write(7,*)mass(1)/smassyr
        write(7,*)mass(2)/smassyr
        write(7,*)mass(3)/smassyr

        gm = mass(1)+mass(2)
        gmp = gm + mass(3)
        call orbite(gm,vecte,vecth,a,e,inc,capom,omega)
        call orbite(gmp,vectep,vecthp,ap,ep,incp,capomp,omegap)
        capm = 0.d0
        capmp = 0.d0

        n = sqrt(gm/a/a/a)    
        ci = cos(inc)                                            
        si = sin(inc)                                            
        co = cos(capom)                                          
        so = sin(capom)     
        rrot = sqrt(rotx(1)*rotx(1)+roty(1)*roty(1)+rotz(1)*rotz(1))
        rot(1) = rrot/n
        irot(1) = acos((si*(rotx(1)*so-roty(1)*co)+ci*rotz(1))/rrot)
        orot(1) = atan2(co*rotx(1)+so*roty(1),
     &              ci*(so*rotx(1)-co*roty(1))+si*rotz(1))
        rrot = sqrt(rotx(2)*rotx(2)+roty(2)*roty(2)+rotz(2)*rotz(2))
        rot(2) = rrot/n
        irot(2) = acos((si*(rotx(2)*so-roty(2)*co)+ci*rotz(2))/rrot)
        orot(2) = atan2(co*rotx(2)+so*roty(2),
     &              ci*(so*rotx(2)-co*roty(2))+si*rotz(2))
        rpl = sqrt(rplsq(1))
        rgyr = sqrt(inert(1)/mass(1)/rplsq(1))
        klov = krpl5(1)/rplsq(1)/rplsq(1)/rpl
        write(7,*)a,e,inc/dr,capom/dr,omega/dr,capm
        write(7,*)rpl*au,klov,Qlov(1),rgyr
        write(7,*)rot(1),irot(1)/dr,orot(1)/dr
        write(7,*)rot(2),irot(2)/dr,orot(2)/dr
        write(7,*)ap,ep,incp/dr,capomp/dr,omegap/dr,capmp
        write(7,*)time,tstop,dt
        write(7,*)dtout,dtdump,fverb

        write(7,'(a)')dirs
        write(7,'(a)')gname
        write(7,'(a)')outfile
        write(7,'(a)')'append'

	close(unit = 7)
	return
	end    ! io_dump_pl_av.f
c--------------------------------------------------------------------------

        subroutine io_init_av(mass,vecte,vecth,vectep,vecthp,
     &         rotx,roty,rotz,inert,krpl5,Qlov,rplsq,
     &         t0,tstop,dt,dtout,dtdump,fverb,diro,dirs,gname,
     &         outfile,fopenstat)

	include 'swift.inc'

c...  Outputs: 
	real*8 mass(3),rplsq(3)
	real*8 vecte(3),vecth(3),vectep(3),vecthp(3)
        real*8 rotx(3),roty(3),rotz(3)
        real*8 inert(3),krpl5(3),Qlov(3)
	real*8 t0,tstop,dt,dtout,dtdump
	real*8 rmin,rmax,rmaxu,qmin
        logical*2 lclose
	character*(chlen) outfile,fopenstat,diro,dirs,gname,dataname
        integer fverb

c...  Internals
        real*8 a,e,inc,capom,omega,capm
        real*8 ap,ep,incp,capomp,omegap,capmp
        real*8 ci,si,co,so,com,som,gm,h,n,e1(3),e2(3),e3(3)
        real*8 rpl,rot,orot,irot,klov,rgyr
        real*8 co1,si1,ci1,so1

        REAL*8, PARAMETER :: SMASSYR = TWOPI*TWOPI,
     &          DR = 1.7453292519943294d-2,      ! 180/pi 
     &          EPOCH = 2449101.0d0,                      
     &          YEAR = 365.2422d0,                        
     &          AU = 1.495978707d8               ! km     

 999    format(a)
c-----
c...  Executable code 

        read(*,*)mass(1)
        read(*,*)mass(2)
        read(*,*)mass(3)
        mass = mass*smassyr
        rplsq = 0.d0
        rotx = 0.d0
        roty = 0.d0
        rotz = 0.d0
        krpl5 = 0.d0
        inert = 0.d0
        qlov = 0.d0
        read(*,*)a,e,inc,capom,omega,capm
        omega = omega*dr
        capom = capom*dr
        inc = inc*dr
        ci = cos(inc)
        si = sin(inc)
        co = cos(capom)
        so = sin(capom)
        com = cos(omega)
        som = sin(omega)
        gm = mass(1)+mass(2)
        h = sqrt(a*gm*(1.d0-e*e))
        n = sqrt(gm/a/a/a)
        e1 = (/ co, so, 0.d0 /)
        e2 = (/ -ci*so, co*ci, si /)
        e3 = (/ so*si, -co*si, ci /)
        vecte = (/ e*(com*co-som*ci*so),
     &                e*(com*so+som*ci*co),e*som*si /)      
        vecth = (/ h*si*so,-h*si*co,h*ci /)
        read(*,*)rpl,klov,qlov(1),rgyr
        rpl = rpl/au
        rplsq(1) = rpl*rpl
        krpl5(1) = klov*rplsq(1)*rplsq(1)*rpl
        inert(1) = rgyr*rgyr*mass(1)*rplsq(1)
        read(*,*)rpl,klov,qlov(2),rgyr
        rpl = rpl/au
        rplsq(2) = rpl*rpl
        krpl5(2) = klov*rplsq(2)*rplsq(2)*rpl
        inert(2) = rgyr*rgyr*mass(2)*rplsq(2)
        read(*,*)rot,irot,orot
        irot = irot*dr
        orot = orot*dr
        rot = rot*n
        ci1 = cos(irot)
        si1 = sin(irot)
        co1 = cos(orot)
        so1 = sin(orot)
        rotx(1) = rot*(so1*si1*e1(1)-co1*si1*e2(1)+ci1*e3(1))
        roty(1) = rot*(so1*si1*e1(2)-co1*si1*e2(2)+ci1*e3(2))
        rotz(1) = rot*(so1*si1*e1(3)-co1*si1*e2(3)+ci1*e3(3))
        read(*,*)rot,irot,orot
        irot = irot*dr
        orot = orot*dr
        rot = rot*n
        ci1 = cos(irot)
        si1 = sin(irot)
        co1 = cos(orot)
        so1 = sin(orot)
        rotx(2) = rot*(so1*si1*e1(1)-co1*si1*e2(1)+ci1*e3(1))
        roty(2) = rot*(so1*si1*e1(2)-co1*si1*e2(2)+ci1*e3(2))
        rotz(2) = rot*(so1*si1*e1(3)-co1*si1*e2(3)+ci1*e3(3))
c.....
        read(*,*)ap,ep,incp,capomp,omegap,capmp
        omegap = omegap*dr
        capomp = capomp*dr
        incp = incp*dr
        ci = cos(incp)
        si = sin(incp)
        co = cos(capomp)
        so = sin(capomp)
        com = cos(omegap)
        som = sin(omegap)
        gm = mass(1)+mass(2)+mass(3)
        h = sqrt(ap*gm*(1.d0-ep*ep))
        vectep = (/ ep*(com*co-som*ci*so),
     &                      ep*(com*so+som*ci*co),ep*som*si /)      
        vecthp = (/ h*si*so,-h*si*co,h*ci /)
c
	read(*,*) t0,tstop,dt
	write(*,*) 't0,tstop,dt : ',t0,tstop,dt
	read(*,*) dtout,dtdump,fverb
	write(*,*) 'dtout,dtdump : ',dtout,dtdump

        read(*,999) dirs
        read(*,999) gname
        diro = trim(dirs)//'/'//gname
        dataname = 'mkdir '//trim(diro)
        call system(dataname)
        read(*,999) outfile
        write(*,*) 'outfile : ', trim(diro)//'/'//outfile

        read(*,999) fopenstat
        
	return
	end     ! io_init_av
c____________________________________________________________________________
c
c
C**************************************************************************
C	    		        ORDRE3
C**************************************************************************

      subroutine ordre3(mass,ap,gm,gmp,nn,nnp,alpha3,alpha4,e,ex2,
     &         exq,exq2,excp,excp2,exqp4,exqp5,exqp6,z,
     &         dhdtx,dhdty,dhdtz,drundtx,drundty,drundtz,
     &         dhpdtx,dhpdty,dhpdtz,drunpdtx,drunpdty,drunpdtz)
       
      include 'swift.inc'

c...  Inputs: 
      real*8 mass(3)
      real*8 ap,gm,gmp,nn,nnp,alpha3,alpha4
      real*8 e,ex2,exq,exq2,excp,excp2,exqp4,exqp5,exqp6
      real*8 z(3,3)

c...  Outputs
      real*8 dhdtx,dhdty,dhdtz,drundtx,drundty,drundtz
      real*8 dhpdtx,dhpdty,dhpdtz,drunpdtx,drunpdty,drunpdtz

c...  Internals:
      real*8 zz0,zz1,zz2,zz3,zz4,zz5,zz6,zz7,zz8,zz9,zz10,zz11,dd
      real*8 ff,gm2,gm3,fac1,fac2,fac3,fac4,z232,z332,z222,z322

c----
c...  Executable code 

c.......................................

      gm2 = gm*gm
      gm3 = gm2*gm

      fac1 = 2.d0+5.d0*ex2
      fac2 = 1.d0+6.d0*ex2
      fac3 = 3.d0+4.d0*ex2
      fac4 = -1.d0+8.d0*ex2
      z232 = z(2,3)*z(2,3)
      z332 = z(3,3)*z(3,3)
      z222 = z(2,2)*z(2,2)
      z322 = z(3,2)*z(3,2)

      zz0 = -z(3,3)*z(1,1)*z(2,3)-z(2,1)*z(3,3)*z(1,3)
     &      -z(3,1)*z(2,3)*z(1,3)

      zz1 = -(17.d0*ex2+11.d0)*z(3,1)+5.d0*fac1*z(3,1)*z232
     &       +15.d0*fac3*z(3,1)*z332
     &       +10.d0*fac1*z(3,3)*z(2,3)*z(2,1)
      zz2 = -(27.d0*ex2+1.d0)*z(2,1)+10.d0*fac2*z(3,3)*z(3,1)*z(2,3)
     &      +105.d0*ex2*z232*z(2,1)+5.d0*fac2*z(2,1)*z332

      dd = mass(3)*(mass(1)-mass(2))*alpha3*excp*e/gm/ap/exqp5
      dhdtx = -2.34375d0*dd*exq2*zz0
      dhdty = 0.234375d0*dd*zz1
      dhdtz = -0.234375d0*dd*zz2 

      zz4 = (24.d0*ex2-1.d0)*z(1,1)+105.d0*ex2*z(2,1)*z(1,3)*z(2,3)
     &     +5.d0*(15.d0*ex2+2.d0)*z(3,1)*z(1,3)*z(3,3)
     &     -5.d0*(3.d0*ex2-1.d0)*z(1,1)*z332

      dd = mass(3)*(mass(1)-mass(2))*alpha4*nn*excp*exq/gm2/exqp5

      drundtx = +0.234375d0*dd*zz2     ! zz3 = zz2
c..............        signe change dans drundtx
      drundty = -0.234375d0*dd*zz4
c..............        signa change dans drundty
      drundtz = +2.34375d0*dd*ex2*zz0  ! zz0 = zz9

      zz6 = -exq2*z(1,2)*z(3,3)*z(3,1)-fac2*z(1,3)*z(3,2)*z(3,1)
     &     -7.d0*z(2,2)*z(2,1)*z(1,3)*ex2-exq2*z(1,1)*z(3,3)*z(3,2)

      zz7 = 70.d0*z222*z(1,3)*ex2-(11.d0+87.d0*ex2)*z(1,3)
     &     +105.d0*z232*ex2*z(1,3)
     &     +10.d0*fac2*z(1,3)*z322
     &     +20.d0*exq2*z(1,2)*z(3,3)*z(3,2)
     &     +15.d0*fac3*z332*z(1,3)

      zz8 = fac4*z(1,2)+5.d0*exq2*z(1,2)*z332
     &     +5.d0*fac1*z(3,2)*z(1,3)*z(3,3)
     &     +35.d0*z(2,2)*ex2*z(2,3)*z(1,3)

      ff = mass(1)*mass(2)*(mass(1)-mass(2))

      dd = gmp*ff*alpha3*excp*e/ap/gm3/exqp5

      dhpdtx = -2.34375d0*dd*zz6
      dhpdty = 0.234375d0*dd*zz7
      dhpdtz = -0.234375d0*dd*zz8

      zz10 = fac4*z(1,1)+5.d0*fac1*z(3,1)*z(1,3)*z(3,3)
     &      +35*ex2*z(2,1)*z(1,3)*z(2,3)+5*exq2*z(1,1)*z332

      dd = ff*alpha3*nnp*e/gm3      
      drunpdtx = 0.234375d0*dd*zz8/exqp4
      drunpdty = -0.234375d0*dd*(4.d0*excp2+1.d0)*zz10/exqp6
      drunpdtz = 2.34375d0*dd*excp2*zz6/exqp6
c                         zz10=z11 car vecte.vecth=0
      end 


C**************************************************************************
C	    		        ORDRE4
C**************************************************************************

      subroutine ordre4(mass,ap,gm,gmp,nn,nnp,alpha4,alpha5,e,ex2,
     &         exq,exq2,exq4,excp,excp2,exqp2,exqp6,exqp7,exqp8,
     &         z,dhdtx,dhdty,dhdtz,drundtx,drundty,drundtz,
     &         dhpdtx,dhpdty,dhpdtz,drunpdtx,drunpdty,drunpdtz)
       
      include 'swift.inc'

c...  Inputs: 
      real*8 mass(3)
      real*8 ap,gm,gmp,nn,nnp,alpha4,alpha5
      real*8 e,ex2,exq,exq2,exq4,excp,excp2,exqp2,exqp6,exqp7,exqp8
      real*8 z(3,3)

c...  Outputs
      real*8 dhdtx,dhdty,dhdtz,drundtx,drundty,drundtz
      real*8 dhpdtx,dhpdty,dhpdtz,drunpdtx,drunpdty,drunpdtz

c...  Internals:
      real*8 zz0,zz1,zz2,zz3,zz4,zz5,zz6,zz7,zz8,zz9,zz10,zz11,dd
      real*8 fac1,fac2,fac3,fac4,fac5,fac6,fac7,fac8,gm2,gm4,gm3
      real*8 facp1,facp2,facp3,facp4,z222,z322,z332,z333,z232,z233
      real*8 z223,facz1,facz2,facz3,facz4,facz5,facz6,facz7,facz8
      real*8 facz9,facz10,facz11,facz12,facz13,facz14,facz15,facz16
      real*8 facz17,facz18,facz19,facz20,facz21,facz22,facz23,facz24
      real*8 facz25,facz26,facz27,facz28,facz29,facz30,facz31,facz32
      real*8 facz33,facz34,fach,facr,fachp,facrp,ex4,dl(3),mu,mup

c----
c...  Executable code 

c.......................................

      dd = mass(1)*mass(1)-mass(1)*mass(2)+mass(2)*mass(2)
      ex4 = ex2*ex2
      gm2 = gm*gm
      gm3 = gm2*gm
      gm4 = gm2*gm2

      fac1 = 1.d0+6.d0*ex2
      fac2 = 1.d0+13.d0*ex2
      fac3 = 1.d0+2.d0*ex2
      fac4 = 1.d0+4.d0*ex2 
      fac5 = 3.d0+4.d0*ex2
      fac6 = 1.d0+8.d0*ex4+12.d0*ex2     
      fac7 = -1.d0+2.d0*ex2+20.d0*ex4
      fac8 = 2.d0+ex2

      facp1 = 5.d0*excp2+2.d0
      facp2 = 17.d0*excp2+76.d0*ex4+309.d0*excp2*ex2
     &            +86.d0*ex2+346.d0*ex4*excp2+6.d0
      facp3 = 24.d0+25.d0*excp2
      facp4 = 2.d0+3.d0*excp2

      z222 = z(2,2)*z(2,2)
      z322 = z(3,2)*z(3,2)
      z332 = z(3,3)*z(3,3)
      z333 = z332*z(3,3)
      Z232 = z(2,3)*z(2,3)
      z233 = z(2,3)*z(2,3)*z(2,3)
      z223 = z222*z(2,2)
      facz1 = z(3,3)*z(2,2) + z(2,3)*z(3,2)
      facz2 = z(3,3)*z(1,2)*z(3,2)*z(2,3)
      facz3 = z(3,3)*z(2,2)*z(3,2)*z(1,3)
      facz4 = z(3,3)*z(3,2)*z(2,3)*z(2,2)
      facz5 = z(3,3)*z(3,2)*z(2,3)*z(2,1)
      facz6 = z(3,3)*z(3,1)*z(2,2)*z(2,3)
      facz7 = z333*z(2,3)
      facz8 = z333*z(1,3)
      facz9 = z333*z(3,3)
      facz11 = z333*z(3,2)
      facz12 = z333*z(3,1)
      facz10 = z233*z(2,3)
      facz13 = z233*z(1,3)
      facz14 = z233*z(2,2)
      facz15 = z233*z(2,1)
      facz16 = z233*z(3,3)
      facz17 = z322*z(3,2)*z(3,3)
      facz18 = z223*z(2,3)
      facz19 = z(3,3)*z232*z(1,3)     
      facz20 = z332*z(3,2)*z(1,2)
      facz21 = z322*z(3,3)*z(1,3)
      facz22 = z222*z(1,3)*z(2,3)
      facz23 = z332*z(2,3)*z(1,3)
      facz24 = z222*z232
      facz25 = z332*z232
      facz26 = z322*z332
      facz27 = z332*z(2,3)*z(2,2)
      facz28 = z322*z(2,3)*z(1,3)
      facz29 = z332*z(2,3)*z(2,1)
      facz30 = z(3,3)*z(3,1)*z322
      facz31 = z222*z(2,1)*z(2,3)
      facz32 = z(3,3)*z(3,1)*z232
      facz33 = z(2,2)*z(2,1)*z232
      facz34 = z332*z(3,1)*z(3,2)

      zz0 = -42.d0*excp2*fac1*z(3,3)*z(3,2)*facz1
     &     -294.d0*excp2*ex2*z(2,3)*z(2,2)*facz1
     &     +6.d0*excp2*fac2*z(3,2)*z(2,2)
     &     +3.d0*(6.d0+17.d0*excp2)*fac2*z(3,3)*z(2,3)
     &     -21.d0*facp1*fac1*facz7-147.d0*ex2*facp1*facz16

      zz1 = -294.d0*excp2*ex2*fac3*z(2,2)*z(1,3)*facz1
     &     -147.d0*ex2*facp1*fac3*facz19-42.d0*excp2*exq2*fac1*facz20
     &     -42.d0*excp2*(19.d0*ex2+1.d0+22.d0*ex4)*facz21
     &     -21.d0*facp1*fac6*facz8-6.d0*excp2*fac7*z(3,2)*z(1,2)
     &     +3.d0*facp2*z(3,3)*z(1,3)

      zz2 = -84.d0*excp2*ex2*facz22-14.d0*excp2*(exq2*facz2+fac3*facz3)
     &     -7.d0*facp1*fac3*facz23-2.d0*excp2*fac8*z(2,2)*z(1,2)
     &     +(2.d0+22.d0*ex2+excp2+95.d0*excp2*ex2)*z(2,3)*z(1,3)
     &     -42.d0*excp2*ex2*facz28-21.d0*ex2*facp1*facz13

      zz4 = 30.d0+138.d0*ex2+45.d0*excp2*ex2
     &     -3.d0*excp2+588.d0*excp2*fac4*facz4+1764.d0*ex2*excp2*facz24
     &     +42.d0*facp1*fac5*facz9+441.d0*ex2*facp1*facz10
     &     +147.d0*facp1*fac4*facz25+84.d0*excp2*(ex2+1.d0)*z222
     &     +12.d0*excp2*(20.d0*ex2+1.d0)*z322
     &     -21.d0*(2.d0+44.d0*ex2+excp2+106.d0*excp2*ex2)*z232
     &     -3.d0*(300.d0*excp2*ex2+211.d0*excp2+86.d0+152.d0*ex2)*z332
     &     +168.d0*excp2*fac5*facz26

      zz6 = -294.d0*excp2*ex2*fac3*z(2,2)*z(3,2)*facz1
     &     -42.d0*excp2*fac6*facz17-882.d0*ex4*excp2*facz18
     &     -147.d0*ex2*facp4*fac3*z(3,3)*z(2,3)*facz1
     &     +21.d0*ex2*facp1*(1.d0+11.d0*ex2)*z(2,3)*z(2,2)
     &     +3.d0*facp1*(38.d0*ex4+3.d0+43.d0*ex2)*z(3,3)*z(3,2)
     &     -21.d0*facp4*fac6*facz11-441.d0*ex4*facp4*facz14

      zz7 = -147.d0*ex2*facp4*fac3*facz29
     &     -147.d0*ex2*(2.d0+7.d0*excp2)*fac3*facz32
     &     -294.d0*excp2*ex2*fac3*z(3,1)*z(2,2)*facz1
     &     -21.d0*facp1*fac6*facz12
     &     +21.d0*ex2*(excp2+53.d0*excp2*ex2
     &                         +2.d0+22.d0*ex2)*z(2,3)*z(2,1)
     &     +3.d0*facp2*z(3,3)*z(3,1)-42.d0*excp2*fac6*facz30
     &     -882.d0*ex4*excp2*facz31-441.d0*ex4*facp1*facz15

      zz8 = -147.d0*ex4*facz33-49.d0*ex2*fac3*facz6
     &     -7.d0*fac6*facz34-7.d0*ex2*fac8*z(2,2)*z(2,1)
     &     -fac7*z(3,2)*z(3,1)-49.d0*ex2*fac3*facz5

      zz10 = 48.d0+1056.d0*ex4+408.d0*ex2-30.d0*excp2*ex2+57.d0*excp2
     &     +225.d0*ex4*excp2+1176.d0*ex2*facp1*fac3*facz4
     &     +21.d0*facp3*fac6*facz9+294.d0*ex2*facp3*fac3*facz25
     &     +1764.d0*ex4*facp1*facz24
     &     +84.d0*ex2*facp1*fac8*z222+12*facp1*fac7*z322
     &     -42*ex2*(265.d0*excp2*ex2+5.d0*excp2+16+260.d0*ex2)*z232
     &     -6.d0*(832*ex4+1024.d0*ex2+750.d0*ex4*excp2
     &        +1055.d0*excp2*ex2+76+85.d0*excp2)*z332
     &     +441.d0*ex4*facp3*facz10+84.d0*facp1*fac6*facz26

      fach = 0.05859375d0*mass(3)*dd*alpha4/gm2/ap/exqp7
      facr = 0.05859375d0*mass(3)*dd*alpha5*nn*e*exq/gm3/exqp7
      fachp = 0.05859375d0*gmp*mass(1)*mass(2)*dd*alpha4/ap/gm4
      facrp = 0.05859375d0*mass(1)*mass(2)*dd*alpha4*nnp*excp/gm4

      dhdtx = fach*zz0*exq2
      dhdty = -fach*zz1
      dhdtz = 21.d0*fach*ex2*zz2
      drundtx = -21.d0*facr*zz2
      drundty = facr*zz4
      drundtz = -facr*zz0

      dhpdtx = fachp*zz6/exqp7
      dhpdty = -fachp*zz7/exqp7
      dhpdtz = 6.d0*fachp*excp2*zz8/exqp7
      drunpdtx = -6.d0*facrp*zz8/exqp6
      drunpdty = 0.25d0*facrp*zz10/exqp8
      drunpdtz = -facrp*zz6/exqp8

c      if (e.gt.0.6) then
c        dl = (/ dhpdtx*z(1,1)+dhpdty*z(1,2)+dhpdtz*z(1,3),
c     &              dhpdtx*z(2,1)+dhpdty*z(2,2)+dhpdtz*z(2,3),
c     &              dhpdtx*z(3,1)+dhpdty*z(3,2)+dhpdtz*z(3,3) /)
c        mu = mass(1)*mass(2)/gm
c        mup = mass(3)*gm/gmp 
c        print*,'h ', sngl(mu*dhdtx),sngl(mu*dhdty),sngl(mu*dhdtz)
c        print*,'hp',sngl(mup*dl)
c        print*,'sm',sngl(mu*dhdtx+mup*dl(1)),sngl(mu*dhdty+mup*dl(2)),
c     &                          sngl(mu*dhdtz+mup*dl(3))
c      end if

c      if (e.gt.0.6) then
c        print*,'z',z
c        print*,'zz0',sngl(fac1/exq2/exqp2),sngl(fac2/exq2),
c     &               sngl(fac3/exqp2),sngl(fac4)
c      end if

      end

C**************************************************************************
C	    		        ORDRE5
C**************************************************************************

      subroutine ordre5(mass,ap,gm,gmp,nn,nnp,alpha4,alpha5,e,ex2,
     &         exq,exq2,exq4,excp,excp2,exqp2,exqp6,exqp7,exqp8,
     &         z,dhdtx,dhdty,dhdtz,drundtx,drundty,drundtz,
     &         dhpdtx,dhpdty,dhpdtz,drunpdtx,drunpdty,drunpdtz)
       
      include 'swift.inc'

c...  Inputs: 
      real*8 mass(3)
      real*8 ap,gm,gmp,nn,nnp,alpha4,alpha5
      real*8 e,ex2,exq,exq2,exq4,excp,excp2,exqp2,exqp6,exqp7,exqp8
      real*8 z(3,3)

c...  Outputs
      real*8 dhdtx,dhdty,dhdtz,drundtx,drundty,drundtz
      real*8 dhpdtx,dhpdty,dhpdtz,drunpdtx,drunpdty,drunpdtz

c...  Internals:
      real*8 zz0,zz1,zz2,zz3,zz4,zz5,zz6,zz7,zz8,zz9,zz10,zz11,dd
      real*8 fac1,fac2,fac3,fac4,fac5,fac6,fac7,fac8,fac9,gm2,gm4,gm3
      real*8 fac10,fac11,fac12,fac13,fac14,fac15,fac16,fac17,fac18
      real*8 fac19,fac20,fac21,fac22,fac23,fac24,fac25
      real*8 facp1,facp2,facp3,facp4,facp5,facp6,facp7,facp8,facp9
      real*8 facp10,facp11,facp12,facp13
      real*8 z222,z322,z332,z333,z232,z233,z334,z234,z323,z224,z324
      real*8 z223,facz1,facz2,facz3,facz4,facz5,facz6,facz7,facz8
      real*8 facz9,facz10,facz11,facz12,facz13,facz14,facz15,facz16
      real*8 facz17,facz18,facz19,facz20,facz21,facz22,facz23,facz24
      real*8 facz25,facz26,facz27,facz28,facz29,facz30,facz31,facz32
      real*8 facz33,fach,facr,fachp,facrp,ex4,dl(3),mu,mup,excp4
      real*8 facz34,facz35,facz36,facz37,facz38,facz39,facz40
      real*8 facz41,facz42,facz43,facz44,facz45,facz46,facz47
      real*8 facz48,facz49,facz50
      real*8 facm1,facm2,facm3,facm4,facm5,facm6,facm7,facm8
      real*8 facm9,facm10,facm11,facm12,facm13,facm14,facm15,facm16
      real*8 facm17,facm18,facm19,facm20,facm21,facm22,facm23,facm24
      real*8 facm25,facm26,facm27,facm28,facm29,facm30,facm31,facm32
      real*8 facm33,facm34,facm35,facm36,facm37,facm38,facm39


c----
c...  Executable code 

c.......................................

      dd = mass(1)*mass(1)-mass(1)*mass(2)+mass(2)*mass(2)
      ex4 = ex2*ex2
      excp4 = excp2*excp2
      gm2 = gm*gm
      gm3 = gm2*gm
      gm4 = gm2*gm2

      fac1 = 2.d0*ex2+1.d0
      fac2 = 7.d0*ex2-4.d0
      fac3 = -1.d0+7.d0*ex2
      fac4 = 17.d0*ex2+14.d0*ex4+2.d0
      fac5 = 7.d0*ex2+4.d0 
      fac6 = 5.d0+8.d0*ex4+20.d0*ex2
      fac7 = 124.d0*ex2+91.d0*ex4+16
      fac8 = 36.d0*ex4+3.d0+38.d0*ex2
      fac9 = 2.d0*ex2-1.d0+32.d0*ex4
      fac10 = 8.d0+3.d0*ex2
      fac11 = 3.d0+8.d0*ex2
      fac12 = 16.d0*ex2+16.d0*ex4+1.d0
      fac13 = 12.d0*ex2-1.d0+40.d0*ex4
      fac14 = 5.d0*ex2-1.d0
      fac15 = 8.d0+5.d0*ex2
      fac16 = 15.d0*ex2+2.d0
      fac17 = 18.d0+25.d0*ex2
      fac18 = 25.d0*ex4+78.d0*ex2+2.d0
      fac19 = 6.d0+5.d0*ex2
      fac20 = 5.d0+6.d0*ex2
      fac21 = 1.d0+8.d0*ex2+2.d0*ex4
      fac22 = 7.d0*ex4+4.d0+22.d0*ex2
      fac23 = 2.d0+9.d0*ex2
      fac24 = 2.d0+5.d0*ex4+26.d0*ex2
      fac25 = -1.d0+8.d0*ex4+4.d0*ex2

      facp1 = 15.d0*excp2+8.d0
      facp2 = 8.d0+3.d0*excp2
      facp3 = 8.d0+5.d0*excp2
      facp4 = 8.d0+7.d0*excp2
      facp5 = 24.d0+29.d0*excp2
      facp6 = 3.d0*excp2+8.d0
      facp7 = 7.d0*excp2+8.d0
      facp8 = excp2+1.d0
      facp9 = 8.d0+13.d0*excp2
      facp10 = 5.d0*excp2+4.d0
      facp11 = 2.d0*excp2+1.d0
      facp12 = 18.d0*excp4+8.d0+73.d0*excp2
      facp13 = 42.d0*excp4+8.d0+85.d0*excp2

      facm1 = 8.d0+153.d0*excp2*ex2-9.d0*excp2+40.d0*ex2
      facm2 = 8.d0+40.d0*ex2+23.d0*excp2*ex2+excp2
      facm3 = -excp2+22.d0*excp2*ex2+32.d0*ex2-8.d0
      facm4 = 29.d0*excp2+72.d0*ex2+67.d0*excp2*ex2+24.d0
      facm5 = 2526.d0*excp2*ex2+2133.d0*excp2*ex4-232.d0-712.d0*ex4
     &         -1168.d0*ex2+225.d0*excp2
      facm6 = 265.d0*excp2*ex2+110.d0*excp2*ex4+72.d0+65.d0*excp2
     &         +312.d0*ex2+144.d0*ex4
      facm7 = 16.d0+280.d0*ex4+232.d0*ex2-18.d0*excp2+99.d0*excp2*ex2
     &         +315.d0*excp2*ex4
      facm8 = 161.d0*excp2*ex4+101.d0*excp2*ex2+2.d0*excp2+232.d0*ex2
     &         +280.d0*ex4+16.d0
      facm9 = 272.d0*ex2+8.d0+1832.d0*ex4+50.d0*excp2*ex2
     &         +1535.d0*excp2*ex4-excp2
      facm10 = 8.d0+200.d0*ex2+320.d0*ex4+excp2+184.d0*excp2*ex4
     &          +79.d0*excp2*ex2
      facm11 = 8.d0+335.d0*excp2*ex2+200.d0*ex2+320.d0*ex4+8.d0*excp2
     &          +680.d0*excp2*ex4
      facm12 = 24.d0+435.d0*excp2*ex2+5.d0*excp2+504.d0*ex2
      facm13 = 8.d0+excp2-120.d0*ex2+160.d0*ex4-69.d0*excp2*ex2
     &          +110.d0*excp2*ex4
      facm14 = 8.d0+1600.d0*ex4-192.d0*ex2+1480.d0*excp2*ex4
     &          -60.d0*excp2*ex2-excp2
      facm15 = 32.d0+1278.d0*excp2*ex2+32.d0*excp2+625.d0*excp2*ex4
     &          +528.d0*ex2+280.d0*ex4
      facm16 = 96.d0+345.d0*excp2*ex2+156.d0*excp2+280.d0*ex2
      facm17 = 360.d0*ex2+48.d0+295.d0*excp2*ex2+10.d0*excp2
      facm18 = 528.d0*ex2+280.d0*ex4+438.d0*excp2*ex2+265.d0*excp2*ex4
     &          +32.d0+32.d0*excp2
      facm19 = 22.d0*excp2*ex2+32.d0*ex2-excp2-8.d0
      facm20 = 10.d0*excp2*ex2+32.d0*ex2-7.d0*excp2-8.d0
      facm21 = -8.d0+excp2+10.d0*excp2*ex2+8.d0*ex2
      facm22 = -20.d0*excp2+9.d0*excp2*ex2-64.d0-24.d0*ex2
      facm23 = 8.d0+7.d0*excp2-112.d0*ex4+100.d0*excp2*ex4
     &          -160.d0*ex2+58.d0*excp2*ex2
      facm24 = 48.d0*ex4+192.d0*ex2+32.d0*excp2*ex4+116.d0*excp2*ex2
     &          +17.d0*excp2+24.d0
      facm25 = -38.d0*excp2-32.d0+27.d0*excp2*ex2-12.d0*ex2
      facm26 = -232.d0+5438.d0*excp2*ex2+593.d0*excp2*ex4+1304.d0*ex4
     &          +4208.d0*ex2-223.d0*excp2
      facm27 = 600.d0*ex2-72.d0+555.d0*excp2*ex2-115.d0*excp2
      facm28 = 10.d0*excp2*ex2+16.d0*ex2-4.d0-7.d0*excp2
      facm29 = 10.d0*excp2*ex4-50.d0*excp2*ex2+7.d0*excp2-80.d0*ex2
     &          -56.d0*ex4+4.d0
      facm30 = 32.d0+201.d0*excp2*ex4+1002.d0*excp2*ex2+84.d0*excp2
     &          +56.d0*ex4+176.d0*ex2
      facm31 = 8.d0+27.d0*excp2+320.d0*ex4+72.d0*excp2*ex4
     &          -132.d0*excp2*ex2-64.d0*ex2
      facm32 = 72.d0*ex2+16.d0+423.d0*excp2*ex2+6.d0*excp2
      facm33 = 30.d0*excp2*ex2+32.d0*ex2-21.d0*excp2-8.d0
      facm34 = 322.d0*excp2*ex2+132.d0*excp4*ex2+32.d0*ex2-67.d0*excp2
     &          -8.d0-6.d0*excp4
      facm35 = 8.d0-6.d0*excp4+320.d0*ex4-64.d0*ex2+3448.d0*excp2*ex4
     &          -572.d0*excp2*ex2+61.d0*excp2-120.d0*excp4*ex2
     &          +1776.d0*excp4*ex4
      facm36 = 192.d0*excp4+32.d0+352.d0*excp2+318.d0*excp4*ex4
     &          +56.d0*ex4+176.d0*ex2+1846.d0*excp2*ex2
     &          +876.d0*excp4*ex2+607.d0*excp2*ex4
      facm37 = 32.d0+312.d0*excp4+412.d0*excp2+655.d0*excp2*ex2
     &          +414.d0*excp4*ex2+56.d0*ex2
      facm38 = 354.d0*excp4*ex2+20.d0*excp4+753.d0*excp2*ex2+72.d0*ex2
     &          +16.d0+138.d0*excp2
      facm39 = 32.d0+192.d0*excp4+750.d0*excp4*ex4+56.d0*ex4
     &          +2556.d0*excp4*ex2+352.d0*excp2+2686.d0*excp2*ex2
     &          +176.d0*ex2+823.d0*excp2*ex4

      z222 = z(2,2)*z(2,2)
      z322 = z(3,2)*z(3,2)
      z332 = z(3,3)*z(3,3)
      z232 = z(2,3)*z(2,3)
      z333 = z332*z(3,3)
      z233 = z232*z(2,3)
      z223 = z222*z(2,2)
      z323 = z322*z(3,2)
      z334 = z332*z332
      z234 = z232*z232
      z224 = z222*z222
      z324 = z322*z322


      facz1 = z(3,3)*z(3,1)*z(1,3)*z(3,2)*z(2,2)
      facz2 = z(3,3)*z(3,2)*z(3,1)*z(2,2)*z(2,3)
      facz3 = z(3,1)*z(2,2)*z(3,2)*z(2,3)*z(1,3)
      facz4 = z(3,1)*z(2,2)*z(2,3)*z(3,3)*z(1,3)
      facz5 = z(3,3)*z(3,2)*z(2,3)*z(2,2)*z(1,3)
      facz25 = z(3,3)*z(3,1)*z(2,3)*z(3,2)*z(1,2)
      facz39 = z(3,3)*z(3,2)*z(2,3)*z(2,1)*z(1,3)
      facz6 = z(1,3)*z233*z(3,1)
      facz7 = z(3,3)*z233*z(2,1)
      facz8 = z(3,3)*z(3,1)*z233
      facz9 = z(1,3)*z(2,1)*z233
      facz10 = z233*z(2,2)*z(1,3)
      facz11 = z333*z(1,3)*z(2,1)
      facz12 = z333*z(1,1)*z(2,3)
      facz13 = z333*z(2,3)*z(2,1)
      facz14 = z333*z(2,3)*z(3,1)
      facz15 = z333*z(3,1)*z(1,3)
      facz16 = z333*z(3,1)*z(1,2)
      facz17 = z(1,1)*z333*z(3,2)
      facz18 = z333*z(3,2)*z(1,2)
      facz19 = z(3,2)*z333*z(1,3)
      facz20 = z(1,3)*z(3,1)*z323
      facz21 = z(1,2)*z(3,3)*z323
      facz22 = z(3,3)*z323*z(1,3)
      facz23 = z223*z(2,1)*z(1,3)
      facz24 = z223*z(1,3)*z(2,3)
      facz26 = z(2,3)*z(1,3)*z(3,1)*z322
      facz27 = z332*z(2,3)*z(3,1)*z(1,3)
      facz28 = z232*z(1,3)*z(3,3)*z(2,1)
      facz29 = z332*z(1,3)*z(2,3)*z(2,1)
      facz30 = z332*z(3,2)*z(3,1)*z(1,3)
      facz31 = z332*z(3,1)*z(3,2)*z(2,2)
      facz32 = z222*z(3,1)*z(2,3)*z(1,3)
      facz33 = z(3,3)*z(2,3)*z(2,1)*z222
      facz34 = z322*z(2,3)*z(3,3)*z(3,1)
      facz35 = z(3,1)*z222*z(3,3)*z(1,3)
      facz36 = z222*z(1,3)*z(2,3)*z(2,1)
      facz37 = z(3,1)*z232*z(3,3)*z(1,3)
      facz38 = z(1,3)*z322*z(3,3)*z(3,1)
      facz40 = z322*z(3,1)*z(1,2)*z(3,3)
      facz41 = z(3,1)*z232*z(3,2)*z(1,3)
      facz42 = z(3,2)*z(3,1)*z(1,3)*z222
      facz43 = z232*z(2,2)*z(1,3)*z(2,1)
      facz44 = z222*z(3,3)*z(3,2)*z(1,2)
      facz45 = z(3,2)*z232*z(3,3)*z(1,3)
      facz46 = z322*z(2,2)*z(1,3)*z(2,3)
      facz47 = z332*z(2,2)*z(1,3)*z(2,3)
      facz48 = z(3,3)*z(3,2)*z(1,3)*z222
      facz49 = z332*z(3,1)*z(3,2)*z(1,2)
      facz50 = z(2,3)*z(2,1)*z(1,3)*z222

      zz0 = -9.d0*ex2*facp1*facz6-3.d0*facp2*fac1*facz11
     &      -72.d0*z(3,3)*excp2*ex2*z222*z(1,3)*z(2,1)
     &      +facm1*z(3,1)*z(2,3)*z(1,3)
     &      +4.d0*excp2*fac2*z(3,1)*z(2,2)*z(1,2)
     &      +facm2*z(2,1)*z(3,3)*z(1,3)-facm3*z(3,3)*z(1,1)*z(2,3)
     &      -36.d0*excp2*fac1*facz1-108.d0*ex2*excp2*facz32
     &      -36.d0*excp2*facz25*exq2-36.d0*ex2*facp3*facz28
     &      -3.d0*facp2*facz12*exq2-12.d0*excp2*fac3*facz26
     &      -3.d0*facm4*facz27

      zz1 = -168.d0*excp2*fac4*z322*z(3,1)*z232
     &      +504.d0*excp2*ex2*fac5*facz33-facm5*z(3,1)
     &      -42.d0*facm6*z(3,1)*z332+105.d0*facp4*fac6*z334*z(3,1)
     &      -14.d0*facm7*z(3,1)*z232+1008.d0*excp2*fac4*facz2
     &      +420.d0*excp2*fac6*z332*z(3,1)*z322
     &      +42.d0*facp5*fac4*z332*z(3,1)*z232
     &      +28.d0*excp2*fac7*z(3,1)*z222+84.d0*excp2*fac8*z(3,1)*z322
     &      +84.d0*facp6*fac4*facz13-28.d0*facm8*z(3,3)*z(2,3)*z(2,1)
     &      +756.d0*excp2*ex2*fac5*z222*z(3,1)*z232
     &      +63.d0*ex2*facp1*fac5*z(3,1)*z234
     &      +252.d0*ex2*facp3*fac5*facz7

      zz2 = 28.d0*excp2*fac9*z(3,1)*z(3,2)*z(2,2)+facm9*z(2,1)
     &      -14.d0*facm10*z(2,1)*z332
     &      -28.d0*facm11*z(3,3)*z(3,1)*z(2,3)
     &      +3465.d0*ex4*facp7*z234*z(2,1)
     &      +84.d0*excp2*ex2*fac10*z222*z(2,1)
     &      +504.d0*excp2*ex2*fac11*z(3,1)*z232*z(3,2)*z(2,2)
     &      +2016d0*excp2*ex2*fac11*z(3,1)*z222*z(3,3)*z(2,3)
     &      +252.d0*excp2*fac12*facz31
     &      +168.d0*excp2*fac12*facz34+21.d0*facp2*fac12*z(2,1)*z334
     &      +378.d0*ex2*facp2*fac11*z332*z232*z(2,1)
     &      +672.d0*facp8*fac12*facz14-42.d0*ex2*facm12*z232*z(2,1)
     &      +13860.d0*z222*z(2,1)*excp2*ex4*z232
     &      +252.d0*ex2*facp9*fac11*facz8

      zz4 = -14.d0*facm13*z(1,1)*z332
     &      +28.d0*excp2*fac13*z(3,1)*z(3,2)*z(1,2)+facm14*z(1,1)
     &      -252.d0*excp2*fac14*facz49*exq2
     &      -21.d0*facp2*fac14*z(1,1)*z334*exq2
     &      -7.d0*facm15*z(3,1)*z(3,3)*z(1,3)
     &      +84.d0*excp2*ex2*fac15*z(2,1)*z(2,2)*z(1,2)
     &      +756.d0*excp2*ex2*fac15*facz35
     &      +13860.d0*ex4*excp2*facz36+63.d0*ex2*facm16*facz37
     &      -21.d0*ex2*facm17*z(1,3)*z(2,3)*z(2,1)
     &      +756.d0*excp2*ex2*fac16*facz3+63.d0*ex2*facp2*fac17*facz29
     &      +84.d0*excp2*fac18*facz38+3465.d0*ex4*facp7*facz9
     &      +21.d0*facm18*facz15


      zz6 = -9.d0*ex2*facp6*fac19*facz39-6.d0*excp2*fac12*facz20
     &      -facm3*z(1,1)*z(3,3)*z(3,2)*exq2
     &      -facm20*z(1,2)*z(3,3)*z(3,1)*exq2
     &      -24.d0*excp2*fac1*facz40*exq2-18.d0*ex2*facm21*facz41
     &      -36.d0*excp2*ex2*fac11*facz42
     &      +3.d0*ex2*facm22*z(1,3)*z(2,2)*z(2,1)
     &      +facm23*z(1,3)*z(3,2)*z(3,1)-3.d0*facp3*exq4*facz16
     &      -9.d0*ex2*facp3*fac19*facz4-198.d0*ex4*facz23*excp2
     &      -3.d0*facp6*exq4*facz17-99.d0*ex4*facp3*facz43
     &      -72.d0*excp2*ex2*z(3,3)*z222*z(3,1)*z(1,2)*exq2
     &      -3.d0*facm24*facz30

      zz7 = 1008.d0*ex2*facp10*fac19*facz5+336.d0*facp10*exq4*facz18
     &      +3465.d0*ex4*facp7*z(1,3)*z234
     &      +168.d0*excp2*fac12*z(1,3)*z324
     &      +5544.d0*ex4*z224*z(1,3)*excp2
     &      +672.d0*excp2*fac1*facz21*exq2
     &      -168.d0*ex2*facm25*z(1,3)*z222
     &      +105.d0*facp7*fac6*z334*z(1,3)-facm26*z(1,3)
     &      -42.d0*ex2*facm27*z232*z(1,3)
     &      +1008.d0*excp2*ex2*fac11*z322*z(1,3)*z222
     &      +112.d0*facm28*z(1,2)*z(3,3)*z(3,2)*exq2
     &      -56.d0*facm29*z(1,3)*z322
     &      +630.d0*ex2*facp7*fac20*z332*z232*z(1,3)
     &      +5544.d0*ex4*facp10*z232*z(1,3)*z222
     &      -1008.d0*ex2*facp10*z322*z232*z(1,3)*exq2
     &      -42.d0*facm6*z332*z(1,3)
     &      +504.d0*facp10*fac21*z332*z(1,3)*z322
     &      +2016.d0*excp2*ex2*facz44*exq2

      zz8 = 21.d0*facp1*fac22*facz19-7.d0*facm30*z(3,3)*z(3,2)*z(1,3)
     &      +63.d0*ex2*facp1*fac5*facz45
     &      +84.d0*excp2*ex2*fac10*z222*z(1,2)+facm31*z(1,2)
     &      +756.d0*excp2*ex2*fac23*facz46
     &      +63.d0*ex2*facp1*fac19*facz47
     &      +756.d0*excp2*ex2*fac10*facz48
     &      -21.d0*ex2*facm32*z(2,2)*z(1,3)*z(2,3)
     &      +14.d0*facm33*z(1,2)*z332*exq2+252.d0*excp2*fac24*facz22
     &      +21.d0*facp1*exq4*z334*z(1,2)+693.d0*ex4*facp1*facz10
     &      +8316.d0*ex4*excp2*facz24+84.d0*excp2*fac25*z322*z(1,2)
     &      +756.d0*excp2*exq4*z322*z332*z(1,2)

      zz10 = 14.d0*facm34*z(1,1)*z332*exq2
     &      +84.d0*excp2*facp11*fac25*z(3,1)*z(3,2)*z(1,2)
     &      +facm35*z(1,1)+756.d0*excp2*facp11*exq4*facz49
     &      -7.d0*facm39*z(3,1)*z(3,3)*z(1,3)
     &      +84.d0*excp2*ex2*facp11*fac10*z(2,1)*z(2,2)*z(1,2)
     &      +756.d0*excp2*ex2*facp11*fac10*facz35
     &      +8316.d0*excp2*ex4*facp11*facz50
     &      +21.d0*facp12*exq4*z(1,1)*z334+21.d0*facm36*facz15
     &      +63.d0*ex2*facm37*facz37
     &      -21.d0*ex2*facm38*z(1,3)*z(2,3)*z(2,1)
     &      +693.d0*ex4*facp13*facz9
     &      +756.d0*excp2*ex2*facp11*fac23*facz3
     &      +63.d0*ex2*facp12*fac19*facz29
     &      +252.d0*excp2*facp11*fac24*facz38



      fach = 0.05859375d0*mass(3)*dd*alpha4/gm2/ap/exqp7
      facr = 0.05859375d0*mass(3)*dd*alpha5*nn*e*exq/gm3/exqp7
      fachp = 0.05859375d0*gmp*mass(1)*mass(2)*dd*alpha4/ap/gm4
      facrp = 0.05859375d0*mass(1)*mass(2)*dd*alpha4*nnp*excp/gm4

      dhdtx = fach*zz0*exq2
      dhdty = -fach*zz1
      dhdtz = 21.d0*fach*ex2*zz2
      drundtx = -21.d0*facr*zz2
      drundty = facr*zz4
      drundtz = -facr*zz0

      dhpdtx = fachp*zz6/exqp7
      dhpdty = -fachp*zz7/exqp7
      dhpdtz = 6.d0*fachp*excp2*zz8/exqp7
      drunpdtx = -6.d0*facrp*zz8/exqp6
      drunpdty = 0.25d0*facrp*zz10/exqp8
      drunpdtz = -facrp*zz6/exqp8

c      if (e.gt.0.6) then
c        dl = (/ dhpdtx*z(1,1)+dhpdty*z(1,2)+dhpdtz*z(1,3),
c     &              dhpdtx*z(2,1)+dhpdty*z(2,2)+dhpdtz*z(2,3),
c     &              dhpdtx*z(3,1)+dhpdty*z(3,2)+dhpdtz*z(3,3) /)
c        mu = mass(1)*mass(2)/gm
c        mup = mass(3)*gm/gmp 
c        print*,'h ', sngl(mu*dhdtx),sngl(mu*dhdty),sngl(mu*dhdtz)
c        print*,'hp',sngl(mup*dl)
c        print*,'sm',sngl(mu*dhdtx+mup*dl(1)),sngl(mu*dhdty+mup*dl(2)),
c     &                          sngl(mu*dhdtz+mup*dl(3))
c      end if

c      if (e.gt.0.6) then
c        print*,'z',z
c        print*,'zz0',sngl(fac1/exq2/exqp2),sngl(fac2/exq2),
c     &               sngl(fac3/exqp2),sngl(fac4)
c      end if

      end

C
C
C---------------------------------------------------------------------------
C		Adaptation du pas du RUNGE-RUTTA Cash-Karp
C---------------------------------------------------------------------------
C
	SUBROUTINE RKUTTAQS(mass,inert,krpl5,Qlov,X,DPAS,Y,EPS)

	IMPLICIT NONE	

	INTEGER*4	N,		! La dimension de l'espace.
     &                  NST,            ! L'ordre de la methode
     &			K,I,J		! Indices
	PARAMETER(	N=18,NST=4       )

        real*8 mass(3),inert(3),krpl5(3),Qlov(3)
        real*8          ff(10),mu,gm
        REAL*8		YSAV(N),	! La solution non incrementee.
     &			Y(N),		! La solution a incrementer
     &			DXY(N),		! Les derivees en un point.
     &			YSCAL(N),	! Le Vecteur de test pour Delta
     &                  YERR(N),        ! L'erreur
     &			X,DX,		! La variable (le temps)
     &			XSAV,		! La variable non modifiee.
     &			DPAS,H,		! Le pas.
     &			ERRMAX,		! Erreur maximale
     &			EPS		! Tolerance.
        REAL*8, PARAMETER ::  
     &                  PSHRK=-1./NST,  ! Puissances de modification du pas
     &                  PGROW=-1./(NST+1.)     
        REAL*8, PARAMETER :: SAFETY = 0.9d0
        REAL*8, PARAMETER :: ERRCON = 1.89d-4

        LOGICAL :: HDID                 ! Vrai si un pas a ete effectue

      gm = mass(1)+mass(2)
      ff(1) = (1.d0+mass(2)/mass(1))*krpl5(1)
      ff(2) = (1.d0+mass(1)/mass(2))*krpl5(2)
      ff(3) = krpl5(1)*mass(2)/mass(1)/Qlov(1)
      ff(4) = krpl5(2)*mass(1)/mass(2)/Qlov(2)
      mu =  mass(2)*mass(1)/gm
      ff(5) = mu/gm  
      ff(6) = -mu/inert(1)
      ff(8) = -mu/inert(2)
      ff(7) = gm
      ff(9) = ff(1)*mass(2)
      ff(10) = ff(2)*mass(1)

        CALL derivs_tid(mass,ff,y,dxy)
        YSCAL=ABS(Y)+ABS(DPAS*DXY)+1d-30
        YSAV=Y
	XSAV=X
	H=DPAS
        CALL RKUTTACK(mass,ff,H,X,Y,YERR)
	ERRMAX=0.
	DO I=1,N
          ERRMAX=MAX(ERRMAX,ABS(YERR(I)/YSCAL(I)))
	END DO
	ERRMAX=ERRMAX/EPS
	IF (ERRMAX.GT.1.) THEN
	  DPAS=0.9*DPAS*(ERRMAX**PSHRK)
          Y=YSAV
          X=XSAV
          HDID=.FALSE.
          DX=0.
          print*,dpas
        ELSE
	  IF (ERRMAX.GT.ERRCON) THEN
            DPAS=0.9*DPAS*(ERRMAX**PGROW)
          ELSE
            DPAS=DPAS*4.
          END IF
          HDID=.TRUE.
        END IF
        END SUBROUTINE RKUTTAQS

C
C----------------------------------------------------------------------------
C		Procedure RUNGE-KUTTA Cash-Karp
C----------------------------------------------------------------------------
C
        SUBROUTINE RKUTTACK(mass,ff,DPAS,X,Y,YERR) 

	IMPLICIT NONE	

        REAL*8          DPAS		! Le pas.
        real*8          ff(10),mass(3)
        REAL*8          X               ! La variable
        REAL*8          Y(18)            ! La solution a incrementer
        REAL*8          YERR(18)         ! L'erreur
     	INTEGER*4       I		! Indice
        REAL*8 ::	K(1:6,18),       ! Coefficients Ki du Runge-Kutta
     &                  XST,            ! Point courant
     &			YST(18), 	! La solution courante.
     &			DXY(18)		! Les derivees en un point.
        REAL*8, DIMENSION(6), PARAMETER ::
     &    A=(/ 0.,1./5.,3./10.,3./5.,1.,7./8. /)
        REAL*8, DIMENSION(6), PARAMETER ::
     &    C=(/ 37./378.,0.,250./621.,125./594.,0.,512./1771. /)
        REAL*8, DIMENSION(6), PARAMETER ::
     &    DC=(/ C(1)-2825./27648.,C(2),C(3)-18575./48384.,
     &          C(4)-13525./55296.,C(5)-277./14336.,C(6)-1./4. /)
        REAL*8, DIMENSION(6,5), PARAMETER ::
     &    B=RESHAPE((/ 0.,1./5.,3./40.,3./10.,-11./54.,1631./55296.,
     &                 0.,0.,9./40.,-9./10.,5./2.,175./512., 
     &                 0.,0.,0.,6./5.,-70./27.,575./13824.,
     &                 0.,0.,0.,0.,35./27.,44275./110592.,
     &                 0.,0.,0.,0.,0.,253./4096. /),(/ 6,5 /))
c        K(0,:)=0.
	DO I=1,6
	  XST=X+DPAS*A(I)
          IF (I.gt.1) THEN
            YST=Y+MATMUL(B(I,1:I-1),K(1:I-1,:))
          ELSE
            YST=Y
          END IF
          CALL derivs_tid(mass,ff,yst,dxy)
          K(I,:)=DPAS*DXY
	END DO
        Y=Y+MATMUL(C,K(1:6,:))
        YERR=MATMUL(DC,K(1:6,:))
	X=X+DPAS
	END

