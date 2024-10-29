c**********************************************************************
c		      SWIFT_HJS_TMRP.F
c**********************************************************************
c  This code is a realization of Hierarchical Jacobi Symplectic
c  n-body mapping method for hierarchical stellar systems
c  (Beust 2003), with Poincare time regularization
c         Source : Mikkola, 1997, CeMDA 67, 147
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
c Date:    Sep. 27, 2006
c Last revision:
c Remark : based on swift_hjs.f

     
	include 'swift.inc'

	real*8 xjt(NTPMAX),yjt(NTPMAX),zjt(NTPMAX)
	real*8 vxjt(NTPMAX),vyjt(NTPMAX),vzjt(NTPMAX)

	real*8 mass(NPLMAX),coef(NPLMAX)
	real*8 xj(NPLMAX),yj(NPLMAX),zj(NPLMAX)
	real*8 vxj(NPLMAX),vyj(NPLMAX),vzj(NPLMAX)

	integer istat(NTPMAX,NSTAT),i1st,istati(NSTAT)
        integer oloc(NPLMAX,NPLMAX),oloct(NPLMAX,NTPMAX)
	integer nbod,ntp,nleft
	integer iflgchk,iub,iuj,iud,iue,i,j
        real*8 rstat(NTPMAX,NSTATR),rstati(NSTATR)
        real*8 umat(NPLMAX,NPLMAX),mat(NPLMAX,NPLMAX)
        real*8 umatp(NPLMAX,NTPMAX),matp(NPLMAX,NTPMAX)
        real*8 mu(NPLMAX),eta(NPLMAX),etatp(NTPMAX)

	real*8 t0,tstop,h,dt,dtout,dtdump
	real*8 t,tout,tdump,tfrac,eoff,p0

	real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
        logical*2 lclose 

	character*(chlen) outfile,inparfile,inplfile,intpfile,
     &            fopenstat,diro,dirs,gname,dataname,genfile

        integer ialpha
        real*8 a,e,inc,capom,omega,capm
        logical ok

c-----
c...    Executable code 

c...    print version number
        call util_version

c Get data for generic file
        write(*,*) 'Enter name of generic file'
        read(*,999) genfile 

c Get data for the run and the test particles
	write(*,*) 'Enter name of parameter data file : '
	read(*,999) inparfile
	call io_init_param(inparfile,t0,tstop,h,dtout,dtdump,
     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,diro,dirs,
     &         gname,outfile,fopenstat)

c Prompt and read name of massive bodies data file
	write(*,*) ' '
	write(*,*) 'Enter name of massive bodies data file : '
	read(*,999) inplfile
999 	format(a)
	call io_init_pl_hjs_tmrp(inplfile,lclose,iflgchk,nbod,oloc,
     &       coef,mass,eta,mu,mat,umat,xj,yj,zj,vxj,vyj,vzj,rplsq)

c Get data for the run and the test paqrticles
	write(*,*) 'Enter name of test particle data file : '
	read(*,999) intpfile
        call io_init_tp_hjs(intpfile,nbod,ntp,oloc,oloct,mass,umat,
     &         etatp,matp,umatp,xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,rstat)

c Copy initial files input into work directory
        if ((fopenstat(1:3).eq.'new')
     &       .or.(fopenstat(1:3).eq.'NEW')) then
          dataname = 'cp '//trim(inparfile)//' '//trim(diro)
          call system(dataname)
          dataname = 'cp '//trim(inplfile)//' '//trim(diro)
          call system(dataname)
          dataname = 'cp '//trim(intpfile)//' '//trim(diro)
          call system(dataname)
          dataname = 'cp '//trim(genfile)//' '//trim(diro)
          call system(dataname)
          dataname = 'cp matpass.dat '//trim(diro)
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
        if(btest(iflgchk,0))  then ! bit 0 is set
           call io_write_frame_hjs(t0,nbod,ntp,oloc,matp,umat,
     &           mass,eta,mu,xj,yj,zj,vxj,vyj,vzj,etatp,xjt,yjt,zjt,
     &           vxjt,vyjt,vzjt,istat,trim(diro)//'/'//outfile,
     &           iub,fopenstat)
        endif
        if(btest(iflgchk,1))  then ! bit 1 is set
           call io_write_frame_r_hjs(t0,nbod,ntp,oloc,matp,umat,
     &           mass,eta,mu,xj,yj,zj,vxj,vyj,vzj,etatp,xjt,yjt,zjt,
     &           vxjt,vyjt,vzjt,istat,trim(diro)//'/'//outfile,
     &           iub,fopenstat)
        endif
        if(btest(iflgchk,2))  then    ! bit 2 is set
           eoff = 0.0d0
           call anal_energy_write_hjs(t,nbod,umat,mass,xj,yj,zj,
     &                         vxj,vyj,vzj,iue,fopenstat,diro,eoff)
        endif
c        if(btest(iflgchk,3))  then    ! bit 3 is set
c           call anal_jacobi_write(t0,nbod,ntp,mass,xh,yh,zh,vxh,
c     &        vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,2,iuj,fopenstat)
c        endif            No Jacobi intagral in HJS...

c...    must initize discard io routine
        if(btest(iflgchk,4))  then ! bit 4 is set
           call io_discard_write(0,t,nbod,ntp,xj,yj,zj,vxj,vyj,
     &          vzj,xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,rstat,iud,
     &          trim(diro)//'/'//'discard.out',fopenstat,nleft)
        endif

Compute the energy -> p0

        call anal_energy_hjs(nbod,umat,mass,xj,yj,zj,
     &                                       vxj,vyj,vzj,eoff,p0)
        p0 = -p0


        nleft = ntp
	i1st = 0

c***************here's the big loop *************************************
        write(*,*) ' ************** MAIN LOOP ****************** '

          ok = .true.
c	  do while ( (t .le. tstop) .and. 
c     &       ((ntp.eq.0).or.(nleft.gt.0)) )
          do while (ok)

c             call step_dkd_hjs(1,i1st,t,nbod,ntp,oloc,oloct,mat,umat,
             call step_kdk_hjs_tmrp(i1st,nbod,ntp,oloc,oloct,coef,
     &     mat,umat,matp,umatp,mass,eta,mu,etatp,p0,xj,yj,zj,
     &     vxj,vyj,vzj,xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,rstat,h,dt)	

	     t = t + dt
cccccccccccccccccccccc

             if (t.ge.tout) then 
             print*,'t=',t

             do j=2,nbod
               call orbel_xv2el(xj(j),yj(j),zj(j),vxj(j),vyj(j),
     &           vzj(j),eta(j)+mu(j),ialpha,a,e,inc,capom,
     &           omega,capm)
               print*,j,sngl(a),sngl(e)
!     &            ,sngl(inc),sngl(capm*180./PI)
               ok = ok.and.(e.lt.1.0d0)
             end do
             end if
             ok = (ok.and.(t.le.tstop))

ccccccccccccccccccc




             if (btest(iflgchk,4))  then ! bit 4 is set
               do i=1,ntp 
                 if (istat(i,1).eq.0) THEN
                   do j=1,NSTAT
                     istati(j) = istat(i,j)
                     rstati(j) = rstat(i,j)
                   end do
                   call discard_hjs(t,dt,nbod,i,mass,xj,yj,zj,vxj,
     &     vyj,vzj,xjt(i),yjt(i),zjt(i),vxjt(i),vyjt(i),vzjt(i),
     &     umatp(1,i),etatp(i),rmin,rmax,rmaxu,qmin,lclose,rplsq,
     &     istati,rstati)
                   do j=1,NSTAT
                     istat(i,j) = istati(j)
                     rstat(i,j) = rstati(j)
                   end do
                 end if
               end do
               call io_discard_write(1,t,nbod,ntp,xj,yj,zj,vxj,vyj,
     &               vzj,xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,rstat,iud,
     &               trim(diro)//'/'//'discard.out',fopenstat,nleft)
             else
                nleft = ntp
             endif

c if it is time, output orb. elements, 
	  if (t .ge. tout) then 
            if (btest(iflgchk,0))  then    ! bit 0 is set
              call io_write_frame_hjs(t,nbod,ntp,oloc,matp,umat,
     &           mass,eta,mu,xj,yj,zj,vxj,vyj,vzj,etatp,xjt,yjt,zjt,
     &           vxjt,vyjt,vzjt,istat,trim(diro)//'/'//outfile,
     &           iub,fopenstat)
            endif
            if (btest(iflgchk,1))  then ! bit 1 is set
              call io_write_frame_r_hjs(t,nbod,ntp,oloc,matp,umat,
     &           mass,eta,mu,xj,yj,zj,vxj,vyj,vzj,etatp,xjt,yjt,zjt,
     &           vxjt,vyjt,vzjt,istat,trim(diro)//'/'//outfile,
     &           iub,fopenstat)
            endif

	    tout = tout + dtout
	  endif

c If it is time, do a dump
          if(t.ge.tdump) then

             tfrac = (t-t0)/(tstop-t0)
             write(*,998) t,tfrac,nleft
 998         format(' Time = ',1p1e12.5,': fraction done = ',0pf5.3,
     &            ': Number of active tp =',i4)
             call io_dump_pl_hjs_tmrp(trim(diro)//'/'//'dump_pl.dat',
     &                    nbod,oloc,coef,mass,umat,
     &                    xj,yj,zj,vxj,vyj,vzj,lclose,iflgchk,rplsq)
             call io_dump_tp_hjs(trim(diro)//'/'//'dump_tp.dat',
     &          nbod,ntp,matp,xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,rstat)
	     call io_dump_param(trim(diro)//'/'//'dump_param.dat',
     &           t,tstop,dt,dtout,
     &           dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)
	     tdump = tdump + dtdump

             if (btest(iflgchk,2))  then    ! bit 2 is set
               call anal_energy_write_hjs(t,nbod,umat,mass,xj,yj,zj,
     &                           vxj,vyj,vzj,iue,fopenstat,diro,eoff)
             endif
c        if(btest(iflgchk,3))  then    ! bit 3 is set
c           call anal_jacobi_write(t0,nbod,ntp,mass,xh,yh,zh,vxh,
c     &        vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,2,iuj,fopenstat)
c        endif            No Jacobi intagral in HJS...

	  endif

	enddo
c********** end of the big loop from time 't0' to time 'tstop'

c Do a final dump for possible resumption later 

        call io_dump_pl_hjs_tmrp(trim(diro)//'/'//'dump_pl.dat',
     &                 nbod,oloc,coef,mass,umat,
     &                 xj,yj,zj,vxj,vyj,vzj,lclose,iflgchk,rplsq)
        call io_dump_tp_hjs(trim(diro)//'/'//'dump_tp.dat',
     &                 nbod,ntp,matp,
     &                 xjt,yjt,zjt,vxjt,vyjt,vzjt,istat,rstat)
	call io_dump_param(trim(diro)//'/'//'dump_param.dat',
     &                 t,tstop,dt,dtout,
     &         dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)

        call util_exit(0)
        end    ! swift_hjs.f
c---------------------------------------------------------------------
