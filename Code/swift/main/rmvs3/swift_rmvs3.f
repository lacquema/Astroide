c**********************************************************************
c		      SWIFT_RMVS3.F
c**********************************************************************
c
c                 INCLUDES CLOSE ENCOUNTERS
c                 To run, need 3 input files. The code prompts for
c                 the file names, but examples are :
c
c                   parameter file like       param.in
c		    planet file like          pl.in
c                   test particle file like   tp.in
c
c Authors:  Hal Levison \& Martin Duncan
c Date:    8/25/94
c Last revision: Feb 14 2023 H Beust (save mvsfile)

     
	include '../../sub/swift.inc'

	real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

	real*8 mass(NPLMAX),j2rp2,j4rp4
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

	integer istat(NTPMAX,NSTAT),i1st
	integer nbod,ntp,nleft
	integer iflgchk,iub,iuj,iud,iue
        real*8 rstat(NTPMAX,NSTATR)

	real*8 t0,tstop,dt,dtout,dtdump
	real*8 t,tout,tdump,tfrac,eoff

	real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
        logical*2 lclose 

	character*(chlen) outfile,inparfile,inplfile,intpfile,
     &       fopenstat,diro,dirs,gname,genfile,dataname
        integer*8 iloop,noutloop,ndumploop


c-----
c...    Executable code 

c...    print version number
        call util_version

c Get data for the generic file
        write(*,*) 'Enter name of generic file : '
        read(*,999) genfile

c Get data for the run and the test particles
	write(*,*) 'Enter name of parameter data file : '
	read(*,999) inparfile
        call io_init_param_hb(inparfile,t0,tstop,dt,dtout,dtdump,
     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,diro,dirs,
     &         gname,outfile,fopenstat)

        write(*,*) 'iflgchk=',iflgchk

c Prompt and read name of planet data file
	write(*,*) ' '
	write(*,*) 'Enter name of planet data file : '
	read(*,999) inplfile
999 	format(a)
	call io_init_pl(inplfile,lclose,iflgchk,nbod,mass,xh,yh,zh,
     &       vxh,vyh,vzh,rplsq,j2rp2,j4rp4)

c Get data for the run and the test particles
	write(*,*) 'Enter name of test particle data file : '
	read(*,999) intpfile
	call io_init_tp(intpfile,ntp,xht,yht,zht,vxht,vyht,
     &               vzht,istat,rstat)

c       Save input file for further continuation (H Beust Feb 14, 2023)
        call io_write_mvsfile(diro,genfile)
c
	
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
        end if


c Initialize initial time and times for first output and first dump
	t = t0
c	tout = t0 + dtout
c	tdump = t0 + dtdump

c New HBeust 28 mar. 2011
        iloop = 0_8
        noutloop = nint(dtout/dt,8)
        ndumploop = nint(dtdump/dt,8)
c

        iub = 20
        iuj = 30
        iud = 40
        iue = 60


c...    Do the initial io write
        if(btest(iflgchk,0))  then ! bit 0 is set
           write(*,*) 'bit 1 is set'
           call io_write_frame(t0,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
     &         xht,yht,zht,vxht,vyht,vzht,istat,
     &         trim(diro)//'/'//outfile,iub,fopenstat)
        endif
        
        if(btest(iflgchk,1))  then ! bit 1 is set
           write(*,*) 'bit 2 is set'
           call io_write_frame_r(t0,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
     &         xht,yht,zht,vxht,vyht,vzht,istat,
     &         trim(diro)//'/'//outfile,iub,fopenstat)
        endif
        if(btest(iflgchk,2))  then    ! bit 2 is set
           write(*,*) 'bit 3 is set'
           eoff = 0.0d0
           call anal_energy_write(t0,nbod,mass,j2rp2,j4rp4,xh,yh,zh,vxh,
     &          vyh,vzh,iue,fopenstat,eoff)
        endif
        if(btest(iflgchk,3))  then    ! bit 3 is set
           write(*,*) 'bit 4 is set'
           call anal_jacobi_write(t0,nbod,ntp,mass,xh,yh,zh,vxh,
     &        vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,2,iuj,fopenstat)
        endif

c...    must initize discard io routine
        if(btest(iflgchk,4))  then ! bit 4 is set
           write(*,*) 'bit 5 is set'
           call io_discard_write(0,t,nbod,ntp,xh,yh,zh,vxh,vyh,
     &          vzh,xht,yht,zht,vxht,vyht,vzht,istat,rstat,iud,
     &           trim(diro)//'/'//'discard.out',fopenstat,nleft)
        endif

        nleft = ntp
	i1st = 0
 
c***************here's the big loop *************************************
        write(*,*) ' ************** MAIN LOOP ****************** '

	  do while ( (t .le. tstop) .and. 
     &       ((ntp.eq.0).or.(nleft.gt.0)) )

 	     call rmvs3_step(i1st,t,nbod,ntp,mass,j2rp2,j4rp4,
     &         xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,
     &         vzht,istat,rstat,dt)

c	     t = t + dt
             iloop = iloop+1_8
             t = t0 + dble(iloop)*dt

             if(btest(iflgchk,4))  then    ! bit 4 is set
                call discard(t,dt,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
     &               xht,yht,zht,vxht,vyht,vzht,rmin,rmax,rmaxu,
     &               qmin,lclose,rplsq,istat,rstat)
                call io_discard_write(1,t,nbod,ntp,xh,yh,zh,vxh,vyh,
     &               vzh,xht,yht,zht,vxht,vyht,vzht,istat,rstat,iud,
     &               trim(diro)//'/'//'discard.out',fopenstat,nleft)
             else
                nleft = ntp
             endif

c if it is time, output orb. elements, 
c	  if(t .ge. tout) then 
          if (mod(iloop,noutloop).eq.0_8) then

             if(btest(iflgchk,0))  then    ! bit 0 is set
                call  io_write_frame(t,nbod,ntp,mass,xh,yh,zh,vxh,
     &               vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,
     &               trim(diro)//'/'//outfile,iub,fopenstat)
              endif
             if(btest(iflgchk,1))  then    ! bit 1 is set
                call  io_write_frame_r(t,nbod,ntp,mass,xh,yh,zh,vxh,
     &               vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,
     &               trim(diro)//'/'//outfile,iub,fopenstat)
              endif

c	    tout = tout + dtout
	  endif

c If it is time, do a dump
c          if(t.ge.tdump) then
          if (mod(iloop,ndumploop).eq.0_8) then

             tfrac = (t-t0)/(tstop-t0)
             write(*,998) t,tfrac,nleft
 998         format(' Time = ',1p1e12.5,': fraction done = ',0pf5.3,
     &            ': Number of active tp =',i6)
	     call io_dump_pl(trim(diro)//'/dump_pl.dat',nbod,mass,
     &           xh,yh,zh,vxh,vyh,vzh,lclose,iflgchk,rplsq,j2rp2,j4rp4)
	     call io_dump_tp(trim(diro)//'/dump_tp.dat',ntp,
     &           xht,yht,zht,vxht,vyht,vzht,istat,rstat)
	     call io_dump_param(trim(diro)//'/'//'dump_param.dat',
     &          t,tstop,dt,dtout,dtdump,iflgchk,rmin,rmax,rmaxu,
     &	        qmin,lclose,dirs,gname,outfile) 
c	     tdump = tdump + dtdump

             if(btest(iflgchk,2))  then    ! bit 2 is set
                call anal_energy_write(t,nbod,mass,j2rp2,j4rp4,
     &               xh,yh,zh,vxh,vyh,vzh,iue,fopenstat,eoff)
             endif
             if(btest(iflgchk,3))  then    ! bit 3 is set
                call anal_jacobi_write(t,nbod,ntp,mass,xh,yh,zh,vxh,
     &               vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,2,
     &               iuj,fopenstat)
             endif

	  endif

	enddo
c********** end of the big loop from time 't0' to time 'tstop'

c Do a final dump for possible resumption later 

	call io_dump_pl(trim(diro)//'/dump_pl.dat',nbod,mass,xh,yh,zh,
     &            vxh,vyh,vzh,lclose,iflgchk,rplsq,j2rp2,j4rp4)
	call io_dump_tp(trim(diro)//'/dump_tp.dat',ntp,xht,yht,zht,
     &              vxht,vyht,vzht,istat,rstat)
	call io_dump_param(trim(diro)//'/'//'dump_param.dat',
     &          t,tstop,dt,dtout,dtdump,iflgchk,rmin,rmax,rmaxu,
     &	        qmin,lclose,dirs,gname,outfile)
        call util_exit(0)
        end    ! swift_rmvs3.f
c---------------------------------------------------------------------
