c*************************************************************************
c                            RMVS3_ELOG_ONE.F
c*************************************************************************
c This subroutine keeps a log of the close encounters using istat (1 tp)
c
c             Input:
c                 icflg         ==> ecounters? = 1 Yes, in inner region only
c                                              =  0 No (integer scalar)  
c                 ienc,iencio,       ==> ienc = 0 if tp j not involved in enc 
c                                   in inner region: = planet# if it is. 
c                                     (integer scalar)
c                 istat           ==>  status of the test particle
c                                      istat(1) = 0 ==> active:  = 1 not
c                                      istat(NSTATP+n-1) is number of
c                                          inner encounters that tp had
c                                          with planet n 
c             Output:
c                 istat           ==>  status of the test particles (updated)
c
c
c Remarks: Based on rmvs3_elog
c Authors:  Hervé Beust
c Date:    Apr 20, 2023
c Last revision: 

      subroutine rmvs3_elog_one(icflg,iencio,ienc,istat)

      include '../swift.inc'
      include '../rmvs/rmvs.inc'

c...  Inputs Only: 
      integer ntp,icflg,ienc

c...  Inputs and Outputs:
      integer istat(NSTAT),iencio

c...  Internals
      integer i,i1st,np

c----
c...  Executable code 

c...  check for new encounters
      if (icflg.eq.1) then
         if ((iencio.eq.0).and.(ienc.ne.0)) then 
               np = NSTATP + ienc - 1
               istat(np) = istat(np) + 1
          endif
      endif

c.... save the current value
      iencio = ienc 

      return
      end                       ! rmvs3_elog_one
c------------------------------------------------------

