C ----------------------------------------------
C STPOUT: terminate output and close open files.
C ----------------------------------------------

      SUBROUTINE symbatr_io_STPOUT
 
	INCLUDE './inc/treedefs.f'
 
        WRITE (ULOG, '(//,1X,''CALCULATION TERMINATED'')')
C       ---------------------
C       Close the open files.
C       ---------------------
        CLOSE(UNIT = UBODO)
        CLOSE(UNIT = UBODF)
        CLOSE(UNIT = ULOG)
      END

C -------------------------------------------------------------
C TERROR: terminate the program as the result of a fatal error,
C closing the output files and writing timing information.
C -------------------------------------------------------------

      SUBROUTINE symbatr_io_TERROR(MSG)
      CHARACTER*(*) MSG

	INCLUDE './inc/treedefs.f'

C       ------------------------------------
C       Write error message to the log file.
C       ------------------------------------
        CALL symbatr_io_OUTERR(MSG)
C       ---------------------------------------------------
C       Output timing data, close files, terminate the run.
C       ---------------------------------------------------
c        CALL OUTCPU
        CALL symbatr_io_STPOUT
	STOP
      END

C -----------------------------------------------------------
C OUTERR: output error messages to the log file and terminal.
C -----------------------------------------------------------

      SUBROUTINE symbatr_io_OUTERR(MSG)
      CHARACTER*(*) MSG

	INCLUDE './inc/treedefs.f'
 
C       ------------------------------------------------------
C       Write the message, surrounded by stars for visibility.
C       ------------------------------------------------------
        WRITE (ULOG, '(/,1X,72(''*''))')
        WRITE (ULOG, '(/,A)') MSG
        WRITE (ULOG, '(/,1X,72(''*''))')
C       ---------------------------
C       Repeat message to terminal.
C	---------------------------
        WRITE (UTERM, '(/,1X,72(''*''))')
        WRITE (UTERM, '(/,A)') MSG
        WRITE (UTERM, '(/,1X,72(''*''))')
      END

