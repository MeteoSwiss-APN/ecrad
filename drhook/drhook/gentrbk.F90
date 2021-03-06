! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!-- Generic traceback calls here

SUBROUTINE GENTRBK_DUMMY
END SUBROUTINE GENTRBK_DUMMY

SUBROUTINE INTEL_TRBK()
#ifdef __INTEL_COMPILER
USE IFCORE
#endif
USE MPL_DATA_MODULE, ONLY : MPL_RANK
CHARACTER*80 MESSAGE
LOGICAL :: DONE_TRACEBACK = .FALSE.
INTEGER :: MYPROC,MYTHREAD

#ifdef _OPENMP
INTEGER,EXTERNAL :: OMP_GET_THREAD_NUM
#endif

IF(DONE_TRACEBACK) THEN
  WRITE(0,*) "INTEL_TRBK already called"
  RETURN
ENDIF

MYPROC=MPL_RANK
#ifdef _OPENMP
MYTHREAD=OMP_GET_THREAD_NUM() + 1
#else
MYTHREAD=1
#endif

#ifndef BOM
  WRITE(MESSAGE,'(A,I4,A,I2,A)') &
  &           "Process ",MYPROC," thread ",MYTHREAD, &
  &           " calling tracebackqq from intel_trbk()"
#ifdef __INTEL_COMPILER
! Commented out, as TRACEBACKQQ is strangely missing from IFCORE module
! on leap42 ifort version 18.1
!  CALL TRACEBACKQQ(MESSAGE, USER_EXIT_CODE=-1)
#endif
#endif
#ifdef LINUX
  WRITE(0,*) "Process ",MYPROC," thread ",MYTHREAD, &
 &           " calling linux_trbk from intel_trbk()"
  CALL LINUX_TRBK() ! See ifsaux/utilities/linuxtrbk.c
#endif
DONE_TRACEBACK=.TRUE.
END SUBROUTINE INTEL_TRBK

#ifndef VPP
SUBROUTINE ERRTRA
END SUBROUTINE ERRTRA
#endif

#ifdef NECSX
SUBROUTINE NECSX_TRBK(CDMESS)
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: CDMESS
CALL MESPUT(CDMESS, LEN(CDMESS), 1)
CALL DBX_TRBK()
END SUBROUTINE NECSX_TRBK

SUBROUTINE NECSX_TRBK_FL(CDMESS, CDFILENAME, KLINENO)
USE PARKIND1  ,ONLY : JPIM
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: CDMESS
CHARACTER(LEN=*), INTENT(IN) :: CDFILENAME
INTEGER(KIND=JPIM), INTENT(IN) :: KLINENO
CHARACTER(LEN=LEN(CDMESS)+LEN(CDFILENAME)+30) CLOCAL
WRITE(CLOCAL,'(a," at ",a,":",i6.6)') TRIM(CDMESS),TRIM(CDFILENAME),KLINENO
CALL NECSX_TRBK(TRIM(CLOCAL))
END SUBROUTINE NECSX_TRBK_FL
#endif
