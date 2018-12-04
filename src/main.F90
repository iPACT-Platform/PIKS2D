module Fortran_Sleep
   use, intrinsic :: iso_c_binding, only: c_int

   implicit none

   interface

      !  should be unsigned int ... not available in Fortran
      !  OK until highest bit gets set.
      function FortSleep (seconds)  bind ( C, name="sleep" )
          import
          integer (c_int) :: FortSleep
          integer (c_int), intent (in), VALUE :: seconds
      end function FortSleep

   end interface

end module Fortran_Sleep

PROGRAM main
!Common Variables
USE parameters
USE flow
USE MPIParams
USE solver
use, intrinsic :: iso_c_binding, only: c_int
use Fortran_Sleep

IMPLICIT NONE
include "mpif.h"

! Local variables
INTEGER :: MPI_ERR, MPI_PROVIDED
INTEGER(c_int) :: sleep
double precision :: startTime, endTime
integer :: kni, str_start_i ! used to parse token


! Initialize MPI environment
CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, MPI_PROVIDED, MPI_ERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, MPI_ERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc, MPI_ERR)

! initialize (read) user input parmeters
call initParams
! check if read ok?
if (proc == master) then
  call printParams
endif

! setup discrete velocity grid
CALL setupVelocityGrid
! read and create virtual CPU grid
CALL setupVirtualProcessGrid
! setup global and local grid system
CALL setupPhysicalGrid
! allocate flow data array and initialize
CALL setupFlow
! save nodeCounts, switch off here
!call saveNodeCounts

! wait for debuger to attcah
!sleep = FortSleep(20)
!write(*,*) sleep

startTime = MPI_Wtime()
! Main iteration loop
! for each Kn in allKn

str_start_i = 1 ! used to parse KniStr

do  kni = 1, nKn
    Kn = allKn(kni)
    mu = dsqrt(PI)/2.0d0/Kn

    ! now get the KniStr, the results for each Kn will be saved in its
    ! own directory
    call str_parse_next_token(allKnStr, str_start_i, KniStr)
    dataSaveDir = "Kn"//KniStr
    if (nKn .eq. 1) then
        ! if only one Kndsen number, will store in current directory
        dataSaveDir = "."
    endif

    ! set error
    error = 1.D0

    do iStep = 1, maxStep
    ! Save data if required
        CALL iterate
        IF ( MOD(iStep,chkConvergeStep) == 0 ) CALL chkConverge
        if ( MOD(iStep,saveStep) == 0 )  then
            call saveFlowField(saveFormat)
        endif
        if ( error <= eps ) then
            call saveFlowField(saveFormat)
            if(proc ==  master) then
                print*, "Converged, exit"
            endif
            exit
        endif
    enddo

    if(proc == master) then 
        print*, "Reaches maximum of steps, end of ", KniStr
    endif
    call saveFlowField(saveFormat)

enddo ! end of each Kn

! calculating wall time
endTime = MPI_Wtime()
if(proc==master) then
    write(*,'(A,ES11.2, A, I6, A, ES15.6)') "Walltime= ", endTime - startTime, &
    ", Steps= ", iStep, ", K= ", permeability
endif

! Free memory, close MPI environment and end program
CALL memFree
CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_ERR)
CALL MPI_FINALIZE(MPI_ERR)

END PROGRAM

! Minh
