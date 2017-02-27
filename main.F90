PROGRAM main
!Common Variables
USE flow
USE MPIParams
USE MPI ! this is system mpi
USE solver
IMPLICIT NONE

! Local variables
INTEGER :: MPI_ERR, MPI_PROVIDED

! Initialize MPI environment
CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, MPI_PROVIDED, MPI_ERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, MPI_ERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc, MPI_ERR)

! setup discrete velocity grid
CALL setupVelocityGrid
! read and create virtual CPU grid
CALL setupVirtualProcessGrid
! setup global and local grid system
CALL setupPhysicalGrid
! allocate flow data array and initialize
CALL setupFlow

! set error
error = 1.D0
! Main iteration loop
DO iStep = 1, MaxStep
! Save data if required
    CALL iterate
    !if(proc==master) PRINT*, "STEP: ", iStep
    IF ( MOD(iStep,interval) == 0 ) CALL chkConverge
    IF ( MOD(iStep,interval) == 0 ) CALL saveFlowField
! Test flow field convergence
    IF ( error <= eps ) then
        EXIT
    ENDIF
END DO

! Save final data
CALL saveFlowField

! Free memory, close MPI environment and end program
CALL memFree
CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_ERR)
CALL MPI_FINALIZE(MPI_ERR)

END PROGRAM
