!----------------------------------------------------------------------
!> @brief MPI related parameters and information exchange buffer arrays
!> @details
!! @param nprocs  : total number of processors
!! @param proc    : processor ID
!! @param vproc   : processor ID in virtual grid
!! @param mpi_dim : virtual grid partition scheme (1->stripes, 2->boxes, 3->cubes)
!----------------------------------------------------------------------
MODULE mpiParams
IMPLICIT NONE
SAVE


! Constant tags used in the mpi exchanges
INTEGER, PARAMETER :: TAG1 = 1, TAG2 = 2, TAG3 = 3, TAG4 = 4

! Communication parameters
INTEGER :: nprocs, proc, vproc
INTEGER, parameter :: mpi_xdim = 1
INTEGER, parameter :: mpi_ydim = 1
INTEGER :: east, west, north, south, MPI_COMM_VGRID
INTEGER, PARAMETER :: master  = 0
INTEGER, PARAMETER :: mpi_dim = 2


! Information exchange buffers (x direction)
double precision, ALLOCATABLE, DIMENSION(:) :: f1_west_snd, f1_east_snd
double precision, ALLOCATABLE, DIMENSION(:) :: f1_west_rcv, f1_east_rcv

! Information exchange buffers (y direction)
double precision, ALLOCATABLE, DIMENSION(:) :: f1_south_snd, f1_north_snd
double precision, ALLOCATABLE, DIMENSION(:) :: f1_south_rcv, f1_north_rcv

contains
!-------------------------------------------------------------------------------
! Subroutine : setupVirtualProcessGrid
!-------------------------------------------------------------------------------
!> @file
!! Creates a virtual cartesian grid (2D) linking all MPI-allocated processors.
!> @details
!! Creates a virtual cartesian grid (2D) linking all MPI-allocated processors.
!! The domain limits corresponding to each processor [(xl,xu),(yl,yu)] are
!! assigned and the rank of each processor and its neigbors in all streaming
!! directions are obtained. By default the processor network is periodic and
!! reordered for faster communication
!! @note
!! Some MPI implementations do not reorder correctly. This may affect the
!! performance but not the results.
!! @warning
!! If you use OpenMPI (v1.2.5) you may have to comment out 'USE MPI' and use
!! instead 'INCLUDE "mpif.h"' because a current bug in OpenMPI prevents the
!! command MPI_CART_CREATE() from being recognized.
    subroutine setupVirtualProcessGrid    
        !Variables to be set
        use physicalGrid, only : xl, xlg, xmax, xmin, xu, xug, yl, ylg, ymax, ymin, yu, yug, ghostLayers
        use MPI
        IMPLICIT NONE
        
        !Local variables
        INTEGER :: complete, direction, partial, shift
        INTEGER :: xeast, xwest, ynorth, ysouth, MPI_ERR
        INTEGER, DIMENSION(1:mpi_dim) :: dims, mpi_coords
        LOGICAL, DIMENSION(1:mpi_dim) :: periodic
        LOGICAL :: reorder
        
        !Initialize data for domain partitioning. Defaults are:
        !Partitioning is periodic in all dimensions (periodic = .true.)
        !CPUs are reordered in the grid for proximity (reorder = .true.)
        dims(1)     = mpi_xdim
        dims(2)     = mpi_ydim
        periodic(1) = .true.
        periodic(2) = .false.
        reorder     = .true.
        
        !Create the new virtual connectivity grid
        CALL MPI_CART_CREATE(MPI_COMM_WORLD, mpi_dim, dims, periodic, reorder, MPI_COMM_VGRID, MPI_ERR)
        
        write(*,*) "I'm ",  proc

        !Get this processor ID within the virtual grid
        CALL MPI_COMM_RANK(MPI_COMM_VGRID, vproc, MPI_ERR)
        !write(*,*) "vproc = ",  vproc

        CALL MPI_CART_COORDS(MPI_COMM_VGRID, vproc, mpi_dim, mpi_coords, MPI_ERR)
        !PRINT*, "After first mpi_cart_coords", proc
        
        !------- Compute the limits [(xl,xu),(yl,yu)] assigned to this processor ------
        !Partitioning in the x direction
        complete = (xmax - xmin) / dims(1)
        partial  = (xmax - xmin) - complete*dims(1)
        IF(mpi_coords(1) + 1 <= partial) THEN
          xl = xmin + (complete + 1)*mpi_coords(1)
          xu = xmin + (complete + 1)*(mpi_coords(1) + 1) - 1
        ELSE
          xl = xmin + complete*mpi_coords(1) + partial
          xu = xmin + complete*(mpi_coords(1) + 1) + partial - 1
        END IF
        IF(MOD(mpi_coords(1) + 1,dims(1)) == 0) xu = xu + 1
        
        !Partitioning in the y direction
        complete = (ymax - ymin) / dims(2)
        partial  = (ymax - ymin) - complete*dims(2)
        IF(mpi_coords(2) + 1 <= partial) THEN
          yl = ymin + (complete + 1)*mpi_coords(2)
          yu = ymin + (complete + 1)*(mpi_coords(2) + 1) - 1
        ELSE
          yl = ymin + complete*mpi_coords(2) + partial
          yu = ymin + complete*(mpi_coords(2) + 1) + partial - 1
        END IF
        IF(MOD(mpi_coords(2) + 1,dims(2)) == 0) yu = yu + 1
        
        !Ghost layers
        xlg = xl - ghostLayers
        xug = xu + ghostLayers
        ylg = yl - ghostLayers
        yug = yu + ghostLayers

        !------- Determine neighbours of this processor -------------------------------
        !MPI_CART counts dimensions using 0-based arithmetic so that
        !direction = 0 -> x  |  direction = 1 -> y
        !Ranks of neighbours of this processor in the x and y directions
        shift     = 1
        direction = 0
        CALL MPI_CART_SHIFT(MPI_COMM_VGRID, direction, shift, west, east, MPI_ERR)
        direction = 1
        CALL MPI_CART_SHIFT(MPI_COMM_VGRID, direction, shift, south, north, MPI_ERR)
        ! west/east and sourth/north has been assigned above
        
        !Coordinates of neighbours of this processor in the x and y directions
        !CALL MPI_CART_COORDS(MPI_COMM_VGRID, east, mpi_dim, mpi_coords, MPI_ERR)
        !xeast = mpi_coords(1)
        !PRINT*, "After 2 mpi_cart_coords", proc
        !CALL MPI_CART_COORDS(MPI_COMM_VGRID, west, mpi_dim, mpi_coords, MPI_ERR)
        !xwest = mpi_coords(1)
        !PRINT*, "After 3 mpi_cart_coords", proc
        !CALL MPI_CART_COORDS(MPI_COMM_VGRID, north, mpi_dim, mpi_coords, MPI_ERR)
        !ynorth = mpi_coords(2)
        !PRINT*, "After 4 mpi_cart_coords", proc
        !CALL MPI_CART_COORDS(MPI_COMM_VGRID, south, mpi_dim, mpi_coords, MPI_ERR)
        !ysouth = mpi_coords(2)
        !RINT*, "After 5 mpi_cart_coords", proc
        
    end subroutine setupVirtualProcessGrid
end module mpiParams
