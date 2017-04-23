module parameters
use physicalGrid
use velocityGrid
use mpiParams
use flow
use solver
implicit none

namelist /physicalNml/ imageFileName, Nx, Ny, wallExtOrder
namelist /velocityNml/ Nc_fundamental, halfRange
namelist /mpiNml/ mpi_xdim, mpi_ydim
namelist /solverNml/ maxStep, chkConvergeStep, saveStep, eps
namelist /flowNml/ Kn, pressDrop, accom
! file units
integer, parameter :: PARAFILE = 10

contains 
    subroutine initParams
        integer :: ios
        ! read file called "para.in" using namelist of Fortran 90
        open(unit=PARAFILE,file='para.in',status='old',iostat=ios)
        if (ios /= 0) then
            print*,'ERROR: could not open namelist file'
            stop
        end if
        ! read data into the declared namelist
        read(UNIT=PARAFILE,NML=physicalNml,IOSTAT=ios)
        read(UNIT=PARAFILE,NML=velocityNml,IOSTAT=ios)
        read(UNIT=PARAFILE,NML=mpiNml,IOSTAT=ios)
        read(UNIT=PARAFILE,NML=solverNml,IOSTAT=ios)
        read(UNIT=PARAFILE,NML=flowNml,IOSTAT=ios) 
        if (ios /= 0) then
            print*,'ERROR: could not read example namelist'
            stop
        else 
            close(PARAFILE)
        end if
    end subroutine initParams


    subroutine printParams
        ! print parameters
        print*, "========== Parameters ================"
        print*, "imageFileName = ", imageFileName
        print*, "Nx = ", Nx
        print*, "Ny = ", Ny
        print*, "wallExtOrder = ", wallExtOrder
        print*, "Nc_fundamental = ", Nc_fundamental
        print*, "halfRange = ", halfRange
        print*, "mpi_xdim = ", mpi_xdim
        print*, "mpi_ydim = ", mpi_ydim
        print*, "maxStep = ", maxStep
        print*, "chkConvergeStep = ", chkConvergeStep
        print*, "saveStep = ", saveStep
        print*, "eps = ", eps
        print*, "Kn = ", Kn
        print*, "pressDrop = ", pressDrop
        print*, "accom = ", accom
        print*, "======================================"
    end subroutine

end module parameters
