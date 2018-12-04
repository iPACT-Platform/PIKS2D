module parameters
use physicalGrid
use velocityGrid
use mpiParams
use flow
use solver
use string_utility_module
implicit none

namelist /physicalNml/ imageFileName, Nx, Ny, xPadRatio, wallExtOrder
namelist /velocityNml/ Nc_fundamental, halfRange
namelist /mpiNml/ mpi_xdim, mpi_ydim, block_repx, block_repy
namelist /solverNml/ maxStep, chkConvergeStep, saveStep, eps, saveFormat
namelist /flowNml/ allKnStr, pressDrop, accom
! file units
integer, parameter :: PARAFILE = 10

contains 
    subroutine initParams
        integer :: ios
        integer :: i0, ii, i
        character(kind=CK,len=:), allocatable :: KniStr
        

        ! set default nml variables
        block_repx = 1
        block_repy = 1
        saveFormat = 1 ! default saving format is vti
        xPadRatio = 0.d0

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
            print*,'ERROR: could not read namelist, may be format is wrong'
            stop
        else 
            close(PARAFILE)
        end if

        ! map varialbes, for weak scaling study
        Nx_base = Nx
        Ny_base = Ny
        Nx = block_repx * Nx
        Ny = block_repy * Ny

        nKn = str_count_tokens(allKnStr)

        ! allocate allKn
        allocate(allKn(nKn))
        ! convert from allKnStr to allKn
        call str_parse_all_token_to_real(trim(allKnStr), allKn)

        ! create directories for the serial of Knudsen numbers
        ii = 1
        do i = 1, nKn
           call str_parse_next_token(allKnStr, ii, KniStr)
           ! only master rank create dir for data output
           if(proc==master) then
               call system("mkdir "//"Kn"//KniStr)
           endif
           ! if dir exist, loop will continuue
        end do

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
        print*, "block_repx = ", block_repx
        print*, "block_repy = ", block_repy
        print*, "maxStep = ", maxStep
        print*, "chkConvergeStep = ", chkConvergeStep
        print*, "saveStep = ", saveStep
        print*, "eps = ", eps
        print*, "saveFormat = ", saveFormat
        print*, "allKnStr = ", allKnStr
        print*, "nKn = ", nKn
        print*, "pressDrop = ", pressDrop
        print*, "accom = ", accom
        print*, "======================================"
    end subroutine

end module parameters
