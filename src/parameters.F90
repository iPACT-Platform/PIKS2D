!/------------------------------------------------------------\
!                                                             |
!  $$$$$$$\ $$$$$$\ $$\   $$\  $$$$$$\   $$$$$$\  $$$$$$$\    |
!  $$  __$$\\_$$  _|$$ | $$  |$$  __$$\ $$  __$$\ $$  __$$\   |
!  $$ |  $$ | $$ |  $$ |$$  / $$ /  \__|\__/  $$ |$$ |  $$ |  |
!  $$$$$$$  | $$ |  $$$$$  /  \$$$$$$\   $$$$$$  |$$ |  $$ |  |
!  $$  ____/  $$ |  $$  $$<    \____$$\ $$  ____/ $$ |  $$ |  |
!  $$ |       $$ |  $$ |\$$\  $$\   $$ |$$ |      $$ |  $$ |  |
!  $$ |     $$$$$$\ $$ | \$$\ \$$$$$$  |$$$$$$$$\ $$$$$$$  |  |
!  \__|     \______|\__|  \__| \______/ \________|\_______/   |
!                                                             |
!       Parallel Image-based Kinetic Solver (2D)              |
!\------------------------------------------------------------/

! The MIT License (MIT)

! Copyright (c) 2017-2020 
!   Lianhua Zhu <zhulianhua121@gmail.com> 
!   and Minh-Tuan Ho <minhtuanho.vn@gmail.com>
!   and Yonghao Zhang <y.h.zhang168@gmail.com>

! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:

! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

!-------------------------------------------------------------------------------
! module    : parameters
!-------------------------------------------------------------------------------
! This is a module reading simulation parameters specified by users
! in "para.in" file for 2D DVM solver.
! For details:
!
! [1]   M.T. Ho, L. Zhu, L. Wu, P. Wang, Z. Guo, Z.-H. Li, Y. Zhang
!       "A multi-level parallel solver for rarefied gas flows in porous media"
! 		Computer Physics Communications, 234 (2019), pp. 14-25
!
!! @param imageFileName : file name of the digital porous image
!! @param Nx, Ny   : number of pixels in x, y direction
!! @param wallExtOrder  : order of extrapolation scheme for incomming f at solid 
!! @param fluidLayer 	: number of fluid layers will be added at (each) inlet, outlet
!!						  allows periodic BC available (see page 21 of Ref.[1])
!! @param Ref_L 		: reference length defined by Nx/Ny 
!! @param Nc_fundamental: half of the number of the discrete velocities in one axis
!!						  The total number of discrete velocities is Nv=(2*Nc_fundamental)^3
!! @param halfRange 	: half-range or full-range Gauss-Hermit is chosen
!!						  halfRange=T means half-range, halfRange=F means full-range
!! @param mpi_xdim, mpi_ydim, mpi_zdim : number of MPI (processes) subdomains in x,y,z direction, Section 5.2 of Ref.[1]
!! @param maxStep 		: maximum number of iterations
!! @param chkConvergeStep : iteration interval for checking convergence criteria Eq.(18) of Ref.[1]	
!! @param saveStep 		: iteration interval for exporting output files	
!! @param eps 			: convergence tolerance Eq.(18) of Ref.[1]	
!! @param saveLast 		: whether data from the last iteration is save or not
!! @param saveFormat 	: flow field format (1) VTI (2) Tecplot (3) VTK	
!! @param Kn		 	: Knudsen number defined in Eq.(1) of [1]	
!! @param pressDrop		: pressure drop applied at inlet/outlet Eq.(10) of [1]
!! @param accom		 	: tangential momentum accommodation coefficient (TMAC)  Eq.(8) of [1]				
!-------------------------------------------------------------------------------

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
        read(unit=PARAFILE,nml=physicalNml,iostat=ios)
        read(unit=PARAFILE,nml=velocityNml,iostat=ios)
        read(unit=PARAFILE,nml=mpiNml,iostat=ios)
        read(unit=PARAFILE,nml=solverNml,iostat=ios)
        read(unit=PARAFILE,nml=flowNml,iostat=ios) 

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
