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
! Program    : main
!-------------------------------------------------------------------------------
! This is the main file of 2D DVM parallel solver. 
! For details:
!
! [1]   M.T. Ho, L. Zhu, L. Wu, P. Wang, Z. Guo, Z.-H. Li, Y. Zhang
!       "A multi-level parallel solver for rarefied gas flows in porous media"
! 		Computer Physics Communications, 234 (2019), pp. 14-25
!
! A digital image (2D array of binary data) of porous medium will be read.
! Linearized BGK kinetic model of Boltzmann equation will be solved to find
! gas apparent permeability as a function of Knudsen number.
!-------------------------------------------------------------------------------

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
