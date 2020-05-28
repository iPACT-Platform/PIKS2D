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
! module    : flow
!-------------------------------------------------------------------------------
! This is a module for the flow parameters of 2D DVM parallel solver. 
! For details:
!
! [1]   M.T. Ho, L. Zhu, L. Wu, P. Wang, Z. Guo, Z.-H. Li, Y. Zhang
!       "A multi-level parallel solver for rarefied gas flows in porous media"
! 		Comput. Phys. Commun., 234 (2019), pp. 14-25
!
!	Initial conditions for velocity distribution function, velocity are set.
!-------------------------------------------------------------------------------

module flow
use velocityGrid, only: PI
implicit none
save

integer, parameter :: CK = selected_char_kind('DEFAULT')
! flow parameters, to be read from NML: flowNml
character(len=200) :: allKnStr
double precision :: pressDrop, accom
double precision :: mu


double precision, dimension(:), allocatable :: allKn
! Kn : Knudsen number defined in Eq.(1) of [1]
double precision :: Kn
character(kind=CK,len=:), allocatable :: KniStr
character(kind=CK,len=:), allocatable :: dataSaveDir
integer :: nKn
! f1(spatial_id,velocity_id,sweeppath_id) : velocity distribution function in Eq.(5) of [1]
double precision, DIMENSION(:,:), ALLOCATABLE :: f1
! Rho : number density in Eq.(6) of [1]
! Ux, Uy : two component of velocity vector U in Eq.(6) of [1]
double precision, DIMENSION(:), ALLOCATABLE :: Rho, Ux, Uy
double precision :: mass

contains
    subroutine setupFlow     
        use physicalGrid, only: Ntotal, Nxtotal, Nytotal, ghostLayers
        use velocityGrid, only: Nc
        use mpiParams, only: f1_west_snd, f1_east_snd, f1_west_rcv, f1_east_rcv, &
                             f1_south_snd, f1_north_snd, f1_south_rcv, f1_north_rcv
        implicit none

        allocate(f1(Ntotal,Nc))
        allocate(Rho(Ntotal), Ux(Ntotal), Uy(Ntotal))

        f1=0.d0
        Rho = 0.d0
        Ux = 0.d0
        Uy = 0.d0

        mass = 1.d0

        ! Buffer for exchange data of flull domain of the velocity
        allocate( f1_west_snd(Nc*Nytotal*ghostLayers) ) 
        allocate( f1_west_rcv(Nc*Nytotal*ghostLayers) ) 
        allocate( f1_east_snd(Nc*Nytotal*ghostLayers) ) 
        allocate( f1_east_rcv(Nc*Nytotal*ghostLayers) ) 

        allocate( f1_north_snd(Nc*Nxtotal*ghostLayers) ) 
        allocate( f1_north_rcv(Nc*Nxtotal*ghostLayers) )
        allocate( f1_south_snd(Nc*Nxtotal*ghostLayers) ) 
        allocate( f1_south_rcv(Nc*Nxtotal*ghostLayers) )

        ! Allways initialize allocated arrays
        f1_west_rcv = 0.d0
        f1_east_rcv = 0.d0
        f1_north_rcv = 0.d0
        f1_south_rcv = 0.d0
        f1_west_snd = 0.d0
        f1_east_snd = 0.d0
        f1_north_snd = 0.d0
        f1_south_snd = 0.d0

    end subroutine setupFlow
end module flow
