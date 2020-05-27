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
! module    : solver
!-------------------------------------------------------------------------------
! This is a module for each iteration solving numerically the govering Equation (4) of Ref. [1]. 
! For details:
!
! [1]   M.T. Ho, L. Zhu, L. Wu, P. Wang, Z. Guo, Z.-H. Li, Y. Zhang
!       "A multi-level parallel solver for rarefied gas flows in porous media"
! 		Computer Physics Communications, 234 (2019), pp. 14-25
!
!	For fluid nodes, there are 8 groups of molecular velocity and hence 8 paths of sweep. 
!   Boundary conditions on the solid surface are calculated.
! 	Macroscopic parameters for fluid nodes such as density number, velocity are calculated.
!	OpenMP for loop accelaration and  buffer pack/unpack are implemented.
!   See Section 2.3, 3.1, 3.2, Algorithm 1, Figure 1 of Ref.[1]
!
!
!!  f(spatial_id,velocity_id,sweeppath_id)  : velocity distribution function in Eq.(5) of [1], main unknow to be solved
!!  fEq: equilibrium distribution function in Eq.(3) of [1]
!!  cx, cy: discrete velocity components in x, y directions 
!!  dir1,dir2,..,dir4 : 4 paths of sweep
!!  coef1,coef2,..,coef4: arrays storing weighting of fluid nodes in upwind schemes
!!  RhoWall:  number density on the solid surface in Eq.(9) of Ref.[1]
!
!
!       ---------------------------------
!				symmetric BC
!       ---------------------------------
!       ||i                           o||                   y
!       ||n                           u||                   /\
!       ||-                           t||                   |
!       ||l                           l||                   |
!       ||e                           e||                   |
!       ||t                           t||                   |
!       ---------------------------------                   |----------->x
!				symmetric BC
!       ---------------------------------
!
!
!
!                       II      |       I
!                               |
!                               |
!               ----------------------------------  molecular velocity group
!                               |
!                               |
!                       III     |       IV
!
! Periodic & pressure drop on inlet&outlet
! Symmetry on lateral wall
!-------------------------------------------------------------------------------

module solver
use flow
use velocityGrid
use physicalGrid
use mpiParams

implicit none

! to be read from NML: solverNml
double precision :: eps
integer :: maxStep
integer :: chkConvergeStep
integer :: saveStep
integer :: saveFormat ! 1 (default) for vti, 2 for tecplot, 3 for vtk

integer :: iStep
double precision :: error
double precision :: permeability

contains
    subroutine iterate
        implicit none
        include "mpif.h"
        integer :: k, l, i, j, shiftll, shiftuu
        integer :: MPI_ERR
        integer :: MPI_REQ_X(4), MPI_REQ_Y(4)
        integer :: MPI_STAT(MPI_status_SIZE,4)
        integer :: xsize, ysize
        double precision :: feq, RhoWall

        xsize = Nytotal*Nc*ghostLayers
        ysize = Nxtotal*Nc*ghostLayers

        MPI_REQ_X = MPI_REQUEST_NULL
        MPI_REQ_Y = MPI_REQUEST_NULL

!$OMP PARALLEL &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(l, i, j, k, fEq, RhoWall)

        ! Start Recieving
!$OMP SINGLE
        call MPI_IRECV( f1_east_rcv, xsize, MPI_DOUBLE_PRECISION, east,  tag1, &
                        MPI_COMM_VGRID, MPI_REQ_X(1), MPI_ERR )
        call MPI_IRECV( f1_west_rcv, xsize, MPI_DOUBLE_PRECISION, west,  tag2, &
                        MPI_COMM_VGRID, MPI_REQ_X(2), MPI_ERR )
        call MPI_IRECV( f1_north_rcv, ysize, MPI_DOUBLE_PRECISION, north,  tag3, &
                        MPI_COMM_VGRID, MPI_REQ_Y(1), MPI_ERR )
        call MPI_IRECV( f1_south_rcv, ysize, MPI_DOUBLE_PRECISION, south,  tag4, &
                        MPI_COMM_VGRID, MPI_REQ_Y(2), MPI_ERR )     

        ! Start Sending
        call MPI_ISEND( f1_west_snd, xsize, MPI_DOUBLE_PRECISION, west, tag1, &
                        MPI_COMM_VGRID, MPI_REQ_X(3), MPI_ERR )
        call MPI_ISEND( f1_east_snd, xsize, MPI_DOUBLE_PRECISION, east, tag2, &
                        MPI_COMM_VGRID, MPI_REQ_X(4), MPI_ERR )
        call MPI_ISEND( f1_south_snd, ysize, MPI_DOUBLE_PRECISION, south, tag3, &
                        MPI_COMM_VGRID, MPI_REQ_Y(3), MPI_ERR )
        call MPI_ISEND( f1_north_snd, ysize, MPI_DOUBLE_PRECISION, north, tag4, &
                        MPI_COMM_VGRID, MPI_REQ_Y(4), MPI_ERR )       
!$OMP END SINGLE NOWAIT

!$OMP do SCHEDULE(RUNTIME)
        ! sweep direction 1
        do l=1,Nc/4
            do i=1,Nstencil1
                k=dir1(i)
                ! Switch from reflected in x-direction (default) to in y-direction
                ! only for group of velocity overlapped by two reflected direction
                if (image(k-Nxtotal)==WallXpYp) then
                    f1(k-Nxtotal,l)=f1(k-Nxtotal,l+Nc/2)
                endif
    
                fEq=w(l)*(Rho(k)+2.d0*(cx(l)*Ux(k)+cy(l)*Uy(k)))
    
                f1(k,l)=(mu*(fEq-0.5d0*f1(k,l)) &           
                &        + cx(l)*coefI(i,2)*f1(k-1,l) &
                &        + cx(l)*coefI(i,3)*f1(k-2,l) &
                &        + cy(l)*coefI(i,5)*f1(k-Nxtotal,l) &
                &        + cy(l)*coefI(i,6)*f1(k-2*Nxtotal,l) &
                & )/(0.5d0*mu+cx(l)*coefI(i,1)+cy(l)*coefI(i,4))   
            enddo
        enddo
!$OMP END do NOWAIT

!$OMP do SCHEDULE(RUNTIME)
        ! sweep direction 2
        do l=Nc/4+1,Nc/2
            do i=1,Nstencil2
                k=dir2(i)
                ! Switch from reflected in x-direction (default) to in y-direction
                ! only for group of velocity overlapped by two reflected direction
                if (image(k-Nxtotal)==WallXnYp) then
                    f1(k-Nxtotal,l)=f1(k-Nxtotal,l+Nc/2)
                endif
    
                fEq=w(l)*(Rho(k)+2.d0*(cx(l)*Ux(k)+cy(l)*Uy(k)))
                f1(k,l)=(mu*(fEq-0.5d0*f1(k,l)) &
                ! f1(k,l)=(mu*(fEq) &    
                ! f1(k,l)=(mu*(fEq-f1(k,l)) &        
                &        + cx(l)*coefII(i,2)*f1(k+1,l) &
                &        + cx(l)*coefII(i,3)*f1(k+2,l) &
                &        + cy(l)*coefII(i,5)*f1(k-Nxtotal,l) &
                &        + cy(l)*coefII(i,6)*f1(k-2*Nxtotal,l) &
                & )/(0.5d0*mu+cx(l)*coefII(i,1)+cy(l)*coefII(i,4))
                ! & )/(mu+cx(l)*coefII(i,1)+cy(l)*coefII(i,4))   
                ! & )/(cx(l)*coefII(i,1)+cy(l)*coefII(i,4))          
            enddo
        enddo
!$OMP END do NOWAIT

!$OMP do SCHEDULE(RUNTIME)
        ! sweep direction 3
        do l=Nc/2+1,Nc*3/4
            do i=1,Nstencil3
                k=dir3(i)
                ! Switch from reflected in x-direction (default) to in y-direction
                ! only for group of velocity overlapped by two reflected direction
                if (image(k+Nxtotal)==WallXnYn) then
                    f1(k+Nxtotal,l)=f1(k+Nxtotal,l-Nc/2)
                endif

                fEq=w(l)*(Rho(k)+2.d0*(cx(l)*Ux(k)+cy(l)*Uy(k)))
                f1(k,l)=(mu*(fEq-0.5d0*f1(k,l)) &
                ! f1(k,l)=(mu*(fEq) &    
                ! f1(k,l)=(mu*(fEq-f1(k,l)) &        
                &        + cx(l)*coefIII(i,2)*f1(k+1,l) &
                &        + cx(l)*coefIII(i,3)*f1(k+2,l) &
                &        + cy(l)*coefIII(i,5)*f1(k+Nxtotal,l) &
                &        + cy(l)*coefIII(i,6)*f1(k+2*Nxtotal,l) &
                & )/(0.5d0*mu+cx(l)*coefIII(i,1)+cy(l)*coefIII(i,4))
                ! & )/(mu+cx(l)*coefIII(i,1)+cy(l)*coefIII(i,4)) 
                ! & )/(cx(l)*coefIII(i,1)+cy(l)*coefIII(i,4))            
            enddo
        enddo
!$OMP END do NOWAIT

!$OMP do SCHEDULE(RUNTIME)
        ! sweep direction 4
        do l=Nc*3/4+1,Nc
            do i=1,Nstencil4
                k=dir4(i)
                ! Switch from reflected in x-direction (default) to in y-direction
                ! only for group of velocity overlapped by two reflected direction
                if (image(k+Nxtotal)==WallXpYn) then
                    f1(k+Nxtotal,l)=f1(k+Nxtotal,l-Nc/2)
                endif
    
                fEq=w(l)*(Rho(k)+2.d0*(cx(l)*Ux(k)+cy(l)*Uy(k)))
                f1(k,l)=(mu*(fEq-0.5d0*f1(k,l)) &
                ! f1(k,l)=(mu*(fEq) &    
                ! f1(k,l)=(mu*(fEq-f1(k,l)) &        
                &        + cx(l)*coefIV(i,2)*f1(k-1,l) &
                &        + cx(l)*coefIV(i,3)*f1(k-2,l) &
                &        + cy(l)*coefIV(i,5)*f1(k+Nxtotal,l) &
                &        + cy(l)*coefIV(i,6)*f1(k+2*Nxtotal,l) &
                & )/(0.5d0*mu+cx(l)*coefIV(i,1)+cy(l)*coefIV(i,4))
                !& )/(mu+cx(l)*coefIV(i,1)+cy(l)*coefIV(i,4))   
                !& )/(cx(l)*coefIV(i,1)+cy(l)*coefIV(i,4))  
            enddo
        enddo
!$OMP END do 

        ! Wait until send and recv done
!$OMP SINGLE
        call MPI_WAITALL(4, MPI_REQ_X, MPI_STAT, MPI_ERR)
        call MPI_WAITALL(4, MPI_REQ_Y, MPI_STAT, MPI_ERR)
!$OMP END SINGLE 

!-------------------------------------------------------------------
!> Pack/unpack boundary data
!-------------------------------------------------------------------
!$OMP DO
        do j = 1, Nytotal
            do i = 1, ghostLayers
                do l = 1, Nc ! dir 1, 2, 3 and 4
                    f1_west_snd((j-1)*ghostLayers*Nc + (i-1)*Nc + l) &
                    &  = f1((j-1)*Nxtotal + i+ghostLayers, l)
                    f1((j-1)*Nxtotal + i+Nxsub+ghostLayers, l) = &
                    &   f1_east_rcv((j-1)*ghostLayers*Nc + (i-1)*Nc + l)

                    f1_east_snd((j-1)*ghostLayers*Nc + (i-1)*Nc + l) &
                    &  = f1((j-1)*Nxtotal + i+Nxsub, l)
                    f1((j-1)*Nxtotal + i, l) = &
                    &   f1_west_rcv((j-1)*ghostLayers*Nc + (i-1)*Nc + l)                   
                enddo
            enddo
        enddo
!$OMP END do NOWAIT

!$OMP DO
        ! pack&unpack south&north buffer
        do j = 1, ghostLayers
            do i = 1, Nxtotal
                do l = 1, Nc ! dir 1, 2, 3 and 4
                    f1_north_snd((j-1)*Nxtotal*Nc + (i-1)*Nc + l) &
                    &  = f1((Nysub+j-1)*Nxtotal + i, l)
                    f1((j-1)*Nxtotal + i, l) = &
                    &   f1_south_rcv((j-1)*Nxtotal*Nc + (i-1)*Nc + l)
                    f1_south_snd((j-1)*Nxtotal*Nc + (i-1)*Nc + l) &
                    &  = f1((j-1+ghostLayers)*Nxtotal + i, l)
                    f1((Nysub+ghostLayers+j-1)*Nxtotal + i,l) = &
                    &   f1_north_rcv((j-1)*Nxtotal*Nc + (i-1)*Nc + l)
                enddo
            enddo
        enddo
!$OMP enddo


!-------------------------------------------------------------------
!> Processing wall nodes
!-------------------------------------------------------------------
!$OMP do SCHEDULE(RUNTIME)
        do i=1,nWall
            k=vecWall(i)
            RhoWall=0.d0
            select case (image(k))
                case (WallXp)
                    do l=Nc/4+1,3*Nc/4
                        !f1(k,l)=2.d0*f1(k+1,l)-f1(k+2,l)
                        !f1(k,l)=1.d0*f1(k+1,l)
                        f1(k,l)=extCoef(i,1)*f1(k+1,l) + extCoef(i,2)*f1(k+2,l)
                        RhoWall=RhoWall-cx(l)*f1(k,l)
                    enddo
                    RhoWall=RhoWall/DiffFlux
                    do l=1,Nc/4
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeX(l))
                    enddo
                    do l=3*Nc/4+1,Nc
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeX(l))
                    enddo
                case (WallXn)
                    do l=1,Nc/4
                        !f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                        !f1(k,l)=1.d0*f1(k-1,l)
                        f1(k,l)=extCoef(i,1)*f1(k-1,l) + extCoef(i,2)*f1(k-2,l)
                        RhoWall=RhoWall+cx(l)*f1(k,l)
                    enddo
                    do l=3*Nc/4+1,Nc
                        !f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                        !f1(k,l)=1.d0*f1(k-1,l)
                        f1(k,l)=extCoef(i,1)*f1(k-1,l) + extCoef(i,2)*f1(k-2,l)
                        RhoWall=RhoWall+cx(l)*f1(k,l)
                    enddo
                    RhoWall=RhoWall/DiffFlux
                    do l=Nc/4+1,3*Nc/4
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeX(l))
                    enddo
                case (WallYp)
                    do l=Nc/2+1,Nc
                        !f1(k,l)=2.d0*f1(k+Nxtotal,l)-f1(k+2*Nxtotal,l)
                        !f1(k,l)=1.d0*f1(k+Nxtotal,l)
                        f1(k,l)=extCoef(i,1)*f1(k+Nxtotal,l) + extCoef(i,2)*f1(k+2*Nxtotal,l)
                        RhoWall=RhoWall-cy(l)*f1(k,l)
                    enddo
                    RhoWall=RhoWall/DiffFlux
                    do l=1,Nc/2
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeY(l))
                    enddo
                case (WallYn)
                    do l=1,Nc/2
                        !f1(k,l)=2.d0*f1(k-Nxtotal,l)-f1(k-2*Nxtotal,l)
                        !f1(k,l)=1.d0*f1(k-Nxtotal,l)
                        f1(k,l)=extCoef(i,1)*f1(k-Nxtotal,l) + extCoef(i,2)*f1(k-2*Nxtotal,l)
                        RhoWall=RhoWall+cy(l)*f1(k,l)
                    enddo
                    RhoWall=RhoWall/DiffFlux
                    do l=Nc/2+1,Nc
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeY(l))
                    enddo
                !=======================================================================
                !     Boundary condition on the corner wall
                !=======================================================================
                case (WallXpYp)
                    !Calculate default reflected f1 (x-direction)
                    do l=Nc/4+1,3*Nc/4
                        !f1(k,l)=2.d0*f1(k+1,l)-f1(k+2,l)
                        !f1(k,l)=1.d0*f1(k+1,l)
                        f1(k,l)=extCoef(i,1)*f1(k+1,l) + extCoef(i,2)*f1(k+2,l)
                        RhoWall=RhoWall-cx(l)*f1(k,l)
                    enddo
                    RhoWall=RhoWall/DiffFlux
                    do l=1,Nc/4
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeX(l))
                    enddo
                    do l=3*Nc/4+1,Nc
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeX(l))
                    enddo
                    !Calculate and store RhoWallY
                    RhoWall=0.d0
                    do l=Nc/2+1,Nc
                        !RhoWall=RhoWall-cy(l)*(2.d0*f1(k+Nxtotal,l)-f1(k+2*Nxtotal,l))
                        !RhoWall=RhoWall-cy(l)*(1.d0*f1(k+Nxtotal,l))
                        RhoWall=RhoWall-cy(l)*(extCoef(i,1)*f1(k+Nxtotal,l) + extCoef(i,2)*f1(k+2*Nxtotal,l)) !STOP
                    enddo
                    RhoWall=RhoWall/DiffFlux
                    ! RhoWallY(k)=RhoWall
                    !Calculate default reflected f1 (y-direction)
                    do l=Nc/4+1,Nc/2
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*(2.d0*f1(k+Nxtotal,oppositeY(l))-f1(k+2*Nxtotal,oppositeY(l)))
                    enddo
                    !Calculate and store the secondary reflected f1 (y-direction)               
                    do l=1,Nc/4
                        f1(k,l+Nc/2)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*(2.d0*f1(k+Nxtotal,oppositeY(l))-f1(k+2*Nxtotal,oppositeY(l)))
                    enddo               
                case (WallXnYp)
                    !Calculate default reflected f1 (x-direction)
                    do l=1,Nc/4
                        !f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                        !f1(k,l)=1.d0*f1(k-1,l)
                        f1(k,l)=extCoef(i,1)*f1(k-1,l) + extCoef(i,2)*f1(k-2,l)
                        RhoWall=RhoWall+cx(l)*f1(k,l)
                    enddo
                    do l=3*Nc/4+1,Nc
                        !f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                        !f1(k,l)=1.d0*f1(k-1,l)
                        f1(k,l)=extCoef(i,1)*f1(k-1,l) + extCoef(i,2)*f1(k-2,l)
                        RhoWall=RhoWall+cx(l)*f1(k,l)
                    enddo
                    RhoWall=RhoWall/DiffFlux
                    do l=Nc/4+1,3*Nc/4
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeX(l))
                    enddo
                    !Calculate and store RhoWallY
                    RhoWall=0.d0
                    do l=Nc/2+1,Nc
                        !RhoWall=RhoWall-cy(l)*(2.d0*f1(k+Nxtotal,l)-f1(k+2*Nxtotal,l))
                        !RhoWall=RhoWall-cy(l)*(1.d0*f1(k+Nxtotal,l))
                        RhoWall=RhoWall-cy(l)*(extCoef(i,1)*f1(k+Nxtotal,l) + extCoef(i,2)*f1(k+2*Nxtotal,l))
                    enddo
                    RhoWall=RhoWall/DiffFlux
                    !RhoWallY(k)=RhoWall
                    !Calculate default reflected f1 (y-direction)
                    do l=1,Nc/4
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*(2.d0*f1(k+Nxtotal,oppositeY(l))-f1(k+2*Nxtotal,oppositeY(l)))
                    enddo
                    !Calculate and store the secondary reflected f1 (y-direction)               
                    do l=Nc/4+1,Nc/2
                        f1(k,l+Nc/2)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*(2.d0*f1(k+Nxtotal,oppositeY(l))-f1(k+2*Nxtotal,oppositeY(l)))
                    enddo               
                case (WallXnYn)
                    !Calculate default reflected f1 (x-direction)
                    do l=1,Nc/4
                        !f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                        !f1(k,l)=1.d0*f1(k-1,l)
                        f1(k,l)=extCoef(i,1)*f1(k-1,l) + extCoef(i,2)*f1(k-2,l)
                        RhoWall=RhoWall+cx(l)*f1(k,l)
                    enddo
                    do l=3*Nc/4+1,Nc
                        !f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                        !f1(k,l)=1.d0*f1(k-1,l)
                        f1(k,l)=extCoef(i,1)*f1(k-1,l) + extCoef(i,2)*f1(k-2,l)
                        RhoWall=RhoWall+cx(l)*f1(k,l)
                    enddo
                    RhoWall=RhoWall/DiffFlux
                    do l=Nc/4+1,3*Nc/4
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeX(l))
                    enddo
                    !Calculate and store RhoWallY
                    RhoWall=0.d0
                    do l=1,Nc/2
                        !RhoWall=RhoWall+cy(l)*(2.d0*f1(k-Nxtotal,l)-f1(k-2*Nxtotal,l))
                        !RhoWall=RhoWall+cy(l)*(1.d0*f1(k-Nxtotal,l))
                        RhoWall=RhoWall+cy(l)*(extCoef(i,1)*f1(k-Nxtotal,l) + extCoef(i,2)*f1(k-2*Nxtotal,l))
                    enddo
                    RhoWall=RhoWall/DiffFlux
                    ! RhoWallY(k)=RhoWall
                    !Calculate default reflected f1 (y-direction)
                    do l=3*Nc/4+1,Nc
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*(2.d0*f1(k-Nxtotal,oppositeY(l))-f1(k-2*Nxtotal,oppositeY(l)))
                    enddo
                    !Calculate and store the secondary reflected f1 (y-direction)               
                    do l=Nc/2+1,3*Nc/4
                        f1(k,l-Nc/2)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*(2.d0*f1(k-Nxtotal,oppositeY(l))-f1(k-2*Nxtotal,oppositeY(l)))
                    enddo               
                case (WallXpYn)
                    !Calculate default reflected f1 (x-direction)
                    do l=Nc/4+1,3*Nc/4
                        !f1(k,l)=2.d0*f1(k+1,l)-f1(k+2,l)
                        !f1(k,l)=1.d0*f1(k+1,l)
                        f1(k,l)=extCoef(i,1)*f1(k+1,l) + extCoef(i,2)*f1(k+2,l)
                        RhoWall=RhoWall-cx(l)*f1(k,l)
                    enddo
                    RhoWall=RhoWall/DiffFlux
                    do l=1,Nc/4
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeX(l))
                    enddo
                    do l=3*Nc/4+1,Nc
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeX(l))
                    enddo
                    !Calculate and store RhoWallY
                    RhoWall=0.d0
                    do l=1,Nc/2
                        !RhoWall=RhoWall+cy(l)*(2.d0*f1(k-Nxtotal,l)-f1(k-2*Nxtotal,l))
                        !RhoWall=RhoWall+cy(l)*(1.d0*f1(k-Nxtotal,l))
                        RhoWall=RhoWall+cy(l)*(extCoef(i,1)*f1(k-Nxtotal,l) + extCoef(i,2)*f1(k-2*Nxtotal,l))                        
                    enddo
                    RhoWall=RhoWall/DiffFlux
                    ! RhoWallY(k)=RhoWall
                    !Calculate default reflected f1 (y-direction)
                    do l=Nc/2+1,3*Nc/4
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*(2.d0*f1(k-Nxtotal,oppositeY(l))-f1(k-2*Nxtotal,oppositeY(l)))
                    enddo
                    !Calculate and store the secondary reflected f1 (y-direction)               
                    do l=3*Nc/4+1,Nc
                        f1(k,l-Nc/2)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*(2.d0*f1(k-Nxtotal,oppositeY(l))-f1(k-2*Nxtotal,oppositeY(l)))
                    enddo                               
            END SELECT
        enddo
!$OMP END do  NOWAIT


!--------------------------------------------------------
!> inlet/outlet
!--------------------------------------------------------
        if(xl==xmin) then ! inlet block (west most processor)
!$OMP DO
            do j = ylg, yug
                i = xl
                k = (j-ylg)*Nxtotal + i-xlg+1
                !inlet
                do l=1,Nc/4
                    !f1(k,l)=f1(k-1+Nxsub,l)+w(l)*PressDrop !lhzhu, need other block's info, so vgrid is perodical in x dir
                    f1(k,l)=f1(k-1,l)+w(l)*PressDrop !lhzhu, need other block's info, so vgrid is perodical in x dir
                enddo   
                do l=3*Nc/4+1,Nc
                    !f1(k,l)=f1(k-1+Nxsub,l)+w(l)*PressDrop ! NOTE, for NprocX=1
                    f1(k,l)=f1(k-1,l)+w(l)*PressDrop ! NOTE, for NprocX=1
                enddo
            enddo
!$OMP END do NOWAIT
        endif

        if(xu==xmax) then ! outlet block (east most processor)
!$OMP DO
            do j = ylg, yug
                i = xu
                k = (j-ylg)*Nxtotal + i-xlg+1
                !outlet
                do l=Nc/4+1,3*Nc/4          
                    !f1(k,l)=f1(k-Nxsub+1,l)-w(l)*PressDrop ! NOTE, for NprocX=1
                    f1(k,l)=f1(k+1,l)-w(l)*PressDrop ! NOTE, for NprocX=1
                enddo
            enddo
!$OMP enddo
        endif

!----------------------------------------------------
!> Symmetric BC
!----------------------------------------------------
        if(yl==ymin) then ! south
!$OMP DO
            do i=xl, xu
                k = ghostLayers*Nxtotal + i-xlg+1
                do l=1,Nc/2
                     f1(k,l)=f1(k,oppositeY(l))
                enddo
            enddo
!$OMP END do NOWAIT 
        endif
        if(yu==ymax) then ! north
!$OMP DO
            do i=xl, xu
                k = (Nysub+ghostLayers-1)*Nxtotal + i-xlg+1
                do l=Nc/2+1,Nc
                     f1(k,l)=f1(k,oppositeY(l))
                enddo
            enddo
!$OMP enddo
        endif

        !----------------------------------------------------
        !> Update Macro
        !----------------------------------------------------
!$OMP DO
        do i=1,Nfluid
            k = mapF(i)
            Rho(k)=0.d0
            Ux(k)=0.d0
            Uy(k)=0.d0
            do l=1,Nc
                Rho(k)=Rho(k)+f1(k,l)
                Ux(k)=Ux(k)+cx(l)*f1(k,l)
                Uy(k)=Uy(k)+cy(l)*f1(k,l)
            enddo
        enddo
!$OMP enddo

!$OMP END PARALLEL

    end subroutine iterate
   
    subroutine chkConverge
        implicit none
        include "mpif.h"

        ! local vars
        integer :: byl, byu, k, j, MPI_ERR
        double precision :: massInner, massNorth, massSouth
        double precision :: massLocal, mass2

        massInner = 0.d0
        massNorth = 0.d0
        massSouth = 0.d0
        
        byl = yl
        byu = yu
        if(yl == ymin) byl = yl + 1 !if most south block
        if(yu == ymax) byu = yu - 1 !if most north block

        !mass2=0.d0
        if(xl == xmin) then !only left most processors
            do j=byl, byu
                k=(j-ylg)*Nxtotal + column+ghostLayers
                massInner=massInner+Ux(k)*ds
            enddo
            !print*, "massInner=",  massInner
            if (yl == ymin) then !only south most processors
                massSouth = 0.5d0*ds*Ux(ghostLayers+column + (yl-ylg)*Nxtotal)
            endif
            if (yu == ymax) then !only north most processors
                massNorth = 0.5d0*ds*Ux(ghostLayers+column + (ghostLayers+Nysub-1)*Nxtotal) !bug here
            endif
            ! debug
            massLocal = (massInner + massSouth + massNorth) * 2.d0 / PressDrop

            ! reduction
            call MPI_ALLREDUCE(massLocal, mass2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                               mpi_comm_inlet, MPI_ERR)

            error=dabs(1.d0-mass2/mass)/(chkConvergeStep)
            mass=mass2
            if (proc == master) then           
                permeability=mass*Kn*dsqrt(4.d0/PI)*((Nx-1)/(1.d0+xPadRatio))/(Ny-1)/2.d0
                !permeability=mass*Kn*dsqrt(4.d0/PI)*(3000.0/1500.0)/2.d0
                write(*,"( 1I10, 3ES15.6)")  iStep,  mass,  permeability, error
                open(22,file='Results.dat', position="append")
                write(22,'(4ES15.6, 1I15)') Kn, mass, permeability, error, iStep
                close(22)
            endif
        endif
        !bcast error so every process in WORLD can stop
        call MPI_BCAST(error, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_VGRID, MPI_ERR)
    end subroutine chkConverge

    subroutine saveFlowField(saveFormat)
        integer, intent(in) :: saveFormat
        ! Save final data
        select case (saveFormat)
            case (1)
                call saveFlowFieldVTI
            case (2)
                call saveFlowFieldTec
            case (3)
                call saveFlowFieldVTK
        end select
    end subroutine saveFlowField

    subroutine saveFlowFieldTec
        integer :: j, i, k
        character(21) fname
        write(fname, '(A, I0.3, A)') 'Field_', proc, '.dat'
        !write(fname, '(A, I0.3, A, I0.5, A)') 'Field_', proc, '_T_', iStep, '.dat'

        open(20,file=dataSaveDir//"/"//fname, status='REPLACE')
        write(20,*) ' TITLE=" Field"'
        write(20,*) 'StrandID='//itoa(iStep)//',SolutionTime='//itoa(iStep)
        write(20,*) ' VARIABLES=x,y,flag,Rho,Ux,Uy'
        write(20,'(A,I0.3,A,I4,A,I4,A)') ' ZONE T=proc', proc, ', I=', Nxsub,', J=', Nysub,', F=POINT'

        do j=yl,yu
            do i=xl,xu
                k= (j-ylg)*Nxtotal + i-xlg+1
                if (image(k)==fluid) then
                    write(20,'(3I10,3ES15.6)') i, j, 0, Rho(k)+1.d0, Ux(k), Uy(k)
                else
                    write(20,'(3I10,3ES15.6)') i, j, 1, 0.d0, 0.d0, 0.d0
                endif   
            enddo
        enddo
        close(20)
    end subroutine saveFlowFieldTec

    subroutine saveFlowFieldVTK
        implicit none
        integer :: i, j, k, l, MPI_ERR, IO_ERR
        character(13) fname
        integer :: zl = 1
        integer :: zu = 1
        integer :: Nzsub = 1

        write(fname, '(A, I0.3, A)') 'Field_', proc, '.vtk'
        open(unit = 12, file = dataSaveDir//"/"//fname, status = "REPLACE", position = "APPEND", &
          iostat = IO_ERR)
        if ( IO_ERR == 0 ) then
             write(12,'(A)')"# vtk DataFile Version 2.0"
             write(12,'(A)')"DVM MPI"
             write(12,'(A)')"ASCII"
             write(12,*)
             write(12,'(A)')"DATASET STRUCTURED_POINTS"
             write(12,*)"DIMENSIONS",Nxsub,Nysub,Nzsub
             write(12,*)"ORIGIN",xl,yl,zl
             write(12,*)"SPACING",1,1,1
             write(12,*)
             write(12,*)"POINT_DATA",Nxsub*Nysub*Nzsub
             write(12,*)
             write(12,'(A)')"SCALARS Rho double"
             write(12,'(A)')"LOOKUP_TABLE default"
             do j = yl, yu
                 do i = xl, xu
                     l= (j-ylg)*Nxtotal + i-xlg+1
                     if (image(l)==fluid) then
                         write(12,'(ES15.6)') Rho(l)+1.d0
                     else
                         write(12,'(ES15.6)') 0.d0
                     endif   
                 enddo
             enddo
             write(12,*)
             write(12,'(A)')"VECTORS Velocity double"
             do j = yl, yu
                 do i = xl, xu
                     l= (j-ylg)*Nxtotal + i-xlg+1
                     if (image(l)==fluid) then
                         write(12,'(3ES15.6)') Ux(l), Uy(l), 0.d0
                     else
                         write(12,'(3ES15.6)') 0.d0, 0.d0, 0.d0
                     endif  
                 enddo
            enddo
            close(unit = 12)
        else
            call memFree
            call MPI_FINALIZE(MPI_ERR)
            stop "Error: Unable to open output vtk file."
        endif

        return
    end subroutine saveFlowFieldVTK

    subroutine saveFlowFieldVTI
        integer :: i, j, k, l, MPI_ERR, IO_ERR
        character(13) fname
        character(10) pfname
        integer :: exl, exu, eyl, eyu, ezl, ezu
        integer :: wexl, wexu, weyl, weyu, wezl, wezu
        integer :: zl, zu, zmin, zmax

        zl = 1
        zu = 1
        zmin = 1
        zmax = 1
        ! whole extend
        wexl = xmin - 1
        wexu = xmax
        weyl = ymin - 1
        weyu = ymax
        wezl = zmin - 1
        wezu = zmax
        ! local extend
        exl = xl - 1
        exu = xu
        eyl = yl - 1
        eyu = yu
        ezl = zl - 1
        ezu = zu

        write(fname, '(A, I0.3, A)') 'Field_', proc, '.vti'
        open(unit = 13, file = dataSaveDir//"/"//fname, status = "REPLACE", position = "APPEND", &
          iostat = IO_ERR)
        IF ( IO_ERR == 0 ) then
            write(13,'(A)') '<?xml version="1.0"?>'
            write(13,'(A)') '<VTKFile type="ImageData">'
            write(13,'(A, 6I8, A)') '<ImageData WholeExtent="', exl, exu, & 
                eyl, eyu, ezl, ezu, ' " Origin="0 0 0" Spacing="1 1 1">'
            write(13, '(A, 6I8, A)') '<Piece Extent="', exl, exu, eyl, eyu, &
                ezl, ezu, '">'
            write(13, '(A)') '<CellData Scalars="flag Rho" Vectors ="U">'
            write(13, '(A)') '<DataArray type="Int32" Name="flag" format="ascii">'

            do j = yl, yu
                do i = xl, xu
                    l= (j-ylg)*Nxtotal + i-xlg+1
                    if (image(l)==fluid) then
                        write(13, '(1I4)') 0 ! fluid flag = 0
                    else
                        write(13, '(1I4)') 1
                    endif
                enddo
            enddo

            write(13, '(A)') '</DataArray>'
            write(13, '(A)') '<DataArray type="Float32" Name="Rho" format="ascii">'
            do j = yl, yu
                do i = xl, xu
                    l= (j-ylg)*Nxtotal + i-xlg+1
                    if (image(l)==fluid) then
                        write(13,'(ES15.6)') Rho(l)+1.d0
                    else
                        write(13,'(ES15.6)') 0.d0
                    endif   
                enddo
            enddo
            write(13, '(A)') '</DataArray>'

            write(13, '(A)') '<DataArray type="Float32" Name="U" format="ascii" &
                & NumberOfComponents="3">'
            do j = yl, yu
               do i = xl, xu
                  l= (j-ylg)*Nxtotal + i-xlg+1
                  if (image(l)==fluid) then
                      write(13,'(3ES15.6)') Ux(l), Uy(l), 0.d0
                  else
                      write(13,'(3ES15.6)') 0.d0, 0.d0, 0.d0
                  endif  
               enddo
            enddo
            write(13, '(A)') '</DataArray>'
            write(13, '(A)') '</CellData>'
            write(13, '(A)') '</Piece>'
            write(13, '(A)') '</ImageData>'
            write(13, '(A)') '</VTKFile>'

            close(unit = 13)
        else
            call memFree
            call MPI_FINALIZE(MPI_ERR)
            stop "Error: Unable to open output vti file."
        endif


        if (proc == master) then  
            write(pfname, '(A)') 'Field.pvti'
            open(unit = 14, file = dataSaveDir//"/"//pfname, status = "REPLACE", position = "APPEND", &
              iostat = IO_ERR)
            write(14,'(A)') '<?xml version="1.0"?>'
            write(14,'(A)') '<VTKFile type="PImageData">'
            write(14,'(A, 6I8, A)') '<PImageData WholeExtent="', wexl, wexu, & 
                weyl, weyu, wezl, wezu, ' " Origin="0 0 0" Spacing="1 1 1">'
            write(14,'(A)') '<PCellData Scalars="flag Rho" Vectors ="U">'
            write(14, '(A)') '<DataArray type="Int32" Name="flag" format="ascii"/>'
            write(14, '(A)') '<DataArray type="Float32" Name="Rho" format="ascii"/>'
            write(14, '(A)') '<DataArray type="Float32" Name="U" format="ascii" &
                & NumberOfComponents="3"/>'
            write(14, '(A)') '</PCellData>'
            ! the master process has recorded the extend of the subdomain in sub_ext
            do l=1, nprocs
                exl = sub_ext(1,l) - 1
                exu = sub_ext(2,l)
                eyl = sub_ext(3,l) - 1
                eyu = sub_ext(4,l)
                ezl = 0
                ezu = 1
                write(fname, '(A, I0.3, A)') 'Field_', l-1, '.vti'
                write(14, '(A, 6I8, A)') '<Piece Extent="', exl, exu, eyl, eyu, &
                    ezl, ezu, '" Source="'//fname//'"/>'
            enddo
            write(14, '(A)') '</PImageData>'
            write(14, '(A)') '</VTKFile>'
            close(unit=14)
        endif

    end subroutine saveFlowFieldVTI


    subroutine saveNodeCounts
        implicit none
        include "mpif.h"
        ! array holds fluid node counts and total node counts in each subdomain
        integer, dimension(:), allocatable :: fluidNodeCountAll, totalNodeCountAll
        integer :: fluidNodeCount, totalNodeCount
        integer :: MPI_ERR, IO_ERR, i, j, localid

        allocate(fluidNodeCountAll(nprocs))
        allocate(totalNodeCountAll(nprocs))

        totalNodeCount = Nxsub*Nysub
        fluidNodeCount = 0

        do j=yl,yu
            do i=xl,xu
                localid = (j-ylg)*Nxtotal + i-xlg+1
                if(image(localid) == fluid) then
                    fluidNodeCount = fluidNodeCount + 1
                endif
            enddo
        enddo

        ! do MPI gather put fluidNodeCount to fluidNodeCountAll
        call MPI_GATHER(totalNodeCount, 1, MPI_INTEGER, totalNodeCountAll, 1, &
            MPI_INTEGER, master, MPI_COMM_VGRID, MPI_ERR)
        call MPI_GATHER(fluidNodeCount, 1, MPI_INTEGER, fluidNodeCountAll, 1, &
            MPI_INTEGER, master, MPI_COMM_VGRID, MPI_ERR)
        ! master process write to file
        if(proc == master) then
            open(unit=15, file="nodeCounts", status="replace", iostat=IO_ERR)
            write(15,'(A)') "totalNodeCount fluidNodeCount localPorosity"
            do i=1,nprocs
                write(15,'(2I12, ES15.6)') &
                    totalNodeCountAll(i), fluidNodeCountAll(i), &
                    dble(fluidNodeCountAll(i))/dble(totalNodeCountAll(i))
            enddo
            close(unit=15)
        endif
    end subroutine saveNodeCounts

    subroutine memFree
        deallocate (array2D)
        deallocate (image)
        deallocate (vecWall, dir1, dir2, dir3, dir4)
        deallocate (coefI, coefII, coefIII, coefIV)
        deallocate (f1, Rho, Ux, Uy)
        deallocate (f1_west_snd, f1_west_rcv, f1_east_snd, f1_east_rcv)
        deallocate (f1_north_snd, f1_north_rcv, f1_south_snd, f1_south_rcv)
    end subroutine memFree

    function itoa(i) result(res)
        character(:),allocatable :: res
        integer,intent(in) :: i
        character(range(i)+2) :: tmp
        write(tmp,'(i0)') i
        res = trim(tmp)
    end function
end module solver
