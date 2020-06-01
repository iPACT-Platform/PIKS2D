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

! Copyright (c) 2017-2020 iPACT-Platform

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
! FITNESS .or.A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUT.or..or.COPYRIGHT HOLDERS BE LIABLE .or.ANY CLAIM, DAMAGES.or.OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, .or..or.OTHERWISE, ARISING FROM,
! OUT OF.or.IN CONNECTION WITH THE SOFTWARE.or.THE USE.or.OTHER DEALINGS IN THE
! SOFTWARE.

!-------------------------------------------------------------------------------
! module    : physicalGrid
!-------------------------------------------------------------------------------
! This is a module for physical (spatial) space configurations of 2D DVM parallel solver. 
! For details:
!
! [1]   M.T. Ho, L. Zhu, L. Wu, P. Wang, Z. Guo, Z.-H. Li, Y. Zhang
!       "A multi-level parallel solver for rarefied gas flows in porous media"
! 		Computer Physics Communications, 234 (2019), pp. 14-25
!
!	Surface nodes are categorized into flat, 2-fold corner 
!   so the proper diffusive BC can be applied later (Fig.1 of Ref.[1]). 
!	For fluid nodes, there are 4 groups of molecular velocity and hence 4 paths of sweep. 
!   Weighting of fluid nodes in upwind schemes are store in arrays 
!   (See Section 2.3, 3.1, Figure 1 of Ref.[1]). 
!   At least 2 fluid layers inside each pore is assumed for the input digital image.
!
!! dir1,dir2,..,dir4 : 4 paths of sweep
!! coef1,coef2,..,coef4: arrays storing weighting of fluid nodes in upwind schemes
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

module physicalGrid
implicit none
save

!domain defination, to be read from NML: velocityNml
! Total domain size
integer :: Nx, Ny
double precision :: xPadRatio

! image file name
character(len=30):: imageFileName
! wall extrapolation order, 1, 2, 3(mixed)
integer :: wallExtOrder

!repetition times of the base image, used only for weak scaling study
! default values are both 1.
integer :: block_repx
integer :: block_repy
!the based (readed in Nx,Ny) image size,
! used only for weak scaling study
integer :: Nx_base, Ny_base

!subdomain bounds
integer :: xl, xu, yl, yu
!subdomain with ghostLayers bounds
integer :: xlg, xug, ylg, yug
! ghostLayers
integer, parameter :: ghostLayers = 2
! 2D porous structure image file (ascii)
integer, parameter :: IMAGEFILE = 20

integer :: xmin, xmax, ymin, ymax
integer :: Nxtotal, Nytotal, Nxsub, Nysub, Ntotal
double precision :: ds
integer, parameter :: column=2 ! layer to extract flow rate

! raw and extended flag array
integer, dimension (:,:), allocatable :: array2D, array2Dg
integer, dimension (:), allocatable :: image
! grid point flags
integer, parameter :: fluid = 0, solid = 1, ghost=4
integer, parameter :: WallXpYp = 20, WallYp = 21, WallXnYp = 22, WallXn = 23
integer, parameter :: WallXnYn = 24, WallYn = 25, WallXpYn = 26, WallXp = 27

integer :: Nstencil1, Nstencil2, Nstencil3, Nstencil4
double precision :: RhoWall
double precision :: real_porosity

integer :: Nfluid
integer, allocatable, dimension(:) :: mapF

! number of wall points
integer :: nWall 

! PROC ID
integer :: peid

integer, dimension(:), allocatable :: vecWall
integer, dimension(:), allocatable :: dir1, dir2, dir3, dir4 ! dir1(i) is the ith 
double precision, dimension(:,:), allocatable :: coefI, coefII, coefIII, coefIV
double precision, dimension(:,:), allocatable :: extCoef

contains 
    ! should be called after calling MPIParams::setupVirtualProcessGrid
    ! such that xl, xlg, xu, xug and yl, ylg, yu, yug has been set already
    subroutine setupPhysicalGrid
        implicit none

        ! local vars
        integer :: localid, i, j, l, icount, countVoidP, NneighborFluid
        integer :: bxl, bxu, byl, byu ! bound when counting fluid point for sweeping
        integer :: ii, jj

        !ds = dble(NY)/dble(NX)/(Ny-1)
        ds = (1.0d0 + xPadRatio)/(NX-1)
        print*, "ds = ", ds

        ! set the extend and sizes
        Nxtotal = xug - xlg + 1
        Nytotal = yug - ylg + 1
        Nxsub = xu - xl + 1
        Nysub = yu - yl + 1
        Ntotal = Nxtotal * Nytotal

        allocate(array2D(Nx, Ny)) ! this is global raw geometry data
        !including ghost layer (flag = ghost)
        allocate(array2Dg(xmin-ghostLayers:xmax+ghostLayers, &
                          ymin-ghostLayers:ymax+ghostLayers))
        allocate(image(Ntotal)) ! this is local one

        !switch for debuging
        !read digital image
        array2D = 0 
        open(IMAGEFILE,file=imageFileName,status='OLD')
            do j=1,Ny_base
                read(IMAGEFILE, *) (array2D(i,j), i=1, Nx_base) !NOTE: add extral layer
            enddo
        close(IMAGEFILE)

        ! repeat the base block 
        do j = 0, block_repy-1
            do i = 0, block_repx-1
                do jj = 1, Ny_base
                    do ii = 1, Nx_base
                        array2D(i*Nx_base + ii, j*Ny_base + jj) &
                            = array2D(ii, jj)
                    enddo
                enddo
            enddo
        enddo

        ! set array2g
        array2Dg = ghost ! outer bound
        do j = ymin, ymax
            do i = xmin, xmax
                array2Dg(i,j) = array2D(i,j)
            enddo
        enddo

        ! set array2D
        countVoidP = 0 ! count the void grid points
        do j=1,Ny
            do i=1,Nx
                if (array2D(i,j)==0) then 
                    array2D(i,j)=fluid
                    countVoidP = countVoidP + 1
                else
                    array2D(i,j)=solid
                endif
            enddo
        enddo

        ! calu. poresity
        ! real_porosity=real(countVoidP)/real(Nx*Ny)
        ! print*, real_porosity

        ! set local image
        image = fluid
        bxl = xlg
        byl = ylg
        bxu = xug
        byu = yug
        if(xl == xmin) bxl = xl !if most west block(processor)
        if(xu == xmax) bxu = xu !if most east block(processor)
        if(yl == ymin) byl = yl !if most south block
        if(yu == ymax) byu = yu !if most north block
        do j = byl, byu
            do i = bxl, bxu
               localid = (j-ylg)*Nxtotal + i-xlg+1
               image(localid) = array2D(i,j)
           enddo
        enddo

        !Assign the outerboarder layers (called ghost point in serial program)
        if (ghostLayers>0) then
            do j=ylg,yug
                do i=xlg,xug
                    ! test if boundary processor, here the ghost flag doesn't means 
                    ! the communication boundary, but only means the true global 
                    ! outter boarder.
                    if ((i<xmin).or.(i>xmax).or.(j<ymin).or.(j>ymax)) then 
                        localid = (j-ylg)*Nxtotal + i-xlg+1
                        image(localid) = ghost
                    endif
                enddo
            enddo
        endif

        bxl = xlg
        bxu = xug
        byl = ylg
        byu = yug
        if(xl == xmin) bxl = xl !if most west block(processor)
        if(xu == xmax) bxu = xu !if most east block(processor)
        if(yl == ymin) byl = yl !if most south block
        if(yu == ymax) byu = yu !if most north block

        ! set wall points type based on sournding point type(f/s)
        nWall=0 ! count the wall points, only in the l, u range are counted!
        do j=byl,byu
            do i=bxl,bxu
                !ii = i-xl+1
                !jj = j-yl+1
                localid = (j-ylg)*Nxtotal + i-xlg+1
                if (array2D(i,j)==solid) then
                    NneighborFluid=1 !(1-2-3-4 in D2Q9 corresponding to 2-3-5-7)
                    ! found bug here, array index out of bound of array2D
                    !if (image(localid+1)==fluid)  NneighborFluid=NneighborFluid*2 
                    !if (image(localid+Nxtotal)==fluid)  NneighborFluid=NneighborFluid*3
                    !if (image(localid-1)==fluid)  NneighborFluid=NneighborFluid*5 
                    !if (image(localid-Nxtotal)==fluid)  NneighborFluid=NneighborFluid*7
                    if (array2Dg(i+1, j)==fluid)  NneighborFluid=NneighborFluid*2   !Check neighbor on the East(2)
                    if (array2Dg(i, j+1)==fluid)  NneighborFluid=NneighborFluid*3   !Check neighbor on the North(3)
                    if (array2Dg(i-1, j)==fluid)  NneighborFluid=NneighborFluid*5   !Check neighbor on the West(5)
                    if (array2Dg(i, j-1)==fluid)  NneighborFluid=NneighborFluid*7   !Check neighbor on the South(7)

                    select case (NneighborFluid)
                        case (2)
                            image(localid)=WallXp
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        case (3)
                            image(localid)=WallYp
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        case (5)
                            image(localid)=WallXn
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        case (7)
                            image(localid)=WallYn
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        case (6)
                            image(localid)=WallXpYp
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        case (15)
                            image(localid)=WallXnYp
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        case (35)
                            image(localid)=WallXnYn
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        case (14)
                            image(localid)=WallXpYn
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                    endselect
                endif
            enddo
        enddo

        print*, "nWall = ", nWall



        ! Create wall-type vectors
        !vecWall(walli) mark the global id of the walli'th wall in the image 
        allocate(vecWall(nWall))
        nWall=0
        do j=yl, yu
            do i=xl, xu
                localid = (j-ylg)*Nxtotal + i-xlg+1
                select case (image(localid))
                    case (WallXp)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    case (WallYp)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    case (WallXn)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    case (WallYn)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    case (WallXpYp)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    case (WallXnYp)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    case (WallXnYn)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    case (WallXpYn)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                endselect
            enddo
        enddo

        ! extrapolation coefficients
        allocate(extCoef(nWall,2))
        extCoef(:, 1) =  2.d0
        extCoef(:, 2) = -1.d0
        ! set extCoef, if near communication boundary, use 1st order (mixed)
        do l=1, nWall
            localid = vecWall(l)
            i = xlg -1 + mod(localid, Nxtotal) ! wall location in x 
            j = ylg + localid / Nxtotal     ! wall location in y
            select case (image(localid))
                case (WallXp)
                    if(i == xu-1 .and. image(localid+1) == fluid) then
                        extCoef(l, 1) = 1.d0
                        extCoef(l, 2) = 0.d0
                    endif
                case (WallXn)
                    if(i == xl+1 .and. image(localid-1) == fluid) then
                        extCoef(l, 1) = 1.d0
                        extCoef(l, 2) = 0.d0
                    endif
                case (WallYp)
                    if(j == yu-1 .and. image(localid+Nxtotal) == fluid) then
                        extCoef(l, 1) = 1.d0
                        extCoef(l, 2) = 0.d0
                    endif
                case (WallYn)
                    if(j == yl+1 .and. image(localid-Nxtotal) == fluid) then
                        extCoef(l, 1) = 1.d0
                        extCoef(l, 2) = 0.d0
                    endif
                case (WallXpYp)
                    if(i == xu-1 .and. image(localid+1) == fluid) then
                        extCoef(l, 1) = 1.d0
                        extCoef(l, 2) = 0.d0
                    endif
                    if(j == yu-1 .and. image(localid+Nxtotal) == fluid) then
                        extCoef(l, 1) = 1.d0
                        extCoef(l, 2) = 0.d0
                    endif
                case (WallXnYp)
                    if(i == xl+1 .and. image(localid-1) == fluid) then
                        extCoef(l, 1) = 1.d0
                        extCoef(l, 2) = 0.d0
                    endif
                    if(j == yu-1 .and. image(localid+Nxtotal) == fluid) then
                        extCoef(l, 1) = 1.d0
                        extCoef(l, 2) = 0.d0
                    endif
                case (WallXpYn)
                    if(i == xu-1 .and. image(localid+1) == fluid) then
                        extCoef(l, 1) = 1.d0
                        extCoef(l, 2) = 0.d0
                    endif
                    if(j == yl+1 .and. image(localid-Nxtotal) == fluid) then
                        extCoef(l, 1) = 1.d0
                        extCoef(l, 2) = 0.d0
                    endif
                case (WallXnYn)
                    if(i == xl+1 .and. image(localid-1) == fluid) then
                        extCoef(l, 1) = 1.d0
                        extCoef(l, 2) = 0.d0
                    endif
                    if(j == yl+1 .and. image(localid-Nxtotal) == fluid) then
                        extCoef(l, 1) = 1.d0
                        extCoef(l, 2) = 0.d0
                    endif
            endselect
        enddo

        if(wallExtOrder == 2) then
            ! reset to 2nd order
            extCoef(:, 1) =  2.d0
            extCoef(:, 2) = -1.d0
        elseif (wallExtOrder == 1) then
            ! reset to 1st order
            extCoef(:, 1) =  1.d0
            extCoef(:, 2) =  0.d0
        elseif (wallExtOrder /= 3) then
            print*, "Error: wallExtOrder wroond shoud be [1|2|3]"
        endif

        ! count num of non-solid points
        Nfluid=0
        do j=yl,yu
            do i=xl,xu
                localid = (j-ylg)*Nxtotal + i-xlg+1
                if (image(localid) == fluid) then
                    Nfluid=Nfluid+1
                endif
            enddo
        enddo
        allocate(mapF(Nfluid))

        Nfluid=0
        do j=yl,yu
            do i=xl,xu
                localid = (j-ylg)*Nxtotal + i-xlg+1
                if (image(localid) == fluid) then
                    Nfluid=Nfluid+1
                    mapF(Nfluid) = localid
                endif
            enddo
        enddo

        !Direction 1
        !bound
        bxl = xl
        bxu = xu
        byl = yl
        byu = yu
        if (xl == xmin) bxl = xl + 1 !if most west block(processor)
        if (yl == ymin) byl = yl + 1 !if most south block
        !if(xu == xmax) bxu = xu     !if most east block(processor)
        !if(yu == ymax) byu = yu     !if most north block   
        !count fluid points when sweeping from 1st direction
        Nstencil1=0
        do j=byl,byu
            do i=bxl,bxu
                localid = (j-ylg)*Nxtotal + i-xlg+1
                if (image(localid)==fluid) then
                    Nstencil1=Nstencil1+1
                endif
            enddo
        enddo
        !DEBUG
        print*, "Proc ", peid, "Nstencil1 = ", Nstencil1
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir1(Nstencil1))
        icount=0
        do j=byl,byu
            do i=bxl,bxu
                localid = (j-ylg)*Nxtotal + i-xlg+1
                if (image(localid)==fluid) then
                    icount=icount+1
                    dir1(icount)=localid
                endif
            enddo
        enddo

        !Direction 2
        !bound
        bxl = xl
        bxu = xu
        byl = yl
        byu = yu
        if (xu == xmax) bxu = xu - 1 !if most east block
        if (yl == ymin) byl = yl + 1 !if most south block
        !if(xl == xmin) bxl = xl     !if most west block
        !if(yu == ymax) byu = yu     !if most north block
        !count fluid points when sweeping from 2nd direction
        Nstencil2=0
        do j=byl,byu
            do i=bxu,bxl,-1
                localid = (j-ylg)*Nxtotal + i-xlg+1
                if (image(localid)==fluid) then
                    Nstencil2=Nstencil2+1
                endif
            enddo
        enddo
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir2(Nstencil2))
        icount=0
        do j=byl,byu
            do i=bxu,bxl,-1
                localid = (j-ylg)*Nxtotal + i-xlg+1
                if (image(localid)==fluid) then
                    icount=icount+1
                    dir2(icount)=localid
                endif
            enddo
        enddo

        !Direction 3
        !bound
        bxl = xl
        bxu = xu
        byl = yl
        byu = yu
        if (xu == xmax) bxu = xu - 1 !if most east block
        if (yu == ymax) byu = yu - 1 !if most north block
        !if(xl == xmin) bxl = xl     !if most west block
        !if(yl == ymin) byl = yl     !if most south block
        !count fluid points when sweeping from 3rd direction
        Nstencil3=0
        do j=byu,byl,-1
            do i=bxu,bxl,-1
                localid = (j-ylg)*Nxtotal + i-xlg+1
                if (image(localid)==fluid) then
                    Nstencil3=Nstencil3+1
                endif
            enddo
        enddo
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir3(Nstencil3))
        icount=0
        do j=byu,byl,-1
            do i=bxu,bxl,-1
                localid = (j-ylg)*Nxtotal + i-xlg+1
                if (image(localid)==fluid) then
                    icount=icount+1
                    dir3(icount)=localid
                endif
            enddo
        enddo

        !Direction 4
        !bound
        bxl = xl
        bxu = xu
        byl = yl
        byu = yu
        if (xl == xmin) bxl = xl + 1 !if most west block
        if (yu == ymax) byu = yu - 1 !if most north block
        !if(xu == xmax) bxu = xu     !if most east block
        !if(yl == ymin) byl = yl     !if most south block
        !count fluid points when sweeping from 4th direction
        Nstencil4=0
        do j=byu,byl,-1
            do i=bxl,bxu
                localid = (j-ylg)*Nxtotal + i-xlg+1
                if (image(localid)==fluid) then
                    Nstencil4=Nstencil4+1
                endif
            enddo
        enddo
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir4(Nstencil4))
        icount=0
        do j=byu,byl,-1
            do i=bxl,bxu
                localid = (j-ylg)*Nxtotal + i-xlg+1
                if (image(localid)==fluid) then
                    icount=icount+1
                    dir4(icount)=localid
                endif
            enddo
        enddo

        !construct the deferencial coefficients
        allocate(coefI(Nstencil1,6),coefII(Nstencil2,6),coefIII(Nstencil3,6),coefIV(Nstencil4,6))
        do i=1,Nstencil1
            localid=dir1(i) 
            !2nd order of accuracy
            coefI(i,1)=1.5d0/ds !x0n
            coefI(i,2)=2.d0/ds  !x1n
            coefI(i,3)=-0.5d0/ds   !x2n
            coefI(i,4)=1.5d0/ds !y0n
            coefI(i,5)=2.d0/ds  !y1n
            coefI(i,6)=-0.5d0/ds   !y2n
            if ((image(localid-2)==ghost).or.(image(localid-1)==WallXp).or.(image(localid-1)==WallXpYn) &
            & .or.(image(localid-1)==WallXpYp)) then !1st order of accuracy in x
                coefI(i,1)=1.d0/ds  !x0n
                coefI(i,2)=1.d0/ds  !x1n
                coefI(i,3)=0.d0     !x2n
            endif
            if ((image(localid-2*Nxtotal)==ghost).or.(image(localid-Nxtotal)==WallYp) &
            & .or.(image(localid-Nxtotal)==WallXpYp).or.(image(localid-Nxtotal)==WallXnYp)) then   !1st order of accuracy in y
                coefI(i,4)=1.d0/ds  !y0n
                coefI(i,5)=1.d0/ds  !y1n
                coefI(i,6)=0.d0     !y2n
            endif
        enddo

        do i=1,Nstencil2
            localid=dir2(i)
            !2nd order of accuracy
            coefII(i,1)=-1.5d0/ds !x0n
            coefII(i,2)=-2.d0/ds  !x1n
            coefII(i,3)=0.5d0/ds   !x2n
            coefII(i,4)=1.5d0/ds !y0n
            coefII(i,5)=2.d0/ds  !y1n
            coefII(i,6)=-0.5d0/ds   !y2n
            if ((image(localid+2)==ghost).or.(image(localid+1)==WallXn).or.(image(localid+1)==WallXnYp) &
            & .or.(image(localid+1)==WallXnYn)) then !1st order of accuracy in x
                coefII(i,1)=-1.d0/ds  !x0n
                coefII(i,2)=-1.d0/ds  !x1n
                coefII(i,3)=0.d0     !x2n
            endif
            if ((image(localid-2*Nxtotal)==ghost).or.(image(localid-Nxtotal)==WallYp) &
            & .or.(image(localid-Nxtotal)==WallXpYp).or.(image(localid-Nxtotal)==WallXnYp)) then   !1st order of accuracy in y
                coefII(i,4)=1.d0/ds  !y0n
                coefII(i,5)=1.d0/ds  !y1n
                coefII(i,6)=0.d0     !y2n
            endif
        enddo

        do i=1,Nstencil3
            localid=dir3(i)
            !2nd order of accuracy
            coefIII(i,1)=-1.5d0/ds !x0n
            coefIII(i,2)=-2.d0/ds  !x1n
            coefIII(i,3)=0.5d0/ds   !x2n
            coefIII(i,4)=-1.5d0/ds !y0n
            coefIII(i,5)=-2.d0/ds  !y1n
            coefIII(i,6)=0.5d0/ds   !y2n
            if ((image(localid+2)==ghost).or.(image(localid+1)==WallXn).or.(image(localid+1)==WallXnYp) &
            & .or.(image(localid+1)==WallXnYn)) then !1st order of accuracy in x
                coefIII(i,1)=-1.d0/ds  !x0n
                coefIII(i,2)=-1.d0/ds  !x1n
                coefIII(i,3)=0.d0     !x2n
            endif
            if ((image(localid+2*Nxtotal)==ghost).or.(image(localid+Nxtotal)==WallYn) &
            & .or.(image(localid+Nxtotal)==WallXnYn).or.(image(localid+Nxtotal)==WallXpYn)) then   !1st order of accuracy in y
                coefIII(i,4)=-1.d0/ds  !y0n
                coefIII(i,5)=-1.d0/ds  !y1n
                coefIII(i,6)=0.d0     !y2n
            endif
        enddo

        do i=1,Nstencil4
            localid=dir4(i)
            !2nd order of accuracy
            coefIV(i,1)=1.5d0/ds !x0n
            coefIV(i,2)=2.d0/ds  !x1n
            coefIV(i,3)=-0.5d0/ds   !x2n
            coefIV(i,4)=-1.5d0/ds !y0n
            coefIV(i,5)=-2.d0/ds  !y1n
            coefIV(i,6)=0.5d0/ds   !y2n
            if ((image(localid-2)==ghost).or.(image(localid-1)==WallXp).or.(image(localid-1)==WallXpYn) &
            & .or.(image(localid-1)==WallXpYp)) then !1st order of accuracy in x
                coefIV(i,1)=1.d0/ds  !x0n
                coefIV(i,2)=1.d0/ds  !x1n
                coefIV(i,3)=0.d0     !x2n
            endif
            if ((image(localid+2*Nxtotal)==ghost).or.(image(localid+Nxtotal)==WallYn) &
            & .or.(image(localid+Nxtotal)==WallXnYn).or.(image(localid+Nxtotal)==WallXpYn)) then   !1st order of accuracy in y
                coefIV(i,4)=-1.d0/ds  !y0n
                coefIV(i,5)=-1.d0/ds  !y1n
                coefIV(i,6)=0.d0     !y2n
            endif
        enddo
    end subroutine setupPhysicalGrid
end module physicalGrid
