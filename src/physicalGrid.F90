!=======================================================================
!> @brief Physical space configurations
!=======================================================================
module physicalGrid
implicit none
save

!domain defination, to be read from NML: velocityNml
! Total domain size
integer :: Nx, Ny
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

! number of wall points
integer :: nWall 

! PROC ID
integer :: peid

integer, DIMENSION(:), ALLOCATABLE :: vecWall
integer, DIMENSION(:), ALLOCATABLE :: dir1, dir2, dir3, dir4 ! dir1(i) is the ith 
double precision, DIMENSION(:,:), ALLOCATABLE :: coefI, coefII, coefIII, coefIV
double precision, DIMENSION(:,:), ALLOCATABLE :: extCoef

contains 
    ! should be called after calling MPIParams::setupVirtualProcessGrid
    ! such that xl, xlg, xu, xug and yl, ylg, yu, yug has been set already
    subroutine setupPhysicalGrid
        implicit none

        ! local vars
        integer :: localid, i, j, l, icount, countVoidP, NneighborFluid
        integer :: bxl, bxu, byl, byu ! bound when counting fluid point for sweeping
        integer :: ii, jj

        ds = 0.5d0/(Ny-1)

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
        Open(IMAGEFILE,file=imageFileName,status='OLD')
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
        ! PRINT*, real_porosity

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
        If (ghostLayers>0) then
            Do j=ylg,yug
                Do i=xlg,xug
                    ! test if boundary processor, here the ghost flag doesn't means 
                    ! the communication boundary, but only means the true global 
                    ! outter boarder.
                    If ((i<xmin).OR.(i>xmax).OR.(j<ymin).OR.(j>ymax)) then 
                        localid = (j-ylg)*Nxtotal + i-xlg+1
                        image(localid) = ghost
                    End if
                Enddo
            Enddo
        End if

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
        Do j=byl,byu
            Do i=bxl,bxu
                !ii = i-xl+1
                !jj = j-yl+1
                localid = (j-ylg)*Nxtotal + i-xlg+1
                If (array2D(i,j)==solid) then
                    NneighborFluid=1 !(1-2-3-4 in D2Q9 corresponding to 2-3-5-7)
                    ! found bug here, array index out of bound of array2D
                    !If (image(localid+1)==fluid)  NneighborFluid=NneighborFluid*2 
                    !If (image(localid+Nxtotal)==fluid)  NneighborFluid=NneighborFluid*3
                    !If (image(localid-1)==fluid)  NneighborFluid=NneighborFluid*5 
                    !If (image(localid-Nxtotal)==fluid)  NneighborFluid=NneighborFluid*7
                    If (array2Dg(i+1, j)==fluid)  NneighborFluid=NneighborFluid*2   !Check neighbor on the East(2)
                    If (array2Dg(i, j+1)==fluid)  NneighborFluid=NneighborFluid*3   !Check neighbor on the North(3)
                    If (array2Dg(i-1, j)==fluid)  NneighborFluid=NneighborFluid*5   !Check neighbor on the West(5)
                    If (array2Dg(i, j-1)==fluid)  NneighborFluid=NneighborFluid*7   !Check neighbor on the South(7)

                    SELECT Case (NneighborFluid)
                        CASE (2)
                            image(localid)=WallXp
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        CASE (3)
                            image(localid)=WallYp
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        CASE (5)
                            image(localid)=WallXn
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        CASE (7)
                            image(localid)=WallYn
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        CASE (6)
                            image(localid)=WallXpYp
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        CASE (15)
                            image(localid)=WallXnYp
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        CASE (35)
                            image(localid)=WallXnYn
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        CASE (14)
                            image(localid)=WallXpYn
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                    END SELECT
                Endif
            Enddo
        Enddo

        PRINT*, "nWall = ", nWall



        ! Create wall-type vectors
        !vecWall(walli) mark the global id of the walli'th wall in the image 
        ALLOCATE(vecWall(nWall))
        nWall=0
        Do j=yl, yu
            Do i=xl, xu
                localid = (j-ylg)*Nxtotal + i-xlg+1
                SELECT Case (image(localid))
                    CASE (WallXp)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    CASE (WallYp)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    CASE (WallXn)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    CASE (WallYn)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    CASE (WallXpYp)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    CASE (WallXnYp)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    CASE (WallXnYn)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    CASE (WallXpYn)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                END SELECT
            Enddo
        Enddo

        ! extrapolation coefficients
        ALLOCATE(extCoef(nWall,2))
        extCoef(:, 1) =  2.d0
        extCoef(:, 2) = -1.d0
        ! set extCoef, if near communication boundary, use 1st order (mixed)
        do l=1, nWall
            localid = vecWall(l)
            i = xlg -1 + mod(localid, Nxtotal) ! wall location in x 
            j = ylg + localid / Nxtotal     ! wall location in y
            SELECT Case (image(localid))
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
            end select
        enddo

        if(wallExtOrder == 2) then
            ! reset to 2nd order
            extCoef(:, 1) =  2.d0
            extCoef(:, 2) = -1.d0
        else if (wallExtOrder == 1) then
            ! reset to 1st order
            extCoef(:, 1) =  1.d0
            extCoef(:, 2) =  0.d0
        else if (wallExtOrder /= 3) then
            PRINT*, "Error: wallExtOrder wroond shoud be [1|2|3]"
        endif

        !Direction 1
        !bound
        bxl = xl
        bxu = xu
        byl = yl
        byu = yu
        if(xl == xmin) bxl = xl + 1 !if most west block(processor)
        if(yl == ymin) byl = yl + 1 !if most south block
        !if(xu == xmax) bxu = xu     !if most east block(processor)
        !if(yu == ymax) byu = yu     !if most north block   
        !count fluid points when sweeping from 1st direction
        Nstencil1=0
        Do j=byl,byu
            Do i=bxl,bxu
                localid = (j-ylg)*Nxtotal + i-xlg+1
                If (image(localid)==fluid) then
                    Nstencil1=Nstencil1+1
                End if
            End do
        End do
        !DEBUG
        print*, "Proc ", peid, "Nstencil1 = ", Nstencil1
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir1(Nstencil1))
        icount=0
        Do j=byl,byu
            Do i=bxl,bxu
                localid = (j-ylg)*Nxtotal + i-xlg+1
                If (image(localid)==fluid) then
                    icount=icount+1
                    dir1(icount)=localid
                End if
            End do
        End do

        !Direction 2
        !bound
        bxl = xl
        bxu = xu
        byl = yl
        byu = yu
        if(xu == xmax) bxu = xu - 1 !if most east block
        if(yl == ymin) byl = yl + 1 !if most south block
        !if(xl == xmin) bxl = xl     !if most west block
        !if(yu == ymax) byu = yu     !if most north block
        !count fluid points when sweeping from 2nd direction
        Nstencil2=0
        Do j=byl,byu
            Do i=bxu,bxl,-1
                localid = (j-ylg)*Nxtotal + i-xlg+1
                If (image(localid)==fluid) then
                    Nstencil2=Nstencil2+1
                End if
            End do
        End do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir2(Nstencil2))
        icount=0
        Do j=byl,byu
            Do i=bxu,bxl,-1
                localid = (j-ylg)*Nxtotal + i-xlg+1
                If (image(localid)==fluid) then
                    icount=icount+1
                    dir2(icount)=localid
                End if
            End do
        End do

        !Direction 3
        !bound
        bxl = xl
        bxu = xu
        byl = yl
        byu = yu
        if(xu == xmax) bxu = xu - 1 !if most east block
        if(yu == ymax) byu = yu - 1 !if most north block
        !if(xl == xmin) bxl = xl     !if most west block
        !if(yl == ymin) byl = yl     !if most south block
        !count fluid points when sweeping from 3rd direction
        Nstencil3=0
        Do j=byu,byl,-1
            Do i=bxu,bxl,-1
                localid = (j-ylg)*Nxtotal + i-xlg+1
                If (image(localid)==fluid) then
                    Nstencil3=Nstencil3+1
                End if
            End do
        End do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir3(Nstencil3))
        icount=0
        Do j=byu,byl,-1
            Do i=bxu,bxl,-1
                localid = (j-ylg)*Nxtotal + i-xlg+1
                If (image(localid)==fluid) then
                    icount=icount+1
                    dir3(icount)=localid
                End if
            End do
        End do

        !Direction 4
        !bound
        bxl = xl
        bxu = xu
        byl = yl
        byu = yu
        if(xl == xmin) bxl = xl + 1 !if most west block
        if(yu == ymax) byu = yu - 1 !if most north block
        !if(xu == xmax) bxu = xu     !if most east block
        !if(yl == ymin) byl = yl     !if most south block
        !count fluid points when sweeping from 4th direction
        Nstencil4=0
        Do j=byu,byl,-1
            Do i=bxl,bxu
                localid = (j-ylg)*Nxtotal + i-xlg+1
                If (image(localid)==fluid) then
                    Nstencil4=Nstencil4+1
                End if
            End do
        End do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir4(Nstencil4))
        icount=0
        Do j=byu,byl,-1
            Do i=bxl,bxu
                localid = (j-ylg)*Nxtotal + i-xlg+1
                If (image(localid)==fluid) then
                    icount=icount+1
                    dir4(icount)=localid
                End if
            End do
        End do

        !construct the deferencial coefficients
        ALLOCATE(coefI(Nstencil1,6),coefII(Nstencil2,6),coefIII(Nstencil3,6),coefIV(Nstencil4,6))
        Do i=1,Nstencil1
            localid=dir1(i) 
            !2nd order of accuracy
            coefI(i,1)=1.5d0/ds !x0n
            coefI(i,2)=2.d0/ds  !x1n
            coefI(i,3)=-0.5d0/ds   !x2n
            coefI(i,4)=1.5d0/ds !y0n
            coefI(i,5)=2.d0/ds  !y1n
            coefI(i,6)=-0.5d0/ds   !y2n
            if ((image(localid-2)==ghost).OR.(image(localid-1)==WallXp).OR.(image(localid-1)==WallXpYn) &
            & .OR.(image(localid-1)==WallXpYp)) then !1st order of accuracy in x
                coefI(i,1)=1.d0/ds  !x0n
                coefI(i,2)=1.d0/ds  !x1n
                coefI(i,3)=0.d0     !x2n
            end if
            if ((image(localid-2*Nxtotal)==ghost).OR.(image(localid-Nxtotal)==WallYp) &
            & .OR.(image(localid-Nxtotal)==WallXpYp).OR.(image(localid-Nxtotal)==WallXnYp)) then   !1st order of accuracy in y
                coefI(i,4)=1.d0/ds  !y0n
                coefI(i,5)=1.d0/ds  !y1n
                coefI(i,6)=0.d0     !y2n
            end if
        End do

        Do i=1,Nstencil2
            localid=dir2(i)
            !2nd order of accuracy
            coefII(i,1)=-1.5d0/ds !x0n
            coefII(i,2)=-2.d0/ds  !x1n
            coefII(i,3)=0.5d0/ds   !x2n
            coefII(i,4)=1.5d0/ds !y0n
            coefII(i,5)=2.d0/ds  !y1n
            coefII(i,6)=-0.5d0/ds   !y2n
            if ((image(localid+2)==ghost).OR.(image(localid+1)==WallXn).OR.(image(localid+1)==WallXnYp) &
            & .OR.(image(localid+1)==WallXnYn)) then !1st order of accuracy in x
                coefII(i,1)=-1.d0/ds  !x0n
                coefII(i,2)=-1.d0/ds  !x1n
                coefII(i,3)=0.d0     !x2n
            end if
            if ((image(localid-2*Nxtotal)==ghost).OR.(image(localid-Nxtotal)==WallYp) &
            & .OR.(image(localid-Nxtotal)==WallXpYp).OR.(image(localid-Nxtotal)==WallXnYp)) then   !1st order of accuracy in y
                coefII(i,4)=1.d0/ds  !y0n
                coefII(i,5)=1.d0/ds  !y1n
                coefII(i,6)=0.d0     !y2n
            end if
        End do

        Do i=1,Nstencil3
            localid=dir3(i)
            !2nd order of accuracy
            coefIII(i,1)=-1.5d0/ds !x0n
            coefIII(i,2)=-2.d0/ds  !x1n
            coefIII(i,3)=0.5d0/ds   !x2n
            coefIII(i,4)=-1.5d0/ds !y0n
            coefIII(i,5)=-2.d0/ds  !y1n
            coefIII(i,6)=0.5d0/ds   !y2n
            if ((image(localid+2)==ghost).OR.(image(localid+1)==WallXn).OR.(image(localid+1)==WallXnYp) &
            & .OR.(image(localid+1)==WallXnYn)) then !1st order of accuracy in x
                coefIII(i,1)=-1.d0/ds  !x0n
                coefIII(i,2)=-1.d0/ds  !x1n
                coefIII(i,3)=0.d0     !x2n
            end if
            if ((image(localid+2*Nxtotal)==ghost).OR.(image(localid+Nxtotal)==WallYn) &
            & .OR.(image(localid+Nxtotal)==WallXnYn).OR.(image(localid+Nxtotal)==WallXpYn)) then   !1st order of accuracy in y
                coefIII(i,4)=-1.d0/ds  !y0n
                coefIII(i,5)=-1.d0/ds  !y1n
                coefIII(i,6)=0.d0     !y2n
            end if
        End do

        Do i=1,Nstencil4
            localid=dir4(i)
            !2nd order of accuracy
            coefIV(i,1)=1.5d0/ds !x0n
            coefIV(i,2)=2.d0/ds  !x1n
            coefIV(i,3)=-0.5d0/ds   !x2n
            coefIV(i,4)=-1.5d0/ds !y0n
            coefIV(i,5)=-2.d0/ds  !y1n
            coefIV(i,6)=0.5d0/ds   !y2n
            if ((image(localid-2)==ghost).OR.(image(localid-1)==WallXp).OR.(image(localid-1)==WallXpYn) &
            & .OR.(image(localid-1)==WallXpYp)) then !1st order of accuracy in x
                coefIV(i,1)=1.d0/ds  !x0n
                coefIV(i,2)=1.d0/ds  !x1n
                coefIV(i,3)=0.d0     !x2n
            end if
            if ((image(localid+2*Nxtotal)==ghost).OR.(image(localid+Nxtotal)==WallYn) &
            & .OR.(image(localid+Nxtotal)==WallXnYn).OR.(image(localid+Nxtotal)==WallXpYn)) then   !1st order of accuracy in y
                coefIV(i,4)=-1.d0/ds  !y0n
                coefIV(i,5)=-1.d0/ds  !y1n
                coefIV(i,6)=0.d0     !y2n
            end if
        End do
    end subroutine setupPhysicalGrid
end module physicalGrid
