!=======================================================================
!> @brief Physical space configurations
!=======================================================================
module physicalGrid
IMPLICIT NONE
SAVE

! Domain size

integer :: xl, xu, yl, yu
integer :: xlg, xug, ylg, yug
integer, parameter :: ghostLayers = 2


! NX and NY is the global grid size
integer, parameter :: Nx = 120, Ny = 120
integer, parameter :: xmin = 1
integer, parameter :: xmax = Nx
integer, parameter :: ymin = 1
integer, parameter :: ymax = Ny
integer :: Nxtotal, Nytotal, Nxsub, Nysub, Ntotal
double precision, parameter :: ds =  0.5d0/(Ny-1)
integer, parameter :: column=2 ! layer to extract flow rate

! raw and extended flag array
integer, dimension (:,:), allocatable :: array2D
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

integer, DIMENSION(:), ALLOCATABLE :: vecWall
integer, DIMENSION(:), ALLOCATABLE :: dir1, dir2, dir3, dir4 ! dir1(i) is the ith 
double precision, DIMENSION(:,:), ALLOCATABLE :: coefI, coefII, coefIII, coefIV

contains 
    ! should be called after calling MPIParams::setupVirtualProcessGrid
    ! such that xl, xlg, xu, xug and yl, ylg, yu, yug has been set already
    subroutine setupPhysicalGrid
        implicit none

        ! local vars
        integer :: localid, i, j, icount, countVoidP, NneighborFluid
        integer :: bxl, bxu, byl, byu ! bound when counting fluid point for sweeping

        ! set the extend and sizes
        Nxtotal = xug - xlg + 1
        Nytotal = yug - ylg + 1
        Nxsub = xu - xl + 1
        Nysub = yu - yl + 1
        Ntotal = Nxtotal * Nytotal

        allocate(array2D(Nx, Ny)) ! this is global raw geometry data
        allocate(image(Ntotal))   ! this is local one

        ! read digital image
        !Open(200,file='Processed_2D_Berea.dat',status='OLD')
        !    do j=1,Ny
        !        read(200, *) (array2D(i,j), i=1,Nx)
        !    enddo
        !Close(200)
        array2D=0 !NOTE, for debug
        !for debug
        do j = 1, 10
            array2D(:,j) = 1
        end do
        do j = 110, 120
            array2D(:,j) = 1
        end do

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
        real_porosity=real(countVoidP)/real(Nx*Ny)

        ! set local image
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
        nWall=0 ! count the wall points
        Do j=byl,byu
            Do i=bxl,bxu
                localid = (j-ylg)*Nxtotal + i-xlg+1
                If (array2D(i,j)==solid) then
                    NneighborFluid=1 !(1-2-3-4 in D2Q9 corresponding to 2-3-5-7)
                    If (array2D(i+1, j)==fluid)  NneighborFluid=NneighborFluid*2   !Check neighbor on the East(2)
                    If (array2D(i, j+1)==fluid)  NneighborFluid=NneighborFluid*3   !Check neighbor on the North(3)
                    If (array2D(i-1, j)==fluid)  NneighborFluid=NneighborFluid*5   !Check neighbor on the West(5)
                    If (array2D(i, j-1)==fluid)  NneighborFluid=NneighborFluid*7   !Check neighbor on the South(7)

                    SELECT Case (NneighborFluid)
                        CASE (2)
                            image(localid)=WallXp
                            nWall=nWall+1
                        CASE (3)
                            image(localid)=WallYp
                            nWall=nWall+1
                        CASE (5)
                            image(localid)=WallXn
                            nWall=nWall+1
                        CASE (7)
                            image(localid)=WallYn
                            nWall=nWall+1
                        CASE (6)
                            image(localid)=WallXpYp
                            nWall=nWall+1
                        CASE (15)
                            image(localid)=WallXnYp
                            nWall=nWall+1
                        CASE (35)
                            image(localid)=WallXnYn
                            nWall=nWall+1
                        CASE (14)
                            image(localid)=WallXpYn
                            nWall=nWall+1
                    END SELECT
                Endif
            Enddo
        Enddo


        ! Create wall-type vectors
        !vecWall(walli) mark the global id of the walli'th wall in the image 
        ALLOCATE(vecWall(nWall))
        nWall=0
        Do j=ylg, yug
            Do i=xlg, xug
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

        !Direction 1
        !bound
        bxl = xl
        bxu = xug
        byl = yl
        byu = yug
        if(xl == xmin) bxl = xl + 1 !if most west block(processor)
        if(yl == ymin) byl = yl + 1 !if most south block
        if(xu == xmax) bxu = xu     !if most east block(processor)
        if(yu == ymax) byu = yu     !if most north block   
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
        bxl = xlg
        bxu = xu
        byl = yl
        byu = yug
        if(xu == xmax) bxu = xu - 1 !if most east block
        if(yl == ymin) byl = yl + 1 !if most south block
        if(xl == xmin) bxl = xl     !if most west block
        if(yu == ymax) byu = yu     !if most north block
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
        bxl = xlg
        bxu = xu
        byl = ylg
        byu = yu
        if(xu == xmax) bxu = xu - 1 !if most east block
        if(yu == ymax) byu = yu - 1 !if most north block
        if(xl == xmin) bxl = xl     !if most west block
        if(yl == ymin) byl = yl     !if most south block
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
        bxu = xug
        byl = ylg
        byu = yu
        if(xl == xmin) bxl = xl + 1 !if most west block
        if(yu == ymax) byu = yu - 1 !if most north block
        if(xu == xmax) bxu = xu     !if most east block
        if(yl == ymin) byl = yl     !if most south block
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
