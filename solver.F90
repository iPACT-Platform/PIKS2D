module solver

use flow
use velocityGrid
use physicalGrid
use mpiParams

implicit none

double precision, parameter :: eps=1.d-10
integer, parameter :: maxStep = 100000000
integer, parameter :: interval = 1000
integer :: iStep
double precision :: error

contains
    subroutine iterate
        call sweep
        call BCwall
        call BCinletOutlet
        call BCsymmetry
        call updateMacro
    end subroutine iterate

    subroutine sweep
        use MPI
    	implicit none
        integer :: k, l, i, j
        INTEGER :: MPI_ERR
        INTEGER :: MPI_REQ_X(4), MPI_REQ_Y(4)
        INTEGER :: MPI_STAT(MPI_STATUS_SIZE,4)
        integer :: xsize
        integer :: ysize
        double precision :: feq

        xsize =  Nytotal*Nc/2*ghostLayers
        ysize = Nxtotal*Nc/2*ghostLayers
        ! Start Recieving
        CALL MPI_IRECV( f1_east_rcv, xsize, MPI_DOUBLE_PRECISION, east,  TAG1, &
                        MPI_COMM_VGRID, MPI_REQ_X(1), MPI_ERR )
        CALL MPI_IRECV( f1_west_rcv, xsize, MPI_DOUBLE_PRECISION, west,  TAG2, &
                        MPI_COMM_VGRID, MPI_REQ_X(2), MPI_ERR )
        CALL MPI_IRECV( f1_north_rcv, ysize, MPI_DOUBLE_PRECISION, north,  TAG3, &
                        MPI_COMM_VGRID, MPI_REQ_Y(1), MPI_ERR )
        CALL MPI_IRECV( f1_south_rcv, ysize, MPI_DOUBLE_PRECISION, south,  TAG4, &
                        MPI_COMM_VGRID, MPI_REQ_Y(2), MPI_ERR )     

        ! Start Sending
        CALL MPI_ISEND( f1_west_snd, xsize, MPI_DOUBLE_PRECISION, west, TAG1, &
                        MPI_COMM_VGRID, MPI_REQ_X(3), MPI_ERR )
        CALL MPI_ISEND( f1_east_snd, xsize, MPI_DOUBLE_PRECISION, east, TAG2, &
                        MPI_COMM_VGRID, MPI_REQ_X(4), MPI_ERR )
        CALL MPI_ISEND( f1_south_snd, ysize, MPI_DOUBLE_PRECISION, south, TAG3, &
                        MPI_COMM_VGRID, MPI_REQ_Y(3), MPI_ERR )
        CALL MPI_ISEND( f1_north_snd, ysize, MPI_DOUBLE_PRECISION, north, TAG4, &
                        MPI_COMM_VGRID, MPI_REQ_Y(4), MPI_ERR )

        ! sweep direction 1
        Do l=1,Nc/4
            Do i=1,Nstencil1
                k=dir1(i)
                ! Switch from reflected in x-direction (default) to in y-direction
                ! only for group of velocity overlapped by two reflected direction
                If (image(k-Nxtotal)==WallXpYp) then
                    f1(k-Nxtotal,l)=f1(k-Nxtotal,l+Nc/2)
                End if
    
                fEq=w(l)*(Rho(k)+2.d0*(cx(l)*Ux(k)+cy(l)*Uy(k)))
    
                f1(k,l)=(mu*(fEq-0.5d0*f1(k,l)) &           
                &        + cx(l)*coefI(i,2)*f1(k-1,l) &
                &        + cx(l)*coefI(i,3)*f1(k-2,l) &
                &        + cy(l)*coefI(i,5)*f1(k-Nxtotal,l) &
                &        + cy(l)*coefI(i,6)*f1(k-2*Nxtotal,l) &
                & )/(0.5d0*mu+cx(l)*coefI(i,1)+cy(l)*coefI(i,4))   
            End do
   	    End do

        ! sweep direction 2
        Do l=Nc/4+1,Nc/2
            Do i=1,Nstencil2
                k=dir2(i)
                ! Switch from reflected in x-direction (default) to in y-direction
                ! only for group of velocity overlapped by two reflected direction
                If (image(k-Nxtotal)==WallXnYp) then
                    f1(k-Nxtotal,l)=f1(k-Nxtotal,l+Nc/2)
                End if
    
                fEq=w(l)*(Rho(k)+2.d0*(cx(l)*Ux(k)+cy(l)*Uy(k)))
                f1(k,l)=(mu*(fEq-0.5d0*f1(k,l)) &
!                f1(k,l)=(mu*(fEq) &    
!                f1(k,l)=(mu*(fEq-f1(k,l)) &        
                &        + cx(l)*coefII(i,2)*f1(k+1,l) &
                &        + cx(l)*coefII(i,3)*f1(k+2,l) &
                &        + cy(l)*coefII(i,5)*f1(k-Nxtotal,l) &
                &        + cy(l)*coefII(i,6)*f1(k-2*Nxtotal,l) &
                & )/(0.5d0*mu+cx(l)*coefII(i,1)+cy(l)*coefII(i,4))
!                & )/(mu+cx(l)*coefII(i,1)+cy(l)*coefII(i,4))   
!                & )/(cx(l)*coefII(i,1)+cy(l)*coefII(i,4))          
            End do
        End do

        ! sweep direction 3
        Do l=Nc/2+1,Nc*3/4
            Do i=1,Nstencil3
                k=dir3(i)
                ! Switch from reflected in x-direction (default) to in y-direction
                ! only for group of velocity overlapped by two reflected direction
                If (image(k+Nxtotal)==WallXnYn) then
                    f1(k+Nxtotal,l)=f1(k+Nxtotal,l-Nc/2)
                End if

                fEq=w(l)*(Rho(k)+2.d0*(cx(l)*Ux(k)+cy(l)*Uy(k)))
                f1(k,l)=(mu*(fEq-0.5d0*f1(k,l)) &
!                f1(k,l)=(mu*(fEq) &    
!                f1(k,l)=(mu*(fEq-f1(k,l)) &        
                &        + cx(l)*coefIII(i,2)*f1(k+1,l) &
                &        + cx(l)*coefIII(i,3)*f1(k+2,l) &
                &        + cy(l)*coefIII(i,5)*f1(k+Nxtotal,l) &
                &        + cy(l)*coefIII(i,6)*f1(k+2*Nxtotal,l) &
                & )/(0.5d0*mu+cx(l)*coefIII(i,1)+cy(l)*coefIII(i,4))
!                & )/(mu+cx(l)*coefIII(i,1)+cy(l)*coefIII(i,4)) 
!                & )/(cx(l)*coefIII(i,1)+cy(l)*coefIII(i,4))            
            End do
        End do

        ! sweep direction 4
        Do l=Nc*3/4+1,Nc
            Do i=1,Nstencil4
                k=dir4(i)
                ! Switch from reflected in x-direction (default) to in y-direction
                ! only for group of velocity overlapped by two reflected direction
                If (image(k+Nxtotal)==WallXpYn) then
                    f1(k+Nxtotal,l)=f1(k+Nxtotal,l-Nc/2)
                End if
    
                fEq=w(l)*(Rho(k)+2.d0*(cx(l)*Ux(k)+cy(l)*Uy(k)))
                f1(k,l)=(mu*(fEq-0.5d0*f1(k,l)) &
!                f1(k,l)=(mu*(fEq) &    
!                f1(k,l)=(mu*(fEq-f1(k,l)) &        
                &        + cx(l)*coefIV(i,2)*f1(k-1,l) &
                &        + cx(l)*coefIV(i,3)*f1(k-2,l) &
                &        + cy(l)*coefIV(i,5)*f1(k+Nxtotal,l) &
                &        + cy(l)*coefIV(i,6)*f1(k+2*Nxtotal,l) &
                & )/(0.5d0*mu+cx(l)*coefIV(i,1)+cy(l)*coefIV(i,4))
!                & )/(mu+cx(l)*coefIV(i,1)+cy(l)*coefIV(i,4))   
!                & )/(cx(l)*coefIV(i,1)+cy(l)*coefIV(i,4))  
            End do
        End do

        ! Wait until send and recv done
        CALL MPI_WAITALL(4, MPI_REQ_X, MPI_STAT, MPI_ERR)
        CALL MPI_WAITALL(4, MPI_REQ_Y, MPI_STAT, MPI_ERR)

        ! pack&unpack west&east buffer
        do j = 1, Nytotal
        	do i = 1, ghostLayers
        	    do l = Nc/4+1, Nc*3/4 ! dir 2 and 3
                    f1_west_snd((j-1)*ghostLayers*Nc/2 + (i-1)*Nc/2 + l-Nc/4) &
                    &  = f1((j-1)*Nxtotal + i, l)
                    f1((j-1)*Nxtotal + i+Nxsub+ghostLayers, l) = &
                    &   f1_east_rcv((j-1)*ghostLayers*Nc/2 + (i-1)*Nc/2 + l-Nc/4)
                enddo
                do l = 1, Nc/4      !dir 1
                    f1_east_snd((j-1)*ghostLayers*Nc/2 + (i-1)*Nc/2 + l) &
                    &  = f1((j-1)*Nxtotal + i+Nxsub+ghostLayers, l)
                    f1((j-1)*Nxtotal + i, l) = &
                    &   f1_west_rcv((j-1)*ghostLayers*Nc/2 + (i-1)*Nc/2 + l)                   
                enddo
                do l = Nc*3/4+1, Nc ! dir 4
                    f1_east_snd((j-1)*ghostLayers*Nc/2 + (i-1)*Nc/2 + l-Nc/2) &
                    &  = f1((j-1)*Nxtotal + i+Nxsub+ghostLayers, l)
                    f1((j-1)*Nxtotal + i, l) = &
                    &   f1_west_rcv((j-1)*ghostLayers*Nc/2 + (i-1)*Nc/2 + l-Nc/2)                   
                enddo
            enddo
        enddo
           
        ! pack&unpack south&north buffer
        do j = 1, ghostLayers
            do i = 1, Nxtotal
                do l = 1, Nc/2 ! dir 1, 2
                    f1_north_snd((j-1)*Nxtotal*Nc/2 + (i-1)*Nc/2 + l) &
                    &  = f1((Nysub+ghostLayers+j-1)*Nxtotal + i, l)
                    f1((j-1)*Nxtotal + i, l) = &
                    &   f1_south_rcv((j-1)*Nxtotal*Nc/2 + (i-1)*Nc/2 + l)
                enddo
                do l = Nc/2+1, Nc ! dir 3, 4
                    f1_south_snd((j-1)*Nxtotal*Nc/2 + (i-1)*Nc/2 + l-Nc/2) &
                    &  = f1((j-1)*Nxtotal + i, l)
                    f1((Nysub+ghostLayers+j-1)*Nxtotal + i,l) = &
                    &   f1_north_rcv((j-1)*Nxtotal*Nc/2 + (i-1)*Nc/2 + l-Nc/2)
                enddo
            enddo
        enddo
    end subroutine sweep

    subroutine BCwall
        implicit none
        integer :: i, j, k, l
        double precision :: RhoWall

        Do i=1,nWall
            k=vecWall(i)
            RhoWall=0.d0
            SELECT CASE (image(k))
                CASE (WallXp)
                    Do l=Nc/4+1,3*Nc/4
                        f1(k,l)=2.d0*f1(k+1,l)-f1(k+2,l)
                        RhoWall=RhoWall-cx(l)*f1(k,l)
                    Enddo
                    RhoWall=RhoWall/DiffFlux
                    Do l=1,Nc/4
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeX(l))
                    Enddo
                    Do l=3*Nc/4+1,Nc
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeX(l))
                    Enddo
                CASE (WallXn)
                    Do l=1,Nc/4
                        f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                        RhoWall=RhoWall+cx(l)*f1(k,l)
                    Enddo
                    Do l=3*Nc/4+1,Nc
                        f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                        RhoWall=RhoWall+cx(l)*f1(k,l)
                    Enddo
                    RhoWall=RhoWall/DiffFlux
                    Do l=Nc/4+1,3*Nc/4
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeX(l))
                    Enddo
                CASE (WallYp)
                    Do l=Nc/2+1,Nc
                        f1(k,l)=2.d0*f1(k+Nxtotal,l)-f1(k+2*Nxtotal,l)
                        RhoWall=RhoWall-cy(l)*f1(k,l)
                    Enddo
                    RhoWall=RhoWall/DiffFlux
                    Do l=1,Nc/2
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeY(l))
                    Enddo
                CASE (WallYn)
                    Do l=1,Nc/2
                        f1(k,l)=2.d0*f1(k-Nxtotal,l)-f1(k-2*Nxtotal,l)
                        RhoWall=RhoWall+cy(l)*f1(k,l)
                    Enddo
                    RhoWall=RhoWall/DiffFlux
                    Do l=Nc/2+1,Nc
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeY(l))
                    Enddo
                !=======================================================================
                !     Boundary condition on the corner wall
                !=======================================================================
                CASE (WallXpYp)
                    !Calculate default reflected f1 (x-direction)
                    Do l=Nc/4+1,3*Nc/4
                        f1(k,l)=2.d0*f1(k+1,l)-f1(k+2,l)
                        RhoWall=RhoWall-cx(l)*f1(k,l)
                    Enddo
                    RhoWall=RhoWall/DiffFlux
                    Do l=1,Nc/4
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeX(l))
                    Enddo
                    Do l=3*Nc/4+1,Nc
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeX(l))
                    Enddo
                    !Calculate and store RhoWallY
                    RhoWall=0.d0
                    Do l=Nc/2+1,Nc
                        RhoWall=RhoWall-cy(l)*(2.d0*f1(k+Nxtotal,l)-f1(k+2*Nxtotal,l))
                    Enddo
                    RhoWall=RhoWall/DiffFlux
!                    RhoWallY(k)=RhoWall
                    !Calculate default reflected f1 (y-direction)
                    Do l=Nc/4+1,Nc/2
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*(2.d0*f1(k+Nxtotal,oppositeY(l))-f1(k+2*Nxtotal,oppositeY(l)))
                    End do
                    !Calculate and store the secondary reflected f1 (y-direction)               
                    Do l=1,Nc/4
                        f1(k,l+Nc/2)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*(2.d0*f1(k+Nxtotal,oppositeY(l))-f1(k+2*Nxtotal,oppositeY(l)))
                    Enddo               
                CASE (WallXnYp)
                    !Calculate default reflected f1 (x-direction)
                    Do l=1,Nc/4
                        f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                        RhoWall=RhoWall+cx(l)*f1(k,l)
                    Enddo
                    Do l=3*Nc/4+1,Nc
                        f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                        RhoWall=RhoWall+cx(l)*f1(k,l)
                    Enddo
                    RhoWall=RhoWall/DiffFlux
                    Do l=Nc/4+1,3*Nc/4
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeX(l))
                    Enddo
                    !Calculate and store RhoWallY
                    RhoWall=0.d0
                    Do l=Nc/2+1,Nc
                        RhoWall=RhoWall-cy(l)*(2.d0*f1(k+Nxtotal,l)-f1(k+2*Nxtotal,l))
                    Enddo
                    RhoWall=RhoWall/DiffFlux
!                    RhoWallY(k)=RhoWall
                    !Calculate default reflected f1 (y-direction)
                    Do l=1,Nc/4
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*(2.d0*f1(k+Nxtotal,oppositeY(l))-f1(k+2*Nxtotal,oppositeY(l)))
                    End do
                    !Calculate and store the secondary reflected f1 (y-direction)               
                    Do l=Nc/4+1,Nc/2
                        f1(k,l+Nc/2)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*(2.d0*f1(k+Nxtotal,oppositeY(l))-f1(k+2*Nxtotal,oppositeY(l)))
                    Enddo               
                CASE (WallXnYn)
                    !Calculate default reflected f1 (x-direction)
                    Do l=1,Nc/4
                        f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                        RhoWall=RhoWall+cx(l)*f1(k,l)
                    Enddo
                    Do l=3*Nc/4+1,Nc
                        f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                        RhoWall=RhoWall+cx(l)*f1(k,l)
                    Enddo
                    RhoWall=RhoWall/DiffFlux
                    Do l=Nc/4+1,3*Nc/4
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeX(l))
                    Enddo
                    !Calculate and store RhoWallY
                    RhoWall=0.d0
                    Do l=1,Nc/2
                        RhoWall=RhoWall+cy(l)*(2.d0*f1(k-Nxtotal,l)-f1(k-2*Nxtotal,l))
                    Enddo
                    RhoWall=RhoWall/DiffFlux
!                    RhoWallY(k)=RhoWall
                    !Calculate default reflected f1 (y-direction)
                    Do l=3*Nc/4+1,Nc
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*(2.d0*f1(k-Nxtotal,oppositeY(l))-f1(k-2*Nxtotal,oppositeY(l)))
                    Enddo
                    !Calculate and store the secondary reflected f1 (y-direction)               
                    Do l=Nc/2+1,3*Nc/4
                        f1(k,l-Nc/2)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*(2.d0*f1(k-Nxtotal,oppositeY(l))-f1(k-2*Nxtotal,oppositeY(l)))
                    Enddo               
                CASE (WallXpYn)
                    !Calculate default reflected f1 (x-direction)
                    Do l=Nc/4+1,3*Nc/4
                        f1(k,l)=2.d0*f1(k+1,l)-f1(k+2,l)
                        RhoWall=RhoWall-cx(l)*f1(k,l)
                    Enddo
                    RhoWall=RhoWall/DiffFlux
                    Do l=1,Nc/4
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeX(l))
                    Enddo
                    Do l=3*Nc/4+1,Nc
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*f1(k,oppositeX(l))
                    Enddo
                    !Calculate and store RhoWallY
                    RhoWall=0.d0
                    Do l=1,Nc/2
                        RhoWall=RhoWall+cy(l)*(2.d0*f1(k-Nxtotal,l)-f1(k-2*Nxtotal,l))
                    Enddo
                    RhoWall=RhoWall/DiffFlux
!                    RhoWallY(k)=RhoWall
                    !Calculate default reflected f1 (y-direction)
                    Do l=Nc/2+1,3*Nc/4
                        f1(k,l)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*(2.d0*f1(k-Nxtotal,oppositeY(l))-f1(k-2*Nxtotal,oppositeY(l)))
                    Enddo
                    !Calculate and store the secondary reflected f1 (y-direction)               
                    Do l=3*Nc/4+1,Nc
                        f1(k,l-Nc/2)=accom*w(l)*RhoWall &
                        &       + (1.d0-accom)*(2.d0*f1(k-Nxtotal,oppositeY(l))-f1(k-2*Nxtotal,oppositeY(l)))
                    Enddo                               
            END SELECT
        End do
    end subroutine

    subroutine BCinletOutlet
        implicit none
        integer :: i, j, k, l
        if(xl == xmin) then ! inlet block (west most processor)
            Do j = yl, yu
                i = xl
                k = (j-ylg)*Nxtotal + i-xlg+1
                !inlet
                Do l=1,Nc/4
                    f1(k,l)=f1(k-1,l)+w(l)*PressDrop !lhzhu, need other block's info, so vgrid is perodical in x dir
                Enddo   
                Do l=3*Nc/4+1,Nc
                    f1(k,l)=f1(k-1,l)+w(l)*PressDrop
                Enddo
            End do
        endif

        if(xu == xmax) then ! outlet block (east most processor)
            Do j = yl, yu
                i = xu
                k = (j-ylg)*Nxtotal + i-xlg+1
                !outlet
                Do l=Nc/4+1,3*Nc/4          
                    f1(k,l)=f1(k+1,l)-w(l)*PressDrop
                Enddo
            Enddo            
        endif
    end subroutine BCinletOutlet
    !
    subroutine BCsymmetry
        implicit none
        integer :: i, j, k, l
        Do i=xl, xu
            !bottom plane
            k = ghostLayers*Nxtotal + i-xlg+1
            Do l=1,Nc/2
                 f1(k,l)=f1(k,oppositeY(l))
            End do
            !top plane  
            k = (Nysub+ghostLayers-1)*Nxtotal + i-xlg+1
            Do l=Nc/2+1,Nc
                 f1(k,l)=f1(k,oppositeY(l))
            End do
        End do
    end subroutine BCsymmetry

    subroutine updateMacro
    implicit none
    integer :: k, l, j
    	Do k=1,Ntotal
            Rho(k)=0.d0
            Ux(k)=0.d0
            Uy(k)=0.d0
            Do l=1,Nc
                Rho(k)=Rho(k)+f1(k,l)
                Ux(k)=Ux(k)+cx(l)*f1(k,l)
                Uy(k)=Uy(k)+cy(l)*f1(k,l)
            End do
        End do
    end subroutine updateMacro
    
    subroutine chkConverge
        use MPI
        implicit none

        ! local vars
        double precision :: permeability
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
        endif

        if(yl == ymin) then !only south most processors
            massSouth = 0.5d0*ds*Ux(ghostLayers+column + (yl-ylg)*Nxtotal)
        endif

        if (yu == ymax) then !only north most processors
            massNorth = 0.5d0*ds*Ux(ghostLayers+column + (ghostLayers+Nysub)*Nxtotal)
        endif
        
        massLocal = (massInner + massSouth + massNorth) * 2.d0 / PressDrop

        ! reduction
        call MPI_REDUCE(massLocal, mass2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        master, MPI_COMM_VGRID, MPI_ERR)

        if (proc == master) then 
            error=dabs(1.d0-mass2/mass)/(interval)
            mass=mass2
            permeability=mass*Kn*sqrt(4.d0/PI)*(Nx-1)/(Ny-1)/2    
            write(*,"( 1I10, 3ES15.6)")  iStep,  mass,  permeability, error
            open(22,file='Results.dat', position="append")
            write(22,'(4ES15.6, 1I15)') Kn, mass, permeability, error, iStep
            close(22)
        endif
    end subroutine chkConverge

    subroutine saveFlowField
        integer :: j, i, k
        character(13) fname
        write(fname, '(A, I0.3, A)') 'Field_', proc, '.dat'
        open(20,file=fname,STATUS="REPLACE")
        write(20,*) ' TITLE=" Field"'
        write(20,*) ' VARIABLES=x,y,Rho,Ux,Uy'
        write(20,*) ' ZONE T=final, I=', Nxsub,', J=', Nysub,', F=POINT'

        Do j=yl,yu
            Do i=xl,xu
                k= (j-ylg)*Nxtotal + i-xlg+1
                If (image(k)==fluid) then
                    write(20,'(2I10,3ES15.6)') i, j, Rho(k)+1.d0, Ux(k), Uy(k)
                else
                    write(20,'(2I10,3ES15.6)') i, j, 0.d0, 0.d0, 0.d0
                Endif   
            Enddo
        Enddo
        close(20)
    end subroutine saveFlowField

    subroutine memFree
        deallocate (array2D)
        deallocate (image)
        deallocate (vecWall, dir1, dir2, dir3, dir4)
        deallocate (coefI, coefII, coefIII, coefIV)
        deallocate (f1, Rho, Ux, Uy)
        deallocate (f1_west_snd, f1_west_rcv, f1_east_snd, f1_east_rcv)
        deallocate (f1_north_snd, f1_north_rcv, f1_south_snd, f1_south_rcv)
    end subroutine memFree
end module solver