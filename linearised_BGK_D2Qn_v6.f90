! The ghost flag has to be redefined
!
!  This program calculates the pressure driven gas flow throught a
!  two-dimensional channel with obstacle using the linearised BGK kinetic equation
!  fix-point method
!  Modified to adapt tomography data
!   Correct specular reflection on the corner point. Error due to not store the incoming f for the second-direction
!       ---------------------------------
!       ---------------------------------
!       ||i                           o||                   y
!       ||n                           u||                   /\
!       ||-                           t||                   |
!       ||l                           l||                   |
!       ||e                           e||                   |
!       ||t                           t||                   |
!       ---------------------------------                   |----------->x
!       ---------------------------------
!
!
!
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
PROGRAM linearised_BGK_D2Qn
use omp_lib
implicit none
!***********************************************************************
!   Constant and variable declarations
!***********************************************************************

double precision, parameter :: accom=1.d0 ! accommodation coefficient
double precision, parameter :: porosity=0.75d0
double precision, parameter :: eps=1.d-10
double precision :: real_porosity
!=======================================================================
!       Constant and indexed, temporary variables
!=======================================================================
double precision,parameter :: PI=datan(1.d0)*4.d0
!double precision,parameter :: PI=3.1416d0
integer :: i,j,k    !index variable for physical space
integer :: l,m,n  !index variable for molecular velocity space
logical :: fileexist
!=======================================================================
!       Physical space configurations
!=======================================================================
!integer, parameter :: Ny = 401 !number of point half height channel, including two end
!integer, parameter :: Nx = (Ny-1)*2+1  !length of the channel
!integer, parameter :: Ny=42
!integer, parameter :: Nx=100
!integer, parameter :: Ny=856
!integer, parameter :: Nx=1066
integer, parameter :: Ny=(428-1)*2+1
integer, parameter :: Nx=(533-1)*2+1

integer, Dimension (Nx,Ny) ::array2D
integer, parameter :: ghostLayer=2 !number of ghost layer at each boundary
!integer, parameter :: obstR=(Nx-1)/4     !square cylinder with porosity=0.75
!integer, parameter :: obstR=0     !square cylinder with porosity=0.75
integer, parameter :: obstR = ceiling((Nx-1)*dsqrt((1.d0-porosity)/PI))     !circular cylinder
!integer, parameter :: obstR = ceiling( dsqrt(1-porosity)/2.0*(Nx-1) )       !square cylinder
integer, parameter :: obstX = ghostLayer+Ny-2
integer, parameter :: obstY = ghostLayer-1
integer, parameter :: Nxmin=1-ghostLayer,Nxmax=Nx+ghostLayer
integer, parameter :: Nymin=1-ghostLayer,Nymax=Ny+ghostLayer
integer, parameter :: Nxtotal=Nxmax-Nxmin+1
integer, parameter :: Nytotal=Nymax-Nymin+1
integer, parameter :: Ntotal=Nxtotal*Nytotal

integer, parameter:: fluid = 0, solid = 1, ghost=4
integer, parameter:: WallXpYp = 20, WallYp = 21, WallXnYp = 22, WallXn = 23
integer, parameter:: WallXnYn = 24, WallYn = 25, WallXpYn = 26, WallXp = 27
integer :: NneighborSolid, NneighborFluid, Nstencil, icount, iteration1
integer :: Nstencil1, Nstencil2, Nstencil3, Nstencil4
double precision :: RhoWall
!=======================================================================
!       Flow configuration configurations
!=======================================================================

integer, parameter :: column=2 ! layer to extract flow rate
double precision,parameter :: PressDrop=1.0d-3 !Pressure drop
integer, parameter :: iteration_max=100,iteration_min=200000,interval=1000
!double precision :: ds = 1.d0/(Nx-1) !uniform grid spacing in physical space
double precision :: ds = 0.5d0/(Ny-1) !uniform grid spacing in physical space
!double precision :: dt  !uniform grid spacing in  time space
double precision :: Kn,mass,mass2,permeability
double precision :: DiffFlux, error, fEq, mu

integer, parameter :: nKn= 1!28 number of Kn case
!integer, parameter :: nKn=6! number of Kn case
integer :: iKn
double precision, parameter, DIMENSION(1:nKn):: seriesKn=(/1.0d-1/)
!=======================================================================
!       Molecular velocity space configurations
!=======================================================================
integer, parameter :: Nc_fundamental=2 !Number of fundamental molecular velocity
integer, parameter :: Nc=(2*Nc_fundamental)**2 !Number of moleculer velocity in 2D-Gaussian Hermite
integer, parameter :: Vmax=5
integer, parameter :: power_law=3
!integer, parameter :: Nphi=100 !Number of discrete point in angular of polar coordinate mode(Nphi,4)=0
!integer, parameter :: Nc=Nc_fundamental*Nphi !Number of molecular velocity in polar coordinate
double precision, DIMENSION (1:Nc_fundamental) :: xi,weight1D !abscissae and weighting Hermite quadrature
double precision, DIMENSION (1:Nc) :: cx, cy, w !molecular velocity and weighting

integer, DIMENSION (1:Nc) :: oppositeX, oppositeY! specular wall's normal vector in X, Y direction

!double precision :: phi
!=======================================================================
!       Microscopic and macroscopic parameters
!=======================================================================

integer :: nWall
integer, DIMENSION(:), ALLOCATABLE :: vecWall
integer, DIMENSION(:), ALLOCATABLE :: dir1, dir2, dir3, dir4 ! dir1(i) is the ith 
double precision, DIMENSION(:,:), ALLOCATABLE :: coefI, coefII, coefIII, coefIV

integer*4:: startday(3), starttime(3),endday(3), endtime(3) !Real time

!integer(kind=int64) :: wall_time
integer :: wall_time

integer, DIMENSION(:), ALLOCATABLE :: image
double precision, DIMENSION(:,:), ALLOCATABLE :: f1
double precision, DIMENSION(:), ALLOCATABLE :: Rho, Ux, Uy
ALLOCATE(image(Ntotal))
ALLOCATE(f1(Ntotal,Nc))
ALLOCATE(Rho(Ntotal), Ux(Ntotal), Uy(Ntotal))

!Gaussian Hermite 4th order (Half-range)
!xi(1) = 0.3001939310608394d0         !fundamental abscissae
!xi(2) = 0.1252421045333717d1
!weight1D(1) = 0.6405291796843786d0/dsqrt(PI)    !fundamental weighting
!weight1D(2) = 0.2456977457683793d0/dsqrt(PI)
xi(1) = dsqrt(3.d0-dsqrt(6.d0)) /dsqrt(2.d0)         !fundamental abscissae
xi(2) = dsqrt(3.d0+dsqrt(6.d0)) /dsqrt(2.d0)
weight1D(1) = (3.d0+dsqrt(6.d0))/12.d0    !fundamental weighting
weight1D(2) = (3.d0-dsqrt(6.d0))/12.d0

!------------------------------------------------------------------------
!           Two-dimensional Hermite quadrature: tensor production formulae
!           Mapping from  2D array to 1D array
!------------------------------------------------------------------------
Do i=1,4    ! index for molecular velocity group I,II,III,IV
    Do m=1,Nc_fundamental
        Do n=1,Nc_fundamental
            l=n+(m-1)*Nc_fundamental+(i-1)*Nc_fundamental**2
            SELECT CASE (i)
                CASE (1)
                    cx(l) = xi(n)
                    cy(l) = xi(m)
                    oppositeY(l)=l+3*Nc/4   ! specular wall's normal vector in Y direction
                    oppositeX(l)=l+Nc/4     ! specular wall's normal vector in X direction
                CASE (2)
                    cx(l) = -xi(n)
                    cy(l) = xi(m)
                    oppositeY(l)=l+Nc/4
                    oppositeX(l)=l-Nc/4
                CASE (3)
                    cx(l) = -xi(n)
                    cy(l) = -xi(m)
                    oppositeY(l)=l-Nc/4
                    oppositeX(l)=l+Nc/4
                CASE (4)
                    cx(l) = xi(n)
                    cy(l) = -xi(m)
                    oppositeY(l)=l-3*Nc/4
                    oppositeX(l)=l-Nc/4
            END SELECT
            w(l) = weight1D(n)*weight1D(m)
        End do
    End do
End do

!***********************************************************************
!   Initialisation
!***********************************************************************
    inquire(file="Results.dat",exist=fileexist)
    if (.not.fileexist) then
        open(22,file="Results.dat",STATUS='NEW')
        write(22,*) ' TITLE=" Results"'
        write(22,*) ' VARIABLES=Kn,mass,permeability,error,iteration,wall_time'
        write(22,*) ' ZONE T="final"'
        close(22)
    end if

!=======================================================================
!       Physical image of flow+solid+wall domain
!=======================================================================
image = fluid
!------------------------------------------------------------------------
!           Prefill ghosts layer by solid
!------------------------------------------------------------------------
If (ghostLayer>0) then
    Do j=Nymin,Nymax
        Do i=Nxmin,Nxmax
            If ((i<1).OR.(i>Nx).OR.(j<1).OR.(j>Ny)) then
                k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
                image(k) = ghost
            End if
        Enddo
    Enddo
End if
!------------------------------------------------------------------------
!           Read from tomography file
!------------------------------------------------------------------------
!Open(200,file='Processed_2D_Tomography.dat',status='OLD')
!Open(200,file='Processed_2D_Berea.dat',status='OLD')
Open(200,file='Processed_2x_2D_Berea.dat',status='OLD')
!Open(200,file='cylinder.dat',status='OLD')
    Do j=1,Ny
        read(200, *) (array2D(i,j), i=1,Nx)
    Enddo
Close(200)
!array2D = 0
!do j = 1, 2
!    array2D(:,j) = 1
!end do
!do j = Ny-1, Ny
!    array2D(:,j) = 1
!end do

Do j=1, Ny
    Do i=1,Nx
        If (array2D(i,j)==0) then 
            array2D(i,j)=fluid
        End if
        If (array2D(i,j)==1) then 
            array2D(i,j)=solid
        End if
        k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
        image(k)= array2D(i,j)
    Enddo
Enddo

!------------------------------------------------------------------------
!           Calculating real porosity  
!------------------------------------------------------------------------
l=0
Do j=1,Ny
    Do i=1,Nx
        k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
        If (image(k)==solid) then  
            l=l+1   
        endif   
    Enddo
Enddo   

real_porosity=1.d0-real(l)/real(Nx*Ny)
!------------------------------------------------------------------------
!           Assign automatically the wall layer between fluid and solid 
!           What does it do?
!             Identify the wall node types based ont the node's neighbour node types (solid/fluid) ?
!------------------------------------------------------------------------
nWall=0
Do j=1,Ny
    Do i=1,Nx
        k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
        If (image(k)==solid) then
            NneighborFluid=1 !(1-2-3-4 in D2Q9 corresponding to 2-3-5-7)
            If (image(k+1)==fluid)  NneighborFluid=NneighborFluid*2   !Check neighbor on the East(2)
            If (image(k+Nxtotal)==fluid)  NneighborFluid=NneighborFluid*3   !Check neighbor on the North(3)
            If (image(k-1)==fluid)  NneighborFluid=NneighborFluid*5   !Check neighbor on the West(5)
            If (image(k-Nxtotal)==fluid)  NneighborFluid=NneighborFluid*7   !Check neighbor on the South(7)

            SELECT Case (NneighborFluid)
                CASE (2)
                    image(k)=WallXp
                    nWall=nWall+1
                CASE (3)
                    image(k)=WallYp
                    nWall=nWall+1
                CASE (5)
                    image(k)=WallXn
                    nWall=nWall+1
                CASE (7)
                    image(k)=WallYn
                    nWall=nWall+1
                CASE (6)
                    image(k)=WallXpYp
                    nWall=nWall+1
                CASE (15)
                    image(k)=WallXnYp
                    nWall=nWall+1
                CASE (35)
                    image(k)=WallXnYn
                    nWall=nWall+1
                CASE (14)
                    image(k)=WallXpYn
                    nWall=nWall+1
            END SELECT
        Endif
    Enddo
Enddo

!------------------------------------------------------------------------
!           Create wall-type vectors
!           vecWall(walli) mark the global id of the walli'th wall in the image 
!------------------------------------------------------------------------
ALLOCATE(vecWall(nWall))
nWall=0
Do j=1,Ny
    Do i=1,Nx
        k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
        SELECT Case (image(k))
            CASE (WallXp)
                nWall=nWall+1
                vecWall(nWall)=k
            CASE (WallYp)
                nWall=nWall+1
                vecWall(nWall)=k
            CASE (WallXn)
                nWall=nWall+1
                vecWall(nWall)=k
            CASE (WallYn)
                nWall=nWall+1
                vecWall(nWall)=k
            CASE (WallXpYp)
                nWall=nWall+1
                vecWall(nWall)=k
            CASE (WallXnYp)
                nWall=nWall+1
                vecWall(nWall)=k
            CASE (WallXnYn)
                nWall=nWall+1
                vecWall(nWall)=k
            CASE (WallXpYn)
                nWall=nWall+1
                vecWall(nWall)=k
        END SELECT
    Enddo
Enddo

!------------------------------------------------------------------------
!           Assign the ghost layers
!------------------------------------------------------------------------
If (ghostLayer>0) then
    Do j=Nymin,Nymax
        Do i=Nxmin,Nxmax
            If ((i<1).OR.(i>Nx).OR.(j<1).OR.(j>Ny)) then
                k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
                image(k) = ghost
            End if
        Enddo
    Enddo
End if

    open(11,file='Map.dat',STATUS="REPLACE")
        write(11,*) 'Map with porosity =', real_porosity
!        write(11,*) 'Nxtotal',Nxtotal
        Do j=Nymax,Nymin,-1
            Do i=Nxmin,Nxmax
                k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
                if (image(k)==fluid) then
                    write(11,"(A)",Advance='no') "   "
!                   write(11,"(I2,A)",Advance='no') image(k)," "
                else
                    write(11,"(I2,A)",Advance='no') image(k)," "
                end if
            End do
            write(11,*)
        Enddo
    close(11)
!=======================================================================
!       Create advancing-point vectors
!=======================================================================

icount=0
Do j=2,Ny
    Do i=2,Nx
        k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
        If (image(k)==fluid) then
            icount=icount+1
        End if
    End do
End do
Nstencil1=icount !Number of fluid points in the numerical stencil

icount=0
Do j=2,Ny
    Do i=Nx-1,1,-1
        k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
        If (image(k)==fluid) then
            icount=icount+1
        End if
    End do
End do
Nstencil2=icount !Number of fluid points in the numerical stencil

icount=0
Do j=Ny-1,1,-1
    Do i=Nx-1,1,-1
        k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
        If (image(k)==fluid) then
            icount=icount+1
        End if
    End do
End do
Nstencil3=icount !Number of fluid points in the numerical stencil

icount=0
Do j=Ny-1,1,-1
    Do i=2,Nx
        k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
        If (image(k)==fluid) then
            icount=icount+1
        End if
    End do
End do
Nstencil4=icount !Number of fluid points in the numerical stencil

ALLOCATE(dir1(Nstencil1),dir2(Nstencil2),dir3(Nstencil3),dir4(Nstencil4))


icount=0
Do j=2,Ny
    Do i=2,Nx
        k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
        If (image(k)==fluid) then
            icount=icount+1
            dir1(icount)=k
        End if
    End do
End do

icount=0
Do j=2,Ny
    Do i=Nx-1,1,-1
        k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
        If (image(k)==fluid) then
            icount=icount+1
            dir2(icount)=k
        End if
    End do
End do

icount=0
Do j=Ny-1,1,-1
    Do i=Nx-1,1,-1
        k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
        If (image(k)==fluid) then
            icount=icount+1
            dir3(icount)=k
        End if
    End do
End do

icount=0
Do j=Ny-1,1,-1
    Do i=2,Nx
        k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
        If (image(k)==fluid) then
            icount=icount+1
            dir4(icount)=k
        End if
    End do
End do

!=======================================================================
!       Create matrix of coefficients used in numerical stencil
!       coefI(point, weight_x0n & weight_x1n & weight_x2n, weight_y0n & weight_y1n & weight_y2n)
!       coefII(point, weight_x0p & weight_x1p & weight_x2p, weight_y0n & weight_y1n & weight_y2n)
!       coefIII(point, weight_x0p & weight_x1p & weight_x2p, weight_y0p & weight_y1p & weight_y2p)
!       coefIV(point, weight_x0n & weight_x1n & weight_x2n, weight_y0p & weight_y1p & weight_y2p)
!       coefI(i,1)-coefI(i,2)-coefI(i,3)=0, coefI(i,4)-coefI(i,5)-coefI(i,6)=0
!=======================================================================
ALLOCATE(coefI(Nstencil1,6),coefII(Nstencil2,6),coefIII(Nstencil3,6),coefIV(Nstencil4,6))
Do i=1,Nstencil1
    k=dir1(i) 
    !2nd order of accuracy
    coefI(i,1)=1.5d0/ds !x0n
    coefI(i,2)=2.d0/ds  !x1n
    coefI(i,3)=-0.5d0/ds   !x2n
    coefI(i,4)=1.5d0/ds !y0n
    coefI(i,5)=2.d0/ds  !y1n
    coefI(i,6)=-0.5d0/ds   !y2n
    if ((image(k-2)==ghost).OR.(image(k-1)==WallXp).OR.(image(k-1)==WallXpYn) &
    & .OR.(image(k-1)==WallXpYp)) then !1st order of accuracy in x
        coefI(i,1)=1.d0/ds  !x0n
        coefI(i,2)=1.d0/ds  !x1n
        coefI(i,3)=0.d0     !x2n
    end if
    if ((image(k-2*Nxtotal)==ghost).OR.(image(k-Nxtotal)==WallYp) &
    & .OR.(image(k-Nxtotal)==WallXpYp).OR.(image(k-Nxtotal)==WallXnYp)) then   !1st order of accuracy in y
        coefI(i,4)=1.d0/ds  !y0n
        coefI(i,5)=1.d0/ds  !y1n
        coefI(i,6)=0.d0     !y2n
    end if
End do

Do i=1,Nstencil2
    k=dir2(i)
    !2nd order of accuracy
    coefII(i,1)=-1.5d0/ds !x0n
    coefII(i,2)=-2.d0/ds  !x1n
    coefII(i,3)=0.5d0/ds   !x2n
    coefII(i,4)=1.5d0/ds !y0n
    coefII(i,5)=2.d0/ds  !y1n
    coefII(i,6)=-0.5d0/ds   !y2n
    if ((image(k+2)==ghost).OR.(image(k+1)==WallXn).OR.(image(k+1)==WallXnYp) &
    & .OR.(image(k+1)==WallXnYn)) then !1st order of accuracy in x
        coefII(i,1)=-1.d0/ds  !x0n
        coefII(i,2)=-1.d0/ds  !x1n
        coefII(i,3)=0.d0     !x2n
    end if
    if ((image(k-2*Nxtotal)==ghost).OR.(image(k-Nxtotal)==WallYp) &
    & .OR.(image(k-Nxtotal)==WallXpYp).OR.(image(k-Nxtotal)==WallXnYp)) then   !1st order of accuracy in y
        coefII(i,4)=1.d0/ds  !y0n
        coefII(i,5)=1.d0/ds  !y1n
        coefII(i,6)=0.d0     !y2n
    end if
End do

Do i=1,Nstencil3
    k=dir3(i)
    !2nd order of accuracy
    coefIII(i,1)=-1.5d0/ds !x0n
    coefIII(i,2)=-2.d0/ds  !x1n
    coefIII(i,3)=0.5d0/ds   !x2n
    coefIII(i,4)=-1.5d0/ds !y0n
    coefIII(i,5)=-2.d0/ds  !y1n
    coefIII(i,6)=0.5d0/ds   !y2n
    if ((image(k+2)==ghost).OR.(image(k+1)==WallXn).OR.(image(k+1)==WallXnYp) &
    & .OR.(image(k+1)==WallXnYn)) then !1st order of accuracy in x
        coefIII(i,1)=-1.d0/ds  !x0n
        coefIII(i,2)=-1.d0/ds  !x1n
        coefIII(i,3)=0.d0     !x2n
    end if
    if ((image(k+2*Nxtotal)==ghost).OR.(image(k+Nxtotal)==WallYn) &
    & .OR.(image(k+Nxtotal)==WallXnYn).OR.(image(k+Nxtotal)==WallXpYn)) then   !1st order of accuracy in y
        coefIII(i,4)=-1.d0/ds  !y0n
        coefIII(i,5)=-1.d0/ds  !y1n
        coefIII(i,6)=0.d0     !y2n
    end if
End do

Do i=1,Nstencil4
    k=dir4(i)
    !2nd order of accuracy
    coefIV(i,1)=1.5d0/ds !x0n
    coefIV(i,2)=2.d0/ds  !x1n
    coefIV(i,3)=-0.5d0/ds   !x2n
    coefIV(i,4)=-1.5d0/ds !y0n
    coefIV(i,5)=-2.d0/ds  !y1n
    coefIV(i,6)=0.5d0/ds   !y2n
    if ((image(k-2)==ghost).OR.(image(k-1)==WallXp).OR.(image(k-1)==WallXpYn) &
    & .OR.(image(k-1)==WallXpYp)) then !1st order of accuracy in x
        coefIV(i,1)=1.d0/ds  !x0n
        coefIV(i,2)=1.d0/ds  !x1n
        coefIV(i,3)=0.d0     !x2n
    end if
    if ((image(k+2*Nxtotal)==ghost).OR.(image(k+Nxtotal)==WallYn) &
    & .OR.(image(k+Nxtotal)==WallXnYn).OR.(image(k+Nxtotal)==WallXpYn)) then   !1st order of accuracy in y
        coefIV(i,4)=-1.d0/ds  !y0n
        coefIV(i,5)=-1.d0/ds  !y1n
        coefIV(i,6)=0.d0     !y2n
    end if
End do
!=======================================================================
!     Constant diffusive reflected flux at wall boundary
!=======================================================================
    DiffFlux=0.d0
    Do l=1,Nc/2
        DiffFlux=DiffFlux+cy(l)*w(l)
    Enddo
    
!***********************************************************************
!   Print out the input parameter
!***********************************************************************
    open(10,file='Test.dat',STATUS="REPLACE")
        write(10,*) 'Velocity space D2Q16: l, cx, cy, w, oppositeX, oppositeY'
        Do l=1,Nc
            write(10,*) l, cx(l), cy(l), w(l), oppositeX(l), oppositeY(l)
        Enddo
    close(10)   
!=======================================================================
!       Initial condition
!=======================================================================
f1=0.d0
! Do j=1,Ny
     ! Do i=1,Nx
         ! k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
         ! if (image(k)/=solid) then
             ! Rho(k)=PressDrop*(i/2.d0-Nx)/Nx
             ! RhoWallY(k)=Rho(k)
         ! end if
         ! f1(k,:)=w(:)*Rho(k) !Check
     ! Enddo
! Enddo


Rho=0.d0 ! For linearised fomula
!RhoWallY=Rho
Ux=0.d0
Uy=0.d0

Do iKn=1,nKn
!Do iKn=10,10

Kn=seriesKn(iKn)
write(*,*) 'Kn = ', Kn
mu=dsqrt(PI)/2.0d0/Kn 
error=1.d0
mass=1.d0
iteration1=0


!***********************************************************************
!   Main procedure
!***********************************************************************
!Do WHILE ((error>eps).OR.(iteration1<iteration_min))
!Do WHILE ((iteration1<iteration_min))
    call itime(starttime)  ! starttime(1)=hour, (2)=minute, (3)=second

Do WHILE ((error>eps))
!=======================================================================
!     Reset summational variables
!=======================================================================

!=======================================================================
!     Collision and streaming
!=======================================================================
!------------------------------------------------------------------------
!           In the 1st group of direction cx>0 & cy>0
!------------------------------------------------------------------------
!$OMP PARALLEL & 
!$OMP DEFAULT(NONE) &
!$OMP PRIVATE(i,j,k,l) &
!$OMP PRIVATE(fEq,RhoWall) &


!!$OMP PRIVATE() &
!$OMP SHARED(nWall,vecWall) &
!$OMP SHARED(cx,cy,w,mu) &
!$OMP SHARED(Nstencil1,Nstencil2,Nstencil3,Nstencil4) &
!$OMP SHARED(dir1,dir2,dir3,dir4) &
!$OMP SHARED(coefI,coefII,coefIII,coefIV) &
!$OMP SHARED(image) &
!$OMP SHARED(f1) &
!$OMP SHARED(oppositeX,oppositeY) &
!$OMP SHARED(Rho,Ux,Uy,Kn,ds,DiffFlux) 

!!$OMP SHARED() &
!!$OMP SHARED()

!PRINT*, "Before sweep Ftest1 =", f1((ghostLayer+10)*Nxtotal + ghostLayer+1+1,1)
!$OMP DO SCHEDULE(STATIC) 
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
!$OMP END DO NOWAIT 
!PRINT*, "After  sweep Ftest1 =", f1((ghostLayer+10)*Nxtotal + ghostLayer+1+1,1)
!------------------------------------------------------------------------
!           In the 2nd group of direction cx<0 & cy>0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC) 
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
!            f1(k,l)=(mu*(fEq) &    
!            f1(k,l)=(mu*(fEq-f1(k,l)) &        
            &        + cx(l)*coefII(i,2)*f1(k+1,l) &
            &        + cx(l)*coefII(i,3)*f1(k+2,l) &
            &        + cy(l)*coefII(i,5)*f1(k-Nxtotal,l) &
            &        + cy(l)*coefII(i,6)*f1(k-2*Nxtotal,l) &
            & )/(0.5d0*mu+cx(l)*coefII(i,1)+cy(l)*coefII(i,4))
!            & )/(mu+cx(l)*coefII(i,1)+cy(l)*coefII(i,4))   
!            & )/(cx(l)*coefII(i,1)+cy(l)*coefII(i,4))          
        End do
    End do
!$OMP END DO NOWAIT 
!------------------------------------------------------------------------
!           In the 3rd group of direction cx<0 & cy<0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
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
!            f1(k,l)=(mu*(fEq) &    
!            f1(k,l)=(mu*(fEq-f1(k,l)) &        
            &        + cx(l)*coefIII(i,2)*f1(k+1,l) &
            &        + cx(l)*coefIII(i,3)*f1(k+2,l) &
            &        + cy(l)*coefIII(i,5)*f1(k+Nxtotal,l) &
            &        + cy(l)*coefIII(i,6)*f1(k+2*Nxtotal,l) &
            & )/(0.5d0*mu+cx(l)*coefIII(i,1)+cy(l)*coefIII(i,4))
!            & )/(mu+cx(l)*coefIII(i,1)+cy(l)*coefIII(i,4)) 
!            & )/(cx(l)*coefIII(i,1)+cy(l)*coefIII(i,4))            
        End do
    End do
!$OMP END DO NOWAIT 
!------------------------------------------------------------------------
!           In the 4th group of direction cx>0 & cy<0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
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
!            f1(k,l)=(mu*(fEq) &    
!            f1(k,l)=(mu*(fEq-f1(k,l)) &        
            &        + cx(l)*coefIV(i,2)*f1(k-1,l) &
            &        + cx(l)*coefIV(i,3)*f1(k-2,l) &
            &        + cy(l)*coefIV(i,5)*f1(k+Nxtotal,l) &
            &        + cy(l)*coefIV(i,6)*f1(k+2*Nxtotal,l) &
            & )/(0.5d0*mu+cx(l)*coefIV(i,1)+cy(l)*coefIV(i,4))
!            & )/(mu+cx(l)*coefIV(i,1)+cy(l)*coefIV(i,4))   
!            & )/(cx(l)*coefIV(i,1)+cy(l)*coefIV(i,4))  
        End do
    End do
!$OMP END DO    
!=======================================================================
!     Boundary condition on the flat wall
!=======================================================================
!$OMP DO SCHEDULE(STATIC)
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
!                RhoWallY(k)=RhoWall
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
!                RhoWallY(k)=RhoWall
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
!                RhoWallY(k)=RhoWall
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
!                RhoWallY(k)=RhoWall
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
!$OMP END DO  



!=======================================================================
!     Boundary condition on inlet & outlet
!=======================================================================
!$OMP DO SCHEDULE(STATIC) 
!   i=1
!    Do j=1,Ny
!       k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
!        !inlet
!        f1(k,1:Nc/4)=f1(k+Nx-1,1:Nc/4)+w(1:Nc/4)*PressDrop
!        f1(k,3*Nc/4+1:Nc)=f1(k+Nx-1,3*Nc/4+1:Nc)+w(3*Nc/4+1:Nc)*PressDrop
!        !outlet
!        f1(k+Nx-1,Nc/4+1:3*Nc/4)=f1(k,Nc/4+1:3*Nc/4)-w(Nc/4+1:3*Nc/4)*PressDrop
!    End do

    Do j=1,Ny
        i=1 
        k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
        !inlet
        Do l=1,Nc/4
            f1(k,l)=f1(k+Nx-1,l)+w(l)*PressDrop
        Enddo   
        Do l=3*Nc/4+1,Nc
            f1(k,l)=f1(k+Nx-1,l)+w(l)*PressDrop
        Enddo
            !outlet
        Do l=Nc/4+1,3*Nc/4          
            f1(k+Nx-1,l)=f1(k,l)-w(l)*PressDrop
        End do
    End do  
!$OMP END DO   

!PRINT*, "Ftest1 =", f1((ghostLayer+10)*Nxtotal + ghostLayer+1,1)
!PRINT*, "Ftest5 =", f1(ghostLayer*Nxtotal + ghostLayer+1,5)
!PRINT*, "Ftest9 =", f1(ghostLayer*Nxtotal + ghostLayer+1,9)
!PRINT*, "Ftest13=", f1((ghostLayer+10)*Nxtotal + ghostLayer+1,13)

!=======================================================================
!     Boundary condition on symmetric planes
!=======================================================================
!$OMP DO SCHEDULE(STATIC) 
!   j=0
!     Do i=1,Nx
!         !bottom plane
!       k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
!         Do l=1,Nc/2
!             f1(k,l)=f1(k+2*Nxtotal,oppositeY(l))
!         End do
!         !top plane    
!         Do l=Nc/2+1,Nc
!             f1(k+(Ny+1)*Nxtotal,l)=f1(k+(Ny-1)*Nxtotal,oppositeY(l))
!       End do
!     End do

    Do i=1,Nx
        !bottom plane
        j=1 
        k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
 !       Do l=Nc/2+1,Nc
 !            f1(k,l)=2.d0*f1(k+Nxtotal,l)-f1(k+2*Nxtotal,l)
 !      End do      
        Do l=1,Nc/2
             f1(k,l)=f1(k,oppositeY(l))
        End do
        !top plane  
        j=Ny 
        k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal 
!        Do l=1,Nc/2
!             f1(k,l)=2.d0*f1(k-Nxtotal,l)-f1(k-2*Nxtotal,l)
!       End do  
        Do l=Nc/2+1,Nc
             f1(k,l)=f1(k,oppositeY(l))
        End do
    End do
!$OMP END DO 
!=======================================================================
!     Macroscopic parameter evaluation
!=======================================================================
!------------------------------------------------------------------------
!           Rho, Ux, Uy
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC) 
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
!$OMP END DO 
!$OMP END PARALLEL  
!Debug

!PRINT*, Rho(ghostLayer*Nxtotal + ghostLayer+1)
!PRINT*, "cx16 =", cx(5)
!PRINT*, "cx16 =", cy(5)

!------------------------------------------------------------------------
!           Flow rate
!------------------------------------------------------------------------
    iteration1=iteration1+1
    If (mod(iteration1,interval)==0) then
        mass2=0.d0
        Do j=2,Ny-1
            k=(column+ghostLayer)+(j+ghostLayer-1)*Nxtotal
            mass2=mass2+Ux(k)*ds
        Enddo
        mass2=(mass2+0.5d0*ds*(Ux(column+ghostLayer+ghostLayer*Nxtotal) &
        &      +Ux(column+ghostLayer+(ghostLayer+Ny-1)*Nxtotal)))*2.d0/PressDrop
        error=dabs(1.d0-mass2/mass)/(interval)
        mass=mass2
        permeability=mass*seriesKn(iKn)*sqrt(4.d0/pi)*(Nx-1)/(Ny-1)/2

!------------------------------------------------------------------------
!           At the wall nodes
!------------------------------------------------------------------------
        write(*,"( 1I10, 3ES15.6)")  iteration1,  mass,  permeability, error
    endif
Enddo
!***********************************************************************
!   Print out results
!***********************************************************************
    call itime(endtime)     ! endtime(1)=hour, (2)=minute, (3)=second
    wall_time = (endtime(1)-starttime(1))*3600 + (endtime(2)-starttime(2))*60 + (endtime(3)-starttime(3))   
    
    open(22,file='Results.dat', position="append")
        write(22,'(4ES15.6, 1I15, 1I20)') Kn, mass, permeability,error,iteration1,wall_time 
    close(22)
    
    
    open(20,file='Field.dat',STATUS="REPLACE")
        write(20,*) ' TITLE=" Field"'
        write(20,*) ' VARIABLES=x,y,Rho,Ux,Uy'
        write(20,*)'ZONE T=final, I=',Nx,', J=',Ny,', F=POINT'
        Do j=1,Ny
            Do i=1,Nx
                k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
                
                If (image(k)==fluid) then
                    write(20,'(2I10,3ES15.6)') i, j, Rho(k)+1.d0, Ux(k), Uy(k)
                else
                    write(20,'(2I10,3ES15.6)') i, j, 0.d0, 0.d0, 0.d0
                Endif   
            Enddo
        Enddo
    close(20)

    ! open(10,file='Test.dat', STATUS='OLD', position="append")
        ! write(10,*) 'Flow rate = ', mass
    ! close(10)



Enddo

DEALLOCATE(vecWall)
DEALLOCATE(dir1, dir2, dir3, dir4)
DEALLOCATE(coefI, coefII, coefIII, coefIV)
DEALLOCATE(f1)
DEALLOCATE(Rho, Ux, Uy)
DEALLOCATE(image)

END PROGRAM linearised_BGK_D2Qn
