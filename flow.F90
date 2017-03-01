!----------------------------------------------------------------------
!> @brief Flow field array and physical simualtion parameters
!----------------------------------------------------------------------
module flow
use velocityGrid, only: PI
implicit none
save

double precision, parameter :: Kn = 9d-5
double precision, parameter :: mu = dsqrt(PI)/2.0d0/Kn
double precision, parameter :: PressDrop=1.0d-3
double precision, parameter :: accom = 1.d0

double precision, DIMENSION(:,:), ALLOCATABLE :: f1
double precision, DIMENSION(:), ALLOCATABLE :: Rho, Ux, Uy
double precision :: mass

contains
    subroutine setupFlow     
        use physicalGrid, only: Ntotal, Nxtotal, Nytotal, ghostLayers
        use velocityGrid, only: Nc
        use mpiParams, only: f1_west_snd, f1_east_snd, f1_west_rcv, f1_east_rcv, &
                             f1_south_snd, f1_north_snd, f1_south_rcv, f1_north_rcv
        implicit none

        ALLOCATE(f1(Ntotal,Nc))
        ALLOCATE(Rho(Ntotal), Ux(Ntotal), Uy(Ntotal))

        f1=0.d0
        Rho = 0.d0
        Ux = 0.d0
        Uy = 0.d0

        mass = 1.d0

        ! Buffer for exchange data only half domain of the velocity
        allocate( f1_west_snd(Nc/2*Nytotal*ghostLayers) ) 
        allocate( f1_west_rcv(Nc/2*Nytotal*ghostLayers) ) 
        allocate( f1_east_snd(Nc/2*Nytotal*ghostLayers) ) 
        allocate( f1_east_rcv(Nc/2*Nytotal*ghostLayers) ) 

        allocate( f1_north_snd(Nc/2*Nxtotal*ghostLayers) ) 
        allocate( f1_north_rcv(Nc/2*Nxtotal*ghostLayers) )
        allocate( f1_south_snd(Nc/2*Nxtotal*ghostLayers) ) 
        allocate( f1_south_rcv(Nc/2*Nxtotal*ghostLayers) )

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
