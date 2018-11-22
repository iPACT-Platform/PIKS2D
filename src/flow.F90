!----------------------------------------------------------------------
!> @brief Flow field array and physical simualtion parameters
!----------------------------------------------------------------------
module flow
use velocityGrid, only: PI
implicit none
save

! flow parameters, to be read from NML: flowNml
character(len=200) :: allKnStr
double precision :: pressDrop, accom
double precision :: mu


double precision, dimension(:), allocatable :: allKn
double precision :: Kn
integer :: nKn
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
