!=======================================================================
!> @brief Physical space configurations
!=======================================================================
module velocityGrid
use gaussHermite

implicit none

!Number of fundamental molecular velocity, to be read from NML: velocityNml
integer :: Nc_fundamental

!Number of moleculer velocity in 2D-Gaussian Hermite
integer :: Nc
!abscissae and weighting Hermite quadrature, dimension(Nc_fundamental)
double precision, dimension (:), allocatable :: xi, weight1D 
!molecular velocity and weighting, dimension(Nc)
double precision, dimension (:), allocatable :: cx, cy, w
!specular wall's normal vector in X, Y direction, dimension(Nc)
integer, dimension (:), allocatable :: oppositeX, oppositeY
!constant PI
double precision, parameter :: PI=datan(1.d0)*4.d0

!half range flux of discrete velocity grid 
double precision :: DiffFlux

contains
    subroutine setupVelocityGrid
        implicit none
        integer :: l, m, n, i

        ! Nc_fundamental has been initialized from the NML
        allocate(xi(Nc_fundamental))
        allocate(weight1D(Nc_fundamental))

        if ( Nc_fundamental == 2) then
            xi = xi2
            weight1D = wi2
        else if (Nc_fundamental == 10) then
            xi = xi10
            weight1D = wi10
        else
            print*, "Nc_fundamental xi/wi for ", Nc_fundamental, &
             "has not been provided"
            deallocate(xi)
            deallocate(weight1D)
            stop
        endif

        Nc=(2*Nc_fundamental)**2

        allocate(cx(Nc), cy(Nc), w(Nc))
        allocate(oppositeX(Nc), oppositeY(Nc))
        
        Do i=1,4    ! index for molecular velocity group I,II,III,IV
             Do m=1,Nc_fundamental
                Do n=1,Nc_fundamental
                    l=n+(m-1)*Nc_fundamental+(i-1)*Nc_fundamental**2
                    SELECT CASE (i)
                        CASE (1)
                            cx(l) = xi(n)
                            cy(l) = xi(m)
                            oppositeY(l)=l+3*Nc/4
                            oppositeX(l)=l+Nc/4
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

        DiffFlux=0.d0
        Do l=1,Nc/2
            DiffFlux=DiffFlux+cy(l)*w(l)
        Enddo

    end subroutine setupVelocityGrid

end module velocityGrid
