!=======================================================================
!> @brief Physical space configurations
!=======================================================================
module velocityGrid
integer, parameter :: Nc_fundamental=3 !Number of fundamental molecular velocity
integer, parameter :: Nc=(2*Nc_fundamental)**2 !Number of moleculer velocity in 2D-Gaussian Hermite
integer, parameter :: Vmax=5
integer, parameter :: power_law=3
!integer, parameter :: Nphi=100 !Number of discrete point in angular of polar coordinate mode(Nphi,4)=0
!integer, parameter :: Nc=Nc_fundamental*Nphi !Number of molecular velocity in polar coordinate
double precision, DIMENSION (1:Nc_fundamental) :: xi, weight1D !abscissae and weighting Hermite quadrature
double precision, DIMENSION (1:Nc) :: cx, cy, w !molecular velocity and weighting
integer, DIMENSION (1:Nc) :: oppositeX, oppositeY! specular wall's normal vector in X, Y direction
double precision, parameter :: PI=datan(1.d0)*4.d0
double precision :: DiffFlux

contains
    subroutine setupVelocityGrid
        integer :: l, m, n
        !Gaussian Hermite 4th order (Half-range)
        xi(1) = 0.3001939310608394d0         !fundamental abscissae
        xi(2) = 0.1252421045333717d1
        weight1D(1) = 0.6405291796843786d0/dsqrt(PI)    !fundamental weighting
        weight1D(2) = 0.2456977457683793d0/dsqrt(PI)

        !xi(1) = 0.4360774119276165086792d0         !fundamental abscissae
        !xi(2) = 1.335849074013696949715d0
        !xi(3) = 2.350604973674492222834d0 
        !weight1D(1) = 0.724629595224392524092d0/dsqrt(PI)    !fundamental weighting
        !weight1D(2) = 0.1570673203228566439163d0/dsqrt(PI)
        !weight1D(3) = 0.00453000990550884564086d0/dsqrt(PI)

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

        DiffFlux=0.d0
        Do l=1,Nc/2
            DiffFlux=DiffFlux+cy(l)*w(l)
        Enddo

        !open(10,file='Test.dat',STATUS="REPLACE")
        !write(10,*) 'Velocity space D2Q16: l, cx, cy, w, oppositeX, oppositeY'
        !Do l=1,Nc
        !    write(10,*) l, cx(l), cy(l), w(l), oppositeX(l), oppositeY(l)
        !Enddo
        !close(10)
    end subroutine setupVelocityGrid

end module velocityGrid
