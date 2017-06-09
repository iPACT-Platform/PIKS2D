program main
    implicit none
    integer, parameter :: NY = 401
    integer, parameter :: NX = 801
    double precision, parameter :: porosity = 0.75d0
    double precision, parameter :: PI = 4.d0*datan(1.d0)
    integer, parameter :: obstR = ceiling((NX-1)*dsqrt((1.d0-porosity)/PI))     !circular cylinder
    double precision, parameter :: obstX = (NX+1)/2.d0
    double precision, parameter :: obstY = 1.d0
    integer, dimension(NX,NY) :: flag
    integer :: openstatus
    integer :: i,j
    double precision :: l2
    integer :: solidCount

    flag = 0
    solidCount = 0
    print*, "Radius : ", obstR
    do j=1,NY
        do i=1,NX
            l2 = (j-obstY)**2 + (i-obstX)**2
            if (l2 < (obstR+1.d-6)**2 ) then
                flag(i,j) = 1
                solidCount = solidCount + 1
            endif
        enddo
    enddo
    print*, "solidCount : ", solidCount
    print*, "Solid region ratio : " , dble(solidCount)/(NX*NY)

    open(unit=11, file = 'singleCylinder.dat', status='replace', action='write', iostat=openstatus)

    if(openstatus /= 0) then
        print *, 'Could not open flag file'
        stop
    endif

    do j=1, NY
        write(11,'(801I2)') (flag(i,j), i=1,NX)
    enddo
    close(11)

end program
