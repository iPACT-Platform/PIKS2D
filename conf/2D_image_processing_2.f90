! Read the input binary-image 
! Remove the improper solid points for the intended method of simulation
! Export the binary-image formatted output file 
! Scale up and export the scale-up binary-image 
! Scale keep the ratio of cell (area) while the 1st version keep the ratio of node
Program D2Qn_image_processing_2
implicit none
!=======================================================================
!       Physical space configurations
!=======================================================================
integer, parameter :: Ny=1500
integer, parameter :: Nx=3000

integer, parameter :: ghostLayer=1 !number of ghost layer at each boundary
integer, parameter :: Nxmin=1-ghostLayer,Nxmax=Nx+ghostLayer
integer, parameter :: Nymin=1-ghostLayer,Nymax=Ny+ghostLayer
integer, parameter :: Nxtotal=Nxmax-Nxmin+1
integer, parameter :: Nytotal=Nymax-Nymin+1
integer, parameter :: Ntotal=Nxtotal*Nytotal

integer, parameter:: fluid = 0, solid = 1, inlet = 30, outlet = 31, ghost=4
integer, parameter:: wallNE = 20, wallN = 21, wallNW = 22, wallW = 23
integer, parameter:: wallSW = 24, wallS = 25, wallSE = 26, wallE = 27
integer, DIMENSION(:,:), ALLOCATABLE :: array2D,array2D_2x,array2D_3x,array2D_4x
integer, DIMENSION(:), ALLOCATABLE :: image

integer :: i,j,k    !index variable for physical space
integer :: l,m,n  !index variable for molecular velocity space
integer :: mod_x,mod_y  
integer	:: icount,icount2


character(LEN=30)::input_name,output_name

ALLOCATE(array2D(Nx,Ny))	
ALLOCATE(image(Ntotal))	
ALLOCATE(array2D_2x((Nx-1)*2+1,(Ny-1)*2+1))	
ALLOCATE(array2D_3x((Nx-1)*3+1,(Ny-1)*3+1))
ALLOCATE(array2D_4x((Nx-1)*4+1,(Ny-1)*4+1))		

!=======================================================================
!           Read from tomography file  
!=======================================================================
input_name='var2circle.dat'
output_name='Processed_'//TRIM(ADJUSTL(input_name))
Open(200,file=input_name,status='OLD')
	Do j=1,Ny
		read(200, *) (array2D(i,j), i=1,Nx)
	Enddo
Close(200)

icount=0
Do j=1, Ny
	Do i=1,Nx
		If (array2D(i,j)==0) then 
			array2D(i,j)=fluid
			icount=icount+1
		End if
		If (array2D(i,j)==1) then 
			array2D(i,j)=solid
		End if
		l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
		image(l)= array2D(i,j)
	Enddo
Enddo 
open(10,file='Test.dat', STATUS='REPLACE')
	Write(10,*) 'Porosity of the original binary image (node-base): ', Dble(icount)/Nx/Ny

icount=0
	 Do j=1,(Ny-1)
		Do i=1, (Nx-1)
			If ((array2D(i,j)==fluid).OR.(array2D(i+1,j)==fluid).OR.&
			& (array2D(i,j+1)==fluid).OR.(array2D(i+1,j+1)==fluid)) then
 				icount=icount+1
			Endif				
		Enddo
	 Enddo	 
Write(10,*) 'Porosity of the original binary image (cell-base): ', Dble(icount)/((Nx-1))/((Ny-1))
Write(10,*) 
!=======================================================================
!       Remove improper solid nodes
!	No 1-layer thin of fluid in the pore
!	No 1-layer thin of solid (include the diagonal)in the rock
!=======================================================================
!------------------------------------------------------------------------
!           Fill ghosts layer by ghost
!------------------------------------------------------------------------
Do j=Nymin,Nymax
    Do i=Nxmin,Nxmax
        If ((i<1).OR.(i>Nx).OR.(j<1).OR.(j>Ny)) then
            l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
            image(l) = ghost
        End if
    Enddo
Enddo 
!------------------------------------------------------------------------
!           Drill to obtain at least two layer of fluid in pore
!			Convert 1-layer thin of solid to fluid
!------------------------------------------------------------------------
Do j=1,Ny
    Do i=1,Nx
		l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal

        If (image(l)==solid) then
			If(((image(l-1)/=solid).AND.(image(l+1)/=solid)).OR.((image(l-Nxtotal)/=solid).AND.(image(l+Nxtotal)/=solid))) then
				image(l)=fluid
            End if
        End if
		
        If (image(l)==fluid) then
			If(((image(l-1)/=fluid).AND.(image(l+1)/=fluid)).OR.((image(l-Nxtotal)/=fluid).AND.(image(l+Nxtotal)/=fluid))) then
				If (j==Ny) then
					If (i==Nx) then
						image(l-1)=fluid
						image(l-Nxtotal)=fluid
						image(l-1-Nxtotal)=fluid
					else
						image(l+1)=fluid
						image(l-Nxtotal)=fluid
						image(l+1-Nxtotal)=fluid
					End	If
				else
					If (i==Nx) then
						image(l-1)=fluid
						image(l+Nxtotal)=fluid
						image(l-1+Nxtotal)=fluid
					else
						image(l+1)=fluid
						image(l+Nxtotal)=fluid
						image(l+1+Nxtotal)=fluid
					End	If					
				End if
            End if
        End if		
		
		If((image(l)==solid).AND.(image(l+1+Nxtotal)==solid).AND.(image(l-1-Nxtotal)==solid) &
		& .AND.(image(l-1+Nxtotal)==fluid).AND.(image(l+1-Nxtotal)==fluid)) then
			image(l)=fluid
			image(l-1)=fluid
			image(l+1)=fluid
			image(l-Nxtotal)=fluid
			image(l+Nxtotal)=fluid	
		End If

		If((image(l)==solid).AND.(image(l-1+Nxtotal)==solid).AND.(image(l+1-Nxtotal)==solid) &
		& .AND.(image(l-1-Nxtotal)==fluid).AND.(image(l+1+Nxtotal)==fluid)) then
			image(l)=fluid
			image(l-1)=fluid
			image(l+1)=fluid
			image(l-Nxtotal)=fluid
			image(l+Nxtotal)=fluid	
		End If
		
    Enddo
Enddo 

Do j=1,Ny
    Do i=1,Nx
		l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal

        If (image(l)==solid) then
			If(((image(l-1)/=solid).AND.(image(l+1)/=solid)).OR.((image(l-Nxtotal)/=solid).AND.(image(l+Nxtotal)/=solid))) then
				image(l)=fluid
            End if
        End if
		
        If (image(l)==fluid) then
			If(((image(l-1)/=fluid).AND.(image(l+1)/=fluid)).OR.((image(l-Nxtotal)/=fluid).AND.(image(l+Nxtotal)/=fluid))) then
				If (j==Ny) then
					If (i==Nx) then
						image(l-1)=fluid
						image(l-Nxtotal)=fluid
						image(l-1-Nxtotal)=fluid
					else
						image(l+1)=fluid
						image(l-Nxtotal)=fluid
						image(l+1-Nxtotal)=fluid
					End	If
				else
					If (i==Nx) then
						image(l-1)=fluid
						image(l+Nxtotal)=fluid
						image(l-1+Nxtotal)=fluid
					else
						image(l+1)=fluid
						image(l+Nxtotal)=fluid
						image(l+1+Nxtotal)=fluid
					End	If					
				End if
            End if
        End if		

		If((image(l)==solid).AND.(image(l+1+Nxtotal)==solid).AND.(image(l-1-Nxtotal)==solid) &
		& .AND.(image(l-1+Nxtotal)==fluid).AND.(image(l+1-Nxtotal)==fluid)) then
			image(l)=fluid
			image(l-1)=fluid
			image(l+1)=fluid
			image(l-Nxtotal)=fluid
			image(l+Nxtotal)=fluid	
		End If

		If((image(l)==solid).AND.(image(l-1+Nxtotal)==solid).AND.(image(l+1-Nxtotal)==solid) &
		& .AND.(image(l-1-Nxtotal)==fluid).AND.(image(l+1+Nxtotal)==fluid)) then
			image(l)=fluid
			image(l-1)=fluid
			image(l+1)=fluid
			image(l-Nxtotal)=fluid
			image(l+Nxtotal)=fluid	
		End If
		
    Enddo
Enddo 

Do j=1,Ny
    Do i=1,Nx
		l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal

        If (image(l)==solid) then
			If(((image(l-1)/=solid).AND.(image(l+1)/=solid)).OR.((image(l-Nxtotal)/=solid).AND.(image(l+Nxtotal)/=solid))) then
				image(l)=fluid
            End if
        End if
		
        If (image(l)==fluid) then
			If(((image(l-1)/=fluid).AND.(image(l+1)/=fluid)).OR.((image(l-Nxtotal)/=fluid).AND.(image(l+Nxtotal)/=fluid))) then
				If (j==Ny) then
					If (i==Nx) then
						image(l-1)=fluid
						image(l-Nxtotal)=fluid
						image(l-1-Nxtotal)=fluid
					else
						image(l+1)=fluid
						image(l-Nxtotal)=fluid
						image(l+1-Nxtotal)=fluid
					End	If
				else
					If (i==Nx) then
						image(l-1)=fluid
						image(l+Nxtotal)=fluid
						image(l-1+Nxtotal)=fluid
					else
						image(l+1)=fluid
						image(l+Nxtotal)=fluid
						image(l+1+Nxtotal)=fluid
					End	If					
				End if
            End if
        End if		

		If((image(l)==solid).AND.(image(l+1+Nxtotal)==solid).AND.(image(l-1-Nxtotal)==solid) &
		& .AND.(image(l-1+Nxtotal)==fluid).AND.(image(l+1-Nxtotal)==fluid)) then
			image(l)=fluid
			image(l-1)=fluid
			image(l+1)=fluid
			image(l-Nxtotal)=fluid
			image(l+Nxtotal)=fluid	
		End If

		If((image(l)==solid).AND.(image(l-1+Nxtotal)==solid).AND.(image(l+1-Nxtotal)==solid) &
		& .AND.(image(l-1-Nxtotal)==fluid).AND.(image(l+1+Nxtotal)==fluid)) then
			image(l)=fluid
			image(l-1)=fluid
			image(l+1)=fluid
			image(l-Nxtotal)=fluid
			image(l+Nxtotal)=fluid	
		End If
		
    Enddo
Enddo 

Do j=1,Ny
    Do i=1,Nx
		l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal

        If (image(l)==solid) then
			If(((image(l-1)/=solid).AND.(image(l+1)/=solid)).OR.((image(l-Nxtotal)/=solid).AND.(image(l+Nxtotal)/=solid))) then
				image(l)=fluid
            End if
        End if
		
        If (image(l)==fluid) then
			If(((image(l-1)/=fluid).AND.(image(l+1)/=fluid)).OR.((image(l-Nxtotal)/=fluid).AND.(image(l+Nxtotal)/=fluid))) then
				If (j==Ny) then
					If (i==Nx) then
						image(l-1)=fluid
						image(l-Nxtotal)=fluid
						image(l-1-Nxtotal)=fluid
					else
						image(l+1)=fluid
						image(l-Nxtotal)=fluid
						image(l+1-Nxtotal)=fluid
					End	If
				else
					If (i==Nx) then
						image(l-1)=fluid
						image(l+Nxtotal)=fluid
						image(l-1+Nxtotal)=fluid
					else
						image(l+1)=fluid
						image(l+Nxtotal)=fluid
						image(l+1+Nxtotal)=fluid
					End	If					
				End if
            End if
        End if		

		If((image(l)==solid).AND.(image(l+1+Nxtotal)==solid).AND.(image(l-1-Nxtotal)==solid) &
		& .AND.(image(l-1+Nxtotal)==fluid).AND.(image(l+1-Nxtotal)==fluid)) then
			image(l)=fluid
			image(l-1)=fluid
			image(l+1)=fluid
			image(l-Nxtotal)=fluid
			image(l+Nxtotal)=fluid	
		End If

		If((image(l)==solid).AND.(image(l-1+Nxtotal)==solid).AND.(image(l+1-Nxtotal)==solid) &
		& .AND.(image(l-1-Nxtotal)==fluid).AND.(image(l+1+Nxtotal)==fluid)) then
			image(l)=fluid
			image(l-1)=fluid
			image(l+1)=fluid
			image(l-Nxtotal)=fluid
			image(l+Nxtotal)=fluid	
		End If
		
    Enddo
Enddo 
!=======================================================================
!    Print out the binary image
!=======================================================================
!------------------------------------------------------------------------
!           Map of solid part
!------------------------------------------------------------------------
    open(11,file='Map.dat',STATUS="REPLACE")
        write(11,*) 'Map'
        Do j=Nymax,Nymin,-1
			Do i=Nxmin,Nxmax
				l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
				if (image(l)==fluid) then
					write(11,"(A)",Advance='no') "   "
				else
					write(11,"(I2,A)",Advance='no') image(l)," "
				end if
			End do
			write(11,*)
        Enddo
    close(11)
!------------------------------------------------------------------------
!           Binary image
!------------------------------------------------------------------------
icount=0
Do j=1, Ny
	Do i=1,Nx
		l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
		If (image(l)==fluid) then 
			array2D(i,j)=0
			icount=icount+1
		End if
		If (image(l)==solid) then 
			array2D(i,j)=1
		End if
	Enddo
Enddo 
Write(10,*) 'Number-of-node base: '	
Write(10,*) 'Porosity of the processed 1x binary image: ', Dble(icount)/Nx/Ny

! Open(200,file=output_name,status='New')
	! Do j=1,Ny
		! write(200, *) (array2D(i,j), i=1,Nx)
	! Enddo
! Close(200)

 Open(200,file=output_name,STATUS="REPLACE")
	 Do j=1,Ny
		Do i=1, Nx
			write(200, "(I2)",Advance='no') array2D(i,j)
		Enddo
		write(200,*)
	 Enddo
 Close(200)
!=======================================================================
!   Scale-up by 2 times the binary image
!=======================================================================
array2D_2x=ghost
!------------------------------------------------------------------------
!           1x-based image
!------------------------------------------------------------------------
	 Do j=1,Ny
		Do i=1, Nx
			array2D_2x(1+(i-1)*2,1+(j-1)*2)=array2D(i,j)
		Enddo
	 Enddo
!------------------------------------------------------------------------
!           Filled points between based points
!------------------------------------------------------------------------ 
icount=0
	 Do j=1,(Ny-1)*2+1
		Do i=1, (Nx-1)*2+1
			l=1+(i-1)/2
			m=1+(j-1)/2
			mod_x=mod(i-1,2)
			mod_y=mod(j-1,2)
			If ((mod_x/=0).AND.(mod_y==0)) then
				If ((array2D(l,m)==solid).AND.(array2D(l+1,m)==solid)) then
					array2D_2x(i,j)=solid
				else
					array2D_2x(i,j)=fluid
				end if						
			End if
			
			If ((mod_x==0).AND.(mod_y/=0)) then
				If ((array2D(l,m)==solid).AND.(array2D(l,m+1)==solid)) then
					array2D_2x(i,j)=solid
				else
					array2D_2x(i,j)=fluid
				end if						
			End if	

			If ((mod_x/=0).AND.(mod_y/=0).AND.(mod_x==mod_y)) then
				If (((array2D(l,m)==solid).AND.(array2D(l+1,m+1)==solid))&
				&.AND.& !1st and 2nd version switch by AND/OR
				& ((array2D(l,m+1)==solid).AND.(array2D(l+1,m)==solid))) then 
					array2D_2x(i,j)=solid
				else
					array2D_2x(i,j)=fluid
				end if						
			End if				

			If (array2D_2x(i,j)==fluid) icount=icount+1				
		Enddo
	 Enddo
 	 
Write(10,*) 'Porosity of the processed 2x binary image: ', Dble(icount)/((Nx-1)*2+1)/((Ny-1)*2+1)	 
!------------------------------------------------------------------------
!           Map of solid part
!------------------------------------------------------------------------
    open(11,file='Map_2x.dat',STATUS="REPLACE")
        write(11,*) 'Map'
        Do j=(Ny-1)*2+1,1,-1
			Do i=1,(Nx-1)*2+1
				if (array2D_2x(i,j)==fluid) then
					write(11,"(A)",Advance='no') "   "
				else
					write(11,"(I2,A)",Advance='no') array2D_2x(i,j)," "
				end if
			End do
			write(11,*)
        Enddo
    close(11)
!------------------------------------------------------------------------
!           Binary image
!------------------------------------------------------------------------	
output_name='Processed_2x_'//TRIM(ADJUSTL(input_name))	
	Open(200,file=output_name,STATUS="REPLACE")
		Do j=1,(Ny-1)*2+1
			Do i=1, (Nx-1)*2+1
				write(200, "(I2)",Advance='no') array2D_2x(i,j)
			Enddo
			write(200,*)
		Enddo
	Close(200)
!=======================================================================
!   Scale-up by 3 times the binary image
!=======================================================================
array2D_3x=ghost
!------------------------------------------------------------------------
!           1x-based image
!------------------------------------------------------------------------
	 Do j=1,Ny
		Do i=1, Nx
			array2D_3x(1+(i-1)*3,1+(j-1)*3)=array2D(i,j)
		Enddo
	 Enddo
!------------------------------------------------------------------------
!           Filled points between based points
!------------------------------------------------------------------------ 
icount=0
	 ! Do j=1,(Ny-1)*3+1
		! Do i=1, (Nx-1)*3+1
			! l=1+(i-1)/3
			! m=1+(j-1)/3
			! mod_x=mod(i-1,3)
			! mod_y=mod(j-1,3)
			! If ((mod_x/=0).AND.(mod_y==0)) then
				! If ((array2D(l,m)==solid).AND.(array2D(l+1,m)==solid)) then
					! array2D_3x(i,j)=solid
				! else
					! array2D_3x(i,j)=fluid
				! end if						
			! End if
			
			! If ((mod_x==0).AND.(mod_y/=0)) then
				! If ((array2D(l,m)==solid).AND.(array2D(l,m+1)==solid)) then
					! array2D_3x(i,j)=solid
				! else
					! array2D_3x(i,j)=fluid
				! end if						
			! End if	

			! If ((mod_x/=0).AND.(mod_y/=0).AND.(mod_x==mod_y)) then
				! If ((array2D(l,m)==solid).AND.(array2D(l+1,m+1)==solid)) then
					! array2D_3x(i,j)=solid
				! end if						
			! End if				
			! If ((mod_x/=0).AND.(mod_y/=0).AND.((mod_x+mod_y)==3)) then
				! If ((array2D(l+1,m)==solid).AND.(array2D(l,m+1)==solid)) then
					! array2D_3x(i,j)=solid				
				! end if						
			! End if	
			
		! Enddo
	 ! Enddo
	 
	 ! Do j=1,(Ny-1)*3+1
		! Do i=1, (Nx-1)*3+1
			! l=1+(i-1)/3
			! m=1+(j-1)/3
			! mod_x=mod(i-1,3)
			! mod_y=mod(j-1,3)

			! If ((mod_x/=0).AND.(mod_y/=0).AND.(mod_x/=mod_y).AND.((mod_x+mod_y)/=3)) then
				! icount2=0
				! If (array2D_3x(i-1,j)==solid) icount2=icount2+1
				! If (array2D_3x(i+1,j)==solid) icount2=icount2+1
				! If (array2D_3x(i,j-1)==solid) icount2=icount2+1
				! If (array2D_3x(i,j+1)==solid) icount2=icount2+1
				! If (icount2>2) array2D_3x(i,j)=solid
			! End if				
			
		! Enddo
	 ! Enddo	 
	 
	 ! Do j=1,(Ny-1)*3+1
		! Do i=1, (Nx-1)*3+1
			! l=1+(i-1)/3
			! m=1+(j-1)/3
			! mod_x=mod(i-1,3)
			! mod_y=mod(j-1,3)

			! If (array2D_3x(i,j)==ghost) then
				! If ((array2D_3x(i-1,j)/=fluid).AND.(array2D_3x(i+1,j)/=fluid)&
				! &.OR.&  
				! &(array2D_3x(i,j-1)/=fluid).AND.(array2D_3x(i,j+1)/=fluid)) then
					! array2D_3x(i,j)=solid
				! else
					! array2D_3x(i,j)=fluid
				! Endif	
			! End if				

			! If (array2D_3x(i,j)==fluid) icount=icount+1				
		! Enddo
	 ! Enddo	 

	 Do j=1,(Ny-1)*3+1
		Do i=1, (Nx-1)*3+1
			l=1+(i-1)/3
			m=1+(j-1)/3
			mod_x=mod(i-1,3)
			mod_y=mod(j-1,3)
			If ((mod_x/=0).AND.(mod_y==0)) then
				If ((array2D(l,m)==solid).AND.(array2D(l+1,m)==solid)) then
					array2D_3x(i,j)=solid
				else
					array2D_3x(i,j)=fluid
				end if						
			End if
			
			If ((mod_x==0).AND.(mod_y/=0)) then
				If ((array2D(l,m)==solid).AND.(array2D(l,m+1)==solid)) then
					array2D_3x(i,j)=solid
				else
					array2D_3x(i,j)=fluid
				end if						
			End if	

			If ((mod_x/=0).AND.(mod_y/=0)) then
				If (((array2D(l,m)==solid).AND.(array2D(l+1,m+1)==solid))&
				&.AND.& 
				& ((array2D(l,m+1)==solid).AND.(array2D(l+1,m)==solid))) then 
					array2D_3x(i,j)=solid
				else
					array2D_3x(i,j)=fluid
				end if						
			End if				

			If (array2D_3x(i,j)==fluid) icount=icount+1		
		Enddo
	 Enddo
	 
	 
Write(10,*) 'Porosity of the processed 3x binary image: ', Dble(icount)/((Nx-1)*3+1)/((Ny-1)*3+1)	 
!------------------------------------------------------------------------
!           Map of solid part
!------------------------------------------------------------------------
    open(11,file='Map_3x.dat',STATUS="REPLACE")
        write(11,*) 'Map'
        Do j=(Ny-1)*3+1,1,-1
			Do i=1,(Nx-1)*3+1
				if (array2D_3x(i,j)==fluid) then
					write(11,"(A)",Advance='no') "   "
				else
					write(11,"(I2,A)",Advance='no') array2D_3x(i,j)," "
				end if
			End do
			write(11,*)
        Enddo
    close(11)
!------------------------------------------------------------------------
!           Binary image
!------------------------------------------------------------------------	
output_name='Processed_3x_'//TRIM(ADJUSTL(input_name))	
	Open(200,file=output_name,STATUS="REPLACE")
		Do j=1,(Ny-1)*3+1
			Do i=1, (Nx-1)*3+1
				write(200, "(I2)",Advance='no') array2D_3x(i,j)
			Enddo
			write(200,*)
		Enddo
	Close(200) 
!=======================================================================
!   Scale-up by 4 times the binary image
!=======================================================================
array2D_4x=ghost
!------------------------------------------------------------------------
!           1x-based image
!------------------------------------------------------------------------
	 Do j=1,Ny
		Do i=1, Nx
			array2D_4x(1+(i-1)*4,1+(j-1)*4)=array2D(i,j)
		Enddo
	 Enddo
!------------------------------------------------------------------------
!           Filled points between based points
!------------------------------------------------------------------------ 
icount=0
	 ! Do j=1,(Ny-1)*4+1
		! Do i=1, (Nx-1)*4+1
			! l=1+(i-1)/4
			! m=1+(j-1)/4
			! mod_x=mod(i-1,4)
			! mod_y=mod(j-1,4)
			! If ((mod_x/=0).AND.(mod_y==0)) then
				! If ((array2D(l,m)==solid).AND.(array2D(l+1,m)==solid)) then
					! array2D_4x(i,j)=solid
				! else
					! array2D_4x(i,j)=fluid
				! end if						
			! End if
			
			! If ((mod_x==0).AND.(mod_y/=0)) then
				! If ((array2D(l,m)==solid).AND.(array2D(l,m+1)==solid)) then
					! array2D_4x(i,j)=solid
				! else
					! array2D_4x(i,j)=fluid
				! end if						
			! End if	

			! If ((mod_x/=0).AND.(mod_y/=0).AND.(mod_x==mod_y)) then
				! If ((array2D(l,m)==solid).AND.(array2D(l+1,m+1)==solid)) then
					! array2D_4x(i,j)=solid
				! end if						
			! End if				
			! If ((mod_x/=0).AND.(mod_y/=0).AND.((mod_x+mod_y)==4)) then
				! If ((array2D(l+1,m)==solid).AND.(array2D(l,m+1)==solid)) then
					! array2D_4x(i,j)=solid				
				! end if						
			! End if	
			
		! Enddo
	 ! Enddo
	 
	 ! Do j=1,(Ny-1)*4+1
		! Do i=1, (Nx-1)*4+1
			! l=1+(i-1)/4
			! m=1+(j-1)/4
			! mod_x=mod(i-1,4)
			! mod_y=mod(j-1,4)

			! If ((mod_x/=0).AND.(mod_y/=0).AND.(mod_x/=mod_y).AND.((mod_x+mod_y)/=4)) then
				! icount2=0
				! If (array2D_4x(i-1,j)==solid) icount2=icount2+1
				! If (array2D_4x(i+1,j)==solid) icount2=icount2+1
				! If (array2D_4x(i,j-1)==solid) icount2=icount2+1
				! If (array2D_4x(i,j+1)==solid) icount2=icount2+1
				! If (icount2>2) array2D_4x(i,j)=solid
			! End if				
			
		! Enddo
	 ! Enddo	 
	 
	 ! Do j=1,(Ny-1)*4+1
		! Do i=1, (Nx-1)*4+1
			! l=1+(i-1)/4
			! m=1+(j-1)/4
			! mod_x=mod(i-1,4)
			! mod_y=mod(j-1,4)

			! If (array2D_4x(i,j)==ghost) then
				! If ((array2D_4x(i-1,j)/=fluid).AND.(array2D_4x(i+1,j)/=fluid)&
				! &.AND.&  !1st and 2nd version switch by AND/OR
				! &(array2D_4x(i,j-1)/=fluid).AND.(array2D_4x(i,j+1)/=fluid)) then
					! array2D_4x(i,j)=solid
				! else
					! array2D_4x(i,j)=fluid
				! Endif	
			! End if				

			! If (array2D_4x(i,j)==fluid) icount=icount+1				
		! Enddo
	 ! Enddo	 
	 
	 Do j=1,(Ny-1)*4+1
		Do i=1, (Nx-1)*4+1
			l=1+(i-1)/4
			m=1+(j-1)/4
			mod_x=mod(i-1,4)
			mod_y=mod(j-1,4)
			If ((mod_x/=0).AND.(mod_y==0)) then
				If ((array2D(l,m)==solid).AND.(array2D(l+1,m)==solid)) then
					array2D_4x(i,j)=solid
				else
					array2D_4x(i,j)=fluid
				end if						
			End if
			
			If ((mod_x==0).AND.(mod_y/=0)) then
				If ((array2D(l,m)==solid).AND.(array2D(l,m+1)==solid)) then
					array2D_4x(i,j)=solid
				else
					array2D_4x(i,j)=fluid
				end if						
			End if	

			If ((mod_x/=0).AND.(mod_y/=0)) then
				If (((array2D(l,m)==solid).AND.(array2D(l+1,m+1)==solid))&
				&.AND.& 
				& ((array2D(l,m+1)==solid).AND.(array2D(l+1,m)==solid))) then 
					array2D_4x(i,j)=solid
				else
					array2D_4x(i,j)=fluid
				end if							
			End if				

			If (array2D_4x(i,j)==fluid) icount=icount+1				
		Enddo
	 Enddo
	 

	 
Write(10,*) 'Porosity of the processed 4x binary image: ', Dble(icount)/((Nx-1)*4+1)/((Ny-1)*4+1)	 
!------------------------------------------------------------------------
!           Map of solid part
!------------------------------------------------------------------------
    open(11,file='Map_4x.dat',STATUS="REPLACE")
        write(11,*) 'Map'
        Do j=(Ny-1)*4+1,1,-1
			Do i=1,(Nx-1)*4+1
				if (array2D_4x(i,j)==fluid) then
					write(11,"(A)",Advance='no') "   "
				else
					write(11,"(I2,A)",Advance='no') array2D_4x(i,j)," "
				end if
			End do
			write(11,*)
        Enddo
    close(11)
!------------------------------------------------------------------------
!           Binary image
!------------------------------------------------------------------------	
output_name='Processed_4x_'//TRIM(ADJUSTL(input_name))	
	Open(200,file=output_name,STATUS="REPLACE")
		Do j=1,(Ny-1)*4+1
			Do i=1, (Nx-1)*4+1
				write(200, "(I2)",Advance='no') array2D_4x(i,j)
			Enddo
			write(200,*)
		Enddo
	Close(200) 	


!=======================================================================
!   Porosity of the binary image
!=======================================================================
Write(10,*) 
Write(10,*) 'Number-of-cell base: '

icount=0
	 Do j=1,(Ny-1)
		Do i=1, (Nx-1)
			If ((array2D(i,j)==fluid).OR.(array2D(i+1,j)==fluid).OR.&
			& (array2D(i,j+1)==fluid).OR.(array2D(i+1,j+1)==fluid)) then
 				icount=icount+1
			Endif				
		Enddo
	 Enddo	 
Write(10,*) 'Porosity of the processed 1x binary image: ', Dble(icount)/((Nx-1))/((Ny-1))

icount=0
	 Do j=1,(Ny-1)*2
		Do i=1, (Nx-1)*2
			If ((array2D_2x(i,j)==fluid).OR.(array2D_2x(i+1,j)==fluid).OR.&
			& (array2D_2x(i,j+1)==fluid).OR.(array2D_2x(i+1,j+1)==fluid)) then
 				icount=icount+1
			Endif				
		Enddo
	 Enddo	 
Write(10,*) 'Porosity of the processed 2x binary image: ', Dble(icount)/((Nx-1)*2)/((Ny-1)*2)

icount=0
	 Do j=1,(Ny-1)*3
		Do i=1, (Nx-1)*3
			If ((array2D_3x(i,j)==fluid).OR.(array2D_3x(i+1,j)==fluid).OR.&
			& (array2D_3x(i,j+1)==fluid).OR.(array2D_3x(i+1,j+1)==fluid)) then
 				icount=icount+1
			Endif				
		Enddo
	 Enddo	 
Write(10,*) 'Porosity of the processed 3x binary image: ', Dble(icount)/((Nx-1)*3)/((Ny-1)*3)

icount=0
	 Do j=1,(Ny-1)*4
		Do i=1, (Nx-1)*4
			If ((array2D_4x(i,j)==fluid).OR.(array2D_4x(i+1,j)==fluid).OR.&
			& (array2D_4x(i,j+1)==fluid).OR.(array2D_4x(i+1,j+1)==fluid)) then
 				icount=icount+1
			Endif				
		Enddo
	 Enddo	 
Write(10,*) 'Porosity of the processed 4x binary image: ', Dble(icount)/((Nx-1)*4)/((Ny-1)*4)



		
	close (10)
	
DEALLOCATE(array2D,array2D_2x,array2D_3x,array2D_4x)	
End Program D2Qn_image_processing_2
