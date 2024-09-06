subroutine histo(aray,nvals)

!.. This subroutine takes the each value in an array and uses the floor function
!.. to take the value to the closest lower integer. Then it counts the frequency
!.. of each integer value and prints the result in two columns of the fort.32 file.
!.. The first column is the integer value and the second column is the frequency of that
!.. integer value. This fort.32 file can be used to create a histogram or frequency plot
!.. in your favorite plotting software. A GNUplot script is provide in the repository 
!.. for your convience.

implicit none

integer, intent(in) :: nvals !.. input of the number of values in array

double precision , intent(in) :: aray(nvals) !.. input of name of array

integer :: i , j , k , rng , numm , idx , maxi , mini

double precision :: top , bottom

integer, allocatable :: freq(:) , values(:)

top = MAXVAL(aray)
maxi = floor(top)

bottom = MINVAL(aray)
mini = floor(bottom)

rng =  (maxi - mini) + 1

allocate(freq(rng))
allocate(values(rng))

freq = 0.

do j = 1 , rng
	
	values(j) = mini + (j-1) 
	
end do

do i = 1 , nvals

	numm = floor(aray(i))
	
	idx = ((numm - mini) + 1)
	
	freq(idx) = freq(idx) + 1
	
end do

do k = 1 , rng

write(32,*) values(k) , freq(k)

end do

deallocate(freq)
deallocate(values)

end subroutine histo
