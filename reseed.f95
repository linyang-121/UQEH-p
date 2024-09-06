subroutine reseed

!.. This subroutine will generate a new random seed based on the time
!.. and date that the subroutine is called.

implicit none

double precision :: numm , faktor
integer, allocatable :: seed(:)
integer :: i , kount=80 , n=80
character(len=8) :: date
character(len=10) :: time
character(len=80) :: card
character(len=80) :: hodl

call random_seed(size=n)
allocate(seed(kount))

call date_and_time(date,time)

write(card,*) time(5:10)
read(card,*) numm

write(hodl,*) date(1:8)
read(hodl,*) faktor

numm = sqrt(numm)

numm = numm*faktor

do i = 1 , kount
	
	seed(i) = numm + faktor*(i-1)
end do

call random_seed(put = seed)

deallocate(seed)
close(2)

end subroutine reseed
