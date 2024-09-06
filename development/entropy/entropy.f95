#this function calculates the entropy associated with 
#the random sampling used in main code
#may or may not be necessary
#only makes small changes to total electronic energy
function entropy(aray,nvals,temp)

implicit none

real, intent(in) :: temp

integer, intent(in) :: nvals

real, intent(in) :: aray(nvals)

integer :: i , j , w , v , p , rng , numm , idx , maxi , mini , kount , hodl

real :: top , bottom , summ , enti

integer, allocatable :: freq(:) , truefreq(:)

real, parameter :: k=0.000003166811563

real :: entropy !.. output variable

top = MAXVAL(aray)
maxi = floor(top)

bottom = MINVAL(aray)
mini = floor(bottom)

rng =  (maxi - mini) + 1

allocate(freq(rng))

freq = 0

do i = 1 , nvals

	numm = floor(aray(i))
	
	idx = ((numm - mini) + 1)
	
	freq(idx) = freq(idx) + 1
	
end do

kount = 0

do j = 1 , rng

	if (freq(j) > 0) then
	kount = kount + 1
	end if
	
end do

allocate(truefreq(kount))

truefreq = 0

hodl = 0

do v = 1 , rng

	if (freq(v) > 0) then
	hodl = hodl + 1
	truefreq(hodl) = freq(v)
	end if
	
end do

summ = 0

do p = 1 , kount

	summ = summ + real((real(truefreq(p))/real(nvals))*LOG((real(truefreq(p))/real(nvals))))
	
end do

enti = -k*summ

write(*,*) SUM(freq),summ,enti
write(*,*) kount,maxval(truefreq),real(real(maxval(truefreq))/real(nvals))
write(*,*) maxval(freq),MAXLOC(freq),mini+MAXLOC(freq)

entropy = temp*enti

deallocate(freq)
deallocate(truefreq)

end function entropy
