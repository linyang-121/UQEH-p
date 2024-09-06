function periodic(xsub,ysub,zsub,xrel,yrel,zrel,xdim,ydim,zdim,fak,exp,zpot)

implicit none
	
double precision, intent(in) :: xsub , ysub , zsub , xrel , yrel , zrel !.. input variables

double precision, intent(in) :: xdim , ydim , zdim , fak , zpot!.. input variables

integer, intent(in) :: exp

double precision :: twopoints !.. function that calculates distance between two points

double precision :: gap , jumpp , hopp , skipp , xalt , yalt , zalt !.. internal variables

double precision :: etot
	
double precision :: periodic !.. output variable

integer :: i , j , k

etot = 0.

gap = 0.

zcorr: do i = 0 , exp
	skipp = dfloat(i)*zdim
	ycorr: do j = 0 , exp
		hopp =dfloat(j)*ydim
		xcorr: do k = 0 , exp
			if (i==0 .AND. j==0 .AND. k==0) then
			cycle xcorr
			end if
			
			jumpp = dfloat(k)*xdim
			
			xalt = xrel+jumpp
			yalt = yrel+hopp
			zalt = zrel+skipp
			gap = twopoints(xsub,ysub,zsub,xalt,yalt,zalt)
			etot = etot + ((fak/gap) - (fak*zpot))

		end do xcorr
	end do ycorr
end do zcorr

periodic = etot

end function periodic
