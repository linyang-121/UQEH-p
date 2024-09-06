function periodic(xsub,ysub,zsub,xrel,yrel,zrel,xdim,ydim,zdim,fak,expan)

implicit none
	
double precision, intent(in) :: xsub , ysub , zsub , xrel , yrel , zrel !.. input variables

double precision, intent(in) :: xdim , ydim , zdim , fak !.. input variables

integer, intent(in) :: expan

double precision :: twopoints !.. function that calculates distance between two points

double precision, parameter :: pie = 3.14159265359

double precision :: gap , jumpp , hopp , skipp , xalt , yalt , zalt !.. internal variables

double precision :: invjumpp , invhopp , invskipp , invxalt , invyalt , invzalt , invxrel , invyrel , invzrel

double precision :: etot , invx , invy , invz , expo , gsumm , gcosine , faktor , konstant , nu
	
double precision :: periodic !.. output variable

integer :: i , j , k

etot = 0.

gap = 0.

invx = (dfloat(2)*pie)/xdim

invy = (dfloat(2)*pie)/ydim

invz = (dfloat(2)*pie)/zdim

invxrel = (dfloat(2)*pie)/xrel

invyrel = (dfloat(2)*pie)/yrel

invzrel = (dfloat(2)*pie)/zrel

!nu = sqrt((((dfloat(2)*pie)/(0.1))**2+((dfloat(2)*pie)/(0.1))**2+((dfloat(2)*pie)/(0.1))**2))
nu = sqrt((invx)**2+(invy)**2+(invz)**2)

zcorr: do i = 0 , expan
	skipp = dfloat(i)*zdim
	invskipp = dfloat(i)*invz
	ycorr: do j = 0 , expan
		hopp =dfloat(j)*ydim
		invhopp =dfloat(j)*invy
		xcorr: do k = 0 , expan
			if (i==0 .AND. j==0 .AND. k==0) then
			cycle xcorr
			end if
			
			jumpp = dfloat(k)*xdim
			invjumpp = dfloat(k)*invx
			
			xalt = xrel+jumpp
			yalt = yrel+hopp
			zalt = zrel+skipp
			
			invxalt = invxrel+invjumpp
			invyalt = invyrel+invhopp
			invzalt = invzrel+invskipp
			
			gap = twopoints(xsub,ysub,zsub,xalt,yalt,zalt)
			
			expo = exp((-(invxalt**2+invyalt**2+invzalt**2))/(dfloat(4)*nu**2))
			
			gsumm = expo/(invxalt**2+invyalt**2+invzalt**2)
			
			gcosine = cos(((invxalt*xsub)+(invyalt*ysub)+(invzalt*zsub)))
			
			faktor = ((dfloat(4)*pie)/(xdim*ydim*zdim))
			
			konstant = (pie/(xdim*ydim*zdim*(nu**2)))
			
			etot = etot + fak*((erfc(nu*gap)/gap)+(faktor*(gsumm*gcosine-konstant)))

		end do xcorr
	end do ycorr
end do zcorr

periodic = etot

end function periodic
