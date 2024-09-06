function twopoints(x1,y1,z1,x2,y2,z2)

!.. This function finds the distance between to points in 3D
!.. using the distance formula in cartesian corrdinates
!.. i.e. r=sqrt(x^2+y^2+z^2)

implicit none
	
double precision, intent(in) :: x1 , y1 , z1 , x2 , y2 , z2 !.. input variables

double precision :: xdiff , ydiff , zdiff , xdsq , ydsq , zdsq !.. internal variables
	
double precision :: twopoints !.. output variable

!.. calculates difference between each dimension
xdiff = x2 - x1

ydiff = y2 - y1

zdiff = z2 - z1

!.. squared each dimensional difference
xdsq = xdiff**2

ydsq = ydiff**2

zdsq = zdiff**2

!.. square roots the sum of the sqaures to get distance
!.. between the two points
twopoints = sqrt(xdsq+ydsq+zdsq)

end function twopoints
