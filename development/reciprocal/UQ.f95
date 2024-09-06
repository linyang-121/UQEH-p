program UQ

!.. this program uses Hartree atomic units
!.. elementary charge, permittivity, electron rest mass are all 1
!.. this program is limited to one atom type
!.. nuclear radius based on formula R=r_0*(2Z)^(1/3)
!.. makes total kinetic energy positive to align with physical rationale
!.. uses opposite sign convention with coulombic potential to match QE convention
!.. has electrons(:) array allocate and deallocate itself between trials

implicit none

double precision :: dummy , x , y , z , rx , ry , rz , mean , variance , std , min , max

double precision, parameter :: convert=(1./0.529177249) , rsub = 0.00002329545454545455

double precision :: gap , ketot , sumsq=0. , nucradi

double precision :: elec_elec_inter , elec_nuc_inter , nuc_nuc_inter , hammy , peri , ppotential

double precision :: twopoints !.. function that calculates distance between two points

double precision :: periodic

integer :: i , j , k , m , n , p , q , u , v , f , g , b , c , elec , atomkount , atomnum , nbin

integer :: ntrials , expanse

double precision, allocatable :: ke(:) , pos(:,:) , electrons(:,:) , sample(:)

!.. open input files
open(24, FILE='params.txt')
open(25, FILE='ke.txt')
open(26, FILE='pos.txt')

!.. skip comment lines from input files
read(24,*)
read(24,*)
read(24,*)
read(25,*)
read(25,*)
read(25,*)
read(26,*)
read(26,*)
read(26,*)

!.. read input parameters
read(24,*) elec , atomnum , atomkount , x , y , z , ntrials , expanse

!.. defines nuclear radius
nucradi = rsub*(dfloat(2*atomnum)**(1./3.))

allocate(ke(elec))
allocate(pos(atomkount,3))
allocate(sample(ntrials))
allocate(electrons(elec,3))

!.. places electron kinetic energies in array
do i = 1 , (elec)
	read(25,*) ke(i)
end do

!.. places atom positions in array
do j = 1 , atomkount
	read(26,*) dummy , pos(j,1) , pos(j,2) , pos(j,3)
end do

!.. converts all atomic coordinates from angstrom to atomic units
pos = pos*convert

!.. generate random electron positions
!.. treats electrons as infinitesimal point
!.. assumes electrons cannot be within bohr radius of infinitesmial point that defines atomic position

nuc_nuc_inter=0.

!.. calculates potential between nucleus nucleus pairs

focus: do f = 1 , atomkount

relative: do g = 1 , atomkount
	
if (f >= g) then 
peri = periodic(pos(f,1),pos(f,2),pos(f,3),pos(g,1),pos(g,2),pos(g,3),x,y,z,dfloat(atomnum**2),expanse)
nuc_nuc_inter = nuc_nuc_inter + peri
			
cycle relative
end if
peri = periodic(pos(f,1),pos(f,2),pos(f,3),pos(g,1),pos(g,2),pos(g,3),x,y,z,dfloat(atomnum**2),expanse)
		
nuc_nuc_inter = nuc_nuc_inter + peri
		
gap = twopoints(pos(f,1),pos(f,2),pos(f,3),pos(g,1),pos(g,2),pos(g,3))
		
nuc_nuc_inter = nuc_nuc_inter + (dfloat(atomnum**2)/gap)
	
end do relative
	
end do focus

trials: do b = 1 , ntrials

!.. reinitialize parameters
elec_elec_inter=0.
elec_nuc_inter=0.

main: do k = 1  , elec

call reseed
	
electron_position: do
call random_number(rx)
rx=rx*x
call random_number(ry)
ry=ry*y
call random_number(rz)
rz=rz*z
		
check_atoms: do m = 1 , atomkount !..check to make sure electron doesn't overlap atom
		
gap = twopoints(rx,ry,rz,pos(m,1),pos(m,2),pos(m,3))
			
if (gap > nucradi) then 
cycle check_atoms
else
cycle electron_position
end if
end do check_atoms
		
if (k == 1) then !.. if first electron position then generate second electron position
		
electrons(k,1) = rx
electrons(k,2) = ry
electrons(k,3) = rx
exit electron_position
end if
		
check_electrons: do n = 1 , k-1 !.. check to make sure electron doesn't overlap other electrons
		
gap = twopoints(rx,ry,rx,electrons(n,1),electrons(n,2),electrons(n,3))
			
if (gap > 0.) then 
cycle check_electrons
else
cycle electron_position
end if
end do check_electrons
		
electrons(k,1) = rx
electrons(k,2) = ry
electrons(k,3) = rx
exit electron_position
		
end do electron_position

end do main

!.. calculates potential between electron electron pairs

prime: do p = 1 , elec
	
compare: do q = 1 , elec
	
if (p >= q) then
peri = periodic(electrons(p,1),electrons(p,2),electrons(p,3),electrons(q,1),&
electrons(q,2),electrons(q,3),x,y,z,dfloat(1),expanse)
		
elec_elec_inter = elec_elec_inter + peri
			
cycle compare
end if
peri = periodic(electrons(p,1),electrons(p,2),electrons(p,3),electrons(q,1),&
electrons(q,2),electrons(q,3),x,y,z,dfloat(1),expanse)
		
elec_elec_inter = elec_elec_inter + peri
		
gap = twopoints(electrons(p,1),electrons(p,2),electrons(p,3),electrons(q,1),electrons(q,2),electrons(q,3))
		
elec_elec_inter = elec_elec_inter + (dfloat(1)/gap)
	
end do compare
	
end do prime

!.. calculates potential between electron nucleus pairs

alpha: do u = 1 , elec
	
beta: do v = 1 , atomkount
	
peri = periodic(electrons(u,1),electrons(u,2),electrons(u,3),pos(v,1),pos(v,2),pos(v,3),x,y,z,dfloat(atomnum),expanse)
		
elec_nuc_inter = elec_nuc_inter + peri
		
gap = twopoints(electrons(u,1),electrons(u,2),electrons(u,3),pos(v,1),pos(v,2),pos(v,3))
		
elec_nuc_inter = elec_nuc_inter + (dfloat(atomnum)/gap)
	
end do beta
	
end do alpha

!.. calculates sum of electron kinetic energies
ketot = SUM(ke)

!.. calculates value of electronic Hamiltonian of system

ppotential = elec_elec_inter-elec_nuc_inter+nuc_nuc_inter
hammy = ketot+ppotential

sample(b) = hammy

write(100,*) ketot,elec_elec_inter,elec_nuc_inter,nuc_nuc_inter

end do trials

!.. calculate statistics of Hamiltonian sample

mean = SUM(sample)/dfloat(ntrials)

do c = 1 , ntrials

sumsq = sumsq + ((sample(c) - mean)**2)
	
end do

variance = sumsq/dfloat(ntrials)

std = sqrt(variance)

min = MINVAL(sample)
max = MAXVAL(sample)

write(30,*) 'Minimum value for electronic Hamiltonian=',min
write(30,*) 'Maximum value for electronic Hamiltonian=',max
write(30,*) 'Mean value for electronic Hamiltonian=',mean
write(30,*) 'Variance of values for electronic Hamiltonian=',variance
write(30,*) 'Standard Deviation of values for electronic Hamiltonian=',std

call histo(sample,ntrials)	
				
close(24)
close(25)
close(26)

deallocate(ke)
deallocate(pos)
deallocate(sample)
deallocate(electrons)

end program UQ

