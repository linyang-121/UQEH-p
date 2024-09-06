program UQ

!.. This program will generate random points for electrons within the dimension
!.. given by the "params.txt" input file. Based on these randomly generated electron
!.. positions and the nuclear positions given by "pos.txt" it will generate a value for
!.. the total electronic energy of an isolated system using the electronic Hamiltonian for
!.. a monoatomic system of the size defined by the input file "params.txt".

!.. The input file "pos.txt" must be formatted so that there are three commented lines at
!.. the top of the file and each line thereafter contains one point in the format 'A x y z'.
!.. Where A is the atomic number.
!.. Make sure that all the points are positive values to ensure proper sampling of the
!.. isolated system.

!.. This program will randomly generate electron positions the amount of times
!.. defined by the variable "ntrials" as defined in the input file "params.txt".
!.. The maximum value, minimum value, mean value, variance, and standard deviation
!.. of the corresponding electronic energies will be place in the output file "fort.30".

!.. The output file "fort.32" will have the electronic energies organized so that a histogram
!.. or freqency plot may be generated for quick analysis. (see "histo.f95" for details)

!.. The kinetic energy for each electron is defined by the input file "ke.txt" and can be generated
!.. however the user may wish. So long as you leave three commented lines at the top of the input file
!.. and one electron kinetic energy value per line after. This makes sure the main program reads the
!.. input file correctly. The Mathematica notebook "calculate_electron_distro.nb" is the method provided to
!.. quickly generate the required electron kinetic energies based on the Fermi-Dirac or Maxwell-Boltzmann
!.. distribution. Whichever your system requires.

!.. This program uses Hartree atomic units which means the
!.. elementary charge, coulomb constant, and electron rest mass are all 1.

!.. This program is limited to one atom type.

!.. The nuclear radius based on formula R=r_0*(2Z)^(1/3).
!.. This program assumes the atomic nuclei have the same
!.. number of protons and neutrons.

!.. This program makes total kinetic energy positive to align with physical rationale.

!.. Uses opposite sign convention with coulombic potential to match QE convention
!.. coulombic potential sign convention can be swapped to match system you are investigating
!.. so long as you are consistent.

implicit none

double precision :: dummy , x , y , z , rx , ry , rz , hammy , mean , variance , std , min , max

double precision, parameter :: convert=(1./0.529177249) , rsub = 0.00002329545454545455

double precision :: gap !.. temporarily stores output of twopoints function for analysis

double precision :: ketot !.. stores sum of electron kinetic energies

double precision :: sumsq=0. !.. stores sum of squares for variance calculation

double precision :: nucradi !.. store value of atomic nuclei radius

double precision :: elec_elec_inter=0. !.. Stores the sum of electron/electron coulomb interaction

double precision :: elec_nuc_inter=0. !.. Stores the sum of electorn/nuclear coulomb interaction

double precision :: nuc_nuc_inter=0. !.. Stores the sum of nuclear/nuclear coulomb interaction

double precision :: twopoints !.. function that calculates distance between two points

integer :: i , j , k , m , n , p , q , u , v , f , g , b , c !.. each of these are used as an index for a different loop

integer :: elec , atomkount , atomnum !.. stores the electron count, nuceli count, and atomic number

integer :: ntrials !.. stores the number of times electron positions are sampled.

double precision, allocatable :: ke(:) !.. stores electron kinetic energies

double precision, allocatable :: pos(:,:) !.. stores positions of atomic nuceli

double precision, allocatable :: electrons(:,:) !.. stores poitions of electrons

double precision, allocatable :: sample(:) !.. stores value of total electronic energies

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
read(24,*) elec , atomnum , atomkount , x , y , z , ntrials

!.. defines nuclear radius
nucradi = rsub*(real(2*atomnum)**(1./3.))

allocate(ke(elec))
allocate(pos(atomkount,3))
allocate(sample(ntrials))

!.. places electron kinetic energies in array
do i = 1 , elec
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

trials: do b = 1 , ntrials

allocate(electrons(elec,3))

!.. reinitialize parameters
elec_elec_inter=0.
elec_nuc_inter=0.
nuc_nuc_inter=0.


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
			cycle compare
		end if
		
		gap = twopoints(electrons(p,1),electrons(p,2),electrons(p,3),electrons(q,1),electrons(q,2),electrons(q,3))
		
		elec_elec_inter = elec_elec_inter + (1./gap)
	
	end do compare
	
end do prime

!.. calculates potential between electron nucleus pairs

alpha: do u = 1 , elec
	
	beta: do v = 1 , atomkount
		
		gap = twopoints(electrons(u,1),electrons(u,2),electrons(u,3),pos(v,1),pos(v,2),pos(v,3))
		
		elec_nuc_inter = elec_nuc_inter + (real(atomnum)/gap)
	
	end do beta
	
end do alpha

!.. calculates potential between nucleus nucleus pairs

focus: do f = 1 , atomkount

	relative: do g = 1 , atomkount
	
		if (f >= g) then 
			cycle relative
		end if
		
		gap = twopoints(pos(f,1),pos(f,2),pos(f,3),pos(g,1),pos(g,2),pos(g,3))
		
		nuc_nuc_inter = nuc_nuc_inter + ((real(atomnum**2))/gap)
	
	end do relative
	
end do focus

!.. calculates sum of electron kinetic energies
ketot = SUM(ke)
		
!.. calculates value of electronic Hamiltonian of system

!.. Can change the signs on the corresponding terms here
!.. to match system you are comparing to.
!.. Make sure you are careful to be consistant and know
!.. what you are doing.

hammy = ketot-elec_elec_inter+elec_nuc_inter-nuc_nuc_inter

sample(b) = hammy

deallocate(electrons)

end do trials

!.. calculate statistics of Hamiltonian sample

mean = SUM(sample)/real(ntrials)

do c = 1 , ntrials

	sumsq = sumsq + ((sample(c) - mean)**2)
	
end do

variance = sumsq/real(ntrials)

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


end program UQ

