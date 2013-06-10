program tight_binding_model

	!SCALAR VARIABLES
	integer, parameter :: N=4, LWMAX = 1000, output=10
	double precision, parameter :: t=-1.0d0, KB = 8.617d-5, pi = 3.14159265359d0
	character, parameter :: JOBZ='V', UPLO='U'
	complex*16 :: imag = -1.0d0
	integer :: LDA=N, INFO, LWORK, temp
	double precision  B, Nl, Ne, MuGs, high, low, phib, x, u
	double precision RoTar, Ro

	!ARRAY VARIABLES
	double precision   RWORK(3*N-2), E(N)
	complex*16         H0(N,N), WORK(LWMAX)

	!EXTERNAL SUBROUTINES
	external ZHEEV

	!INTRINSIC FUNCTIONS
	intrinsic int, min, sqrt, exp

		
	write(*,*)'Enter magnetic flux'
	read(*,*) phib
	
	x = (phib/dfloat(N))*2.0d0*pi
	imag = sqrt(imag)
	write(*,*) imag

	!CONSTRUCTING THE HAMILTONIAN
	do ix=1, N
		do jx=1, N
			if (((jx-ix) == 1) .or. ((jx-ix) == -(N-1))) then
				H0(ix,jx) = t * exp(imag*x)
			elseif (((jx-ix) == -1) .or. ((jx-ix) == (N-1))) then
				H0(ix,jx) = t * exp(-(imag*x))
			else
				H0(ix,jx) = 0
			endif
		enddo
	enddo

	write(*,*) H0

	LWORK = -1
	call ZHEEV(JOBZ, UPLO, N, H0, LDA, E, WORK, LWORK, RWORK, INFO)
	LWORK = min( LWMAX, int(WORK(1)))

	call ZHEEV(JOBZ, UPLO, N, H0, LDA, E, WORK, LWORK, RWORK, INFO)
	
	if (INFO == 0) then
		write(*,*)'SUCCESS'
	endif
	write (*,*) E
	!WRITING EIGENVALUES TO FILE
	open(unit=output,file='results.txt',action='write',status='replace')
	write(output,*) 'The eigenvalues are:'
	do l=1, n
		write(output,*) E(l)
	enddo

	!GETTING TEMP AND DENSITY FROM USER
	write(*,*)'Enter temperature'
	read(*,*) temp
	write(*,*)'Enter target density'
	read(*,*) RoTar

	!CALCULATING FIRST MuGs AND Ro
	high = E(N)
	low = E(1)	
	MuGs = (high+low) / 2
	B = 1/(KB*temp)
	Ne = 0
	do l=1,N
		Nl = 1/(exp(B*(E(l)-MuGs)) + 1)
		Ne = Ne + (2*Nl)
	enddo

	Ro = Ne/N
	
	!COMPARING Ro TO RoTar
	do while (abs((Ro - RoTar) / RoTar) > 0.01)
		
		!MODIFYING MuGs
		if (Ro > RoTar) then
			high = MuGs
			MuGs = (high + low) / 2	
		else
			low = MuGs
			MuGs = (high + low) / 2
		endif
		write(*,*) MuGs
		
		!MODIFYING Ro		
		Ne = 0
		do l=1,N
			Nl = 1/(exp(B*(E(l)-MuGs)) + 1)
			Ne = Ne + (2*Nl)
		enddo
	
		Ro = Ne/N
		write(*,*) Ro
	enddo

	write(output,*)'The chemical potential is:'
	write(output,*) MuGs
	close(output)
	
	!Calculating total potential energy
	u = 0
	do l=1,N
		Nl = 1/(exp(B*(E(l)-MuGs)) + 1)
		u = E(l)*2*Nl
	enddo

	write(*,*) u
		
	
end program



