program tight_binding_model

	!SCALAR VARIABLES
	integer, parameter :: N=100, LWMAX = 1000
	double precision, parameter :: t=-1.0d0, KB = 8.617d-5, pi = 4.0d0*atan(1.0d0)
	character, parameter :: JOBZ='V', UPLO='U'
	complex*16 :: imag = -1.0d0
	integer :: LDA=N, INFO, LWORK
	double precision  B, Nl, Ne, MuGs, high, low, phib, x, u, RoTar, Ro, temp

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

	!Calling ZHEEV to find eigen energies
	LWORK = -1
	call ZHEEV(JOBZ, UPLO, N, H0, LDA, E, WORK, LWORK, RWORK, INFO)
	LWORK = min( LWMAX, int(WORK(1)))

	call ZHEEV(JOBZ, UPLO, N, H0, LDA, E, WORK, LWORK, RWORK, INFO)
	
	if (INFO == 0) then
		write(*,*)'SUCCESS'
	endif

	!WRITING EIGENVALUES TO FILE
	write(*,*) 'The eigenvalues are:'
	do l=1, n
		write(*,*) E(l)
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
		
		!MODIFYING Ro		
		Ne = 0
		do l=1,N
			Nl = 1/(exp(B*(E(l)-MuGs)) + 1)
			Ne = Ne + (2*Nl)
		enddo
	
		Ro = Ne/N
	enddo

	write(*,*)'The chemical potential is:'
	write(*,*) MuGs
	
	!Calculating total potential energy
	u = 0
	do l=1,N
		Nl = 1/(exp(B*(E(l)-MuGs)) + 1)
		u = E(l)*2*Nl
	enddo

	write(*,*)'Total potential energy:'
	write(*,*) u
		
	
end program



