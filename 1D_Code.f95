program tight_binding_model

	!SCALAR VARIABLES
	integer, parameter :: N=10, LWMAX = 1000, output=10
	double precision, parameter :: t=-1.0, KB = 8.613E-5
	character, parameter :: JOBZ='V', UPLO='U'
	complex :: i = -1
	integer :: LDA=N, INFO, LWORK, temp
	real :: B, Nl, Ne, RoTar, Ro, MuGs, high, low
	
	
	!ARRAY VARIABLES
	double precision   RWORK(3*N-2), E(N)
	complex*16         H0(N,N), WORK(LWMAX)

	!EXTERNAL SUBROUTINES
	external ZHEEV

	!INTRINSIC FUNCTIONS
	intrinsic int, min, sqrt

	i = sqrt(i)

	!CONSTRUCTING THE HAMILTONIAN
	do ix=1, N
		do jx=1, N
			if (abs(jx-ix) == 1) then
				H0(ix,jx) = t
			elseif (abs(jx-ix) == (N-1)) then
				H0(ix,jx) = t
			else
				H0(ix,jx) = 0
			endif
		enddo
	enddo

	LWORK = -1
	call ZHEEV(JOBZ, UPLO, N, H0, LDA, E, WORK, LWORK, RWORK, INFO)
	LWORK = min( LWMAX, int(WORK(1)))

	call ZHEEV(JOBZ, UPLO, N, H0, LDA, E, WORK, LWORK, RWORK, INFO)
	
	if (INFO == 0) then
		write(*,*)'SUCCESS'
	endif

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
	enddo

	write(output,*)'The chemical potential is:'
	write(output,*) MuGs
	close(output)
	
end program



