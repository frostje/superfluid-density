program tight_binding_model

	!SCALAR VARIABLES
	integer, parameter :: N=10, LWMAX = 1000, output=10
	real, parameter :: t=-1.0, KB = 8.613E-5
	character, parameter :: JOBZ='V', UPLO='U'
	integer :: LDA=N, INFO, LWORK, temp
	real B, Nl, Ne, RoTar, Ro, MuGs, high, low
	
	
	!ARRAY VARIABLES
	double precision   RWORK(3*N-2), W(N)
	complex*16         H0(N,N), WORK(LWMAX)

	!EXTERNAL SUBROUTINES
	external ZHEEV

	!INTRINSIC FUNCTIONS
	intrinsic int, min

	!CONSTRUCTING THE HAMILTONIAN
	do i=1, n
		do j=1,n
			if (abs(j-i) == 1) then
				H0(i,j) = t
			else
				H0(i,j) = 0
			endif
		enddo
	enddo

	LWORK = -1
	call ZHEEV(JOBZ, UPLO, N, H0, LDA, W, WORK, LWORK, RWORK, INFO)
	LWORK = min( LWMAX, int(WORK(1)))

	call ZHEEV(JOBZ, UPLO, N, H0, LDA, W, WORK, LWORK, RWORK, INFO)
	
	if (INFO == 0) then
		write(*,*)'SUCCESS'
	endif

	!WRITING EIGENVALUES TO FILE
	open(unit=output,file='results.txt',action='write',status='replace')
	write(output,*) 'The eigenvalues are:'
	do i=1,n
		write(output,*) W(i)
	enddo

	!GETTING TEMP AND DENSITY FROM USER
	write(*,*)'Enter temperature'
	read(*,*) temp
	write(*,*)'Enter target density'
	read(*,*) RoTar

	!CALCULATING FIRST MuGs AND Ro
	high = W(N)
	low = W(1)	
	MuGs = (high+low) / 2
	B = 1/(KB*temp)
	Ne = 0
	do i=1,N
		Nl = 1/(exp(B*(W(i)-MuGs)) + 1)
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
		do i=1,N
			Nl = 1/(exp(B*(W(i)-MuGs)) + 1)
			Ne = Ne + (2*Nl)
		enddo
	
		Ro = Ne/N
	enddo

	write(output,*)'The chemical potential is:'
	write(output,*) MuGs
	close(output)
	
end program



