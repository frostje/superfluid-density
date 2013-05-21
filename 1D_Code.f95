program tight_binding_model

!SCALAR VARIABLES
	integer, parameter :: N=10, LWMAX = 1000
	real, parameter :: t=-1.0
	character, parameter :: JOBZ='V', UPLO='U'
	integer :: LDA=N, INFO, LWORK
	
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
	
	
end program



