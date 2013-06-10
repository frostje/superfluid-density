do iy=1,Ny
	do ix=1,Nx
		
		!Right
		index = (iy-1)*Nx + ix
		if (ix == Nx) then
			ixr = 1
		else
			ixr = ix + 1
		endif
		indexp = (iy-1)*Nx + ixr

		!Left
		if (ix == 1) then
			ixl = Nx
		else
			ixl = ix - 1
		endif
		
		!Up
		if (iy == Ny) then
			iyu = 1
		else
			iyu = iy + 1
		endif

		!Down
		if (iy == 1) then
			iyd = Ny
		else
			iyd = iy -1
		endif












do iy=1,Ny
	do ix=1,Nx
		do jy=1,Ny
			do jx=1,Nx
				ipos = (iy-1)*Nx + ix
				fpos = (jy-1)*Nx + jx
			

