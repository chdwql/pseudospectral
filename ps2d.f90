module ps2d_global

  USE, intrinsic :: iso_c_binding
	include 'fftw3.f03'
	integer :: nx, nz
	type(C_PTR), save :: pfx, pfz, pbx, pbz
	real(C_DOUBLE), save :: pi = 3.141592653589793238462643d0
	real(C_DOUBLE), save :: dx ,dz, dt ,fx, fz
	complex(C_DOUBLE_COMPLEX), save :: cmplx_unit
	complex(C_DOUBLE_COMPLEX), allocatable, save :: kx(:), kz(:)

end module ps2d_global

module ps2d_subs

USE ps2d_global
CONTAINS
!	create dft plan.
	subroutine create_plan()
		implicit none
		real(C_DOUBLE) :: col(nx), row(nz)
		complex(C_DOUBLE_COMPLEX) :: col_f(nx), row_f(nz)
!
		pfx = fftw_plan_dft_r2c_1d(nx, col, col_f, FFTW_ESTIMATE)
		pfz = fftw_plan_dft_r2c_1d(nz, row, row_f, FFTW_ESTIMATE)
		pbx = fftw_plan_dft_c2r_1d(nx, col_f, col, FFTW_ESTIMATE)
		pbz = fftw_plan_dft_c2r_1d(nz, row_f, row, FFTW_ESTIMATE)
	end subroutine create_plan
!	calculate kx, kz
	subroutine cal_kx_kz()
		implicit none
		integer :: ix, iz
		real(C_DOUBLE) :: temp1, temp2
		temp1 = fx * 2.d0 / nx
		temp2 = fz * 2.d0 / nz
		do ix = 1,nx/2
			kx(ix) = (ix * temp1) * cmplx_unit
			kx(ix+nx/2) = (ix * temp1 - fx) * cmplx_unit
		end do

		do iz = 1,nz/2
			kz(iz) = (iz * temp2) * cmplx_unit
			kz(iz+nz/2) = (iz * temp1 - fz) * cmplx_unit
		end do
	end subroutine cal_kx_kz
!	destroy dft plan.
	subroutine destroy_plan()
		implicit none
!
		call fftw_destroy_plan(pfx)
		call fftw_destroy_plan(pfz)
		call fftw_destroy_plan(pbx)
		call fftw_destroy_plan(pbz)
	end subroutine destroy_plan

!	spatial derivative along x.
	subroutine spatial_deriv_x(in, out)
		implicit none
		integer :: ix, iz
		real(C_DOUBLE) :: in(nx, nz), out(nx, nz), col(nx)
		complex(C_DOUBLE_COMPLEX) :: col_f(nx)
!
		do iz = 1, nz
			col = in(:, iz)
			call fftw_execute_dft_r2c(pfx, col, col_f)
			col_f = col_f * kx
			call fftw_execute_dft_c2r(pbx, col_f, col)
			out(:, iz) = col / nx
		end do
	end subroutine spatial_deriv_x

!	spatial derivative along z.
	subroutine spatial_deriv_z(in, out)
		implicit none
		integer :: ix, iz
		real(C_DOUBLE) :: in(nx, nz), out(nx, nz), row(nz)
		complex(C_DOUBLE_COMPLEX) :: row_f(nz)

		do ix = 1, nx
			row = in(ix, :)
			call fftw_execute_dft_r2c(pfz, row, row_f)
			row_f = row_f * kz
			call fftw_execute_dft_c2r(pbz, row_f, row)
			out(ix, :) = row / nz
		end do
	end subroutine spatial_deriv_z

end module ps2d_subs
!pseudo spectral 2d.
program ps2d_aniso
!
	USE ps2d_subs
	implicit none
!
	character(80) :: snapshot
	integer :: it, nsteps, snap
	integer :: xsource, zsource
	integer :: ix, iz, flag
	real(C_DOUBLE) :: a, t0, f0, factor, t, force_x, force_z, source_term
	real(8) :: xlo, xhi, zlo, zhi, ulo, uhi, t_start, t_finish
	real(C_DOUBLE), allocatable, dimension(:, :) :: lambda, mu, rho, cs, cp
	real(C_DOUBLE), allocatable, dimension(:, :) :: vx, vz, sigmaxx, sigmazz, sigmaxz
	real(C_DOUBLE), allocatable, dimension(:, :) :: vx_dx, vx_dz, vz_dx, vz_dz
	real(C_DOUBLE), allocatable, dimension(:, :) :: sigmaxx_dx, sigmazz_dz, sigmaxz_dx, sigmaxz_dz

! model I from Becache, Fauqueux and Joly, which is stable
	real(C_DOUBLE), parameter :: scale_aniso = 1.d10
	real(C_DOUBLE), parameter :: c11 = 4.d0 * scale_aniso
	real(C_DOUBLE), parameter :: c12 = 3.8d0 * scale_aniso
	real(C_DOUBLE), parameter :: c22 = 20.d0 * scale_aniso
	real(C_DOUBLE), parameter :: c33 = 2.d0 * scale_aniso

!	t_start
	call cpu_time(t_start)

!	nx, nz
	nx = 512
	nz = 512

!	allocate array.
	allocate(lambda(nx, nz), mu(nx, nz), rho(nx, nz), cs(nx, nz), cp(nx, nz))
	allocate(vx(nx, nz), vz(nx, nz), sigmaxx(nx, nz), sigmazz(nx, nz), sigmaxz(nx, nz))
	allocate(vx_dx(nx, nz), vx_dz(nx, nz), vz_dx(nx, nz), vz_dz(nx, nz))
	allocate(sigmaxx_dx(nx, nz), sigmazz_dz(nx, nz), sigmaxz_dx(nx, nz), sigmaxz_dz(nx, nz))
	allocate(kx(nx), kz(nz))

!	initialize.
	flag = 1
	snap = 100
	dt = 5.d-4
	f0 = 25.d0
	t0 = 1.2d0/f0
	factor = 1.d7
	xsource = 256
	zsource = 256
	nsteps = 1000
!	cs = 1.3d3
!	cp = 2.d3
	rho = 4.d3
!	mu = rho * cs *cs
!	lambda = rho * (cp * cp - 2.d0 * cs * cs)
	dx = 5
	dz = 5
	xlo = 0.0d0
	xhi = (nx-1)*dx
	zlo = -(nz-1)*dz
	zhi = 0.0d0
	fx = pi / dx
	fz = pi / dz
	cmplx_unit = (0, 1)
	vx = 0.d0
	vz = 0.d0
	sigmaxx = 0.d0
	sigmazz = 0.d0
	sigmaxz = 0.d0

!	create plan.
	call create_plan

!	calculate kx, kz
	call cal_kx_kz

!	cycle time.
	do it = 1, nsteps

!		sigmaxz
		call spatial_deriv_z(vx, vx_dz)
		call spatial_deriv_x(vz, vz_dx)
		sigmaxz = sigmaxz + c33 * (vx_dz + vz_dx) * dt

!		sigmaxx
		call spatial_deriv_x(vx, vx_dx)
		call spatial_deriv_z(vz, vz_dz)
		sigmaxx = sigmaxx + (c11 * vx_dx + c12 * vz_dz) * dt

!		sigmazz
		sigmazz = sigmazz + (c12 * vx_dx + c22 * vz_dz) * dt

!		vx
		call spatial_deriv_x(sigmaxx, sigmaxx_dx)
		call spatial_deriv_z(sigmaxz, sigmaxz_dz)
		vx = vx + (sigmaxx_dx + sigmaxz_dz) / rho * dt

!		vz
		call spatial_deriv_x(sigmaxz, sigmaxz_dx)
		call spatial_deriv_z(sigmazz, sigmazz_dz)
		vz = vz + (sigmaxz_dx + sigmazz_dz) / rho * dt

!		add source
		a = pi*pi * f0*f0
		t = (it-1) * dt
		source_term = -factor * (1.d0 - 2.d0*a*(t-t0)**2)*exp(-a*(t-t0)**2)
		force_x = sin(0.d0) * source_term
		force_z = cos(0.d0) * source_term
		vx(xsource, zsource) = vx(xsource, zsource) + force_x * dt / rho(1,1)
		vz(xsource, zsource) = vz(xsource, zsource) + force_z * dt / rho(1,1)

!		output snapshot.
		if(mod((it-1), snap)==0 .and. flag==1) then
			write(*, *) 'time steps: ', int(dt*(it-1)*1.e3), 'ms'
			write(snapshot, "('snapshot',i6.6,'.grd')") int(dt*(it-1)*1.e3)
			ulo=minval(vx(1:nx, 1:nz))
			uhi=maxval(vx(1:nx, 1:nz))
			open(15, file='snapshots\'//snapshot, status='unknown')
			write(15, '(4a)')'DSAA'
			write(15, *)nx, nz
			write(15, *)xlo, xhi
			write(15, *)zlo, zhi
			write(15, *)ulo, uhi
			do iz = nz, 1, -1
				write(15, *)(vx(ix,iz), ix=1, nx)
			end do
			close(15)
		end if
	end do

!	clean
	call destroy_plan
	deallocate(lambda, mu, rho, cs, cp)
	deallocate(vx, vz, sigmaxx, sigmazz, sigmaxz)
	deallocate(vx_dx, vx_dz, vz_dx, vz_dz)
	deallocate(sigmaxx_dx, sigmazz_dz, sigmaxz_dx, sigmaxz_dz)

!
	call cpu_time(t_finish)
	print*, 'running time: ', t_finish-t_start
	stop
end program ps2d_aniso
