! Program to estimate the best parameters out of the selected few
!
! ****************************************	YET TO ADD ****************************************
!
! estimate other statistically relevant parameters like mean, standard deviation, 95% confidence interval etc.
!
! *********************************************************************************************
module init
	implicit none
	
	integer, parameter :: n_estimate = 100, n_parameter = 14, n_data = 40    ! n_parameter is actualy "number of parameters + 2"
	character (len=9), parameter :: file_input = 'air20.dat'				! enter the name of the input file here
	
	real(kind=8), parameter :: pi = 4.0*atan(1.0)
	
end module init

! THIS FUNCTION NEEDS TO BE UPDATED WHEN THE EQUIVALENT CIRCUIT IS MODIFIED < updated for wagner circuit : 21 jan >
complex (kind=8) function imp_calc(w, ra, ca, phia, rpa, rel, rc, cc, phic, rpc, l, rn, cn)

	implicit none
	real (kind=8), intent(in) :: w, ra, ca, phia, rpa, rel, rc, cc, phic, rpc, l, rn, cn
	complex (kind=8) :: iw, zra, zca, zpa, zel, zrc, zcc, zpc, zl, zrn, zcn, zn, zqa, zqc, cmplx_temp_a, cmplx_temp_c, &
						& cmplx_tanh_a, cmplx_tanh_c, z_anode, z_cathode, ztot
	
	iw = cmplx(0.0,w)
	zra = cmplx(ra)
	zpa = cmplx(rpa)
	zel = cmplx(rel)
	zrc = cmplx(rc)
	zpc = cmplx(rpc)
	zrn = cmplx(rn)
	
	zca = 1.0/ca*(iw**phia)
	
	zcc = 1.0/cc*(iw**phic)
	
	zcn = 1.0/cn*iw
	
	zl = iw*l
	
	zqa = zra*zca/(zra+zca)
	
	zqc = zrc*zcc/(zrc+zcc)
	
	zn = zrn*zcn/(zrn+zcn)
	
	cmplx_temp_a = sqrt(zpa/zqa)
	
	cmplx_temp_c = sqrt(zpc/zqc)

	cmplx_tanh_a = (exp(2.0*cmplx_temp_a)-1.0)/(exp(2.0*cmplx_temp_a)+1.0)
	
	cmplx_tanh_c = (exp(2.0*cmplx_temp_c)-1.0)/(exp(2.0*cmplx_temp_c)+1.0)
	
	z_anode = sqrt(zqa*zpa)/cmplx_tanh_a
	
	z_cathode = sqrt(zqc*zpc)/cmplx_tanh_c
	
	ztot = z_anode + zel + z_cathode + zl + zn
	
	imp_calc = ztot	

end function imp_calc

! THIS FUNCTION NEEDS TO BE UPDATED WHEN THE EQUIVALENT CIRCUIT IS MODIFIED < updated for wagner model on 21 jan >
subroutine model_data(p, handle)
	use init
	implicit none
	
	real (kind=8), intent(in), dimension(n_parameter) :: p
	integer, intent(in) :: handle
	real (kind=8) :: w, phiw, zreal, zimag, freq
	complex (kind=8) :: imp_calc, ztot
	integer :: i
	
	phiw = 0.5
	
	open (unit=1, file=file_input, status='old', action='read')
		
	do i = 1,n_data
		read (1,*) zreal, zimag, freq
		w = 2.0*pi*freq
		ztot = imp_calc(w, p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8), p(9), p(10), p(11), p(12))
		write (handle,*) dreal(ztot), -1.0*dimag(ztot), freq
	end do
	
	close (1)
	
end subroutine model_data

! THIS PART NEEDS TO BE UPDATED WHEN THE EQUIVALENT CIRCUIT IS MODIFIED < updated for wagner model on 21 jan >
program estimates
	use init
	implicit none
	
	real(kind=8),dimension(n_estimate,n_parameter) :: all_data
	real(kind=8),dimension(n_parameter) :: mean , best
	real(kind=8),dimension(n_parameter-2) :: weight
	integer :: i , j, handle
	
	best = 0
	
	open (unit=100,file="all_trials.txt",status='old',action='read')
	
	do i=1,n_estimate
		read (100,*) all_data(i,1), all_data(i,2), all_data(i,3), all_data(i,4), all_data(i,5), all_data(i,6), &
						&all_data(i,7), all_data(i,8), all_data(i,9), all_data(i,10), all_data(i,11), all_data(i,12),& 
							&all_data(i,13), all_data(i,14)
		if (i == 1 ) then
			do j=1,n_parameter
				best(j) = all_data(i,j)
			end do
		else if (best(n_parameter)<all_data(i,n_parameter)) then
			do j=1,n_parameter
				best(j) = all_data(i,j)
			end do
		end if		
	end do

	close (100)
	
	mean = 0
	weight = 0
	
	do j=1,n_parameter-2
		do i=1,n_estimate
			mean(j) = mean(j) + (all_data(i,j)*all_data(i,n_parameter))
			weight(j) = weight(j) + all_data(i,n_parameter)
		end do
		mean(j) = mean(j)/weight(j)
	end do
	
	do i=1,n_estimate
		mean(n_parameter-1) = mean(n_parameter-1) + all_data(i,n_parameter-1)
		mean(n_parameter) = mean(n_parameter) + all_data(i,n_parameter)
	end do
	mean(n_parameter-1) = mean(n_parameter-1)/n_estimate
	mean(n_parameter) = mean(n_parameter)/n_estimate
	
	handle = 200
	open (unit=handle,file="nyquist_best.txt")
	
	call model_data(best, handle)

	close (handle)
	
	handle = 300
	open (unit=handle,file="nyquist_mean.txt")
	
	call model_data(mean, handle)
	
	close (handle)
	
	handle = 400
	open (unit=handle,file="final_estimates.txt")
	
	write (handle,*) "Mean estimates are :"
	write (handle,*) mean
	write (handle,*) "Best estimates are :"
	write (handle,*) best
	
	close (handle)

end program estimates
