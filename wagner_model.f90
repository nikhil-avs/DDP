! progam to estimate the parameters of wagenr circuit for EIS data <ADD LATER> and validate it.

! ****************************************
! input file : air20.dat
! input file format : real_z (ohm), imag_z (ohm), freq (Hz)
! number of data points : 40
! ****************************************

! ****************************************
! Equivalent circuit : modified Gohr's porous electrode model
! circuit diagram : <ADD LATER>
! objective function : <ADD LATER>
! ****************************************

! ****************************************
! population     : PARAMETER SET [P COLUMNS], SS [SUM OF SQUARES], BC [BOUNDARY CONDITIONS] ; all the individuals considered by the genetic algo ; p+2 columns
! generation_sum : GENERATION NUMBER, SUM OF FITNESS VALUES OF ALL THE INDIVIDUALS IN THE GENERATION ; used to check the improvement in fitness values with increasing generations ; 2 columns
! glookup        : FITNESS, INDIVIDAUL NUMMBER, CUMULATIVE SUM OF FITNESS ARRANGED WITH INCREASING FITNESS, PERCENTAGE OF GENERATION FITNESS CONTRIBUTED BY AN INDIVIDUAL ; useful for roulette wheel selection; 4 columns
! plookup        : 
! lookup         :
! residuals      :
! levenberg      :
! nyquist_data   :
! nyquist_lm     :
! ****************************************

! Uses hybrid algorithm
! Stage 1 : genetic algo
! Stage 2 : levenberg marquardart algorithm
! process repeated over several trials and top 'n_guess' number of individuals are selected as initial guess to lm algorithm
!

module global
	implicit none
	
	integer, parameter :: n_gen = 100, g = 100, n_trial = 100			! number of individuals in the population ; maximum number of total generations including seed generation as generation 1
	integer, parameter :: n_data = 40, n_param = 12, n_guess = 10, n_est = 2	! number of data points ; number of parameters in the model ; number of guess values for lm algo
	integer, parameter :: noutput = 500
	integer, parameter :: dpr = SELECTED_REAL_KIND(12, 60)
	
	character (len=9), parameter :: file_input = 'air20.dat'				! enter the name of the input file here
	
! ****************************************	ENTER LIMITS HERE ****************************************

	real (kind=8), dimension(3), parameter :: ra_limit = (/ 0.0, 0.5, 0.5 /)
	real (kind=8), dimension(3), parameter :: ca_limit = (/ 0.0, 20.0, 20.0 /)
	real (kind=8), dimension(3), parameter :: phia_limit = (/ 0.5, 1.0, 0.5 /)
	real (kind=8), dimension(3), parameter :: rpa_limit = (/ 0.0, 5.0, 5.0 /)
	real (kind=8), dimension(3), parameter :: rel_limit = (/ 0.0, 0.5, 0.5 /)
	real (kind=8), dimension(3), parameter :: rc_limit = (/ 0.0, 0.5, 0.5 /)
	real (kind=8), dimension(3), parameter :: cc_limit = (/ 0.0, 20.0, 20.0 /)
	real (kind=8), dimension(3), parameter :: phic_limit = (/ 0.5, 1.0, 0.5 /)
	real (kind=8), dimension(3), parameter :: rpc_limit = (/ 0.0, 10.0, 10.0 /)
	real (kind=8), dimension(3), parameter :: l_limit = (/ 0.0, 50.0, 50.0 /)
	real (kind=8), dimension(3), parameter :: rn_limit = (/ 0.0, 50.0, 50.0 /)
	real (kind=8), dimension(3), parameter :: cn_limit = (/ 0.0, 20.0, 20.0 /)
	
	
	real (kind=8), dimension(3), parameter :: ra_lm_limit = (/ 0.0, 0.5, 0.5 /)
	real (kind=8), dimension(3), parameter :: ca_lm_limit = (/ 0.0, 200.0, 200.0 /)
	real (kind=8), dimension(3), parameter :: phia_lm_limit = (/ 0.5, 1.0, 0.5 /)
	real (kind=8), dimension(3), parameter :: rpa_lm_limit = (/ 0.0, 5.0, 5.0 /)
	real (kind=8), dimension(3), parameter :: rel_lm_limit = (/ 0.0, 0.5, 0.5 /)
	real (kind=8), dimension(3), parameter :: rc_lm_limit = (/ 0.0, 0.5, 0.5 /)
	real (kind=8), dimension(3), parameter :: cc_lm_limit = (/ 0.0, 200.0, 200.0 /)
	real (kind=8), dimension(3), parameter :: phic_lm_limit = (/ 0.5, 1.0, 0.5 /)
	real (kind=8), dimension(3), parameter :: rpc_lm_limit = (/ 0.0, 10.0, 10.0 /)
	real (kind=8), dimension(3), parameter :: l_lm_limit = (/ 0.0, 50.0, 50.0 /)
	real (kind=8), dimension(3), parameter :: rn_lm_limit = (/ 0.0, 50.0, 50.0 /)
	real (kind=8), dimension(3), parameter :: cn_lm_limit = (/ 0.0, 20.0, 20.0 /)
	
! *********************************************************************************************************

	real (kind=8), parameter :: pi = 4.0*atan(1.0)
	
	REAL (kind=8), SAVE  :: dat(2*n_data), freq(n_data)
	
	type individual													! inividual's data type defined - generalised for any equivalent circuit
		real (kind=8), dimension(n_param) :: p
		real (kind=8) :: ss, fit
		integer :: bc, n_indv, f_rank
	end type individual
	
	type (individual), dimension(n_gen*g) :: populations			! should not be the same name as progam name
	type (individual), dimension(n_guess*g) :: lm_param
	real (kind=8), dimension(n_guess*g,2) :: plookup
	real (kind=8), dimension(n_gen*(g+1),4) :: vlookup, glookup
	real (kind=8), dimension(g) :: gsum
	integer :: curr_ind, curr_gen, curr_par				! the latest generation that is being made or will be made  ;; the number of individuals at any time = curr_ind-1 (or) curr_ind is the number for next child
	
	integer (kind=4), dimension(8) :: values   			! variables related to system time used for random number initialization
	integer (kind=4) :: rand_init
	character(len=10), dimension(3) :: b			    ! length of the character array also has to be defined, if running on GNR
	
	character(len=17), dimension(n_trial) :: file_pop 		! variables used for creating output files for all the trials
	character(len=13), dimension(n_trial) :: file_lookup
	character(len=14), dimension(n_trial) :: file_glookup
	character(len=14), dimension(n_trial) :: file_plookup
	character(len=21), dimension(n_trial) :: file_gen
	character(len=19), dimension(n_trial) :: file_data
	character(len=17), dimension(n_trial) :: file_data_lm	
	character(len=16), dimension(n_trial) :: file_lm
	character(len=3) :: char_t

end module global

! THIS FUNCTION NEEDS TO BE UPDATED WHEN THE EQUIVALENT CIRCUIT IS MODIFIED < updated for wagner circuit : 20 jan >
subroutine lm_func()
	use Levenberg_Marquardt								! module in lm.f90
	use global
	implicit none
	
	integer :: m, n, info, iwa(n_param)
	!integer :: nprob, nfev, njev
	REAL (dpr) :: x(n_param),fvec(2*n_data), tol		! REVISIT THE WORKING OF LM METHOD AND LM.F90
	complex (dpr) :: wagner_imp, cmplx_temp				! function calculating the impedance should be initialized here
	integer :: i, j, k, check, check_bounds_final
	real(kind=8) :: sums, ss
	
	INTERFACE
		SUBROUTINE fcn(m, n, x, fvec, iflag)
			IMPLICIT NONE
			INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
			INTEGER, INTENT(IN)        :: m, n
			REAL (dp), INTENT(IN)      :: x(:)
			REAL (dp), INTENT(IN OUT)  :: fvec(:)
			INTEGER, INTENT(IN OUT)    :: iflag
		END SUBROUTINE fcn
	END INTERFACE
	
	m = 2*n_data						! number of data points
	n = n_param							! number of parameters
	tol = 0.0000001

	write (noutput,*) "generation : ", curr_gen

	do i=1,n_guess
	
		k = glookup((n_gen*(curr_gen-1))+i, 2) 			! get the INDIVIDUAL NUMBER of the best 'n_guess' number of individuals from each generation 
		init_guess: do j=1,n_param						! get the parameter set of the best individuals
					
			x(j) = populations(k)%p(j)	
					
		end do init_guess
		
		ss = sums(x)									! calculate the sum of squares of the selected parameter set
		
		write (noutput,*) " ", "guess values :"
		write (noutput,*) " ", x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), x(9), x(10), x(11), x(12), ss		
		
		CALL lmdif1(fcn, m, n, x, fvec, tol, info, iwa)		! call the function in lm.f90
		
		param_store: do j=1,n_param			! output of the LM method is stored in variable x()
					
			lm_param(curr_par)%p(j) = x(j) 
					
		end do param_store
		
		lm_param(curr_par)%ss = sums(x)
		lm_param(curr_par)%bc = check_bounds_final(x)
		if(lm_param(curr_par)%bc == 0) then
			lm_param(curr_par)%fit = 1.0/(1.0 + lm_param(curr_par)%ss)
		else
			lm_param(curr_par)%fit = 0
		end if
		
		plookup(curr_par,1) = lm_param(curr_par)%fit
		plookup(curr_par,2) = curr_par
		
		write (noutput,*) " ", "final values :"		
		write (noutput,*) " ", "Termination criteria : ", info
		write (noutput,*) " ", x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), x(9), x(10), x(11), x(12), &
								&lm_param(curr_par)%ss, lm_param(curr_par)%bc
		write (noutput,*) " "
		
		curr_par = curr_par + 1
		
	end do
	
end subroutine lm_func

! THIS FUNCTION NEEDS TO BE UPDATED WHEN THE EQUIVALENT CIRCUIT IS MODIFIED < updated for wagner circuit : 20 jan >
SUBROUTINE fcn(m, n, x, fvec, iflag)
	USE global
	IMPLICIT NONE

	REAL (dpr), INTENT(IN)   :: x(:)
	INTEGER, INTENT(IN)      :: m, n
	REAL (dpr), INTENT(OUT)  :: fvec(:)
	INTEGER, INTENT(IN OUT)  :: iflag
	
	complex (dpr)		  :: wagner_imp, cmplx_temp
	INTEGER               :: i
	REAL(kind=8)          :: w
	
!     EVALUATE THE RESIDUALS [ fvec(i) ] HERE
	
	DO  i=1,n_data
	
		w = 2.0*pi*freq(i) 
		
		cmplx_temp = wagner_imp(w, x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), x(9), x(10), x(11), x(12))	! impedance value at parameter set x() at given w
		
		fvec(i) = (dat(i) - dreal(cmplx_temp))/abs(cmplx_temp)		! residual in real part -> real part of data - real part of model impedance value at given w
		
		fvec(i+n_data) = -1.0_dpr*(dat(i+n_data) + dimag(cmplx_temp))/abs(cmplx_temp)		
		
	END DO

	RETURN

END SUBROUTINE fcn

! NO UPDATE REQUIRED < last updated 20 jan >
subroutine load_data()
	use global
	implicit none
	
	integer :: i
	
	open (unit=1, file=file_input, status='old', action='read')
	
	do i=1,n_data
		read (1,*) dat(i), dat(i+n_data), freq(i)
	end do
	
	close(1)	
	
end subroutine load_data

! NO UPDATE REQUIRED < revisit the logic >
subroutine gen_sum()
	use global
	implicit none
	
	integer :: vi
		
	gsum(curr_gen) = 0.0
	do vi = curr_ind - n_gen , curr_ind - 1
		gsum(curr_gen) = gsum(curr_gen) + glookup(vi,1)
	end do
			
end subroutine gen_sum

! NO UPDATE REQUIRED < revisit the logic >
subroutine gsort()
	use global
	implicit none
	
	integer :: vi, vj
	real (kind=8) :: temp_fit, temp_ind
		
	do vi = curr_ind - n_gen , curr_ind-2
		do vj = vi+1,curr_ind-1
			if (glookup(vj,1) > glookup(vi,1)) then
				temp_fit = glookup(vi,1)
				temp_ind = glookup(vi,2)
				glookup(vi,1) = glookup(vj,1)
				glookup(vi,2) = glookup(vj,2)
				glookup(vj,1) = temp_fit
				glookup(vj,2) = temp_ind
			end if
		end do
	end do
	
	glookup(curr_ind - n_gen , 3) = gsum(curr_gen)
	glookup(curr_ind - n_gen , 4) = (glookup(1,1)/gsum(curr_gen))*100
	do vi = curr_ind - n_gen + 1,curr_ind-1
		glookup(vi,3) = glookup(vi-1,3) - glookup(vi-1,1)
		glookup(vi,4) = (glookup(vi,1)/gsum(curr_gen))*100
	end do
			
end subroutine gsort

! NO UPDATE REQUIRED < revisit the logic >
subroutine vsort(fsum)
	use global
	implicit none

	real (kind=8) :: fsum
	
	integer :: vi, vj
	real (kind=8) :: temp_fit, temp_ind
		
	do vi = 1,curr_ind-2
		do vj = vi+1,curr_ind-1
			if (vlookup(vj,1) > vlookup(vi,1)) then
				temp_fit = vlookup(vi,1)
				temp_ind = vlookup(vi,2)
				vlookup(vi,1) = vlookup(vj,1)
				vlookup(vi,2) = vlookup(vj,2)
				vlookup(vj,1) = temp_fit
				vlookup(vj,2) = temp_ind
			end if
		end do
	end do
	fsum = 0.0
	do vi = 1,curr_ind-1
		fsum = fsum + vlookup(vi,1)
	end do
	
	vlookup(1,3) = fsum
	vlookup(1,4) = (vlookup(1,1)/fsum)*100
	do vi = 2,curr_ind-1
		vlookup(vi,3) = vlookup(vi-1,3) - vlookup(vi-1,1)
		vlookup(vi,4) = (vlookup(vi,1)/fsum)*100
	end do
			
end subroutine vsort

! NO UPDATE REQUIRED < revisit the logic >
subroutine psort()
	use global
	implicit none
	
	integer :: vi, vj
	real (kind=8) :: temp_fit, temp_ind
		
	do vi = 1,(n_guess*g)-1
		do vj = vi+1,n_guess*g
			if (plookup(vj,1) > plookup(vi,1)) then
				temp_fit = plookup(vi,1)
				temp_ind = plookup(vi,2)
				plookup(vi,1) = plookup(vj,1)
				plookup(vi,2) = plookup(vj,2)
				plookup(vj,1) = temp_fit
				plookup(vj,2) = temp_ind
			end if
		end do
	end do

end subroutine psort

! THIS FUNCTION NEEDS TO BE UPDATED WHEN THE EQUIVALENT CIRCUIT IS MODIFIED < updated for wagner circuit : 20 jan >
subroutine model_data(p, t)
	use global
	implicit none
	
	real (kind=8), intent(in), dimension(n_param) :: p
	integer, intent(in) :: t
	real (kind=8) :: w, phiw
	complex (kind=8) :: wagner_imp, ztot
	integer :: i
	
	open (unit=600, file=file_data(t))
		
	do i = 1,n_data
		w = 2.0*pi*freq(i)
		ztot = wagner_imp(w, p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8), p(9), p(10), p(11), p(12))
		write (600,*) dreal(ztot), -1.0*dimag(ztot), freq(i)
	end do
	
	close(600)
	
end subroutine model_data

! THIS FUNCTION NEEDS TO BE UPDATED WHEN THE EQUIVALENT CIRCUIT IS MODIFIED < updated for wagner circuit : 20 jan >
subroutine model_data_lm(p, t)
	use global
	implicit none
	
	real (kind=8), intent(in), dimension(n_param) :: p
	integer, intent(in) :: t
	real (kind=8) :: w, phiw
	complex (kind=8) :: wagner_imp, ztot
	integer :: i
	
	phiw = 0.5
	
	open (unit=700, file=file_data_lm(t))
		
	do i = 1,n_data
		w = 2.0*pi*freq(i)
		ztot = wagner_imp(w, p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8), p(9), p(10), p(11), p(12))									! CAN BE GENERALISED 
		write (700,*) dreal(ztot), -1.0*dimag(ztot), freq(i)
	end do
	
	close(700)
	
end subroutine model_data_lm

complex (kind=8) function dhirde_imp(w, r1, r2, c1, phi1, rw, tw, r3, c2, phi2, phiw)

	implicit none
	real (kind=8), intent(in) :: w, r1, r2, c1, phi1, rw, tw, r3, c2, phi2, phiw
	complex (kind=8) :: iw, zr1, zr2, zc1, zr3, zc2, zw, ztot, cmplx_temp, cmplx_tanh
	
	iw = cmplx(0.0,w)
	zr1 = cmplx(r1)
	zr2 = cmplx(r2)
	zr3 = cmplx(r3)
	
	zc1 = c1/(iw**phi1)
	
	zc2 = c2/(iw**phi2)
	
	cmplx_temp = (iw*tw)**phiw

	cmplx_tanh = (exp(2.0*cmplx_temp)-1.0)/(exp(2.0*cmplx_temp)+1.0)
	
	zw = rw*cmplx_tanh/cmplx_temp
	
	ztot = zr1 + (zr3*zc2/(zr3+zc2)) + ((zr2+zw)*zc1/(zr2+zw+zc1))
	
	dhirde_imp = ztot	

end function dhirde_imp

complex (kind=8) function wagner_imp(w, ra, ca, phia, rpa, rel, rc, cc, phic, rpc, l, rn, cn)

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
	
	wagner_imp = ztot	

end function wagner_imp

! THIS FUNCTION NEEDS TO BE UPDATED WHEN THE EQUIVALENT CIRCUIT IS MODIFIED < updated for wagner circuit : 20 jan >
real (kind=8) function sums(p)
	use global
	implicit none
	
	real (kind=8), intent(in), dimension(n_param) :: p
	real (kind=8) :: w, r_real, r_imag ,ss_temp, phiw, z_abs
	complex (kind=8) :: wagner_imp, ztot				! the return type of function should also be mentioned in the function from which u are calling the other function
	integer :: i									! counter variable
	
	if (curr_ind == 1) then
		open (unit=289, file='residuals.txt')
	end if
	
	phiw = 0.5
	ss_temp = 0.0
	i=1
	
	do i=1,n_data
	
		w = 2.0*pi*freq(i)
		ztot = wagner_imp(w, p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8), p(9), p(10), p(11), p(12))									! CAN BE GENERALISED 
		z_abs = abs(ztot)
		r_real = ((dat(i) - dreal(ztot))**2)/(z_abs**2)			! real residual is also verified
		r_imag = ((dat(i+n_data) + dimag(ztot))**2)/(z_abs**2)			! corrected :)
		ss_temp = ss_temp + r_real + r_imag
		if (curr_ind == 1) then
			write (289,*) dat(i), dat(i+n_data), dreal(ztot), dimag(ztot), r_real, r_imag, freq(i)
		end if
		
	end do
	ss_temp = ss_temp/(1.0*n_data)
	
	if (curr_ind == 1) close(289)
	
	sums = ss_temp

end function sums

! NO UPDATE REQUIRED < revisit the logic >
real (kind=8) function mutation(param)
	real (kind=8), intent(in) :: param
	real (kind=8) :: new_param
	
	if (rand() > 0.5 ) then
		new_param = param + (0.05*param)
	else
		new_param = param - (0.05*param)
	end if
	
	mutation = new_param
end function mutation

! NO UPDATE REQUIRED < revisit the logic >
real (kind=8) function line1(param1, param2)
	real (kind=8), intent(in) :: param1, param2
	real (kind=8) :: new_param

	if ( param1 > param2) then
		new_param = param2 + (param1 - param2)/2
	else
		new_param = param1 + (param2 - param1)/2
	end if
	
	line1= new_param
end function line1

! NO UPDATE REQUIRED < revisit the logic >
real (kind=8) function line2(param1, param2)
	real (kind=8), intent(in) :: param1, param2
	real (kind=8) :: new_param

	if (rand() > 0.5) then
		if ( param1 > param2) then
			new_param = param1 + (param1 - param2)/2
		else
			new_param = param2 + (param2 - param1)/2
		end if
	else
		if ( param1 > param2) then
			new_param = param2 - (param1 - param2)/2
		else
			new_param = param1 - (param2 - param1)/2
		end if
	end if

	line2 = new_param
end function line2

! NO UPDATE REQUIRED < revisit the logic >
real (kind=8) function minimum(ssq, ssq_size, ind)
	implicit none
	real (kind=8), dimension(ssq_size), intent(in) :: ssq
	integer, intent(in) :: ssq_size, ind
	real (kind=8) :: ssq_temp
	integer :: i
	
	ssq_temp = ssq(1)
	do i = 2,ind
		if (ssq_temp > ssq(i))	ssq_temp = ssq(i)
	end do
	
	minimum = ssq_temp
	
end function minimum

! THIS FUNCTION NEEDS TO BE UPDATED WHEN THE EQUIVALENT CIRCUIT IS MODIFIED < updated for wagner circuit : 20 jan >
integer function check_bounds(p)
	use global
	implicit none
	real (kind=8), dimension(n_param), intent(in) :: p
	integer :: broken
	
	broken = 0
	
	if (p(1) < ra_limit(1) .or. p(1) > ra_limit(2)) 		broken = broken + 1															! CAN BE GENERALISED 
	if (p(2) < ca_limit(1) .or. p(2) > ca_limit(2)) 		broken = broken + 1
	if (p(3) < phia_limit(1) .or. p(3) > phia_limit(2)) 		broken = broken + 1
	if (p(4) < rpa_limit(1) .or. p(4) > rpa_limit(2)) 		broken = broken + 1
	if (p(5) < rel_limit(1) .or. p(5) > rel_limit(2)) 		broken = broken + 1
	if (p(6) < rc_limit(1) .or. p(6) > rc_limit(2))		broken = broken + 1
	if (p(7) < cc_limit(1) .or. p(7) > cc_limit(2))		broken = broken + 1
	if (p(8) < phic_limit(1) .or. p(8) > phic_limit(2))		broken = broken + 1
	if (p(9) < rpc_limit(1) .or. p(9) > rpc_limit(2))		broken = broken + 1
	if (p(10) < l_limit(1) .or. p(10) > l_limit(2))		broken = broken + 1
	if (p(11) < rn_limit(1) .or. p(11) > rn_limit(2))		broken = broken + 1
	if (p(12) < cn_limit(1) .or. p(12) > cn_limit(2))		broken = broken + 1
			
	check_bounds = broken	
	
end function check_bounds

! THIS FUNCTION NEEDS TO BE UPDATED WHEN THE EQUIVALENT CIRCUIT IS MODIFIED < updated for wagner circuit : 20 jan >
integer function check_bounds_final(p)
	use global
	implicit none
	real (kind=8), dimension(n_param), intent(in) :: p
	integer :: broken
	
	broken = 0
	
	if (p(1) < ra_lm_limit(1) .or. p(1) > ra_lm_limit(2)) 		broken = broken + 1
	if (p(2) < ca_lm_limit(1) .or. p(2) > ca_lm_limit(2)) 		broken = broken + 1														! CAN BE GENERALISED 
	if (p(3) < phia_lm_limit(1) .or. p(3) > phia_lm_limit(2)) 		broken = broken + 1
	if (p(4) < rpa_lm_limit(1) .or. p(4) > rpa_lm_limit(2)) 		broken = broken + 1
	if (p(5) < rel_lm_limit(1) .or. p(5) > rel_lm_limit(2)) 		broken = broken + 1
	if (p(6) < rc_lm_limit(1) .or. p(6) > rc_lm_limit(2))		broken = broken + 1
	if (p(7) < cc_lm_limit(1) .or. p(7) > cc_lm_limit(2))		broken = broken + 1
	if (p(8) < phic_lm_limit(1) .or. p(8) > phic_lm_limit(2))		broken = broken + 1
	if (p(9) < rpc_lm_limit(1) .or. p(9) > rpc_lm_limit(2))		broken = broken + 1
	if (p(10) < l_lm_limit(1) .or. p(10) > l_lm_limit(2))		broken = broken + 1
	if (p(11) < rn_lm_limit(1) .or. p(11) > rn_lm_limit(2))		broken = broken + 1
	if (p(12) < cn_lm_limit(1) .or. p(12) > cn_lm_limit(2))		broken = broken + 1
		
	check_bounds_final = broken	
	
end function check_bounds_final

! NO UPDATE REQUIRED < revisit the logic >
integer function roulette()
	use global
	implicit none
	
	real (kind=8) :: temp_rand
	integer :: temp_select
	integer :: i
	
	temp_rand = rand()*gsum(curr_gen-1)
	
	i=1
	rou_select1: do i = (n_gen*(curr_gen-2)) + 1 , n_gen*(curr_gen-1)
		if (i == n_gen*(curr_gen-1)) then
			temp_select = int(glookup(i,2))
		else if (temp_rand <= glookup(i,3) .and. temp_rand >= glookup(i+1,3)) then
			temp_select = int(glookup(i,2))
			exit rou_select1
		end if
	end do rou_select1
		
	roulette = temp_select
	
end function roulette

! THIS PART NEEDS TO BE UPDATED WHEN THE EQUIVALENT CIRCUIT IS MODIFIED < not yet updated >
program wagner_model	
	use global
	implicit none
		
	integer :: i, j, t, temp_var1, ssq_ind, lm_ind				
	real (kind=8) :: temp_rand, ssq_temp, fsum	
	integer, dimension(2) :: rselect
	
	! it is necessary to mention the data type of a function before using it
	real (kind=8) :: mutation, line1, line2, sums, minimum		! type declaration for functions
	integer :: check_bounds, roulette							! type declaration for functions
	complex (kind=8) :: wagner_imp, ztot				! type declaration for functions
	
	call load_data()	
	
	open (unit=900,file="all_trials.txt")
	
	trials: do t=1,n_trial
	
	print *, "trial : ", t
	
	write(char_t, 1250) t
	1250 format(i0.3)
	
	! to add the path for the output files append "file_path//" to the file_pop etc..

	file_pop(t) = 'population'//char_t//'.txt'
	file_lookup(t) = 'lookup'//char_t//'.txt'
	file_glookup(t) = 'glookup'//char_t//'.txt'
	file_plookup(t) = 'plookup'//char_t//'.txt'
	file_gen(t) = 'generation_sum'//char_t//'.txt'
	file_data(t) = 'nyquist_data'//char_t//'.txt'
	file_data_lm(t) = 'nyquist_lm'//char_t//'.txt'
	file_lm(t)   = 'levenberg'//char_t//'.txt'
	
	open (unit=noutput,file=file_lm(t))
	
	call date_and_time(b(1), b(2), b(3), values)
	rand_init = values(8)*values(7)
	
	call srand(rand_init)
	
	ssq_temp = 100000.0
	ssq_ind = 1
	curr_gen = 1
	curr_ind = 1	
	curr_par = 1
	i=1
	
	seeds: do i=1,n_gen
		
		populations(curr_ind)%p(1) = ra_limit(1)+(rand()*ra_limit(3))																		! CAN BE GENERALISED 
		populations(curr_ind)%p(2) = ca_limit(1)+(rand()*ca_limit(3))
		populations(curr_ind)%p(3) = phia_limit(1)+(rand()*phia_limit(3))
		populations(curr_ind)%p(4) = rpa_limit(1)+(rand()*rpa_limit(3))
		populations(curr_ind)%p(5) = rel_limit(1)+(rand()*rel_limit(3))
		populations(curr_ind)%p(6) = rc_limit(1)+(rand()*rc_limit(3))
		populations(curr_ind)%p(7) = cc_limit(1)+(rand()*cc_limit(3))
		populations(curr_ind)%p(8) = phic_limit(1)+(rand()*phic_limit(3))
		populations(curr_ind)%p(9) = rpc_limit(1)+(rand()*rpc_limit(3))
		populations(curr_ind)%p(10) = l_limit(1)+(rand()*l_limit(3))
		populations(curr_ind)%p(11) = rn_limit(1)+(rand()*rn_limit(3))
		populations(curr_ind)%p(12) = cn_limit(1)+(rand()*cn_limit(3))
		populations(curr_ind)%bc = check_bounds(populations(curr_ind)%p)
		populations(curr_ind)%ss = sums(populations(curr_ind)%p)
		if ( populations(curr_ind)%bc == 0 ) then
			populations(curr_ind)%fit = 1.0/(1.0 + populations(curr_ind)%ss)	
		else
			populations(curr_ind)%fit = 0.0
		end if
		if (ssq_temp > populations(curr_ind)%ss .and. populations(curr_ind)%bc == 0) then
			ssq_temp = populations(curr_ind)%ss
			ssq_ind = curr_ind
		end if
		vlookup(curr_ind,1) = populations(curr_ind)%fit
		vlookup(curr_ind,2) = curr_ind*1.0
		glookup(curr_ind,1) = populations(curr_ind)%fit
		glookup(curr_ind,2) = curr_ind*1.0
		curr_ind = curr_ind + 1
	end do seeds
	call gen_sum()
	call gsort()
	call vsort(fsum)							! sorting and ranking should be done here
	call lm_func()
	curr_gen = curr_gen + 1
	
	generations: do
		print *, "generation : ", curr_gen
		j = 1
		operations: do j = 1,int(0.7*n_gen)					! 40 mutaion operations , 20 cross over operations and 10 line operations per generation
			
			rselect(1) = roulette()
			rselect(2) = roulette()
		
		! ****************************************	MUTATION ****************************************	
			if (j < int(0.4*n_gen)+1) then
				do temp_var1 =1,n_param
						populations(curr_ind)%p(temp_var1) = mutation(populations(rselect(1))%p(temp_var1))		
				end do
				populations(curr_ind)%bc = check_bounds(populations(curr_ind)%p)
				populations(curr_ind)%ss = sums(populations(curr_ind)%p)	
				if ( populations(curr_ind)%bc == 0 ) then
					populations(curr_ind)%fit = 1.0/(1.0 + populations(curr_ind)%ss)	
				else
					populations(curr_ind)%fit = 0.0
				end if
				if (ssq_temp > populations(curr_ind)%ss .and. populations(curr_ind)%bc == 0) then
					ssq_temp = populations(curr_ind)%ss
					ssq_ind = curr_ind
				end if
				vlookup(curr_ind,1) = populations(curr_ind)%fit
				vlookup(curr_ind,2) = curr_ind*1.0
				glookup(curr_ind,1) = populations(curr_ind)%fit
				glookup(curr_ind,2) = curr_ind*1.0
				curr_ind = curr_ind + 1
		! *********************************************************************************************************		
		
		! ****************************************	CROSSOVER ****************************************		
			else if (j > int(0.4*n_gen) .and. temp_rand < int(0.6*n_gen)+1) then
				do temp_var1 = 1,n_param-5																		!crossover operation
					populations(curr_ind)%p(temp_var1) = populations(rselect(1))%p(temp_var1)
				end do
				do temp_var1 = n_param-4,n_param
					populations(curr_ind)%p(temp_var1) = populations(rselect(2))%p(temp_var1)
				end do
				populations(curr_ind)%bc = check_bounds(populations(curr_ind)%p)
				populations(curr_ind)%ss = sums(populations(curr_ind)%p)
				if ( populations(curr_ind)%bc == 0 ) then
					populations(curr_ind)%fit = 1.0/(1.0 + populations(curr_ind)%ss)	
				else
					populations(curr_ind)%fit = 0.0
				end if
				if (ssq_temp > populations(curr_ind)%ss .and. populations(curr_ind)%bc == 0) then
					ssq_temp = populations(curr_ind)%ss
					ssq_ind = curr_ind
				end if
				vlookup(curr_ind,1) = populations(curr_ind)%fit
				vlookup(curr_ind,2) = curr_ind*1.0
				glookup(curr_ind,1) = populations(curr_ind)%fit
				glookup(curr_ind,2) = curr_ind*1.0
				curr_ind = curr_ind + 1
				do temp_var1 = 1,n_param-5
					populations(curr_ind)%p(temp_var1) = populations(rselect(2))%p(temp_var1)
				end do
				do temp_var1 = n_param-4,n_param
					populations(curr_ind)%p(temp_var1) = populations(rselect(1))%p(temp_var1)
				end do
				populations(curr_ind)%bc = check_bounds(populations(curr_ind)%p)
				populations(curr_ind)%ss = sums(populations(curr_ind)%p)	
				if ( populations(curr_ind)%bc == 0 ) then
					populations(curr_ind)%fit = 1.0/(1.0 + populations(curr_ind)%ss)	
				else
					populations(curr_ind)%fit = 0.0
				end if
				if (ssq_temp > populations(curr_ind)%ss .and. populations(curr_ind)%bc == 0) then
					ssq_temp = populations(curr_ind)%ss
					ssq_ind = curr_ind
				end if
				vlookup(curr_ind,1) = populations(curr_ind)%fit
				vlookup(curr_ind,2) = curr_ind*1.0
				glookup(curr_ind,1) = populations(curr_ind)%fit
				glookup(curr_ind,2) = curr_ind*1.0
				curr_ind = curr_ind + 1	
		! *********************************************************************************************************		
		
		! ****************************************	LINE  ****************************************
			else
				do temp_var1 = 1,n_param
					populations(curr_ind)%p(temp_var1) = line1(populations(rselect(1))%p(temp_var1),&
								&populations(rselect(2))%p(temp_var1))
				end do
				populations(curr_ind)%bc = check_bounds(populations(curr_ind)%p)
				populations(curr_ind)%ss = sums(populations(curr_ind)%p)	
				if ( populations(curr_ind)%bc == 0 ) then
					populations(curr_ind)%fit = 1.0/(1.0 + populations(curr_ind)%ss)	
				else
					populations(curr_ind)%fit = 0.0
				end if
				if (ssq_temp > populations(curr_ind)%ss .and. populations(curr_ind)%bc == 0) then
					ssq_temp = populations(curr_ind)%ss
					ssq_ind = curr_ind
				end if
				vlookup(curr_ind,1) = populations(curr_ind)%fit
				vlookup(curr_ind,2) = curr_ind*1.0
				glookup(curr_ind,1) = populations(curr_ind)%fit
				glookup(curr_ind,2) = curr_ind*1.0
				curr_ind = curr_ind + 1
				do temp_var1 = 1,n_param
					populations(curr_ind)%p(temp_var1) = line2(populations(rselect(1))%p(temp_var1),&
										&populations(rselect(2))%p(temp_var1))
				end do
				populations(curr_ind)%bc = check_bounds(populations(curr_ind)%p)
				populations(curr_ind)%ss = sums(populations(curr_ind)%p)	
				if ( populations(curr_ind)%bc == 0 ) then
					populations(curr_ind)%fit = 1.0/(1.0 + populations(curr_ind)%ss)	
				else
					populations(curr_ind)%fit = 0.0
				end if
				if (ssq_temp > populations(curr_ind)%ss .and. populations(curr_ind)%bc == 0) then
					ssq_temp = populations(curr_ind)%ss
					ssq_ind = curr_ind
				end if
				vlookup(curr_ind,1) = populations(curr_ind)%fit
				vlookup(curr_ind,2) = curr_ind*1.0
				glookup(curr_ind,1) = populations(curr_ind)%fit
				glookup(curr_ind,2) = curr_ind*1.0
				curr_ind = curr_ind + 1
			end if
		! *********************************************************************************************************
			
		end do operations
		call gen_sum()				!calling gen_sum before vsort is very important
		call gsort()
		call vsort(fsum)
		call lm_func()
		curr_gen = curr_gen + 1
		
		if (curr_gen > g) exit generations
	
	end do generations
	
	call psort()	
	lm_ind = plookup(1,2)
	
	! ******************************************************************** PRINTING STARTS ************************************

	call model_data(populations(ssq_ind)%p, t)
	
	call model_data_lm(lm_param(lm_ind)%p, t)
	
	open (unit=100, file=file_pop(t))
	open (unit=200,file=file_lookup(t))
	open (unit=300,file=file_glookup(t))
	open (unit=400,file=file_gen(t))
	open (unit=800,file=file_plookup(t))
	
	i=1
	do i=1,(curr_ind-1)
		write (100,*) populations(i)%p(1), populations(i)%p(2), populations(i)%p(3), populations(i)%p(4)&								! CAN BE GENERALISED 
			&, populations(i)%p(5), populations(i)%p(6), populations(i)%p(7), populations(i)%p(8),&
				&populations(i)%p(9), populations(i)%p(10), populations(i)%p(11), populations(i)%p(12), &
				  &populations(i)%ss, populations(i)%bc	
		write (200,*) vlookup(i,1), vlookup(i,2), vlookup(i,3), vlookup(i,4)
		write (300,*) glookup(i,1), glookup(i,2), glookup(i,3), glookup(i,4)
	end do
		
	do i=1,n_guess*g
		write (800,*) plookup(i,1), plookup(i,2)
	end do
	
	do i=1,g
		write (400,*) i, gsum(i)
	end do
	
	do i=1,n_est
		lm_ind = plookup(i,2)
		write (900,*) lm_param(lm_ind)%p(1),lm_param(lm_ind)%p(2),lm_param(lm_ind)%p(3),lm_param(lm_ind)%p(4), &							! CAN BE GENERALISED 
					&lm_param(lm_ind)%p(5), lm_param(lm_ind)%p(6), lm_param(lm_ind)%p(7), lm_param(lm_ind)%p(8),&
						&lm_param(lm_ind)%p(9), lm_param(lm_ind)%p(10), lm_param(lm_ind)%p(11), lm_param(lm_ind)%p(12), & 
						  &lm_param(lm_ind)%ss, lm_param(lm_ind)%fit
	end do
	
	close(100)
	close(200)
	close(300)
	close(400)
	close(noutput)
	close(800)
	
		
	print *, curr_ind, ssq_ind, ssq_temp
	
	print *, populations(ssq_ind)%p(1), populations(ssq_ind)%p(2), populations(ssq_ind)%p(3), populations(ssq_ind)%p(4)&
			&, populations(ssq_ind)%p(5), populations(ssq_ind)%p(6), populations(ssq_ind)%p(7),&
				&populations(ssq_ind)%p(8), populations(ssq_ind)%p(9), populations(ssq_ind)%p(10), populations(ssq_ind)%p(11), & 
				  &populations(ssq_ind)%p(12), populations(ssq_ind)%ss
			
	print *, ' '
	
	! ******************************************************************* PRINTING ENDS ***************************************

	end do trials
	
	close(900)
	
end program wagner_model


