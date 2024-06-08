program main
	use omp_lib
	implicit none
	integer :: dimensio = 500
	double precision :: result = 0.0
	double precision :: t1, t2
	integer :: i

	t1 = omp_get_wtime()

	!$omp parallel do schedule(dynamic) reduction(+:result)
	do i=1,dimensio
		result = result + do_some_computation(i)
	end do
	!$omp end parallel do

	t2 = omp_get_wtime()

	write(*,*)'Computation took ', t2-t1, ' seconds'
	write(*,*)'Result is ', result

	contains
	function do_some_computation (i) result(res)
		double precision :: res
		integer, intent(in) :: i
		integer :: j
		double precision :: t
		t = 0.0
		do j = 1, i*i
			t = t + sin(j*1.) * cos(j*1.)
		end do
		res = t
	end function do_some_computation
end program main
