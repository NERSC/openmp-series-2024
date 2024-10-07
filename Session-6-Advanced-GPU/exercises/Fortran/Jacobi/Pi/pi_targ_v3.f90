          program main
          use omp_lib
          implicit none

          integer :: i, id, block_size, ii
          integer, parameter :: num_steps = 100000000
          real*8 :: x, pi, sum, step
          real*8 :: start_time, run_time

          sum = 0.0

          step = 1.0 / num_steps

          start_time = omp_get_wtime()

!$omp target map(sum)
!$omp teams reduction(+:sum)
  block_size = num_steps/omp_get_num_teams()
!$omp distribute dist_schedule(static, 1) 
  do ii = 0, num_steps-1, block_size
!$omp parallel do private(x) reduction(+:sum)
     do i = ii+1, MIN(ii+block_size, num_steps)
        x = (i - 0.5) * step
        sum = sum + 4.0 / (1.0 + x * x)
     end do
!$omp end parallel do
  end do
!$omp end distribute
!$omp end teams
!$omp end target

          pi = step * sum
          run_time = omp_get_wtime() - start_time
          write(*,100) pi, run_time
100       format('pi is ',f15.8,' in ',f8.3,' secs')

          end program main
