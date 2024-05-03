function func(i)
double precision :: func
integer, intent(in) :: i
	func = sin(i * cos(i*1.))
end function func

program main
use omp_lib
implicit none
integer :: dimensio = 9999999 
double precision :: dMin, dMax
double precision, dimension(:), allocatable :: dArray
integer :: i
double precision :: func,t1,t2

allocate (dArray(dimensio))

t1 = omp_get_wtime()


!$omp parallel
!$omp do
	do i=1, dimensio
		dArray(i) = func(i)
		! compute dMin = min(dMin, dArray(i))
		! compute dMax = max(dMax, dArray(i))
	end do
!$omp end parallel

 deallocate (dArray)

t2 = omp_get_wtime()

write (*,*) 'Computation took ', t2-t1, 'seconds.'
write (*,*) 'dMin is ', dMin
write (*,*) 'dMax is ', dMax

end program main

!
!  OpenMP lecture exercises
!  Copyright (C) 2011 by Christian Terboven <terboven@rz.rwth-aachen.de>
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
!
!
