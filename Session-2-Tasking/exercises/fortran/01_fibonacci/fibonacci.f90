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

recursive integer function fib(n) RESULT(res)
integer, intent(in) :: n
integer :: x, y
if (n < 2) then
	res = n
else	

	x = fib(n-1)

	y = fib(n-2)

	res = x + y
end if
end function fib

program main

use omp_lib
integer :: n, fibonacci
integer :: fib
double precision :: starttime
write (*,*) 'Please insert n and press enter, to calculate fib(n)'
read (*,*) n
starttime = omp_get_wtime()

fibonacci = fib(n)

write (*,*) 'fib(', n, ' ) = ', fibonacci

write (*,*) 'calculation took', omp_get_wtime() - starttime, 'sec'

end program main

