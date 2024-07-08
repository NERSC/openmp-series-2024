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

function f(a)
double precision :: f
double precision, intent(in) :: a 
f = 4.0/(1.0 + a*a)
end function f 


function CalcPi(n)
double precision :: CalcPi
integer, intent(in) :: n
double precision :: fH, fSum, fX, f
integer :: i
fH = 1.0 / n
fSum = 0.0

!$omp simd private(fX) reduction(+:fSum)
do i=0,n
	fX = fH * ( 0.5 + i )
	fSum = fSum + f(fX)	
end do
!$omp end simd
 CalcPi = fH * fSum
end function CalcPi



program main
use omp_lib
implicit none
integer :: n = 150000000
double precision :: fPi25DT = 3.141592653589793238462643
double precision :: fPi, calcpi
double precision :: fTimeStart, fTimeEnd


if ( n <= 0 .or. n > 2147483647) then
	write(*,*) 'given value has to be between 0 and 2147483647'
	stop
end if

fTimeStart = omp_get_wtime()

! The calculation is done here
fPi = calcpi(n)

fTimeEnd = omp_get_wtime()

 write (*,100) fPi, abs(fPi - fPi25DT)
 100 format ('pi is approximately = ', e15.8 ,' , Error =' , e15.8)
 write (*,*)'wall clock time	=', fTimeEnd-fTimeStart
 
 end program main

