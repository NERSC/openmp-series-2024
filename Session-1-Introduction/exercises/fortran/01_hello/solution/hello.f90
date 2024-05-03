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

program main
use omp_lib
implicit none

!$omp parallel
write (*,*) 'Hello from thread ', omp_get_thread_num(), ' of ', omp_get_num_threads()
!$omp end parallel

end program main
