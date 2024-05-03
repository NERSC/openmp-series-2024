/*
 *  OpenMP lecture exercises
 *  Copyright (C) 2011 by Christian Terboven <terboven@rz.rwth-aachen.de>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */

#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


double func(int i)
{
	return sin((double)i) * cos((double)i);
}


int main(int argc, char* argv[])
{
	const int dimension = 9999999;
	int i;
	double dArray[dimension];
	double dMin = 0.0, dMax = 0.0;
	
	double t1 = omp_get_wtime();

	#pragma omp parallel
	{

		#pragma omp for
		for (i = 0; i < dimension; i++)
		{
			dArray[i] = func(i);
			// compute dMin = fmin(dMin, dArray[i])
			// compute dMax = fmax(dMax, dArray[i])
		}

	} // end omp parallel
	
	double t2 = omp_get_wtime();
	printf("Computation took %.3lf seconds.\n", t2 - t1);
	printf("dMin is %.3lf.\n", dMin);
	printf("dMax is %.3lf.\n", dMax);
	
	return 0;
}
