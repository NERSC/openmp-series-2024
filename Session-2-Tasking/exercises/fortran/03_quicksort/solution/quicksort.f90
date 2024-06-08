MODULE quick
!$ USE OMP_LIB

use, intrinsic :: ISO_C_BINDING, only: C_INT

IMPLICIT NONE

interface 
  integer (C_INT) function rand()  BIND(C, name='rand')  
        use, intrinsic :: ISO_C_BINDING  
        implicit none 
  end function rand
end interface 



CONTAINS

SUBROUTINE init (array,  length)
  INTEGER, DIMENSION(:) :: array
  INTEGER,  INTENT(IN)    :: length
  INTEGER :: i
  !
  DO i = 1,  length
    array(i) = MODULO( int(rand()), 1234) + 1; 
!     WRITE (*,*) array(i)
  END DO
END SUBROUTINE init



SUBROUTINE swap(a, b)
  INTEGER :: a, b, temp
  !
  temp=a;
  a=b;
  b=temp;
END SUBROUTINE swap


FUNCTION pivot(array, first, last) RESULT (p)!int. funktion
  INTEGER, DIMENSION(:) :: array
  INTEGER :: first, last, p, pivotElement, i
  !
  p = first;
  pivotElement = array(first);
  DO i = first+1, last
    IF ( array(i) <=  pivotElement ) THEN 
      p = p +1 ;
      CALL swap(array(i),array(p));
    END IF
  END DO
  CALL swap(array(p),array(first));
END FUNCTION pivot


RECURSIVE SUBROUTINE serial_quicksort(array, first, last)
  INTEGER, DIMENSION(:) :: array
  INTEGER :: first, last, pivotElement
  !
  IF (first < last) THEN
    pivotElement = pivot(array,first,last);
    CALL serial_quicksort(array,first,pivotElement-1);
    CALL serial_quicksort(array,pivotElement+1,last);
  END IF
END SUBROUTINE serial_quicksort

RECURSIVE SUBROUTINE quicksort(array, first, last)
  INTEGER, DIMENSION(:) :: array
  INTEGER :: first, last, pivotElement
  !
  IF ((last - first + 1) < 10000) THEN
    CALL serial_quicksort(array, first, last);
  ELSE
    pivotElement = pivot(array,first,last);
    !$omp task default(shared)
    CALL quicksort(array,first,pivotElement-1);
    !$omp end task
    !$omp task default(shared)
    CALL quicksort(array,pivotElement+1,last);
    !$omp end task
    !$omp taskwait
  END IF
END SUBROUTINE quicksort


SUBROUTINE checkFn(array,length) ! bool
  INTEGER, DIMENSION(:) :: array
  INTEGER :: i, length
  !
  DO i=1, length - 1
    IF (array(i) > array(i+1)) THEN
      WRITE (*,*)  "ERR! array[" , i , "] = ",array(i), " > array[" , i+1 , "] = ", array(i+1)
    END IF
  END DO
END SUBROUTINE checkFn

END MODULE quick


PROGRAM main
  !$ USE OMP_LIB
  USE quick
  IMPLICIT NONE
  INTEGER :: length
  DOUBLE PRECISION :: t1,t2 = 0.0
  CHARACTER(LEN=256) :: argument
  INTEGER, DIMENSION(:), ALLOCATABLE :: toBeSorted
  !
  IF ( iargc () .EQ. 0 ) THEN
    length=6000000;
  ELSE
    CALL GETARG(1, argument); READ (argument, *) length
  END IF
  !
  ALLOCATE(toBeSorted(length))
  CALL init(toBeSorted, length);
  !        
  WRITE (*,*)  "Sorting an array of " , length , " elements."
  !$ t1=omp_get_wtime();
  !$omp parallel
  !$omp single
    CALL quicksort(toBeSorted,1,length);
  !$omp end single
  !$omp end parallel
  !$ t2=omp_get_wtime()-t1;
  WRITE (*,*)  "quicksort took " , t2 , " sec. to complete" 
  !
  CALL checkFn(toBeSorted, length);
  DEALLOCATE(toBeSorted);
END PROGRAM main
