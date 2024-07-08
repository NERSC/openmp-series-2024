module VariableDef
    implicit none

    TYPE JacobiData
        integer :: iRows
        integer :: iCols
        integer :: iRowFirst
        integer :: iRowLast
        integer :: iIterMax
        double precision :: fAlpha
        double precision :: fRelax
        double precision :: fTolerance
    
!        /* calculated dx & dy */
        double precision :: fDx
        double precision :: fDy

!       /* pointers to the allocated memory */
        double precision, allocatable :: afU(:,:)
        double precision, allocatable :: afF(:,:)
    
!       /* start and end timestamps */
        double precision :: fTimeStart
        double precision :: fTimeStop

!       /* calculated residual (output jacobi) */
        double precision :: fResidual
!       /* effective interation count (output jacobi) */
        integer :: iIterCount

!       /* calculated error (output error_check) */
        double precision :: fError
    
!       /* MPI-Variables */
        integer :: iMyRank   !/* current process rank (number) */
        integer :: iNumProcs !/* how many processes */

   END TYPE JacobiData

end module VariableDef


module JacobiMod
    use VariableDef
    implicit none 

    contains
   
    subroutine Jacobi(myData)
        implicit none
        !********************************************************************
        ! Subroutine HelmholtzJ                                             *
        ! Solves poisson equation on rectangular grid assuming :            *
        ! (1) Uniform discretization in each direction, and                 *
        ! (2) Dirichlect boundary conditions                                *
        ! Jacobi method is used in this routine                             *
        !                                                                   *
        ! Input : n,m   Number of grid points in the X/Y directions         *
        !         dx,dy Grid spacing in the X/Y directions                  *
        !         alpha Helmholtz eqn. coefficient                          *
        !         omega Relaxation factor                                   *
        !         myData%afF(n,m) Right hand side function                  *
        !         myData%afU(n,m) Dependent variable/Solution               *
        !         tol    Tolerance for iterative solver                     *
        !         maxit  Maximum number of iterations                       *
        !                                                                   *
        ! Output : myData%afU(n,m) - Solution                               *
        !******************************************************************** 
  
        !.. Formal Arguments .. 
        type(JacobiData), intent(inout) :: myData 
         
        !.. Local Scalars .. 
        integer :: i, j, iErr
        double precision :: ax, ay, b, residual, fLRes, tmpResd
         
        !.. Local Arrays .. 
        double precision, allocatable :: uold(:,:)
         
        !.. Intrinsic Functions .. 
        intrinsic DBLE, SQRT

        allocate(uold (0 : myData%iCols -1, 0 : myData%iRows -1))

        ! ... Executable Statements ...
        ! Initialize coefficients 
        
        if (allocated(uold)) then    
            ax = 1.0d0 / (myData%fDx * myData%fDx)      ! X-direction coef 
            ay = 1.0d0 / (myData%fDy * myData%fDy)      ! Y-direction coef
            b = -2.0d0 * (ax + ay) - myData%fAlpha      ! Central coeff  
            residual = 10.0d0 * myData%fTolerance
        
            do while (myData%iIterCount < myData%iIterMax .and. residual > myData%fTolerance)
                residual = 0.0d0

                !$omp parallel
                   ! Copy new solution into old
                   !$omp do private(i)
                   do j = 1, myData%iRows - 2
                       do i = 1, myData%iCols - 2
                           uold(i, j) = myData%afU(i, j)
                       end do
                   end do
                   !$omp end do
                   ! Compute stencil, residual, & update
                   !$omp do private(i, fLRes) reduction(+:residual)
                   do j = myData%iRowFirst + 1, myData%iRowLast - 1
                       do i = 1, myData%iCols - 2
                           ! Evaluate residual 
                           fLRes = (ax * (uold(i-1, j) + uold(i+1, j)) &
                                  + ay * (uold(i, j-1) + uold(i, j+1)) &
                                  + b * uold(i, j) - myData%afF(i, j)) / b
                    
                           ! Update solution 
                           myData%afU(i, j) = uold(i, j) - myData%fRelax * fLRes
                    
                           ! Accumulate residual error
                           residual = residual + fLRes * fLRes
                       end do
                   end do
                   !$omp end do
                !$omp end parallel

                 ! Error check 
                 myData%iIterCount = myData%iIterCount + 1      
                 residual = SQRT(residual) / DBLE(myData%iCols * myData%iRows)
             
            ! End iteration loop 
            end do
            myData%fResidual = residual
            deallocate(uold)
        else
           write (*,*) 'Error: cant allocate memory'
           call Finish(myData)
           stop
        end if
    end subroutine Jacobi

end module JacobiMod


program MAIN
    !***********************************************************************
    ! program to solve a finite difference                                 * 
    ! discretization of Helmholtz equation :                               *
    ! (d2/dx2)u + (d2/dy2)u - alpha u = f                                  *
    ! using Jacobi iterative method.                                       *
    !                                                                      *
    ! Modified: Abdelali Malih,    Aachen University (RWTH), 2007          *
    ! Modified: Sanjiv Shah,       Kuck and Associates, Inc. (KAI), 1998   *
    ! Author  : Joseph Robicheaux, Kuck and Associates, Inc. (KAI), 1998   *
    !                                                                      * 
    ! Directives are used in this code to achieve paralleism.              *
    ! All do loops are parallized with default 'static' scheduling.        *
    !                                                                      * 
    ! Input :  n - grid dimension in x direction                           *
    !          m - grid dimension in y direction                           *
    !          alpha - Helmholtz constant (always greater than 0.0)        *
    !          tol   - error tolerance for iterative solver                *
    !          relax - Successice over relaxation parameter                *
    !          mits  - Maximum iterations for iterative solver             *
    !                                                                      *
    ! On output                                                            *
    !       : u(n,m) - Dependent variable (solutions)                      *
    !       : f(n,m) - Right hand side function                            *
    !***********************************************************************
  
    use VariableDef
    use JacobiMod
    use omp_lib
    implicit none
  
    TYPE(JacobiData) :: myData


!   sets default values or reads from stdin
!    * inits MPI and OpenMP if needed
!    * distribute MPI data, calculate MPI bounds
!    */
    call Init(mydata)

    if ( allocated(myData%afU) .and. allocated(myData%afF) ) then
!        /* matrix init */
        call InitializeMatrix(myData)

!        /* starting timer */
        mydata%fTimeStart = omp_get_wtime()

!        /* running calculations */
        call Jacobi(myData)

!        /* stopping timer */
        mydata%fTimeStop = omp_get_wtime()

!        /* error checking */
        call CheckError(myData)

!        /* print result summary */
        call PrintResults(myData)
    else
        write (*,*) " Memory allocation failed ...\n"
    end if

!    /* cleanup */
    call Finish(myData)
    
end program MAIN

subroutine Init (myData)
    use VariableDef
    implicit none
    type(JacobiData), intent(inout) :: myData
    integer :: iErr, i
!/* default medium */
        myData%iCols      = 2000
        myData%iRows      = 2000
        myData%fAlpha     = 0.8
        myData%fRelax     = 1.0
        myData%fTolerance = 1e-10
        myData%iIterMax   = 50
#ifdef READ_INPUT
        write (*,*) 'Input n - matrix size in x direction: '
        read (5,*) myData%iCols
        write (*,*) 'Input m - matrix size in y direction: '
        read (5,*) myData%iRows
        write (*,*) 'Input alpha - Helmholts constant:'
        read (5,*) myData%fAlpha
        write (*,*) 'Input relax - Successive over-relaxation parameter:'
        read (5,*) myData%fRelax
        write (*,*) 'Input tol - error tolerance for iterative solver:'
        read (5,*) myData%fTolerance
        write (*,*) 'Input mits - Maximum iterations for solver:'
        read (5,*) myData%iIterMax
#elif defined DATA_LARGE
        myData%iCols      = 7000
        myData%iRows      = 7000
        myData%fAlpha     = 0.8
        myData%fRelax     = 1.0
        myData%fTolerance = 1e-12
        myData%iIterMax   = 2

#elif defined DATA_SMALL
        myData%iCols      = 200
        myData%iRows      = 200
        myData%fAlpha     = 0.8
        myData%fRelax     = 1.0
        myData%fTolerance = 1e-7
        myData%iIterMax   = 1000
#endif
        write (*,327) "-> matrix size:", myData%iCols, myData%iRows
        write (*,329) "-> alpha: " , myData%fAlpha
        write (*,329) "-> relax: ", myData%fRelax
        write (*,329) "-> tolerance: ", myData%fTolerance
        write (*,328) "-> #of iterations: ", myData%iIterMax
327     format (A22, I10, ' x ', I10)
328     format (A22, I10)
329     format (A22, F10.6)


!    /* MPI values, set to defaults to avoid data inconsistency */
    myData%iMyRank   = 0
    myData%iNumProcs = 1
    myData%iRowFirst = 0
    myData%iRowLast  = myData%iRows - 1

!    /* memory allocation for serial & omp */
    allocate(myData%afU (0 : myData%iCols -1, 0 : myData%iRows -1))
    allocate(myData%afF (0 : myData%iCols -1, 0 : myData%iRows -1))

!    /* calculate dx and dy */
    myData%fDx = 2.0d0 / DBLE(myData%iCols - 1)
    myData%fDy = 2.0d0 / DBLE(myData%iRows - 1)

    myData%iIterCount = 0

end subroutine Init

subroutine InitializeMatrix (myData)
    !*********************************************************************
    ! Initializes data                                                   *
    ! Assumes exact solution is u(x,y) = (1-x^2)*(1-y^2)                 *
    !                                                                    *
    !*********************************************************************
    use VariableDef
    implicit none

    type(JacobiData), intent(inout) :: myData 
    !.. Local Scalars .. 
    integer :: i, j, xx, yy
    !.. Intrinsic Functions .. 
    intrinsic DBLE
   
    ! Initilize initial condition and RHS
  
    !$omp parallel do private(i, xx, yy)
    do j = myData%iRowFirst, myData%iRowLast
        do i = 0, myData%iCols -1
            xx = INT(-1.0 + myData%fDx*DBLE(i)) ! -1 < x < 1
            yy = INT(-1.0 + myData%fDy*DBLE(j)) ! -1 < y < 1
            myData%afU(i, j) = 0.0d0
            myData%afF(i, j) = - myData%fAlpha * (1.0d0 - DBLE(xx*xx))  &
                * (1.0d0 - DBLE(yy*yy)) - 2.0d0 * (1.0d0 - DBLE(xx*xx)) &
                - 2.0d0 * (1.0d0 - DBLE(yy*yy))
        end do
    end do
    !$omp end parallel do
end subroutine InitializeMatrix

subroutine Finish(myData)
    use VariableDef
    implicit none

    integer :: iErr
    type(JacobiData), intent(inout) :: myData

    deallocate (myData%afU)
    deallocate (myData%afF)

end subroutine Finish

subroutine PrintResults(myData)
    use VariableDef
    implicit none

    type(JacobiData), intent(inout) :: myData

    if (myData%iMyRank == 0) then
        write (*,328) " Number of iterations : ", myData%iIterCount
        write (*,*)   " Residual             : ", myData%fResidual
        write (*,329) " Solution Error       : ", myData%fError
        write (*,329) " Elapsed Time         : ", myData%fTimeStop - &
                                                  myData%fTimeStart
        write (*,329) " MFlops               : ", 0.000013 * DBLE (myData%iIterCount &
               * (myData%iCols - 2) * (myData%iRows - 2))                       &
               / (myData%fTimeStop - myData%fTimeStart)
328     format (A22, I10)
329     format (A22, F15.10)
    end if
end subroutine PrintResults


subroutine CheckError(myData)
    use VariableDef
    implicit none
    
    type(JacobiData), intent(inout) :: myData
    !.. Local Scalars .. 
    integer :: i, j, iErr
    double precision :: error, temp, xx, yy
    !.. Intrinsic Functions ..
    intrinsic DBLE, SQRT
    ! ... Executable Statements ...
    error = 0.0d0

    !$omp parallel do private(i, xx, yy, temp) reduction(+:error)
    do j = myData%iRowFirst, myData%iRowLast
        do i = 0, myData%iCols -1
            xx = -1.0d0 + myData%fDx * DBLE(i)
            yy = -1.0d0 + myData%fDy * DBLE(j)
            temp = myData%afU(i, j) - (1.0d0-xx*xx)*(1.0d0-yy*yy)
            error = error + temp*temp
        end do
    end do
    !$omp end parallel do

    myData%fError = sqrt(error) / DBLE(myData%iCols * myData%iRows)
   
end subroutine CheckError

