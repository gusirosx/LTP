module global
    implicit none
    ! PARAMETERS
    double precision, allocatable, dimension(:,:):: &
    S_exact,            & !Exact solution
    S_calc                !Calculated solution
    double precision, allocatable, dimension(:):: &
    x,                  & !Determine the x positions
    y                     !Determine the y positions
    double precision :: &
    lamb,               & !Eigenvalue 
    theta,              & !Mesh offset angle
    alpha,              &
    dx,                 &
    dy,                 &
    dt,                 &
    TOL,                & !Convergence tolerance
    Pe,                 & !Global Peclet Number
    pi = 3.14159265359d0  !Pi number

    integer ::          &
    N                     !Nodes in x and y direction                        
    Character(len=60):: &
    OPfunc,             & !Select the Function
    OPshem,             & !Select the scheme 
    name
    logical ::          & 
    rms_flag,           & !rms routine flag
    sheme_flag,         & !scheme flag
    angle_flag,         & !scheme flag
    Pe_flag,            & !Pe flag
    file_flag,          &
    HSF_flag,           & !HSF flag
    eigenvalue_flag       !eigenvalue flag
    !====== Error Evaluation Routine ======!
    integer ::          & 
    N_i,                & !Initial mesh
    N_f,                & !Final mesh
    IJ                    !Jump Interval
    !==== EIGENVALUE ERROR EVALUATION ROUTINE ====!
    double precision :: &
    lambda_i,           &!Initial lambda
    lambda_f,           &!Final lambda
    lambda_j             !Jump Interval
    !==== ANGLE ERROR EVALUATION ROUTINE ====!
    double precision :: &
    angle_i,            &!Initial angle
    angle_f,            &!Final angle
    angle_j             !Jump Interval
    !==== PECLET ERROR EVALUATION ROUTINE ====!
    double precision :: &
    Pe_i,               &!Initial angle
    Pe_f,               &!Final angle
    Pe_j                 !Jump Interval

end module global
