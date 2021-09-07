 module routines
    use global
    use functions
    use schemes
    implicit none
    contains
    !================================================================
    subroutine Scheme_routine
        !Scheme evaluation routine
        double precision, allocatable:: S_copy(:,:),schemes(:,:),Sline(:)
        integer :: sop,pos
        integer :: op
        pos = 10; sop = 1
        write(*,*)'============================================================'
        write(*,*)'                                                            '
        write(*,*)'                 Scheme Evaluation Routine                  '
        write(*,*)'                                                            '
        write(*,*)'============================================================'
        allocate(S_copy(N,N),S_exact(N,N),S_calc(N,N),x(N),y(N),schemes(N,9),Sline(N))
        call deltas()   !calculates dx, dy, dt
        call exact_functions(S_exact,S_calc)
        schemes = 0.0d0; Sline = 0.0d0
        schemes(1:N,1) = x(1:N)
        call row_col(S_exact,Sline,pos,sop)
        schemes(1:N,2) = Sline(1:N)
        do op=1,7
            S_copy = S_calc !Copy the elements from S_calc to S_copy
            call ADI(S_copy,op)
            call row_col(S_copy,Sline,pos,sop)
            schemes(1:N,op + 2) = Sline(1:N)
            if(file_flag)then
                name = trim('output/')//trim(OPshem)//trim('_function_')//trim(OPfunc)//trim('.dat')
                call w_fille(S_copy,name,20)
            end if
        end do
        if(file_flag)then
            name = trim('output/')//trim('Exact_solution_')//trim(OPfunc)//trim('.dat')
            call w_fille(S_exact,name,20)
        end if
        name = trim('output/')//trim('Sol_func_')//trim(OPfunc)//trim('.dat')
        call w_fille(schemes,name,20) 
        deallocate(S_copy,S_exact,S_calc,x,y,schemes,Sline)
    end subroutine Scheme_routine
    !================================================================
    subroutine row_col(Sol,Sline,pos,op)
        !Function to extract the solution line/row on the plot
        integer :: op,pos,i
        double precision :: Sol(:,:),Sline(N)
        if(op == 1)then
            do i=1,N
                Sline(i) = Sol(i,pos)
            end do
        else
            do i=1,N
                Sline(i) = Sol(pos,i)
            end do
        end if  
    end subroutine row_col
    !================================================================
    subroutine ERROR_routine
        !Calculates the RMS error
        integer, allocatable, dimension(:):: rms_v
        double precision, allocatable, dimension(:,:):: RMS_e, S_copy
        integer :: i,op,N_a
        N_a = N_F/IJ
        allocate(rms_v(N_a),RMS_e(N_a,8))
        rms_v = (/(i, i=N_i,N_F, IJ)/) + 1
        write(*,*)'============================================================'
        write(*,*)'                                                            '
        write(*,*)'                    Linear Test Problem                     '
        write(*,*)'                      RMS Calculation                       '
        write(*,*)'                                                            '
        write(*,*)'============================================================'
        do i = 1,N_a
            N = rms_v(i)
            write(*,*) 'Grid:',N
            allocate(S_copy(N,N),S_exact(N,N),S_calc(N,N),x(N),y(N))
            call deltas()   !calculates dx, dy, dt
            call exact_functions(S_exact,S_calc)
            RMS_e(i,1) = x(2) - x(1)
            do op = 1,7
                S_copy = S_calc !Copy the elements from S_calc to S_copy
                call ADI(S_copy,op)
                call RMS(RMS_e(i,op+1),S_copy,S_exact)
            end do
            deallocate(S_copy,S_exact,S_calc,x,y)
        end do
        name = trim('output/')//trim('RMS_function_') // trim(OPfunc)//trim('.dat')
        call w_fille(RMS_e,name,20)
        deallocate(rms_v,RMS_e)
    end subroutine ERROR_routine
    !================================================================
    subroutine Eigenvalue_routine
        !Calculates the eigenvalues error
        double precision, allocatable :: RMS_e(:,:), S_copy(:,:),Lambda(:)
        integer :: i,op,N_a
        if(lambda_i<=0.1d0)then
            N_a = int(lambda_f/lambda_j)+1
            allocate(Lambda(N_a),RMS_e(N_a,8))
            Lambda(1) = lambda_i
            Lambda(2) = lambda_j
            do i = 3, N_a
                Lambda(i) = Lambda(i-1) + lambda_j
            end do
        else
            N_a = int((lambda_f-lambda_i)/lambda_j)+1
            allocate(Lambda(N_a),RMS_e(N_a,8))
            Lambda(1) = lambda_i
            do i = 2, N_a
                Lambda(i) = Lambda(i-1) + lambda_j
            end do
        end if
        do i = 1, N_a
            RMS_e(i,1) = Lambda(i)/Pe
        end do
        write(*,*)'============================================================'
        write(*,*)'                                                            '
        write(*,*)'                    Linear Test Problem                     '
        write(*,*)'                   Eigenvalue Calculation                   '
        write(*,*)'                                                            '
        write(*,*)'============================================================'
        allocate(x(N),y(N))
        call deltas()   !calculates dx, dy, dt
        do i = 1,N_a
            allocate(S_copy(N,N),S_exact(N,N),S_calc(N,N))
            lamb = Lambda(i) !Lambda
            write(*,*) 'Lambda:',lamb
            write(*,*) 'Lambda/Pe:',lamb/Pe
            call exact_functions(S_exact,S_calc)
            do op = 1,7
                S_copy = S_calc !Copy the elements from S_calc to S_copy
                call ADI(S_copy,op)
                call RMS(RMS_e(i,op+1),S_copy,S_exact)
            end do
            deallocate(S_copy,S_exact,S_calc)
        end do
        deallocate(x,y)
        name = trim('output/')//trim('Eigenvalue_function_')// trim(OPfunc)//trim('.dat')
        call w_fille(RMS_e,name,20)
        deallocate(Lambda,RMS_e)
    end subroutine Eigenvalue_routine
    !================================================================
    subroutine Angle_routine
        !Calculates the angle error
        double precision, allocatable :: RMS_e(:,:),thetaV(:), S_copy(:,:)
        integer :: i,op,N_a
        N_a = int(angle_f/angle_j)+1
        allocate(thetaV(N_a),RMS_e(N_a,8))
        thetaV(1) = angle_i
        do i = 2, N_a
            thetaV(i) = thetaV(i-1) + angle_j
        end do
        do i = 1, N_a
            RMS_e(i,1) = thetaV(i)
        end do
        write(*,*)'============================================================'
        write(*,*)'                                                            '
        write(*,*)'                    Linear Test Problem                     '
        write(*,*)'                   Angle Calculation                        '
        write(*,*)'                                                            '
        write(*,*)'============================================================'
        allocate(x(N),y(N))
        do i = 1,N_a
            allocate(S_copy(N,N),S_exact(N,N),S_calc(N,N))
            theta = thetaV(i) !angle
            write(*,*) 'Angle:',theta
            call deltas()   !calculates dx, dy, dt
            call exact_functions(S_exact,S_calc)
            do op = 1,7
                S_copy = S_calc !Copy the elements from S_calc to S_copy
                call ADI(S_copy,op)
                call RMS(RMS_e(i,op+1),S_copy,S_exact)
            end do
            deallocate(S_copy,S_exact,S_calc)
        end do
        deallocate(x,y)
        name = trim('output/')//trim('Angle_function_')// trim(OPfunc)//trim('.dat')
        call w_fille(RMS_e,name,20)
        deallocate(thetaV,RMS_e)
    end subroutine Angle_routine
    !================================================================
    subroutine Pe_routine
        !Calculates the Pe error
        double precision, allocatable :: RMS_e(:,:),PeV(:), S_copy(:,:)
        integer :: i,op,N_a
        N_a = int(Pe_f/Pe_j)+1
        allocate(PeV(N_a),RMS_e(N_a,8))
        PeV(1) = Pe_i
        do i = 2, N_a
            PeV(i) = PeV(i-1) + Pe_j
        end do
        do i = 1, N_a
            RMS_e(i,1) = PeV(i)
        end do
        write(*,*)'============================================================'
        write(*,*)'                                                            '
        write(*,*)'                    Linear Test Problem                     '
        write(*,*)'                   Peclet Calculation                        '
        write(*,*)'                                                            '
        write(*,*)'============================================================'
        allocate(x(N),y(N))
        do i = 1,N_a
            allocate(S_copy(N,N),S_exact(N,N),S_calc(N,N))
            Pe = PeV(i) !Peclet
            write(*,*) 'Peclet:',Pe
            call deltas()   !calculates dx, dy, dt
            call exact_functions(S_exact,S_calc)
            do op = 1,7
                S_copy = S_calc !Copy the elements from S_calc to S_copy
                call ADI(S_copy,op)
                call RMS(RMS_e(i,op+1),S_copy,S_exact)
            end do
            deallocate(S_copy,S_exact,S_calc)
        end do
        deallocate(x,y)
        name = trim('output/')//trim('Peclet_function_')// trim(OPfunc)//trim('.dat')
        call w_fille(RMS_e,name,20)
        deallocate(PeV,RMS_e)
    end subroutine Pe_routine
    !================================================================
    subroutine HSF_routine
        !Scheme evaluation routine
        double precision, allocatable:: S_copy(:,:),schemes(:,:),Sline(:)
        integer :: sop,pos
        integer :: op
        pos = 39; sop = 1
        write(*,*)'============================================================'
        write(*,*)'                                                            '
        write(*,*)'               Heaviside Step Function Routine              '
        write(*,*)'                                                            '
        write(*,*)'============================================================'
        allocate(S_copy(N,N),S_exact(N,N),S_calc(N,N),x(N),y(N),schemes(N,9),Sline(N))
        call deltas()   !calculates dx, dy, dt
        call HSF_functions(S_exact,S_calc)
        schemes = 0.0d0; Sline = 0.0d0
        schemes(1:N,1) = x(1:N)
        call row_col(S_exact,Sline,pos,sop)
        schemes(1:N,2) = Sline(1:N)
        do op=1,7
            S_copy = S_calc !Copy the elements from S_calc to S_copy
            call ADI(S_copy,op)
            call row_col(S_copy,Sline,pos,sop)
            schemes(1:N,op + 2) = Sline(1:N)
            if(file_flag)then
                name = trim('output/')//trim(OPshem)//trim('_function_HSF')//trim('.dat')
                call w_fille(S_copy,name,20)
            end if
        end do
        if(file_flag)then
            name = trim('output/')//trim('Exact_solution_')//trim(OPfunc)//trim('.dat')
            call w_fille(S_exact,name,20)
        end if
        name = trim('output/')//trim('Sol_func_HSF')//trim('.dat')
        call w_fille(schemes,name,20) 
        deallocate(S_copy,S_exact,S_calc,x,y,schemes,Sline)
    end subroutine HSF_routine
    !================================================================
    !CDS,FOU,EXP,LOADS,UNI,SOU,QUICK
end module routines

