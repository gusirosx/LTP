 module routines
    use global
    use functions
    use schemes
    implicit none
    contains
    !================================================================
    subroutine Scheme_routine
        !Scheme evaluation routine
        double precision, allocatable, dimension(:,:):: S_copy
        integer :: op 
        write(*,*)'============================================================'
        write(*,*)'                                                            '
        write(*,*)'                 Scheme Evaluation Routine                  '
        write(*,*)'                                                            '
        write(*,*)'============================================================'
        allocate(S_copy(N,N),S_exact(N,N),S_calc(N,N),x(N),y(N))
        call deltas()   !calculates dx, dy, dt
        call exact_functions(S_exact,S_calc)
        do op=1,7
            S_copy = S_calc !Copy the elements from S_calc to S_copy
            call ADI(S_copy,op)
            name = trim('output/')//trim(OPshem)//trim('_function_') // trim(OPfunc)//trim('.dat')
            call w_fille(S_copy,name,20)
        end do
        name = trim('output/')//trim('Exact_solution_') // trim(OPfunc)//trim('.dat')
        call w_fille(S_exact,name,20)        
        deallocate(S_copy,S_exact,S_calc,x,y)
    end subroutine Scheme_routine
    !================================================================
    subroutine ERROR_routine
        !Calculates the RMS error
        integer, allocatable, dimension(:):: rms_v
        double precision, allocatable, dimension(:,:):: RMS_e, S_copy
        integer :: i,op,N_a
        N_a = N_F/IJ
        allocate(rms_v(N_a),RMS_e(N_a,7))
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
            do op = 1,7
                S_copy = S_calc !Copy the elements from S_calc to S_copy
                call ADI(S_copy,op)
                call RMS(RMS_e(i,op),S_copy,S_exact)
            end do
            deallocate(S_copy,S_exact,S_calc,x,y)
        end do
        name = trim('output/')//trim('RMS_function_') // trim(OPfunc)//trim('.dat')
        call w_fille(RMS_e,name,20)
        deallocate(rms_v,RMS_e)
    end subroutine ERROR_routine
    !================================================================
end module routines
