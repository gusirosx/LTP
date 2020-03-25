 module routines
    use global
    use functions
    use schemes
    implicit none
    contains
    !================================================================
    subroutine Scheme_routine
        !Scheme evaluation routine
        integer :: i,op
        write(*,*)'============================================================'
        write(*,*)'                                                            '
        write(*,*)'                 Scheme Evaluation Routine                  '
        write(*,*)'                                                            '
        write(*,*)'============================================================'
        allocate(S_exact(N,N),S_calc(N,N),x(N),y(N))
        call deltas()   !calculates dx, dy, dt
        do op=1,7
            call exact_functions(S_exact,S_calc)
            call ADI(S_calc,op)
            name = trim(OPshem)//trim('_function_') // trim(OPfunc)//trim('.txt')
            call w_fille(S_calc,name,20)
        end do
        deallocate(S_exact,S_calc,x,y)
    end subroutine Scheme_routine
    !================================================================
    subroutine ERROR_routine
        !Calculates the RMS error
        integer, allocatable, dimension(:):: rms_v
        double precision, allocatable, dimension(:,:):: RMS_e
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
        do op=1,7
            do i=1,N_a
                N = rms_v(i)
                write(*,*) 'Grid',N
                allocate(S_exact(N,N),S_calc(N,N),x(N),y(N))
                call deltas()   !calculates dx, dy, dt
                call exact_functions(S_exact,S_calc)
                call ADI(S_calc,op)
                call RMS(RMS_e(i,op),S_calc,S_exact)
                deallocate(S_exact,S_calc,x,y)
            end do
        end do
        name = trim('RMS_function_') // trim(OPfunc)//trim('.txt')
        call w_fille(RMS_e,name,20)
        deallocate(rms_v,RMS_e)
    end subroutine ERROR_routine
    !================================================================
end module routines
