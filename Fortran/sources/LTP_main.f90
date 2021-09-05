program LTP
    use routines
    implicit none
    
    call input_data
    if(rms_flag) call ERROR_routine
    
    if(sheme_flag) call Scheme_routine
    
    if(eigenvalue_flag) call Eigenvalue_routine
    
end program LTP
