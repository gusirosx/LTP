program LTP
    use routines
    implicit none
    
    call input_data
    if(rms_flag) call ERROR_routine
    
    if(sheme_flag) call Scheme_routine
    
    if(eigenvalue_flag) call Eigenvalue_routine
    
    if(angle_flag) call Angle_routine
    
    if(Pe_flag) call Pe_routine
    
    if(HSF_flag) call HSF_routine
    
end program LTP
