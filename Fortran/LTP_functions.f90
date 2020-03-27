module functions
    use global
    use schemes
    implicit none
    contains
    !================================================================
    subroutine norm_function(S)
        !Normalization by the max value of a 2D array
        double precision, dimension(:,:):: S
        integer :: i,j
        real :: F
        F = maxval(S) - minval(S)
        S = S/F
    end subroutine norm_function
    !================================================================
    subroutine deltas()
        integer :: i
        alpha = theta*pi/180d0
        dx = 1.0d0/(N-1.0d0) 
        dy = dx
        !dt = 0.98d0*dx*dx
        dt = 0.98*(4*dx)/(3*Pe*cos(alpha))
        !Determine the x and y positions for solutions A, B, C, D
        x = 0.0d0; y = 0.0d0
        do i = 1, N
            x(i) = -0.5 + (i-1.0d0)*dx
            y(i) = -0.5 + (i-1.0d0)*dy
        end do
    end subroutine deltas
    !================================================================
    subroutine TDMA(a,b,c,d,S,n)
        integer i,n
        double precision,dimension(:):: a,b,c,d,S
        double precision,dimension(n):: cp,dp
        double precision :: m
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
        do i = 2,n
            m = b(i)-cp(i-1)*a(i)
            cp(i) = c(i)/m
            dp(i) = (d(i)-dp(i-1)*a(i))/m
        enddo
        S(n) = dp(n)
        do i = n-1, 1, -1
            S(i) = dp(i)-cp(i)*S(i+1)
        end do
    end subroutine TDMA
    !================================================================
    subroutine ADI(S_sch,op)
        !select which ADI to use
        double precision, allocatable, dimension(:,:):: S_sch
        integer :: op
        if(op.ge.1 .and. op.le.5) then
            call ADI_S(S_sch,op)
        else if(op.ge.6 .and. op.le.7) then 
            call ADI_MOD(S_sch,op)
        end if
    end subroutine ADI
    !================================================================
    subroutine ADI_S(S_sch,op)!Standard ADI
        double precision:: Ax,Bx,Cx,Dpx,Ay,By,Cy,Dpy,VTOL
        double precision, allocatable, dimension(:):: A,B,C,D,TS
        double precision, allocatable, dimension(:,:):: S_a,S_rhs,S_sch
        integer :: i,j,icnt,op
        call scheme_Selection(op,Ax,Bx,Cx,Dpx,Ay,By,Cy,Dpy)
        allocate(S_a(N,N),S_rhs(N,N),A(N-2),B(N-2),C(N-2),D(N-2),TS(N-2))
        S_rhs = 0.0d0
        icnt = 0
        VTOL =0.1d0
        do while (TOL <= VTOL)
            !X Direction
            S_a = S_sch !Copy the elements from S_sch to S_a
            if(op.eq.4 .or. op.eq.5) call RHS_schemes(op,S_a,S_rhs)
            do i = 2,N-1
                do j = 2,N-1
                    if (j.eq.2) then
                        A(j-1) = 0.0d0
                        B(j-1) = Bx
                        C(j-1) = Cx
                        D(j-1) = - Ay*S_a(i-1,j) + Dpx*S_a(i,j) - Cy*S_a(i+1,j) - Ax*S_a(i,j-1) + S_rhs(i,j)
                    else if (j.gt.2 .and. j.lt.N-1) then
                        A(j-1) = Ax
                        B(j-1) = Bx
                        C(j-1) = Cx
                        D(j-1) = - Ay*S_a(i-1,j) + Dpx*S_a(i,j) - Cy*S_a(i+1,j) + S_rhs(i,j)
                    else if (j.eq.N-1) then
                        A(j-1) = Ax
                        B(j-1) = Bx
                        C(j-1) = 0.0d0
                        D(j-1) = - Ay*S_a(i-1,j) + Dpx*S_a(i,j) - Cy*S_a(i+1,j) - Cx*S_a(i,j+1) + S_rhs(i,j)
                    end if
                end do
                call TDMA(A,B,C,D,TS,N-2)
                do j = 2,N-1
                    S_sch(i,j) = TS(j-1) !Stores TDMA values ​​in matrix
                end do
            end do
            !Y Direction
            S_a = S_sch
            if(op.eq.4 .or. op.eq.5) call RHS_schemes(op,S_a,S_rhs)
            do j = 2,N-1
                do i = 2,N-1
                    if (i.eq.2) then
                        A(i-1) = 0.0d0
                        B(i-1) = By
                        C(i-1) = Cy
                        D(i-1) = - Ax*S_a(i,j-1) + Dpy*S_a(i,j) - Cx*S_a(i,j+1) - Ay*S_a(i-1,j) + S_rhs(i,j)
                    else if (i.gt.2 .and. i.lt.N-1) then
                        A(i-1) = Ay
                        B(i-1) = By
                        C(i-1) = Cy
                        D(i-1) = - Ax*S_a(i,j-1) + Dpy*S_a(i,j) - Cx*S_a(i,j+1) + S_rhs(i,j)
                    else if (i.eq.N-1) then
                        A(i-1) = Ay
                        B(i-1) = By
                        C(i-1) = 0.0d0
                        D(i-1) = - Ax*S_a(i,j-1) + Dpy*S_a(i,j) - Cx*S_a(i,j+1) - Cy*S_a(i+1,j) + S_rhs(i,j)
                    end if
                end do
                call TDMA(A,B,C,D,TS,N-2)
                do i = 2,N-1
                    S_sch(i,j) = TS(i-1)
                end do
            end do 
            VTOL = maxval(S_sch - S_a)/dt
            icnt = icnt + 1
            if (icnt.eq.1000000) then
                exit
            end if
        end do
        write(*,*) 'C.R.I:',icnt,trim(OPshem) 
        deallocate(S_a,S_rhs,A,B,C,D,TS)
    end subroutine ADI_S
    !================================================================
    subroutine ADI_MOD(S_sch,op)!Modified ADI for WW and SS terms
        double precision:: Ax1,Bx1,Dpx1,Dpy1,Ay1,By1,Cx,Cy,VTOL, &
                           Ax2,Bx2,Dpx2,Dpy2,Ay2,By2,Yss,Xww
        double precision, allocatable, dimension(:):: A,B,C,D,TS
        double precision, allocatable, dimension(:,:):: S_a,S_sch
        integer :: i,j,icnt,op
        call scheme_Selection(op,Ax1,Bx1,Cx,Dpx1,Ay1,By1,Cy,Dpy1,Ax2,Bx2,Dpx2,Ay2,By2,Dpy2,Xww,Yss)
        allocate(S_a(N,N),A(N-2),B(N-2),C(N-2),D(N-2),TS(N-2))
        icnt = 0
        VTOL =0.1d0
        do while (TOL <= VTOL)
            !X Direction
            S_a = S_sch !Copy the elements from S_sch to S_a
            do i = 2,N-1
                if (i.eq.2) then
                    do j = 2,N-1
                        if (j.eq.2) then
                            A(j-1) = 0.0d0
                            B(j-1) = Bx2
                            C(j-1) = Cx
                            D(j-1) = - Ay2*S_a(i-1,j) + Dpx2*S_a(i,j) - Cy*S_a(i+1,j) - Ax2*S_a(i,j-1)
                        else if (j.gt.2 .and. j.lt.N-1) then
                            A(j-1) = Ax1
                            B(j-1) = Bx1
                            C(j-1) = Cx
                            D(j-1) = - Ay2*S_a(i-1,j) + Dpx2*S_a(i,j) - Cy*S_a(i+1,j) &
                                     + Xww*S_a(i,j-2)
                        else if (j.eq.N-1) then
                            A(j-1) = Ax1
                            B(j-1) = Bx1
                            C(j-1) = 0.0d0
                            D(j-1) = - Ay2*S_a(i-1,j) + Dpx2*S_a(i,j) - Cy*S_a(i+1,j) - Cx*S_a(i,j+1) &
                                     + Xww*S_a(i,j-2)
                        end if
                    end do
                else if (i.gt.2) then
                    do j = 2,N-1
                        if (j.eq.2) then
                            A(j-1) = 0.0d0
                            B(j-1) = Bx2
                            C(j-1) = Cx
                            D(j-1) = - Ay1*S_a(i-1,j) + Dpx1*S_a(i,j) - Cy*S_a(i+1,j) - Ax2*S_a(i,j-1) &
                                     + Yss*S_a(i-2,j)
                        else if (j.gt.2 .and. j.lt.N-1) then
                            A(j-1) = Ax1
                            B(j-1) = Bx1
                            C(j-1) = Cx
                            D(j-1) = - Ay1*S_a(i-1,j) + Dpx1*S_a(i,j) - Cy*S_a(i+1,j) + Xww*S_a(i,j-2) &
                                     + Yss*S_a(i-2,j)
                        else if (j.eq.N-1) then
                            A(j-1) = Ax1
                            B(j-1) = Bx1
                            C(j-1) = 0.0d0
                            D(j-1) = - Ay1*S_a(i-1,j) + Dpx1*S_a(i,j) - Cy*S_a(i+1,j) - Cx*S_a(i,j+1) &
                                     + Xww*S_a(i,j-2) + Yss*S_a(i-2,j)
                        end if
                    end do
                end if
                call TDMA(A,B,C,D,TS,N-2)
                do j = 2,N-1
                    S_sch(i,j) = TS(j-1) !Stores TDMA values ​​in matrix
                end do
            end do
            !Y Direction
            S_a = S_sch
            do j = 2,N-1
                if (j.eq.2) then
                    do i = 2,N-1
                        if (i.eq.2) then
                            A(i-1) = 0.0d0
                            B(i-1) = By2
                            C(i-1) = Cy
                            D(i-1) = - Ax2*S_a(i,j-1) + Dpy2*S_a(i,j) - Cx*S_a(i,j+1) - Ay2*S_a(i-1,j)
                        else if (i.gt.2 .and. i.lt.N-1) then
                            A(i-1) = Ay1
                            B(i-1) = By1
                            C(i-1) = Cy
                            D(i-1) = - Ax2*S_a(i,j-1) + Dpy2*S_a(i,j) - Cx*S_a(i,j+1) &
                                     + Yss*S_a(i-2,j)
                        else if (i.eq.N-1) then
                            A(i-1) = Ay1
                            B(i-1) = By1
                            C(i-1) = 0.0d0
                            D(i-1) = - Ax2*S_a(i,j-1) + Dpy2*S_a(i,j) - Cx*S_a(i,j+1) - Cy*S_a(i+1,j) &
                                     + Yss*S_a(i-2,j)
                        end if
                    end do
                else if (j.gt.2) then
                    do i = 2,N-1
                        if (i.eq.2) then
                            A(i-1) = 0.0d0
                            B(i-1) = By2
                            C(i-1) = Cy
                            D(i-1) = - Ax1*S_a(i,j-1) + Dpy1*S_a(i,j) - Cx*S_a(i,j+1) - Ay2*S_a(i-1,j) &
                                     + Xww*S_a(i,j-2)
                        else if (i.gt.2 .and. i.lt.N-1) then
                            A(i-1) = Ay1
                            B(i-1) = By1
                            C(i-1) = Cy
                            D(i-1) = - Ax1*S_a(i,j-1) + Dpy1*S_a(i,j) - Cx*S_a(i,j+1) &
                                     + Xww*S_a(i,j-2) + Yss*S_a(i-2,j)
                        else if (i.eq.N-1) then
                            A(i-1) = Ay1
                            B(i-1) = By1
                            C(i-1) = 0.0d0
                            D(i-1) = - Ax1*S_a(i,j-1) + Dpy1*S_a(i,j) - Cx*S_a(i,j+1) - Cy*S_a(i+1,j) &
                                     + Xww*S_a(i,j-2) + Yss*S_a(i-2,j)
                        end if
                    end do
                end if
                call TDMA(A,B,C,D,TS,N-2)
                do i = 2,N-1
                    S_sch(i,j) = TS(i-1)
                end do
            end do 
            VTOL = maxval(S_sch - S_a)/dt
            icnt = icnt + 1
            if (icnt.eq.1000000) then
                exit
            end if
        end do
        !write(*,*) trim(OPshem)//trim(' C.R.I:'),icnt 
        write(*,*) 'C.R.I:',icnt,trim(OPshem) 
        deallocate(S_a,A,B,C,D,TS)
    end subroutine ADI_MOD
    !================================================================
    subroutine exact_functions(S_sol,S_c)
        integer :: i,j
        double precision :: C,S,M
        double precision, dimension(:,:):: S_sol,S_c
        !Calculation of the exact solution - TYPE A, B, C, D
        if (OPfunc == 'A') then
            C = (Pe - (Pe**2.0d0 + 4.0d0*lamb**2.0d0)**0.5d0)/2.0d0
        else if (OPfunc == 'B') then
            C = (Pe + (Pe**2.0d0 + 4.0d0*lamb**2.0d0)**0.5d0)/2.0d0
        else if (OPfunc == 'C') then
            C = (Pe - (Pe**2.0d0 - 4.0d0*lamb**2.0d0)**0.5d0)/2.0d0
        else if (OPfunc == 'D') then
            C = (Pe + (Pe**2.0d0 - 4.0d0*lamb**2.0d0)**0.5d0)/2.0d0
        else if (OPfunc == 'CD'.or.OPfunc == 'DC') then
            C = (lamb**2.0d0 - (Pe**2.0d0)*0.25d0)**0.5d0
        end if
        !Loop for functions A and B
        if (OPfunc == 'A'.or.OPfunc == 'B') then
            do i=1,N
                do j=1,N
                    S = x(j)*cos(alpha) + y(i)*sin(alpha)
                    M = y(i)*cos(alpha) - x(j)*sin(alpha)
                    S_sol(i,j) = exp(C*S)*sin(lamb*M)
                end do
            end do 
        !Loop for functions C and D
        else if (OPfunc == 'C'.or.OPfunc == 'D') then
            do i=1,N
                do j=1,N
                    S = x(j)*cos(alpha) + y(i)*sin(alpha)
                    M = y(i)*cos(alpha) - x(j)*sin(alpha)
                    S_sol(i,j) = exp(C*S)*exp(lamb*M)
                end do
            end do 
        !Loop for functions CD
        else if (OPfunc == 'CD') then
            do i=1,N
                do j=1,N
                    S = x(j)*cos(alpha) + y(i)*sin(alpha)
                    M = y(i)*cos(alpha) - x(j)*sin(alpha)
                    S_sol(i,j) = exp(Pe*S*0.5)*sin(C*S)*exp(lamb*M)
                end do
            end do
        !Loop for functions DC
        else if (OPfunc == 'DC') then
            do i=1,N
                do j=1,N
                    S = x(j)*cos(alpha) + y(i)*sin(alpha)
                    M = y(i)*cos(alpha) - x(j)*sin(alpha)
                    S_sol(i,j) = exp(Pe*S*0.5)*cos(C*S)*exp(lamb*M)
                end do
            end do 
        end if
        call norm_function(S_sol)
        !Applies boundary conditions
        S_c = 0.0d0
        do i=1,N
            S_c(1,i) = S_sol(1,i)
            S_c(N,i) = S_sol(N,i)
            S_c(i,1) = S_sol(i,1)
            S_c(i,N) = S_sol(i,N)
        end do 
    end subroutine exact_functions
    !================================================================
    subroutine write_matrix(a)
        double precision, dimension(:,:):: a
        integer :: i,j
        write(*,*)
        do i = lbound(a,1), ubound(a,1)
            write(*,*) (a(i,j), j = lbound(a,2), ubound(a,2))
        end do
    end subroutine write_matrix
    !================================================================
    subroutine w_fille(mat,Name,flagfile)
        double precision, dimension(:,:):: mat
        character(len=*) :: Name
        integer :: flagfile,i,j
        open (unit=flagfile,file=Name,action="write",status="replace")
        do i = lbound(mat,1), ubound(mat,1)
            write(flagfile,*) (mat(i,j), j = lbound(mat,2), ubound(mat,2))
        end do
        close(flagfile)
    end subroutine w_fille
    !================================================================
    subroutine RMS(e_rms,S_aprox,S_exat)
        !Calculating the RMS of a matrix
        double precision, allocatable, dimension(:,:):: S_aprox,S_exat
        double precision :: Ei, e_rms
        integer :: i,j
        e_rms = 0.0d0
        do i = 1,N-1
            do j = 1,N-1
            Ei = (S_aprox(i,j) - S_exat(i,j))**2.0d0
            e_rms = e_rms + Ei
            end do
        end do
        Ei = (e_rms/((N-2)*(N-2)))
        e_rms = sqrt(Ei)
    end subroutine RMS
    !================================================================
    subroutine input_data
        implicit none
        Character(len = 60) :: aux
        integer:: error
        open(unit=60,file="input.ltp",form="formatted",&
                    status="unknown",access = "sequential",&
                    action="read",iostat=error)
       if(error /= 0) then
          print*, "error during opening of input.dat"
          stop
       end if
        read(60,*)aux
        read(60,*)OPfunc
        read(60,*)sheme_flag
        read(60,*)Pe
        read(60,*)lamb
        read(60,*)theta
        read(60,*)TOL
        read(60,*)N
        N = N + 1
        read(60,*)aux
        read(60,*)rms_flag
        read(60,*)N_i
        read(60,*)N_f
        read(60,*)IJ        
        close(60)
    end subroutine input_data
    !================================================================
end module functions

