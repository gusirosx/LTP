module schemes
    use global
    implicit none
    contains
    !================================================================
    subroutine scheme_Selection(op,ax,bx,cx,dpx,ay,by,cy,dpy,axm,bxm,dpxm,aym,bym,dpym,xww,yss)
        !Calculates the coefficients for the selected scheme
        double precision :: ax,bx,cx,dpx,ay,by,cy,dpy
        double precision,optional :: axm,bxm,dpxm,aym,bym,dpym,xww,yss
        integer :: op
        if(op.eq.1) then 
            call CDS(ax,bx,cx,dpx,ay,by,cy,dpy)
        else if(op.eq.2) then 
            call FOU(ax,bx,cx,dpx,ay,by,cy,dpy)
        else if(op.ge.3 .and. op.le.5) then 
            call EXPO(ax,bx,cx,dpx,ay,by,cy,dpy)
        else if(op.eq.6) then 
            call SOU(ax,bx,dpx,ay,by,dpy,axm,bxm,dpxm,aym,bym,dpym,cx,cy,xww,yss)
        else if(op.eq.7) then 
            call QUICK(ax,bx,dpx,ay,by,dpy,axm,bxm,dpxm,aym,bym,dpym,cx,cy,xww,yss)
        end if
    end subroutine scheme_Selection
    !================================================================
    subroutine RHS_schemes(op,Taux,RHS)
        !Calculates RHS of LOADS or UNIFAES
        double precision, dimension(:,:) :: Taux,RHS
        integer :: op
        if(op.eq.4) then
            call LOADS_RHS(Taux,RHS)
        else if(op.eq.5) then
            call UNIFAES_RHS(Taux,RHS)
        end if
    end subroutine RHS_schemes
    !================================================================
    subroutine CDS(ax,bx,cx,dpx,ay,by,cy,dpy)
        !CDS Coeficients
        double precision :: Lx,Ly,L1,L2,ax,bx,cx,dpx,ay,by,cy,dpy
        OPshem = 'CDS'
        Lx  = dt/(2*dx*dx)
        Ly  = dt/(2*dy*dy)
        L1  = (Pe*dt*cos(alpha))/(4*dx)
        L2  = (Pe*dt*sin(alpha))/(4*dy)
        ax  =-(L1 + Lx)
        bx  = (1 + 2*Lx)
        cx  = (L1 - Lx)
        ay  =-(L2 + Ly)
        by  = (1 + 2*Ly)
        cy  = (L2 - Ly)
        dpx = (1 - 2*Ly)
        dpy = (1 - 2*Lx)
    end subroutine CDS
    !================================================================
    subroutine FOU(ax,bx,cx,dpx,ay,by,cy,dpy)
        !FOU Coeficients
        double precision :: Lx,Ly,L1,L2,ax,bx,cx,dpx,ay,by,cy,dpy
        OPshem = 'FOU'
        Lx  = dt/(2*dx*dx)
        Ly  = dt/(2*dy*dy)
        L1  = (Pe*dt*cos(alpha))/(2*dx)
        L2  = (Pe*dt*sin(alpha))/(2*dy)
        ax  =-(L1 + Lx)
        bx  = (1 + L1 + 2*Lx)
        cx  =- Lx
        ay  =-(L2 + Ly)
        by  = (1 + L2 + 2*Ly)
        cy  =- Ly
        dpx = (1 - L2 - 2*Ly)
        dpy = (1 - L1 - 2*Lx)
    end subroutine FOU
    !================================================================
    subroutine SOU(ax1,bx1,dpx1,ay1,by1,dpy1,ax2,bx2,dpx2,ay2,by2,dpy2,cx,cy,xww,yss)
        !SOU Coeficients
        double precision :: Lx,Ly,L1,L2,cx,cy,xww,yss, &
                            ax1,bx1,dpx1,ay1,by1,dpy1, &
                            ax2,bx2,dpx2,ay2,by2,dpy2
        OPshem = 'SOU'
        Lx   = dt/(2*dx*dx)
        Ly   = dt/(2*dy*dy)
        L1   = (Pe*dt*cos(alpha))/(4*dx)
        L2   = (Pe*dt*sin(alpha))/(4*dy)
        cx   = - Lx
        cy   = - Ly
        xww  = - L1
        yss  = - L2
        ax1  =-(4*L1 + Lx)
        bx1  = (1 + 3*L1 + 2*Lx)
        dpx1 = (1 - 3*L2 - 2*Ly)
        ay1  =-(4*L2 + Ly)
        by1  = (1 + 3*L2 + 2*Ly)
        dpy1 = (1 - 3*L1 - 2*Lx)
        ax2  =-(2*L1 + Lx)
        bx2  = (1 + 2*L1 + 2*Lx)
        dpx2 = (1 - 2*L2 - 2*Ly)
        ay2  =-(2*L2 + Ly)
        by2  = (1 + 2*L2 + 2*Ly)
        dpy2 = (1 - 2*L1 - 2*Lx)
    end subroutine SOU
    !================================================================
    subroutine QUICK(ax1,bx1,dpx1,ay1,by1,dpy1,ax2,bx2,dpx2,ay2,by2,dpy2,cx,cy,xww,yss)
        !SOU Coeficients
        double precision :: Lx,Ly,L1,L2,cx,cy,xww,yss, &
                            ax1,bx1,dpx1,ay1,by1,dpy1, &
                            ax2,bx2,dpx2,ay2,by2,dpy2
        OPshem = 'QUICK'
        Lx   = dt/(2*dx*dx)
        Ly   = dt/(2*dy*dy)
        L1   = (Pe*dt*cos(alpha))/(16*dx)
        L2   = (Pe*dt*sin(alpha))/(16*dy)
        cx   = (3*L1 - Lx)
        cy   = (3*L2 - Ly)  
        xww  = - L1
        yss  = - L2
        ax1  =-(7*L1 + Lx)
        bx1  = (1 + 3*L1 + 2*Lx)
        dpx1 = (1 - 3*L2 - 2*Ly)
        ay1  =-(7*L2 + Ly)
        by1  = (1 + 3*L2 + 2*Ly)
        dpy1 = (1 - 3*L1 - 2*Lx)
        ax2  =-(5*L1 + Lx) 
        bx2  = (1 + 2*L1 + 2*Lx)
        dpx2 = (1 - 2*L2 - 2*Ly)
        ay2  =-(5*L2 + Ly)
        by2  = (1 + 2*L2 + 2*Ly)
        dpy2 = (1 - 2*L1 - 2*Lx)
    end subroutine QUICK
    !================================================================
    subroutine EXPO(ax,bx,cx,dpx,ay,by,cy,dpy)
        !EXP Coeficients
        double precision :: L1,L2,L3,L4,ax,bx,cx,dpx,ay,by,cy,dpy
        OPshem = 'EXP'
        L1  = (Pe*dt*cos(alpha))/(2*dx)
        L2  = (Pe*dt*sin(alpha))/(2*dy)
        L3  = exp(Pe*dx*cos(alpha))
        L4  = exp(Pe*dy*sin(alpha))
        ax  =-(L3/(L3 - 1))*L1
        bx  = (1 + ((L3 + 1)/(L3 - 1))*L1)
        cx  =-(1/(L3 - 1))*L1
        ay  =-(L4/(L4 - 1))*L2
        by  = (1 + ((L4 + 1)/(L4 - 1))*L2)
        cy  =-(1/(L4 - 1))*L2
        dpx = (1 - ((L4 + 1)/(L4 - 1))*L2)
        dpy = (1 - ((L3 + 1)/(L3 - 1))*L1)
    end subroutine EXPO
    !================================================================
    subroutine LOADS_RHS(Ta,RHS_LOADS)
        !LOADS source term
        double precision, dimension(:,:):: Ta,RHS_LOADS
        double precision :: Px,Py,Pix,Piy,Xx,Xy,Zx,Zy,Aux,Auy
        integer ::i,j
        OPshem = 'LOADS'
        Px = Pe*dx*cos(alpha)
        Py = Pe*dy*sin(alpha)
        Pix = Px/(exp(Px) - 1.0)
        Piy = Py/(exp(Py) - 1.0)
        Xx = ((Pix - 1)/Px) + 0.5
        Xy = ((Piy - 1)/Py) + 0.5
        Aux = (Xx/(4*dy*dy))*dt
        Auy = (Xy/(4*dx*dx))*dt
        do i = 2,N-1
            do j = 2,N-1
                Zx = Aux*Py*(Ta(i,j+1) - Ta(i,j-1) + Ta(i-1,j-1) - Ta(i-1,j+1))
                Zx = Zx + Aux*Piy*(Ta(i+1,j-1) + Ta(i-1,j-1) - 2*Ta(i,j-1) - Ta(i+1,j+1) - Ta(i-1,j+1) + 2*Ta(i,j+1))
                Zy = Auy*Px*(Ta(i+1,j) - Ta(i-1,j) + Ta(i-1,j-1) - Ta(i+1,j-1))
                Zy = Zy + Auy*Pix*(Ta(i-1,j+1) + Ta(i-1,j-1) - 2*Ta(i-1,j) - Ta(i+1,j+1) - Ta(i+1,j-1) + 2*Ta(i+1,j))
                RHS_LOADS(i,j) = Zx + Zy
            end do
        end do
    end subroutine LOADS_RHS
    !================================================================
    subroutine UNIFAES_RHS(Ta,RHS_UNI)
        !UNIFAES source term
        double precision, dimension(:,:):: Ta,RHS_UNI
        double precision :: Px,Py,Pix,Piy,Xx,Xy,Zx,Zy,Kn,Ks,Kw,Ke, &
                            PIax,PIay,PIsx,PIsy
        integer ::i,j
        OPshem = 'UNIFAES'
        Px = Pe*dx*cos(alpha)
        Py = Pe*dy*sin(alpha)
        Pix = Px/(exp(Px) - 1.0)
        Piy = Py/(exp(Py) - 1.0)
        Xx = ((Pix - 1)/Px) + 0.5
        Xy = ((Piy - 1)/Py) + 0.5
        PIax = Pix/(dx*dx)
        PIsx = PIax + Px/(dx*dx)
        PIay = Piy/(dy*dy)       !PI+
        PIsy = PIay + Py/(dy*dy) !PI-
        do i = 2,N-1
            if (i.eq.2) then
                do j = 2,N-1
                    Kn = ((Ta(i,j) - Ta(i+2,j))*PIay + (Ta(i+1,j) - Ta(i-1,j))*PIsy)
                    Ks = ((3*Ta(i,j)  - 4*Ta(i+1,j) + Ta(i+2,j))*PIay + (4*Ta(i,j)  - 3*Ta(i-1,j) - Ta(i+1,j))*PIsy)
                    if (j.eq.2) then
                        Ke = ((Ta(i,j) - Ta(i,j+2))*PIax + (Ta(i,j+1) - Ta(i,j-1))*PIsx)
                        Kw = ((3*Ta(i,j) - 4*Ta(i,j+1) + Ta(i,j+2))*PIax + (4*Ta(i,j) - 3*Ta(i,j-1) - Ta(i,j+1))*PIsx)
                    else if (j.gt.2 .and. j.lt.N-1) then
                        Kw = ((Ta(i,j-1) - Ta(i,j+1))*PIax + (Ta(i,j) - Ta(i,j-2))*PIsx)
                        Ke = ((Ta(i,j) - Ta(i,j+2))*PIax + (Ta(i,j+1) - Ta(i,j-1))*PIsx)
                    else if (j.eq.N-1) then
                        Kw = ((Ta(i,j-1) - Ta(i,j+1))*PIax + (Ta(i,j) - Ta(i,j-2))*PIsx)
                        Ke = ((4*Ta(i,j) - 3*Ta(i,j+1) - Ta(i,j-1))*PIax + (3*Ta(i,j) - 4*Ta(i,j-1) + Ta(i,j-2))*PIsx)
                    end if
                    RHS_UNI(i,j) = - (Xx*(Ke - Kw) + Xy*(Kn - Ks))*0.25*dt
                end do
            else if (i.gt.2 .and. i.lt.N-1) then
                do j = 2,N-1
                    Ks = ((Ta(i-1,j) - Ta(i+1,j))*PIay + (Ta(i,j) - Ta(i-2,j))*PIsy)
                    Kn = ((Ta(i,j) - Ta(i+2,j))*PIay + (Ta(i+1,j) - Ta(i-1,j))*PIsy)
                    if (j.eq.2) then
                        Ke = ((Ta(i,j) - Ta(i,j+2))*PIax + (Ta(i,j+1) - Ta(i,j-1))*PIsx)
                        Kw = ((3*Ta(i,j) - 4*Ta(i,j+1) + Ta(i,j+2))*PIax + (4*Ta(i,j) - 3*Ta(i,j-1) - Ta(i,j+1))*PIsx)
                    else if (j.gt.2 .and. j.lt.N-1) then
                        Kw = ((Ta(i,j-1) - Ta(i,j+1))*PIax + (Ta(i,j) - Ta(i,j-2))*PIsx)
                        Ke = ((Ta(i,j) - Ta(i,j+2))*PIax + (Ta(i,j+1) - Ta(i,j-1))*PIsx)
                    else if (j.eq.N-1) then
                        Kw = ((Ta(i,j-1) - Ta(i,j+1))*PIax + (Ta(i,j) - Ta(i,j-2))*PIsx)
                        Ke = ((4*Ta(i,j) - 3*Ta(i,j+1) - Ta(i,j-1))*PIax + (3*Ta(i,j) - 4*Ta(i,j-1) + Ta(i,j-2))*PIsx)
                    end if
                    RHS_UNI(i,j) = - (Xx*(Ke - Kw) + Xy*(Kn - Ks))*0.25*dt
                end do
            else if (i.eq.N-1) then
                do j = 2,N-1
                    Ks = ((Ta(i-1,j) - Ta(i+1,j))*PIay + (Ta(i,j) - Ta(i-2,j))*PIsy)
                    Kn = ((4*Ta(i,j) - 3*Ta(i+1,j) - Ta(i-1,j))*PIay + (3*Ta(i,j) - 4*Ta(i-1,j) + Ta(i-2,j))*PIsy)
                    if (j.eq.2) then
                        Ke = ((Ta(i,j) - Ta(i,j+2))*PIax + (Ta(i,j+1) - Ta(i,j-1))*PIsx)
                        Kw = ((3*Ta(i,j) - 4*Ta(i,j+1) + Ta(i,j+2))*PIax + (4*Ta(i,j) - 3*Ta(i,j-1) - Ta(i,j+1))*PIsx)
                    else if (j.gt.2 .and. j.lt.N-1) then
                        Kw = ((Ta(i,j-1) - Ta(i,j+1))*PIax + (Ta(i,j) - Ta(i,j-2))*PIsx)
                        Ke = ((Ta(i,j) - Ta(i,j+2))*PIax + (Ta(i,j+1) - Ta(i,j-1))*PIsx)
                    else if (j.eq.N-1) then
                        Kw = ((Ta(i,j-1) - Ta(i,j+1))*PIax + (Ta(i,j) - Ta(i,j-2))*PIsx)
                        Ke = ((4*Ta(i,j) - 3*Ta(i,j+1) - Ta(i,j-1))*PIax + (3*Ta(i,j) - 4*Ta(i,j-1) + Ta(i,j-2))*PIsx)                        
                    end if
                    RHS_UNI(i,j) = - (Xx*(Ke - Kw) + Xy*(Kn - Ks))*0.25*dt
                end do
            end if
        end do
    end subroutine UNIFAES_RHS
    !================================================================    
end module schemes
