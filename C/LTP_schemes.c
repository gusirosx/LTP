#include <global.h>
//================================================================
void scheme_Selection(int op)
{
    double Lx,Ly,L1,L2,L3,L4;
    if (op == 1)
    {
        //CDS Coeficients
        strcpy(OPshem,"CDS");
        Lx  = dt/(2.0*dx*dx);
        Ly  = dt/(2.0*dy*dy);
        L1  = (Pe*dt*cos(alpha))/(4.0*dx);
        L2  = (Pe*dt*sin(alpha))/(4.0*dy);
        Ax  =-(L1 + Lx);
        Bx  = (1.0 + 2.0*Lx);
        Cx  = (L1 - Lx);
        Ay  =-(L2 + Ly);
        By  = (1.0 + 2.0*Ly);
        Cy  = (L2 - Ly);
        Dpx = (1.0 - 2.0*Ly);
        Dpy = (1.0 - 2.0*Lx);
    }
    else if (op == 2) // Esquema FOU
    {
        //FOU Coeficients
        strcpy(OPshem,"FOU");
        Lx  = dt/(2.0*dx*dx);
        Ly  = dt/(2.0*dy*dy);
        L1  = (Pe*dt*cos(alpha))/(2.0*dx);
        L2  = (Pe*dt*sin(alpha))/(2.0*dy);
        Ax  =-(L1 + Lx);
        Bx  = (1.0 + L1 + 2.0*Lx);
        Cx  =- Lx;
        Ay  =-(L2 + Ly);
        By  = (1.0 + L2 + 2.0*Ly);
        Cy  =- Ly;
        Dpx = (1.0 - L2 - 2.0*Ly);
        Dpy = (1.0 - L1 - 2.0*Lx);
    }    
    else if (op >= 3 && op <= 5)//Exponencial
    {
        //EXP Coeficients
        strcpy(OPshem,"EXP");
        L1  = (Pe*dt*cos(alpha))/(2.0*dx);
        L2  = (Pe*dt*sin(alpha))/(2.0*dy);
        L3  = exp(Pe*dx*cos(alpha));
        L4  = exp(Pe*dy*sin(alpha));
        Ax  =-(L3/(L3 - 1.0))*L1;
        Bx  = (1.0 + ((L3 + 1.0)/(L3 - 1.0))*L1);
        Cx  =-(1.0/(L3 - 1.0))*L1;
        Ay  =-(L4/(L4 - 1.0))*L2;
        By  = (1.0 + ((L4 + 1.0)/(L4 - 1.0))*L2);
        Cy  =-(1.0/(L4 - 1.0))*L2;
        Dpx = (1.0 - ((L4 + 1.0)/(L4 - 1.0))*L2);
        Dpy = (1.0 - ((L3 + 1.0)/(L3 - 1.0))*L1);
    }
    else if(op == 6)
    {
        //SOU Coeficients
        strcpy(OPshem,"SOU");
        Lx   = dt/(2.0*dx*dx);
        Ly   = dt/(2.0*dy*dy);
        L1   = (Pe*dt*cos(alpha))/(4.0*dx);
        L2   = (Pe*dt*sin(alpha))/(4.0*dy);
        Ax   =-(4.0*L1 + Lx);
        Ax2  =-(2.0*L1 + Lx);
        Bx   = (1.0 + 3.0*L1 + 2.0*Lx);
        Bx2  = (1.0 + 2.0*L1 + 2.0*Lx);
        Cx   = - Lx;
        Ay   =-(4.0*L2 + Ly);
        Ay2  =-(2.0*L2 + Ly);
        By   = (1.0 + 3.0*L2 + 2.0*Ly);
        By2  = (1.0 + 2.0*L2 + 2.0*Ly);
        Cy   = - Ly;   
        Xww  = - L1;
        Yss  = - L2;
        Dpx  = (1.0 - 3.0*L2 - 2.0*Ly);
        Dpx2 = (1.0 - 2.0*L2 - 2.0*Ly);
        Dpy  = (1.0 - 3.0*L1 - 2.0*Lx);
        Dpy2 = (1.0 - 2.0*L1 - 2.0*Lx);
    }
    else if(op == 7)
    {
        //QUICK Coeficients
        strcpy(OPshem,"QUICK");
        Lx   = dt/(2*dx*dx);
        Ly   = dt/(2*dy*dy);
        L1   = (Pe*dt*cos(alpha))/(16*dx);
        L2   = (Pe*dt*sin(alpha))/(16*dy);
        Ax   =-(7*L1 + Lx);
        Ax2  =-(5*L1 + Lx);
        Bx   = (1 + 3*L1 + 2*Lx);
        Bx2  = (1 + 2*L1 + 2*Lx);
        Cx   = (3*L1 - Lx);
        Ay   =-(7*L2 + Ly);
        Ay2  =-(5*L2 + Ly);
        By   = (1 + 3*L2 + 2*Ly);
        By2  = (1 + 2*L2 + 2*Ly);
        Cy   = (3*L2 - Ly);   
        Xww  = - L1;
        Yss  = - L2;
        Dpx  = (1 - 3*L2 - 2*Ly);
        Dpx2 = (1 - 2*L2 - 2*Ly);
        Dpy  = (1 - 3*L1 - 2*Lx);
        Dpy2 = (1 - 2*L1 - 2*Lx);
    }
}
//================================================================
void RHS_schemes(int op)
{
    //Calculates RHS of LOADS or UNIFAES
    for (int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            S_rhs[i][j] = 0.0;
        }
    }
    if(op == 4) LOADS_RHS();
    else if(op == 5) UNIFAES_RHS();

}
//================================================================
void LOADS_RHS()
{
    //LOADS source term
    double Px,Py,Pix,Piy,Xx,Xy,Zx,Zy,Aux,Auy;
    strcpy(OPshem,"LOADS");
    Px = Pe*dx*cos(alpha);
    Py = Pe*dy*sin(alpha);
    Pix = Px/(exp(Px) - 1.0);
    Piy = Py/(exp(Py) - 1.0);
    Xx = ((Pix - 1.0)/Px) + 0.5;
    Xy = ((Piy - 1.0)/Py) + 0.5;
    Aux = (Xx/(4.0*dy*dy))*dt;
    Auy = (Xy/(4.0*dx*dx))*dt; 
    for(int i=1;i<N-1;i++)
    {
        for(int j=1;j<N-1;j++)
        {
            Zx = Aux*Py*(S_a[i][j+1] - S_a[i][j-1] + S_a[i-1][j-1] - S_a[i-1][j+1]);
            Zx = Zx + Aux*Piy*(S_a[i+1][j-1] + S_a[i-1][j-1] - 2*S_a[i][j-1] - S_a[i+1][j+1] - S_a[i-1][j+1] + 2*S_a[i][j+1]);
            Zy = Auy*Px*(S_a[i+1][j] - S_a[i-1][j] + S_a[i-1][j-1] - S_a[i+1][j-1]);
            Zy = Zy + Auy*Pix*(S_a[i-1][j+1] + S_a[i-1][j-1] - 2*S_a[i-1][j] - S_a[i+1][j+1] - S_a[i+1][j-1] + 2*S_a[i+1][j]);
            S_rhs[i][j] = Zx + Zy;
        }
    }
}
//================================================================
void UNIFAES_RHS()
{
    int i,j;
    double Px,Py,Pix,Piy,Xx,Xy,Zx,Zy,Kn,Ks,Kw,Ke,PIax,PIay,PIsx,PIsy;
    strcpy(OPshem,"UNIFAES");
    Px = Pe*dx*cos(alpha);
    Py = Pe*dy*sin(alpha);
    Pix = Px/(exp(Px) - 1.0);
    Piy = Py/(exp(Py) - 1.0);
    Xx = ((Pix - 1.0)/Px) + 0.5;
    Xy = ((Piy - 1.0)/Py) + 0.5;
    PIax = Pix/(dx*dx);
    PIsx = PIax + Px/(dx*dx);
    PIay = Piy/(dy*dy);       //PI+
    PIsy = PIay + Py/(dy*dy); //PI- 
    for(i=1;i<N-1;i++)
    {   
        if(i == 1)
        {
            for(j=1;j<(N-1);j++)
            {
                Kn = (S_a[i][j] - S_a[i+2][j])*PIay + (S_a[i+1][j] - S_a[i-1][j])*PIsy;
                Ks = (3.0*S_a[i][j] - 4.0*S_a[i+1][j] + S_a[i+2][j])*PIay + (4.0*S_a[i][j] - 3.0*S_a[i-1][j] - S_a[i+1][j])*PIsy;
                if(j == 1)
                {
                    Ke = (S_a[i][j] - S_a[i][j+2])*PIax + (S_a[i][j+1] - S_a[i][j-1])*PIsx;
                    Kw = (3.0*S_a[i][j] - 4.0*S_a[i][j+1] + S_a[i][j+2])*PIax + (4.0*S_a[i][j] - 3.0*S_a[i][j-1] - S_a[i][j+1])*PIsx;
                }
                else if(j > 1 && j < (N - 2))
                {
                    Kw = (S_a[i][j-1] - S_a[i][j+1])*PIax + (S_a[i][j] - S_a[i][j-2])*PIsx;
                    Ke = (S_a[i][j] - S_a[i][j+2])*PIax + (S_a[i][j+1] - S_a[i][j-1])*PIsx;    
                }
                else if(j == (N - 2))
                {
                    Kw = (S_a[i][j-1] - S_a[i][j+1])*PIax + (S_a[i][j] - S_a[i][j-2])*PIsx;
                    Ke = (4.0*S_a[i][j] - 3.0*S_a[i][j+1] - S_a[i][j-1])*PIax + (3.0*S_a[i][j] - 4.0*S_a[i][j-1] + S_a[i][j-2])*PIsx;
                }
                S_rhs[i][j] = - (Xx*(Ke - Kw) + Xy*(Kn - Ks))*0.25*dt;
            }
        }
        else if(i > 1 && i<(N-2))
        {
            for(j=1;j<(N-1);j++)
            {
                Ks = (S_a[i-1][j] - S_a[i+1][j])*PIay + (S_a[i][j] - S_a[i-2][j])*PIsy;
                Kn = (S_a[i][j] - S_a[i+2][j])*PIay + (S_a[i+1][j] - S_a[i-1][j])*PIsy;
                if(j == 1)
                {
                    Ke = (S_a[i][j] - S_a[i][j+2])*PIax + (S_a[i][j+1] - S_a[i][j-1])*PIsx;
                    Kw = (3.0*S_a[i][j] - 4.0*S_a[i][j+1] + S_a[i][j+2])*PIax + (4.0*S_a[i][j] - 3.0*S_a[i][j-1] - S_a[i][j+1])*PIsx;
                }
                else if(j > 1 && j < (N - 2))
                {
                    Kw = (S_a[i][j-1] - S_a[i][j+1])*PIax + (S_a[i][j] - S_a[i][j-2])*PIsx;
                    Ke = (S_a[i][j] - S_a[i][j+2])*PIax + (S_a[i][j+1] - S_a[i][j-1])*PIsx;    
                }
                else if(j == (N - 2))
                {
                    Kw = (S_a[i][j-1] - S_a[i][j+1])*PIax + (S_a[i][j] - S_a[i][j-2])*PIsx;
                    Ke = (4.0*S_a[i][j] - 3.0*S_a[i][j+1] - S_a[i][j-1])*PIax + (3.0*S_a[i][j] - 4.0*S_a[i][j-1] + S_a[i][j-2])*PIsx;    
                }
                S_rhs[i][j] = - (Xx*(Ke - Kw) + Xy*(Kn - Ks))*0.25*dt;
            }
        }
        else if(i == (N-2))
        {
            for(j=1;j<(N-1);j++)
            {
                Ks = (S_a[i-1][j] - S_a[i+1][j])*PIay + (S_a[i][j] - S_a[i-2][j])*PIsy;
                Kn = (4.0*S_a[i][j] - 3.0*S_a[i+1][j] - S_a[i-1][j])*PIay + (3.0*S_a[i][j] - 4.0*S_a[i-1][j] + S_a[i-2][j])*PIsy;
                if(j == 1)
                {
                    Ke = (S_a[i][j] - S_a[i][j+2])*PIax + (S_a[i][j+1] - S_a[i][j-1])*PIsx;
                    Kw = (3.0*S_a[i][j] - 4.0*S_a[i][j+1] + S_a[i][j+2])*PIax + (4.0*S_a[i][j] - 3.0*S_a[i][j-1] - S_a[i][j+1])*PIsx;
                }
                else if(j > 1 && j < (N - 2))
                {
                    Kw = ((S_a[i][j-1] - S_a[i][j+1])*PIax + (S_a[i][j] - S_a[i][j-2])*PIsx);
                    Ke = ((S_a[i][j] - S_a[i][j+2])*PIax + (S_a[i][j+1] - S_a[i][j-1])*PIsx);
                }
                else if(j == (N - 2))
                {
                    Kw = (S_a[i][j-1] - S_a[i][j+1])*PIax + (S_a[i][j] - S_a[i][j-2])*PIsx;
                    Ke = (4.0*S_a[i][j] - 3.0*S_a[i][j+1] - S_a[i][j-1])*PIax + (3.0*S_a[i][j] - 4.0*S_a[i][j-1] + S_a[i][j-2])*PIsx;
                }
                S_rhs[i][j] = - (Xx*(Ke - Kw) + Xy*(Kn - Ks))*0.25*dt;
            }
        } 
    }
}
//================================================================
