import numpy as np
#=================== Global ====================#
#--------------- INPUT PARAMETERS --------------#
OPfunc     = 'A'                #Function Type: A,B,C,D,CD,DC
sheme_flag = True               #Scheme Standard routine
Pe         = 60.0               #Global Peclet Number
lamb       = 6.0                #Eigenvalue 
theta      = 22.5               #Mesh offset angle
TOL        = 1e-5               #Convergence tolerance
N          = 11                 #Number of Control Volumes in x and y
#----------- ERROR EVALUATION ROUTINE ----------#
rms_flag   = False              #Enables or disables the rms routine
N_i        = 10                 #Initial mesh
N_f        = 100                #Final mesh
IJ         = 10                 #Jump Interval
pi         = np.pi              #Pi number
x = np.zeros(N,dtype = float)   #Determine the x positions
y = np.zeros(N,dtype = float)   #Determine the y positions
S_exact = np.zeros((N,N),dtype = float)#Exact solution
S_calc  = np.zeros((N,N),dtype = float)#Calculated solution
#================== Functions ==================#
#----------------------------------------------------------------
def norm_function(S):
    #Normalization by the max value of a 2D array
    global N
    F = np.max(S) - np.min(S)
    for i in range(0,N):
        for j in range(0,N):
            S[i,j] = S[i,j]/F
    return
#----------------------------------------------------------------
def deltas():
    global dx,dy,dt,alpha,x,y
    alpha = theta*pi/180
    dx = 1.0/(N-1.0) 
    dy = dx
    #dt = 0.98d0*dx*dx
    dt = 0.98*(4*dx)/(3*Pe*np.cos(alpha))
    #Determine the x and y positions for solutions A, B, C, D
    for i in range(0,N):
        x[i] = -0.5 + i*dx
        y[i] = -0.5 + i*dy
#----------------------------------------------------------------
def TDMA(A,B,C,D):
    #Thomas Algorithm
    n = len(B)
    W = np.zeros(n,dtype = float)
    G = np.zeros(n,dtype = float)
    T = np.zeros(n,dtype = float)
    W[0] = C[0]/B[0]
    G[0] = D[0]/B[0]
    for i in range(1,n-1):
        W[i] = C[i]/(B[i] - A[i-1]*W[i-1])
    for i in range(1,n):
        G[i] = (D[i] - A[i-1]*G[i-1])/(B[i] - A[i-1]*W[i-1])
    T[n-1] = G[n-1]
    for i in range(n-1,0,-1):
        T[i-1] = G[i-1] - W[i-1]*T[i]
    return T
#----------------------------------------------------------------
def exact_functions(S_sol,S_c):
    #Calculation of the exact solution - TYPE A, B, C, D
    global Pe,lamb,N,x,y,OPfunc
    if(OPfunc == 'A'):
        C = (Pe - (Pe**2.0 + 4.0*lamb**2.0)**0.5)/2.0
    elif(OPfunc == 'B'):
        C = (Pe + (Pe**2.0 + 4.0*lamb**2.0)**0.5)/2.0
    elif (OPfunc == 'C'):
        C = (Pe - (Pe**2.0 - 4.0*lamb**2.0)**0.5)/2.0
    elif (OPfunc == 'D'):
        C = (Pe + (Pe**2.0 - 4.0*lamb**2.0)**0.5)/2.0
    elif (OPfunc == 'CD'or OPfunc == 'DC'):
        C = (lamb**2.0 - (Pe**2.0)*0.25)**0.5
    #Loop for functions A and B
    if(OPfunc == 'A'or OPfunc == 'B'):
        for i in range(0,N):
            for j in range(0,N):
                S = x[j]*np.cos(alpha) + y[i]*np.sin(alpha)
                M = y[i]*np.cos(alpha) - x[j]*np.sin(alpha)
                S_sol[i,j] = np.exp(C*S)*np.sin(lamb*M)              
    #Loop for functions C and D
    if (OPfunc == 'C' or OPfunc == 'D'):
        for i in range(0,N):
            for j in range(0,N):
                S = x[j]*np.cos(alpha) + y[i]*np.sin(alpha)
                M = y[i]*np.cos(alpha) - x[j]*np.sin(alpha)
                S_sol[i,j] = np.exp(C*S)*np.exp(lamb*M)  
    #Loop for functions CD
    if (OPfunc == 'CD'):
        for i in range(0,N):
            for j in range(0,N):
                S = x[j]*np.cos(alpha) + y[i]*np.sin(alpha)
                M = y[i]*np.cos(alpha) - x[j]*np.sin(alpha)
                S_sol[i,j] = np.exp(Pe*S*0.5)*np.sin(C*S)*np.exp(lamb*M)
    #Loop for functions DC
    if (OPfunc == 'DC'):
        for i in range(0,N):
            for j in range(0,N):
                S = x[j]*np.cos(alpha) + y[i]*np.sin(alpha)
                M = y[i]*np.cos(alpha) - x[j]*np.sin(alpha)
                S_sol[i,j] = np.exp(Pe*S*0.5)*np.cos(C*S)*np.exp(lamb*M)
    norm_function(S_sol)
    #Extract the boundary conditions
    for i in range(0,N):
        S_c[0,i]    = S_sol[0,i]
        S_c[N-1,i]  = S_sol[N-1,i]
        S_c[i,0]   = S_sol[i,0]
        S_c[i,N-1] = S_sol[i,N-1]
    return
#----------------------------------------------------------------
def ADI(S_sch,op):
    if(op >= 1 and op <= 5): 
        ADI_S(S_sch,op)
    elif(op >= 6 and op <= 7): 
        ADI_MOD(S_sch,op)
    return
#----------------------------------------------------------------
def ADI_S(S_sch,op):#Standard ADI
    global N,TOL,OPshem
    B, D, TS = np.zeros((3, N-2),dtype = float)
    A, C     = np.zeros((2, N-3),dtype = float)
    Ax,Bx,Cx,Dpx,Ay,By,Cy,Dpy = scheme_Selection(op)
    S_rhs = np.zeros((N,N),dtype = float)
    for k in range(1,1000000):
        #X Direction
        S_a = S_sch.copy() # Copy the elements from S_sch to S_a
        if(op == 4 or op == 5): S_rhs = RHS_schemes(op,S_a)
        #Direção X
        for i in range(1,N-1):
            for j in range (1,N-1):
                if j == 1 :
                    B[j-1] = Bx
                    C[j-1] = Cx
                    D[j-1] = - Ay*S_a[i-1,j] + Dpx*S_a[i,j] - Cy*S_a[i+1,j] - Ax*S_a[i,j-1] + S_rhs[i,j]
                if j > 1 and j < (N - 2):
                    A[j-2] = Ax
                    B[j-1] = Bx
                    C[j-1] = Cx
                    D[j-1] = - Ay*S_a[i-1,j] + Dpx*S_a[i,j] - Cy*S_a[i+1,j] + S_rhs[i,j]
                if j == (N - 2):
                    A[j-2] = Ax
                    B[j-1] = Bx
                    D[j-1] = - Ay*S_a[i-1,j] + Dpx*S_a[i,j] - Cy*S_a[i+1,j] - Cx*S_a[i,j+1] + S_rhs[i,j]
            TS = TDMA(A,B,C,D)
            for j in range (1,N-1):
                S_sch[i,j] = TS[j-1] #Stores TDMA values ​​in matrix S_sch
        S_a = S_sch.copy()
        if(op == 4 or op == 5): S_rhs = RHS_schemes(op,S_a)
        #Direção Y
        for j in range(1,N-1):
            for i in range (1,N-1):
                if i == 1:
                    B[i-1] = By
                    C[i-1] = Cy
                    D[i-1] = - Ax*S_a[i,j-1] + Dpy*S_a[i,j] - Cx*S_a[i,j+1] - Ay*S_a[i-1,j] + S_rhs[i,j]
                if i > 1 and i < (N - 2):
                    A[i-2] = Ay
                    B[i-1] = By
                    C[i-1] = Cy
                    D[i-1] = - Ax*S_a[i,j-1] + Dpy*S_a[i,j] - Cx*S_a[i,j+1] + S_rhs[i,j]
                if i == (N - 2):
                    A[i-2] = Ay
                    B[i-1] = By
                    D[i-1] = - Ax*S_a[i,j-1] + Dpy*S_a[i,j] - Cx*S_a[i,j+1] - Cy*S_a[i+1,j] + S_rhs[i,j]
            TS = TDMA(A,B,C,D)
            for i in range (1,N-1):
                S_sch[i,j] = TS[i-1]
        VTOL = np.max(abs(S_sch - S_a))/dt
        if (VTOL <= TOL):
            print( 'C.R.I:', k, OPshem)
            break
    return
#----------------------------------------------------------------
def ADI_MOD(S_sch,op):#Modified ADI for WW and SS terms
    global N,TOL,OPshem
    B, D, TS = np.zeros((3, N-2),dtype = float)
    A, C     = np.zeros((2, N-3),dtype = float)
    Ax1,Bx1,Dpx1,Ay1,By1,Dpy1,Ax2,Bx2,Dpx2,Ay2,By2,Dpy2,Cx,Cy,Xww,Yss = scheme_SelectionM(op)
    for k in range(0,1000000):
        #X Direction
        S_a = S_sch.copy() # Copy the elements from S_sch to S_a
        for i in range(1,N-1):
            if i == 1:
                for j in range (1,N-1):
                    if j == 1 :
                        B[j-1] = Bx2
                        C[j-1] = Cx
                        D[j-1] = - Ay2*S_a[i-1,j] + Dpx2*S_a[i,j] - Cy*S_a[i+1,j] - Ax2*S_a[i,j-1]
                    if j > 1 and j < (N - 2):
                        A[j-2] = Ax1
                        B[j-1] = Bx1
                        C[j-1] = Cx
                        D[j-1] = - Ay2*S_a[i-1,j] + Dpx2*S_a[i,j] - Cy*S_a[i+1,j] + Xww*S_a[i,j-2]
                    if j == (N - 2):
                        A[j-2] = Ax1
                        B[j-1] = Bx1
                        D[j-1] = - Ay2*S_a[i-1,j] + Dpx2*S_a[i,j] - Cy*S_a[i+1,j] - Cx*S_a[i,j+1] + Xww*S_a[i,j-2]
            if i > 1:
                for j in range (1,N-1):
                    if j == 1 :
                        B[j-1] = Bx2
                        C[j-1] = Cx
                        D[j-1] = - Ay1*S_a[i-1,j] + Dpx1*S_a[i,j] - Cy*S_a[i+1,j] - Ax2*S_a[i,j-1] + Yss*S_a[i-2,j]
                    if j > 1 and j < (N - 2):
                        A[j-2] = Ax1
                        B[j-1] = Bx1
                        C[j-1] = Cx
                        D[j-1] = - Ay1*S_a[i-1,j] + Dpx1*S_a[i,j] - Cy*S_a[i+1,j] + Xww*S_a[i,j-2] + Yss*S_a[i-2,j]  
                    if j == (N - 2):
                        A[j-2] = Ax1
                        B[j-1] = Bx1
                        D[j-1] = - Ay1*S_a[i-1,j] + Dpx1*S_a[i,j] - Cy*S_a[i+1,j] - Cx*S_a[i,j+1] + Xww*S_a[i,j-2] + Yss*S_a[i-2,j]  
            TS = TDMA(A,B,C,D)
            for j in range (1,N-1):
                S_sch[i,j] = TS[j-1] #Stores TDMA values ​​in matrix S_sch
        S_a = S_sch.copy()
        #Y Direction
        for j in range(1,N-1):
            if j == 1:
                for i in range (1,N-1):
                    if i == 1:
                        B[i-1] = By2
                        C[i-1] = Cy
                        D[i-1] = - Ax2*S_a[i,j-1] + Dpy2*S_a[i,j] - Cx*S_a[i,j+1] - Ay2*S_a[i-1,j]
                    if i > 1 and i < (N - 2):
                        A[i-2] = Ay1
                        B[i-1] = By1
                        C[i-1] = Cy
                        D[i-1] = - Ax2*S_a[i,j-1] + Dpy2*S_a[i,j] - Cx*S_a[i,j+1]+ Yss*S_a[i-2,j]
                    if i == (N - 2):
                        A[i-2] = Ay1
                        B[i-1] = By1
                        D[i-1] = - Ax2*S_a[i,j-1] + Dpy2*S_a[i,j] - Cx*S_a[i,j+1] - Cy*S_a[i+1,j] + Yss*S_a[i-2,j]
            if j > 1:
                for i in range (1,N-1):
                    if i == 1:
                        B[i-1] = By2
                        C[i-1] = Cy
                        D[i-1] = - Ax1*S_a[i,j-1] + Dpy1*S_a[i,j] - Cx*S_a[i,j+1] - Ay2*S_a[i-1,j] + Xww*S_a[i,j-2]
                    if i > 1 and i < (N - 2):
                        A[i-2] = Ay1
                        B[i-1] = By1
                        C[i-1] = Cy
                        D[i-1] = - Ax1*S_a[i,j-1] + Dpy1*S_a[i,j] - Cx*S_a[i,j+1] + Xww*S_a[i,j-2] + Yss*S_a[i-2,j]
                    if i == (N - 2):
                        A[i-2] = Ay1
                        B[i-1] = By1
                        D[i-1] = - Ax1*S_a[i,j-1] + Dpy1*S_a[i,j] - Cx*S_a[i,j+1] - Cy*S_a[i+1,j] + Xww*S_a[i,j-2] + Yss*S_a[i-2,j]
            TS = TDMA(A,B,C,D)
            for i in range (1,N-1):
                S_sch[i,j] = TS[i-1]
        VTOL = np.max(abs(S_sch - S_a))/dt
        if (VTOL <= TOL):
            print( 'C.R.I:', k, OPshem)
            break
    return S_sch
#----------------------------------------------------------------
#=================== Schemes ===================#
def scheme_Selection(op):
    if(op == 1):   ax,bx,cx,dpx,ay,by,cy,dpy = CDS()
    elif(op == 2): ax,bx,cx,dpx,ay,by,cy,dpy = FOU()
    elif(op >= 3 and op <= 5): ax,bx,cx,dpx,ay,by,cy,dpy = EXP()
    #elif(op == 6): SOU(ax,bx,dpx,ay,by,dpy,axm,bxm,dpxm,aym,bym,dpym,cx,cy,xww,yss)
    #elif(op == 6):QUICK(ax,bx,dpx,ay,by,dpy,axm,bxm,dpxm,aym,bym,dpym,cx,cy,xww,yss)    
    return ax,bx,cx,dpx,ay,by,cy,dpy  
#----------------------------------------------------------------
def scheme_SelectionM(op):
    if(op == 6):   ax1,bx1,dpx1,ay1,by1,dpy1,ax2,bx2,dpx2,ay2,by2,dpy2,cx,cy,xww,yss = SOU()
    elif(op == 7): ax1,bx1,dpx1,ay1,by1,dpy1,ax2,bx2,dpx2,ay2,by2,dpy2,cx,cy,xww,yss = QUICK()    
    return ax1,bx1,dpx1,ay1,by1,dpy1,ax2,bx2,dpx2,ay2,by2,dpy2,cx,cy,xww,yss
#----------------------------------------------------------------
def CDS():
    #CDS Coeficients
    global dx,dy,dt,Pe,alpha,OPshem
    OPshem = 'CDS'
    Lx = dt/(2*dx*dx)
    Ly = dt/(2*dy*dy)
    L1 = (Pe*dt*np.cos(alpha))/(4*dx)
    L2 = (Pe*dt*np.sin(alpha))/(4*dy)
    ax =-(L1 + Lx)
    bx = (1 + 2*Lx)
    cx = (L1 - Lx)
    ay =-(L2 + Ly)
    by = (1 + 2*Ly)
    cy = (L2 - Ly)
    dpx = (1 - 2*Ly)
    dpy = (1 - 2*Lx)
    return ax,bx,cx,dpx,ay,by,cy,dpy
#----------------------------------------------------------------
def FOU():
    #FOU Coeficients
    global dx,dy,dt,Pe,alpha,OPshem
    OPshem = 'FOU'
    Lx  = dt/(2*dx*dx)
    Ly  = dt/(2*dy*dy)
    L1  = (Pe*dt*np.cos(alpha))/(2*dx)
    L2  = (Pe*dt*np.sin(alpha))/(2*dy)
    ax  =-(L1 + Lx)
    bx  = (1 + L1 + 2*Lx)
    cx  =- Lx
    ay  =-(L2 + Ly)
    by  = (1 + L2 + 2*Ly)
    cy  =- Ly
    dpx = (1 - L2 - 2*Ly)
    dpy = (1 - L1 - 2*Lx)
    return ax,bx,cx,dpx,ay,by,cy,dpy
#----------------------------------------------------------------
def EXP():
    #EXP Coeficients
    global dx,dy,dt,Pe,alpha,OPshem
    OPshem = 'EXP'
    L1  = (Pe*dt*np.cos(alpha))/(2*dx)
    L2  = (Pe*dt*np.sin(alpha))/(2*dy)
    L3  = np.exp(Pe*dx*np.cos(alpha))
    L4  = np.exp(Pe*dy*np.sin(alpha))
    ax  =-(L3/(L3 - 1))*L1
    bx  = (1 + ((L3 + 1)/(L3 - 1))*L1)
    cx  =-(1/(L3 - 1))*L1
    ay  =-(L4/(L4 - 1))*L2
    by  = (1 + ((L4 + 1)/(L4 - 1))*L2)
    cy  =-(1/(L4 - 1))*L2
    dpx = (1 - ((L4 + 1)/(L4 - 1))*L2)
    dpy = (1 - ((L3 + 1)/(L3 - 1))*L1)
    return ax,bx,cx,dpx,ay,by,cy,dpy
#----------------------------------------------------------------
def QUICK(): 
    #QUICK Coeficients
    global dx,dy,dt,Pe,alpha,OPshem
    OPshem = 'QUICK'
    Lx   = dt/(2*dx*dx)
    Ly   = dt/(2*dy*dy)
    L1   = (Pe*dt*np.cos(alpha))/(16*dx)
    L2   = (Pe*dt*np.sin(alpha))/(16*dy)
    ax1  =-(7*L1 + Lx)
    ax2  =-(5*L1 + Lx)
    bx1  = (1 + 3*L1 + 2*Lx)
    bx2  = (1 + 2*L1 + 2*Lx)
    cx   = (3*L1 - Lx)
    ay1  =-(7*L2 + Ly)
    ay2  =-(5*L2 + Ly)
    by1  = (1 + 3*L2 + 2*Ly)
    by2  = (1 + 2*L2 + 2*Ly)
    cy   = (3*L2 - Ly)   
    xww  = - L1
    yss  = - L2
    dpx1 = (1 - 3*L2 - 2*Ly)
    dpx2 = (1 - 2*L2 - 2*Ly)
    dpy1 = (1 - 3*L1 - 2*Lx)
    dpy2 = (1 - 2*L1 - 2*Lx)
    return ax1,bx1,dpx1,ay1,by1,dpy1,ax2,bx2,dpx2,ay2,by2,dpy2,cx,cy,xww,yss
#----------------------------------------------------------------
def SOU():
    #SOU Coeficients
    global dx,dy,dt,Pe,alpha,OPshem
    OPshem = 'SOU'
    Lx   = dt/(2*dx*dx)
    Ly   = dt/(2*dy*dy)
    L1   = (Pe*dt*np.cos(alpha))/(4*dx)
    L2   = (Pe*dt*np.sin(alpha))/(4*dy)
    ax1  =-(4*L1 + Lx)
    ax2  =-(2*L1 + Lx)
    bx1  = (1 + 3*L1 + 2*Lx)
    bx2  = (1 + 2*L1 + 2*Lx)
    cx   = - Lx
    ay1  =-(4*L2 + Ly)
    ay2  =-(2*L2 + Ly)
    by1  = (1 + 3*L2 + 2*Ly)
    by2  = (1 + 2*L2 + 2*Ly)
    cy   = - Ly   
    xww  = - L1
    yss  = - L2
    dpx1 = (1 - 3*L2 - 2*Ly)
    dpx2 = (1 - 2*L2 - 2*Ly)
    dpy1 = (1 - 3*L1 - 2*Lx)
    dpy2 = (1 - 2*L1 - 2*Lx)
    return ax1,bx1,dpx1,ay1,by1,dpy1,ax2,bx2,dpx2,ay2,by2,dpy2,cx,cy,xww,yss
#----------------------------------------------------------------
def RHS_schemes(op,Taux):
    #Calculates of LOADS or UNIFAES
    if(op == 4):   RHS = LOADS_RHS(Taux)
    elif(op == 5): RHS = UNIFAES_RHS(Taux)  
    return RHS
#----------------------------------------------------------------
def LOADS_RHS(Ta):
    #LOADS source term
    global N,dx,dy,dt,Pe,alpha,OPshem
    OPshem = 'LOADS'
    Px = Pe*dx*np.cos(alpha)
    Py = Pe*dy*np.sin(alpha)
    Pix = Px/(np.exp(Px) - 1.0)
    Piy = Py/(np.exp(Py) - 1.0)
    Xx = ((Pix - 1)/Px) + 0.5
    Xy = ((Piy - 1)/Py) + 0.5
    Aux = (Xx/(4*dy*dy))*dt
    Auy = (Xy/(4*dx*dx))*dt   
    RHS_LOADS = np.zeros((N,N),dtype = float)
    for i in range(1,N-1):
        for j in range (1,N-1):
            Zx = Aux*Py*(Ta[i,j+1] - Ta[i,j-1] + Ta[i-1,j-1] - Ta[i-1,j+1])
            Zx = Zx + Aux*Piy*(Ta[i+1,j-1] + Ta[i-1,j-1] - 2*Ta[i,j-1] - Ta[i+1,j+1] - Ta[i-1,j+1] + 2*Ta[i,j+1])
            Zy = Auy*Px*(Ta[i+1,j] - Ta[i-1,j] + Ta[i-1,j-1] - Ta[i+1,j-1])
            Zy = Zy + Auy*Pix*(Ta[i-1,j+1] + Ta[i-1,j-1] - 2*Ta[i-1,j] - Ta[i+1,j+1] - Ta[i+1,j-1] + 2*Ta[i+1,j])
            RHS_LOADS[i,j] = Zx + Zy
    return RHS_LOADS
#----------------------------------------------------------------
def UNIFAES_RHS(Ta):
    #LOADS source term
    global N,dx,dy,dt,Pe,alpha,OPshem
    OPshem = 'UNIFAES'
    Px = Pe*dx*np.cos(alpha)
    Py = Pe*dy*np.sin(alpha)
    Pix = Px/(np.exp(Px) - 1.0)
    Piy = Py/(np.exp(Py) - 1.0)
    Xx = ((Pix - 1)/Px) + 0.5
    Xy = ((Piy - 1)/Py) + 0.5
    PIax = Pix/(dx*dx)
    PIsx = PIax + Px/(dx*dx)
    PIay = Piy/(dy*dy)       #PI+
    PIsy = PIay + Py/(dy*dy) #PI- 
    RHS_UNI = np.zeros((N,N),dtype = float)
    for i in range(1,N-1):
        if (i == 1):
            for j in range (1,N-1):
                Kn = (Ta[i,j] - Ta[i+2,j])*PIay + (Ta[i+1,j] - Ta[i-1,j])*PIsy
                Ks = (3*Ta[i,j] - 4*Ta[i+1,j] + Ta[i+2,j])*PIay + (4*Ta[i,j] - 3*Ta[i-1,j] - Ta[i+1,j])*PIsy
                if (j == 1):
                    Ke = (Ta[i,j] - Ta[i,j+2])*PIax + (Ta[i,j+1] - Ta[i,j-1])*PIsx
                    Kw = (3*Ta[i,j] - 4*Ta[i,j+1] + Ta[i,j+2])*PIax + (4*Ta[i,j] - 3*Ta[i,j-1] - Ta[i,j+1])*PIsx
                if (j > 1 and j < (N - 2)):
                    Kw = (Ta[i,j-1] - Ta[i,j+1])*PIax + (Ta[i,j] - Ta[i,j-2])*PIsx
                    Ke = (Ta[i,j] - Ta[i,j+2])*PIax + (Ta[i,j+1] - Ta[i,j-1])*PIsx
                if (j == (N - 2)):
                    Kw = (Ta[i,j-1] - Ta[i,j+1])*PIax + (Ta[i,j] - Ta[i,j-2])*PIsx
                    Ke = (4*Ta[i,j] - 3*Ta[i,j+1] - Ta[i,j-1])*PIax + (3*Ta[i,j] - 4*Ta[i,j-1] + Ta[i,j-2])*PIsx
                RHS_UNI[i,j] = - (Xx*(Ke - Kw) + Xy*(Kn - Ks))*0.25*dt
        elif (i > 1 and i < (N - 2)):
            for j in range (1,N-1):
                Ks = (Ta[i-1,j] - Ta[i+1,j])*PIay + (Ta[i,j] - Ta[i-2,j])*PIsy
                Kn = (Ta[i,j] - Ta[i+2,j])*PIay + (Ta[i+1,j] - Ta[i-1,j])*PIsy
                if (j == 1) :
                    Ke = (Ta[i,j] - Ta[i,j+2])*PIax + (Ta[i,j+1] - Ta[i,j-1])*PIsx
                    Kw = (3*Ta[i,j] - 4*Ta[i,j+1] + Ta[i,j+2])*PIax + (4*Ta[i,j] - 3*Ta[i,j-1] - Ta[i,j+1])*PIsx
                if (j > 1 and j < (N - 2)):
                    Kw = (Ta[i,j-1] - Ta[i,j+1])*PIax + (Ta[i,j] - Ta[i,j-2])*PIsx
                    Ke = (Ta[i,j] - Ta[i,j+2])*PIax + (Ta[i,j+1] - Ta[i,j-1])*PIsx
                if (j == (N - 2)):
                    Kw = (Ta[i,j-1] - Ta[i,j+1])*PIax + (Ta[i,j] - Ta[i,j-2])*PIsx
                    Ke = (4*Ta[i,j] - 3*Ta[i,j+1] - Ta[i,j-1])*PIax + (3*Ta[i,j] - 4*Ta[i,j-1] + Ta[i,j-2])*PIsx
                RHS_UNI[i,j] = - (Xx*(Ke - Kw) + Xy*(Kn - Ks))*0.25*dt
        elif (i == (N - 2)):
            for j in range (1,N-1):
                Ks = (Ta[i-1,j] - Ta[i+1,j])*PIay + (Ta[i,j] - Ta[i-2,j])*PIsy
                Kn = (4*Ta[i,j] - 3*Ta[i+1,j] - Ta[i-1,j])*PIay + (3*Ta[i,j] - 4*Ta[i-1,j] + Ta[i-2,j])*PIsy
                if (j == 1):
                    Ke = (Ta[i,j] - Ta[i,j+2])*PIax + (Ta[i,j+1] - Ta[i,j-1])*PIsx
                    Kw = (3*Ta[i,j] - 4*Ta[i,j+1] + Ta[i,j+2])*PIax + (4*Ta[i,j] - 3*Ta[i,j-1] - Ta[i,j+1])*PIsx
                if (j > 1 and j < (N - 2)):
                    Kw = ((Ta[i,j-1] - Ta[i,j+1])*PIax + (Ta[i,j] - Ta[i,j-2])*PIsx)
                    Ke = ((Ta[i,j] - Ta[i,j+2])*PIax + (Ta[i,j+1] - Ta[i,j-1])*PIsx)
                if (j == (N - 2)):
                    Kw = (Ta[i,j-1] - Ta[i,j+1])*PIax + (Ta[i,j] - Ta[i,j-2])*PIsx
                    Ke = (4*Ta[i,j] - 3*Ta[i,j+1] - Ta[i,j-1])*PIax + (3*Ta[i,j] - 4*Ta[i,j-1] + Ta[i,j-2])*PIsx
                RHS_UNI[i,j] = - (Xx*(Ke - Kw) + Xy*(Kn - Ks))*0.25*dt
    return RHS_UNI
#----------------------------------------------------------------
def RMS(S_aprox,S_exat): # Calculating the RMS of a matrix
    global N
    e_rms = 0
    for i in range(0,N):
        for j in range(0,N):
            Ei = (S_aprox[i,j] - S_exat[i,j])**2.0
            e_rms = e_rms + Ei
    Ei = (e_rms/((N-2)*(N-2)))
    e_rms = np.sqrt(Ei)
    return e_rms
#----------------------------------------------------------------
#=================== Routines ==================#
def Scheme_routine():
    global OPfunc,S_exact,S_calc,OPshem
    print('============================================================')
    print('                                                            ')
    print('                 Scheme Evaluation Routine                  ')
    print('                                                            ')
    print('============================================================')
    deltas()   #calculates dx, dy, dt
    exact_functions(S_exact,S_calc)
    for op in range(1,8):
        S_copy = S_calc.copy() # Copy the elements from S_calc to S_copy
        ADI(S_copy,op)
        np.savetxt(OPshem + '_function_' + OPfunc + '.txt',S_copy, fmt="%2.15f", delimiter=" ")
        np.savetxt('Exact_solution_' + OPfunc + '.txt',S_exact, fmt="%2.15f", delimiter=" ")
    return
#----------------------------------------------------------------
def ERROR_routine():
    global N_i,N_f,IF
    N_a = int(N_f/IJ)
    rms_v = np.linspace(N_i,N_f, N_a, dtype='int') + 1
    RMS_e = np.zeros((7,N_a),dtype = float)
    print('============================================================')
    print('                                                            ')
    print('                    Linear Test Problem                     ')
    print('                      RMS Calculation                       ')
    print('                                                            ')
    print('============================================================')
    for i in range(0,N_a):
        N = rms_v[i]
        print( 'Grid:', N)
        x, y = np.zeros((2, N),dtype = float)
        S_exact = np.zeros((N,N),dtype = float)
        S_calc = np.zeros((N,N),dtype = float)
        deltas()   #calculates dx, dy, dt
        exact_functions(S_exact,S_calc)
        for op in range(0,7):
            S_copy = S_calc.copy() # Copy the elements from S_calc to S_copy
            ADI(S_copy,op+1)
            RMS_e[op,i] = RMS(S_copy,S_exact)
    np.savetxt('RMS_function_' + OPfunc + '.txt',RMS_e, fmt="%2.15f", delimiter=" ")
    return
#----------------------------------------------------------------
#===================== Main ====================#
if(rms_flag): ERROR_routine()
    
if(sheme_flag): Scheme_routine()
