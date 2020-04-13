#include <global.h>
//================================================================
void norm_function(int n, double **A)
{
    //Normalization by the max value of a 2D array
    int i,j;
    double max = 0.0, min=0.0, f;
    printf("\n");
    for(i=0; i < n; i++)
    {
        for(j=0; j < n; j++)
        {
            if (A[i][j] > max)
                max = A[i][j];
            if (A[i][j] < min)
                min = A[i][j];
        }
    }
    f = max - min;
    //printf("%lf",f);
    for(i=0; i < n; i++)
    {
        for(j=0; j < n; j++)
        {
            A[i][j] = A[i][j]/f;
        }
    }
}
//================================================================
void deltas(void)
{
    int i,j;
    alpha = theta*pi/180.0;
    dx = 1.0/(N-1.0);
    dy = dx;
    dt = 0.98*(4*dx)/(3*Pe*cos(alpha));
    //Determine the x and y positions for solutions A, B, C, D
    for(i=0; i < N; i++)
    {
        x[i] = -0.5 + i*dx;
        y[i] = -0.5 + i*dy;
    }
}    
//================================================================
void TDMA(int n,double A[],double B[],double C[],double D[],double TT[n])
{
    int i,j;
    double W[n],G[n];
    W[0] = C[0]/B[0];
    G[0] = D[0]/B[0];
    for(i=0;i<n-1;i++)
        W[i] = C[i]/(B[i] - A[i-1]*W[i-1]);
    for(i=0;i<n;i++)
        G[i] = (D[i] - A[i-1]*G[i-1])/(B[i] - A[i-1]*W[i-1]);
    TT[n-1] = G[n-1];
    for(i=n-1;i>=0;i--)
        TT[i-1] = G[i-1] - W[i-1]*TT[i];
}
//================================================================
void exact_functions(double **S_sol, double **S_c)
{
    double C,S,M;
    int i,j;
    //Calculation of the exact solution - TYPE A, B, C, D
    if (strcmp(OPfunc,"A") == 0)
        C = (Pe - pow((Pe*Pe + 4.0*lamb*lamb),0.5))*0.5;
    else if (strcmp(OPfunc,"B") == 0)
        C = (Pe + pow((Pe*Pe + 4.0*lamb*lamb),0.5))*0.5;
    else if (strcmp(OPfunc,"C") == 0)
        C = (Pe - pow((Pe*Pe - 4.0*lamb*lamb),0.5))*0.5;
    else if (strcmp(OPfunc,"D") == 0)
        C = (Pe + pow((Pe*Pe - 4.0*lamb*lamb),0.5))*0.5;
    else if (strcmp(OPfunc,"CD") == 0 || strcmp(OPfunc,"DC") == 0)
        C = sqrt(lamb*lamb - Pe*Pe*0.25);
    //Loop for functions A and B
    if (strcmp(OPfunc,"A") == 0 || strcmp(OPfunc,"B") == 0)
    {
        for(i=0; i < N; i++)
        {
            for(j=0; j < N; j++)
            {
                S = x[j]*cos(alpha) + y[i]*sin(alpha);
                M = y[i]*cos(alpha) - x[j]*sin(alpha);
                S_sol[i][j] = exp(C*S)*sin(lamb*M);
            }
        }
    }
    //Loop for functions C and D
    if (strcmp(OPfunc,"C") == 0 || strcmp(OPfunc,"D") == 0)
    {
        for(i=0; i < N; i++)
        {
            for(j=0; j < N; j++)
            {
                S = x[j]*cos(alpha) + y[i]*sin(alpha);
                M = y[i]*cos(alpha) - x[j]*sin(alpha);
                S_sol[i][j] = exp(C*S)*exp(lamb*M);
            }
        }
    }
    //Loop for functions CD
    if (strcmp(OPfunc,"CD") == 0)
    {
        for(i=0; i < N; i++)
        {
            for(j=0; j < N; j++)
            {
                S = x[j]*cos(alpha) + y[i]*sin(alpha);
                M = y[i]*cos(alpha) - x[j]*sin(alpha);
                S_sol[i][j] = exp(Pe*S*0.5)*sin(C*S)*exp(lamb*M);
            }
        }
    }
    //Loop for functions DC
    if (strcmp(OPfunc,"DC") == 0)
    {
        for(i=0; i < N; i++)
        {
            for(j=0; j < N; j++)
            {
                S = x[j]*cos(alpha) + y[i]*sin(alpha);
                M = y[i]*cos(alpha) - x[j]*sin(alpha);
                S_sol[i][j] = exp(Pe*S*0.5)*cos(C*S)*exp(lamb*M);
            }
        }
    }
    norm_function(N,S_sol);
    for(i=0; i < N; i++)
    {
        for(j=0; j < N; j++)
        {
            S_c[i][j] = 0.0; 
        }
    }
    for(i=0; i < N; i++)
    {
        S_c[0][i]   = S_sol[0][i];
        S_c[N-1][i] = S_sol[N-1][i];
        S_c[i][0]   = S_sol[i][0];
        S_c[i][N-1] = S_sol[i][N-1];
    }
}
//================================================================
void ADI(int op)
{
    //select which ADI to use
    if(op >= 1 && op<= 5)
        ADI_S(op);
    else if(op >= 6 && op <= 7)
        ADI_MOD(op);
    strcpy(name,OPshem);
    strcat(name,"_function_");
    strcat(name,OPfunc);
}
//================================================================
void ADI_S(int op)//Standard ADI
{
    int i,j,l,k,icnt,n = N-2;
    double F,VTOL;
    double A[n-1],B[n],C[n-1],D[n],TS[n];
    S_a   = allocate2D(N,N);
    S_rhs = allocate2D(N,N);
    scheme_Selection(op);    
    for(icnt=0; icnt < 1000000; icnt++)
    {
        //X Direction
        copy2D(N,S_calc,S_a);  //Copy the elements from S_calc to S_a
        RHS_schemes(op);
        for(i=1; i < N-1; i++)
        {
            for(j=1; j < N-1; j++)
            {
                if(j == 1)
                {
                    B[j-1] = Bx;
                    C[j-1] = Cx;
                    D[j-1] = - Ay*S_a[i-1][j] + Dpx*S_a[i][j] - Cy*S_a[i+1][j] - Ax*S_a[i][j-1] + S_rhs[i][j];
                }
                else if(j > 1 && j < (N-2))
                {
                    A[j-2] = Ax;
                    B[j-1] = Bx;
                    C[j-1] = Cx;
                    D[j-1] = - Ay*S_a[i-1][j] + Dpx*S_a[i][j] - Cy*S_a[i+1][j] + S_rhs[i][j];
                }
                else if(j == (N - 2))
                {
                    A[j-2] = Ax;
                    B[j-1] = Bx;
                    D[j-1] = - Ay*S_a[i-1][j] + Dpx*S_a[i][j] - Cy*S_a[i+1][j] - Cx*S_a[i][j+1] + S_rhs[i][j];
                }
            }
            TDMA(n,A,B,C,D,TS);
            for(j=1; j<N-1; j++)
            {
                S_calc[i][j] = TS[j-1]; //Stores TDMA values ​​in matrix
            }
        }
        //Direção Y
        copy2D(N,S_calc,S_a);
        RHS_schemes(op);
        for(j=1; j < N-1; j++)
        {
            for(i=1; i < N-1; i++)
            {
                if(i == 1)
                {
                    B[i-1] = By;
                    C[i-1] = Cy;
                    D[i-1] = - Ax*S_a[i][j-1] + Dpy*S_a[i][j] - Cx*S_a[i][j+1] - Ay*S_a[i-1][j] + S_rhs[i][j];
                }
                else if(i>1 && i<(N-2))
                {
                    A[i-2] = Ay;
                    B[i-1] = By;
                    C[i-1] = Cy;
                    D[i-1] = - Ax*S_a[i][j-1] + Dpy*S_a[i][j] - Cx*S_a[i][j+1] + S_rhs[i][j];
                }
                else if(i == (N - 2))
                {
                    A[i-2] = Ay;
                    B[i-1] = By;
                    D[i-1] = - Ax*S_a[i][j-1] + Dpy*S_a[i][j] - Cx*S_a[i][j+1] - Cy*S_a[i+1][j] + S_rhs[i][j];
                }
            }
            TDMA(n,A,B,C,D,TS);
            for(i=1;i<N-1;i++)
            {
                S_calc[i][j] = TS[i-1];
            }
        }
        VTOL = 0.0;
        for(l=0; l < N ;l++){
            for(k=0; k < N ;k++){
                F = fabs(S_calc[l][k] - S_a[l][k])/dt;
                if (F > VTOL){
                    VTOL = F;
                }
            }
        }
        if (VTOL <= TOL){
            printf("C.R.I: %d %s",icnt,OPshem);
            break;
        }
    }
    deallocate2D(S_a,N);
}
//================================================================
void ADI_MOD(int op)//Modified ADI for WW and SS terms
{
    int i,j,l,k,icnt,n = N-2;
    double F,VTOL;
    double A[n-1],B[n],C[n-1],D[n],TS[n];
    S_a   = allocate2D(N,N);
    S_rhs = allocate2D(N,N);
    scheme_Selection(op);    
    for(icnt=0; icnt < 1000000; icnt++)
    {
        //X Direction
        copy2D(N,S_calc,S_a);//Copy the elements from S_calc to S_a
        for(i=1; i < N-1; i++)
        {
            if(i == 1)
            {
                for(j=1; j < N-1; j++)
                {
                    if(j == 1)
                    {
                        B[j-1] = Bx2;
                        C[j-1] = Cx;
                        D[j-1] = - Ay2*S_a[i-1][j] + Dpx2*S_a[i][j] - Cy*S_a[i+1][j] - Ax2*S_a[i][j-1];
                    }
                    if(j>1 && j<(N-2))
                    {
                        A[j-2] = Ax;
                        B[j-1] = Bx;
                        C[j-1] = Cx;
                        D[j-1] = - Ay2*S_a[i-1][j] + Dpx2*S_a[i][j] - Cy*S_a[i+1][j] + Xww*S_a[i][j-2];
                    }
                    if(j == (N - 2))
                    {
                        A[j-2] = Ax;
                        B[j-1] = Bx;
                        D[j-1] = - Ay2*S_a[i-1][j] + Dpx2*S_a[i][j] - Cy*S_a[i+1][j] - Cx*S_a[i][j+1] + Xww*S_a[i][j-2];
                    }
                }
            }
            else if(i > 1)
            {
                for(j=1; j < N-1; j++)
                {
                    if(j == 1)
                    {
                        B[j-1] = Bx2;
                        C[j-1] = Cx;
                        D[j-1] = - Ay*S_a[i-1][j] + Dpx*S_a[i][j] - Cy*S_a[i+1][j] - Ax2*S_a[i][j-1] + Yss*S_a[i-2][j];
                    }
                    else if(j>1 && j<(N-2))
                    {
                        A[j-2] = Ax;
                        B[j-1] = Bx;
                        C[j-1] = Cx;
                        D[j-1] = - Ay*S_a[i-1][j] + Dpx*S_a[i][j] - Cy*S_a[i+1][j] + Xww*S_a[i][j-2] + Yss*S_a[i-2][j];
                    }
                    else if(j == (N - 2))
                    {
                        A[j-2] = Ax;
                        B[j-1] = Bx;
                        D[j-1] = - Ay*S_a[i-1][j] + Dpx*S_a[i][j] - Cy*S_a[i+1][j] - Cx*S_a[i][j+1] + Xww*S_a[i][j-2] + Yss*S_a[i-2][j];
                    }
                }
            }
            TDMA(n,A,B,C,D,TS);
            for(j=1;j<N-1;j++)
            {
                S_calc[i][j] = TS[j-1]; //Stores TDMA values ​​in matrix
            }
        }
        //Direção Y
        copy2D(N,S_calc,S_a);
        for(j=1; j < N-1; j++)
        {
            if(j == 1)
            {
                for(i=1; i < N-1; i++)
                {
                    if(i == 1)
                    {
                        B[i-1] = By2;
                        C[i-1] = Cy;
                        D[i-1] = - Ax2*S_a[i][j-1] + Dpy2*S_a[i][j] - Cx*S_a[i][j+1] - Ay2*S_a[i-1][j];
                    }
                    else if(i>1 && i<(N-2))
                    {
                        A[i-2] = Ay;
                        B[i-1] = By;
                        C[i-1] = Cy;
                        D[i-1] = - Ax2*S_a[i][j-1] + Dpy2*S_a[i][j] - Cx*S_a[i][j+1]+ Yss*S_a[i-2][j];
                    }
                    else if(i == (N - 2))
                    {
                        A[i-2] = Ay;
                        B[i-1] = By;
                        D[i-1] = - Ax2*S_a[i][j-1] + Dpy2*S_a[i][j] - Cx*S_a[i][j+1] - Cy*S_a[i+1][j] + Yss*S_a[i-2][j];
                    }
                }
            }
            else if(j > 1)
            {
                for(i=1; i < N-1; i++)
                {
                    if(i == 1)
                    {
                        B[i-1] = By2;
                        C[i-1] = Cy;
                        D[i-1] = - Ax*S_a[i][j-1] + Dpy*S_a[i][j] - Cx*S_a[i][j+1] - Ay2*S_a[i-1][j] + Xww*S_a[i][j-2];
                    }
                    else if(i>1 && i<(N-2))
                    {
                        A[i-2] = Ay;
                        B[i-1] = By;
                        C[i-1] = Cy;
                        D[i-1] = - Ax*S_a[i][j-1] + Dpy*S_a[i][j] - Cx*S_a[i][j+1] + Xww*S_a[i][j-2] + Yss*S_a[i-2][j];
                    }
                    else if(i == (N - 2))
                    {
                        A[i-2] = Ay;
                        B[i-1] = By;
                        D[i-1] = - Ax*S_a[i][j-1] + Dpy*S_a[i][j] - Cx*S_a[i][j+1] - Cy*S_a[i+1][j] + Xww*S_a[i][j-2] + Yss*S_a[i-2][j];
                    }
                }
            }
            TDMA(n,A,B,C,D,TS);
            for(i=1;i<N-1;i++)
            {
                S_calc[i][j] = TS[i-1];
            }
        }
        VTOL = 0.0;
        for(l=0; l < N ;l++){
            for(k=0; k < N ;k++){
                F = fabs(S_calc[l][k] - S_a[l][k])/dt;
                if (F > VTOL){
                    VTOL = F;
                }
            }
        }
        if (VTOL <= TOL){
            printf("C.R.I: %d %s",icnt,OPshem);
            break;
        }
    }
    deallocate2D(S_a,N);
}   
//================================================================
double RMS(int N, double **exat, double **aprox)
{
    //Calculating the RMS of a matrix
    double e_rms = 0.0,Ei;
    for(int i=0; i<N-1 ;i++)
    {
        for(int j=0; j<N-1 ;j++)
        {
            Ei = (aprox[i][j]- exat[i][j]);
            e_rms = e_rms + Ei*Ei;
        }
    }
    Ei = (e_rms/((N-2)*(N-2)));
    e_rms = sqrt(Ei);
    return e_rms;
}
//================================================================
void readfile(void)
{
    char line[13][100];
    FILE *fptr = NULL; 
    int i = 0;
    if ((fptr = fopen("LTP.ltp", "r")) == NULL) {
        printf("Error! opening file");
        exit(1);
    }
    while(fgets(line[i],100, fptr))
    {  
        if (*line[i] == '=' ) // Ignore lines starting with "="
            continue;        
        strtok(line[i], " "); // Breaking the string
        i++;
    }
    rms_flag   = false;
    sheme_flag = false;
    strcpy(OPfunc,line[0]);
    if(strcmp(line[1],"true") == 0) sheme_flag = true;
    Pe    = atof(line[2]);
    lamb  = atof(line[3]);
    theta = atof(line[4]);
    TOL   = atof(line[5]);
    N     = atoi(line[6]) + 1;
    if(strcmp(line[7],"true") == 0) rms_flag   = true;
    N_i   = atoi(line[8]);
    N_f   = atoi(line[9]);
    IJ    = atoi(line[10]);    
}
//================================================================
void w_fille(double **mat)
{
    FILE *OUT;
    int i,j;
    strcat(name,".txt");
    OUT = fopen(name,"w");
    for(j=0; j < N; j++)
    {
        for(i=0; i < N; i++)
        {
            if (mat[j][i] > 0.0 || mat[j][i] == 0.0)
                fprintf(OUT,"  %.10lf",mat[j][i]);
            if (mat[j][i] < 0.0)
                fprintf(OUT," %.10lf",mat[j][i]);
        }
        fprintf(OUT,"\n");
    }
    fprintf(OUT,"\n");
    fclose (OUT);
}
//================================================================
void Print1D(int n, float *A)
{
    int i;  
    printf("\n");
    for(i=0; i < n; i++)
    {
        printf("%.1f ",A[i]);
    }
    printf("\n");   
}
//================================================================
void Print2D(int n, double **A)
{
    int i,j;
    printf("\n");
    for(i=0; i < n; i++)
    {
        for(j=0; j < n; j++)
        {
            if (A[i][j] > 0.0 || A[i][j] == 0.0)
                printf("  %3.4lf",A[i][j]);
            if (A[i][j] < 0.0)
                printf(" %3.4lf",A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}
//================================================================
//Allocate a 2D array
double** allocate2D(int n, int m)
{
    double **A = malloc (n * sizeof (double *));
    for (int i = 0; i < n; ++i)
        A[i] = malloc (m * sizeof (double));
    if(A == NULL)                     
    {
        printf("Error! memory not allocated.");
        exit(1);
    }
    return A;
}
//================================================================
//deallocate a 2D array
void deallocate2D(double** arr2D,int rows)
{
    for(int i=0; i<rows;i++)
        free(arr2D[i]);
    free(arr2D);
}  
//================================================================
void copy2D(int n,double** source, double** dest)
{
    for (int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            dest[i][j] = source[i][j];
        }
    }
}
//================================================================
