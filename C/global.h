//=============== Importing libraries ===============//
#define pi 3.14159265359
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
//================ GLOBAL PARAMETERS ================//
    double 
    **S_exact,          // Exact solution
    **S_calc;           // Calculated solution
    float
    *x,                 // Determine the x positions
    *y;                 // Determine the y positions  
 
    double 
    dx,
    dy,                 
    dt,                 
    TOL,                // Convergence tolerance    
    Pe;                 // Global Peclet Number

    float 
    lamb,               // Eigenvalue
    theta,              // Mesh offset angle
    alpha;
    
    int
    N;                  // Nodes in x and y direction
    
    char 
    name[60],           
    OPfunc[2],          // Select the Function
    OPshem[10];         // Select the scheme 
    
    bool 
    rms_flag,           // rms routine flag
    sheme_flag;         // scheme flag

    //====== Error Evaluation Routine ======//
    int
    N_i,                // Initial mesh
    N_f,                // Final mesh
    IJ;                 // Jump Interval
    //======== Auxiliary Variables ========//
    double 
    Ax,Bx,Cx,Dpx,
    Ay,By,Cy,Dpy,
    Ax2,Bx2,Xww,Dpx2,
    Ay2,By2,Yss,Dpy2;
    double 
    **S_a,
    **S_rhs;
//==================== FUNCTIONS ====================//
//Functions
void deltas(void);
void norm_function(int, double**);
void exact_functions(double**, double**);
void TDMA(int,double [*],double [*],double[*],double[*],double[*]);
void ADI(int);
void ADI_S(int);
void ADI_MOD(int);
double RMS(int, double**, double**);
void w_fille(double**);

// schemes
void scheme_Selection(int);
void RHS_schemes(int);
void LOADS_RHS();
void UNIFAES_RHS();

// Routines
void Scheme_routine();
void ERROR_routine();

// C Function
void copy2D(int, double**, double**);
void Print1D(int, float*);
void Print2D(int, double**);
double** allocate2D(int,int);
void deallocate2D(double**,int);
void readfile(void);
