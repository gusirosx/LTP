#include <global.h>
void Scheme_routine()
{
    double **S_copy;
    x = malloc (N*sizeof (float*));
    y = malloc (N*sizeof (float*));
    S_exact = allocate2D(N,N);
    S_calc  = allocate2D(N,N);
    S_copy  = allocate2D(N,N);
    printf("\n============================================================");
    printf("\n                                                            ");
    printf("\n                 Scheme Evaluation Routine                  ");
    printf("\n                                                            ");
    printf("\n============================================================");
    deltas();
    exact_functions(S_exact,S_calc);
    copy2D(N,S_calc,S_copy);
    for(int op=1; op <= 7; op++)
    {
        copy2D(N,S_copy,S_calc); //Copy the elements from S_copy to S_calc
        ADI(op);
        w_fille(S_calc);
    }
    strcpy(name,"Exact_solution_");
    strcat(name,OPfunc);
    w_fille(S_exact);
    free(x); free(y);
    deallocate2D(S_calc,N); deallocate2D(S_exact,N);deallocate2D(S_copy,N);
}
//================================================================
void ERROR_routine()
{
    //Calculates the RMS error
    int i,j,op,N_a;
    N_a = N_f/IJ;
    FILE *OUT;
    int rms_v[N_a];
    double **RMS_e,**S_copy;
    RMS_e = allocate2D(7,N_a);
    printf("\n============================================================");
    printf("\n                                                            ");
    printf("\n                    Linear Test Problem                     ");
    printf("\n                      RMS Calculation                       ");
    printf("\n                                                            ");
    printf("\n============================================================");
    for(i=0; i < N_a; i++){
        rms_v[i] = (i+1)*N_a + 1;
    }
    for(i=0; i < N_a; i++)
    {
        N = rms_v[i];
        printf("\nGrid: %d", N);
        x = malloc (N*sizeof (float*));
        y = malloc (N*sizeof (float*));
        S_exact = allocate2D(N,N);
        S_calc  = allocate2D(N,N); S_copy  = allocate2D(N,N);
        deltas();
        exact_functions(S_exact, S_calc);
        copy2D(N,S_calc,S_copy);
        for(op=0; op < 7; op++)
        {
            copy2D(N,S_copy,S_calc); //Copy the elements from S_copy to S_calc
            ADI(op+1);
            RMS_e[op][i] = RMS(N,S_calc,S_exact);
        }
        free(x); free(y);
        deallocate2D(S_calc,N); deallocate2D(S_exact,N);deallocate2D(S_copy,N);
    }
    strcpy(name,"RMS_function_");
    strcat(name,OPfunc);
    strcat(name,".txt");
    OUT = fopen(name,"w");
    for(i=0; i < N_a; i++)
    {
        for(op=0; op < 7; op++)
        {
            fprintf(OUT,"  %2.10lf",RMS_e[op][i]);
        }
        fprintf(OUT,"\n");
    }
    fclose (OUT);
    deallocate2D(RMS_e,7);
}
