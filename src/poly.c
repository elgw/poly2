#include "poly.h"

void poly_print(FILE * fid, double * P, int n)
{
    fprintf(fid, "P=[");
    for(int kk = 0; kk<n; kk++)
    {
        if(kk > 0)
        {
            fprintf(fid, ", ");
        }
        fprintf(fid, "[%f, %f]", P[2*kk], P[kk*kk+1 ]);
    }
    fprintf(fid, "]\n");
    printf("Area(P) = %f\n", poly_area(P, n));
    printf("Circ(P) = %f\n", poly_circ(P, n));
}

double poly_area(double * P, int n)
{
    // Rourke, Eq. 1.14
    if(n<3)
        return 0;

    double a = (P[2*(n-1)] + P[0])*(P[1]-P[2*(n-1)+1]);

    for(int kk = 0; kk+1<n; kk++)
    {
        a += (P[kk] + P[kk+2])*(P[kk+3]-P[kk+1]);
    }

    return a/2.0;
}

double d2(double * P0, double *P1)
{
    return sqrt( pow(P0[0] - P1[0], 2) + pow(P0[1] - P1[1], 2));
}

double poly_circ(double * P, int n)
{
    if(n<3)
    {
        return 0;
    }
    double c = d2(P+0, P+2*(n-1));
    for(int kk = 0; kk+1 < n; kk++)
    {
        c += d2(P+2*kk, P+2*(kk+1));
    }
    return c;
}
