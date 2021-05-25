#include <stdlib.h>
#include <stdio.h>

#include "poly.h"


void poly_area_ut()
{
    double * P = malloc(12*sizeof(double));
    printf("poly_area(NULL) = %f\n", poly_area(NULL, 0));
    P[0] = 0; P[1] = 0;
    printf("poly_area([0,0]) = %f\n", poly_area(P, 1));
    P[2] = 1; P[3] = 1;
    printf("poly_area([[0,0], [1, 1]]) = %f\n", poly_area(P, 2));
    P[4] = 1;
    P[5] = 0;
    printf("poly_area([[0,0], [1, 1], [1,0]]) = %f\n", poly_area(P, 3));
    poly_print(stdout, P, 3);
    P[0] = 0; P[2] = 1; P[4] = 1; P[6] = 1;
    P[1] = 0; P[3] = 0; P[5] = 1; P[7] = 0;
    poly_print(stdout, P, 4);
    free(P);
}

int main()
{
    poly_area_ut();
    return 0;
}
