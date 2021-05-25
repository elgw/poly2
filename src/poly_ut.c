#include <stdlib.h>
#include <stdio.h>

#include "poly.h"

double * new_square(int * n)
{
    double * P = malloc(8*sizeof(double));
    P[0] = 0; P[2] = 1; P[4] = 1; P[6] = 0;
    P[1] = 0; P[3] = 0; P[5] = 1; P[7] = 1;

    n[0] = 4; // 4 corners
    return P;
}

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
    free(P);
    int n = 0;
    P = new_square(&n);
    poly_print(stdout, P, n);
    free(P);
}

void poly_cbinter_ut()
{
    int n = 0;
    double * P = new_square(&n);
    printf("Square:\n");
    poly_print(stdout, P, n);
    int nI = 0;
    double * I = poly_cbinterp(P, n, 2, &nI);
    poly_to_svg(I, nI, "test2.svg");
    printf("Interpolated square:\n");
    poly_print(stdout, I, nI);
    free(I);

    I = poly_cbinterp(P, n, 8, &nI);
    poly_to_svg(I, nI, "test8.svg");
    printf("Interpolated square:\n");
    poly_print(stdout, I, nI);
    free(I);
    free(P);

    return;
}

int main()
{
    poly_area_ut();
    poly_cbinter_ut();
    return 0;
}
