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

void poly_com_ut()
{
    printf("-- poly_com_ut\n");
    int n = 0;
    double * P = new_square(&n);
    poly_print(stdout, P, n);
    double * C = poly_com(P, n);
    printf("Com: (%f, %f)\n", C[0], C[1]);
    // Translate in x -- com_y does not change
    printf("Translating by (0, 0.2)\n");
    poly_translate(P, n, 0, .2);
    poly_print(stdout, P, n);
    poly_to_svg(P, n, "trans_square.svg");
    free(C);
    C = poly_com(P, n);
    printf("Com: (%f, %f)\n", C[0], C[1]);
    // Translate in y -- com_x does not change
    // Rotate around com -- com does not change
    free(C);
    free(P);
}

void poly_cov_ut()
{
    printf("-- poly_cov_ut\n");
    int n = 0;
    double * P = new_square(&n);
    poly_print(stdout, P, n);
    double * C = poly_cov(P, n);
    printf("  - Translating by (.1, .2)\n");
    poly_translate(P, n, .1, .2);
    C = poly_cov(P, n);
    printf("  - Centering\n");
    double * com = poly_com(P, n);
    poly_translate(P, n, -com[0], -com[1]);
    free(com);
    C = poly_cov(P, n);
    printf("  - Rotating\n");
    poly_rotate(P, n, 0.0, 0.0, 0.2);
    poly_print(stdout, P, n);
    C = poly_cov(P, n);
    poly_to_svg(P, n, "temp.svg");
    free(P);
}

int main()
{
    poly_area_ut();
    poly_cbinter_ut();
    poly_com_ut();
    poly_cov_ut();
    return 0;
}
