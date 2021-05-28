#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "poly.h"

static double timespec_diff(struct timespec* end, struct timespec * start)
{
    double elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}

double * new_square(int * n)
{
    double * P = malloc(8*sizeof(double));
    P[0] = 0; P[2] = 1; P[4] = 1; P[6] = 0;
    P[1] = 0; P[3] = 0; P[5] = 1; P[7] = 1;
    n[0] = 4; // 4 corners
    //poly_reverse(P, 4);
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
    poly_mult(P, n, 100.0, 100.0);
    poly_print(stdout, P, n);
    double * C = poly_cov(P, n);
    printf("  - Translating by (.1, .2)\n");
    free(C);
    poly_translate(P, n, .1, .2);
    C = poly_cov(P, n);
    free(C);
    printf("  - Centering\n");
    double * com = poly_com(P, n);
    poly_translate(P, n, -com[0], -com[1]);
    free(com);
    C = poly_cov(P, n);
    free(C);
    printf("  - Rotating by 0.2\n");
    poly_rotate(P, n, 0.0, 0.0, 0.2);
    poly_print(stdout, P, n);
    poly_orientation(P, n);
    C = poly_cov(P, n);
    free(C);
    poly_mult(P, n, 1.0, 1.01);
    poly_orientation(P, n);
    poly_to_svg(P, n, "temp.svg");
    free(P);
}

void benchmark()
{

    #ifdef NDEBUG
    int N = 1000000;
    #else
    int N = 100;
    #endif
    int V = 32;
    printf("Benchmarking with %d polygons with %d vertices\n", N, V);
    double * P = malloc(V*2*sizeof(double));
    double atotal = 0;
    struct timespec tstart, tend;
    clock_gettime(CLOCK_REALTIME, &tstart);
    //printf("%d\n", kk);
    for(int ll = 0 ; ll<2*V; ll++)
    {
        P[ll] = rand() / (double) RAND_MAX;
    }
    for(int kk = 0; kk < N; kk++)
    {
        poly_props * p = poly_measure(P, V);
        atotal += p->Area;
        poly_props_free(&p);
    }
    printf("Total area: %f\n", atotal);
    free(P);
    clock_gettime(CLOCK_REALTIME, &tend);
    double dt = timespec_diff(&tend, &tstart);
    printf("Took: %f s, i.e. %f polygons / s\n", dt, (double) N / dt);
#ifdef NDEBUG
    FILE * fout = fopen("benchmark.txt", "a");
    time_t now;
    time(&now);
    fprintf(fout, "%s\n", ctime(&now));
    fprintf(fout, "Took: %f s, i.e. %f polygons / s\n", dt, (double) N / dt);
    fclose(fout);
    #endif
}

int main(int argc, char ** argv)
{
    if(argc > 1)
    {
        int n = (argc-1);
        if( n%2 != 0)
        {
            printf("An even number of values has to be give, (x0, y0), (x1, y1), ...\n");
            exit(1);
        }
        n /= 2;
        double * P = malloc(2*n*sizeof(double));
        for(int kk = 1; kk < argc; kk++)
        {
            P[kk-1] = atof(argv[kk]);
        }

        poly_print(stdout, P, n);
        poly_props * props = poly_measure(P, n);
        poly_props_print(stdout, props);
        poly_props_free(&props);
        printf("Writing to argv.svg\n");
        poly_to_svg(P, n, "argv.svg");
        free(P);
        return 0;
    }

    benchmark();

    poly_area_ut();
    poly_cbinter_ut();
    poly_com_ut();
    poly_cov_ut();
    return 0;
}
