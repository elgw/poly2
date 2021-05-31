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

double * new_poly_rand(int n)
{
    // Polygon with n vertices at random location.
    double * P = malloc(2*n*sizeof(double));
    for(int kk = 0; kk<2*n; kk++)
    {
        P[kk] = (double) rand() / (double) RAND_MAX;
    }
    return P;
}

double * new_poly_square(int * n)
{
    double * P = malloc(10*sizeof(double));
    P[0] = 0; P[2] = 1; P[4] = 1; P[6] = 0; P[8] = .1;
    P[1] = 0; P[3] = 0; P[5] = 1; P[7] = 1; P[9] = .4;
    n[0] = 5; // 4 corners
    //poly_reverse(P, 4);
    return P;
}

void poly_hull_ut()
{

    printf(" -> poly_hull\n");

    // Should return null when less than 3 points
    int nH = -1;
    double * H = NULL;
    for(int kk = 0; kk<3; kk++)
    {
        H = poly_hull(NULL, kk, &nH);
        assert(H == NULL);
        assert(nH == 0);
    }

    // Should be able to handle random polygons without crashing
    for(int kk = 0; kk < 10000; kk++)
    {
        int nP = 1+rand() % 100;
        double * P = new_poly_rand(nP);
        //poly_print(stdout, P, nP);
        H = poly_hull(P, nP, &nH);
        free(P);
        if(H != NULL)
        {
            free(H);
        }
    }

    return;
}

void poly_area_ut()
{
    double sum = 0;
    printf(" -> poly_area\n");
    // Should return 0 for points and lines
    double * P = new_poly_rand(2);
    for(int kk = 0; kk<3; kk++)
    {
        assert(poly_area(P, kk) == 0);
    }
    free(P);

    // Test random input;
    for(int kk = 0; kk<10000; kk++)
    {
        int n = rand() % 1000;
        P = new_poly_rand(n);
        sum += poly_area(P, n);
        free(P);
    }

    int n = 0;
    P = new_poly_square(&n);
    assert(poly_area(P, 4) == 1);
    free(P);
    printf("\t%f\n", sum);
}

void poly_cbinter_ut()
{
    printf(" -> poly_cbinter\n");
    int n = 0;
    double * P = new_poly_square(&n);
    printf("\tcreating cbsquare2.svg\n");
    int nI = 0;
    double * I = poly_cbinterp(P, n, 2, &nI);
    poly_to_svg(I, nI, "cbsquare2.svg");
    free(I);

    printf("\tcreating cbsquare8.svg\n");
    I = poly_cbinterp(P, n, 8, &nI);
    poly_to_svg(I, nI, "cbsquare8.svg");
    free(I);
    free(P);

    return;
}

void poly_com_ut()
{
    printf("-- poly_com_ut\n");
    int n = 0;
    double * P = new_poly_square(&n);


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
    double * P = new_poly_square(&n);
    //n = 4; // ignore imperfection
    poly_translate(P, n, 50, -50);
    poly_mult(P, n, 10.0, 20.0);

    char * oname = malloc(50);
    int nrot = 32;
    for(int kk = 0; kk < nrot; kk++)
    {
        sprintf(oname, "rot%02d.svg", kk);
        poly_to_svg(P, n, oname);
        poly_rotate(P, n, 0, 0, 2*M_PI / (double) nrot);
    }
    free(oname);


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
    printf("  - Rotating by 0.4\n");
    poly_rotate(P, n, 0.0, 0.0, 0.4);
    //poly_print(stdout, P, n);
    poly_orientation(P, n);
    printf(" - Multiplying by 1, 1.01\n");
    poly_mult(P, n, 1.0, 1.01);
    C = poly_cov(P, n);
    //printf("C: [%f, %f; %f, %f]\n", C[0], C[1], C[1], C[2]);
    //getchar();
    free(C);
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
    fflush(stdout);
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

void counter_clockwise()
{
    printf(" -> counter_clockwise\n");
    int nP = 0;
    double * P = new_poly_square(&nP);
    poly_print(stdout, P, nP);
    poly_reverse(P, nP);
    poly_print(stdout, P, nP);
    printf("\tWriting counter_clockwise.svg\n");
    poly_to_svg(P, nP, "counter_clockwise.svg");
    free(P);
    return;
}

void non_simple()
{
    printf(" -> non_simple.svg\n");
    double * P = new_poly_rand(100);
    poly_to_svg(P, 100, "non_simple.svg");
    free(P);
    return;
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


    non_simple();
    counter_clockwise();
    poly_hull_ut();
    poly_area_ut();
    poly_cbinter_ut();
    poly_com_ut();
    poly_cov_ut();
    benchmark();


    return 0;
}
