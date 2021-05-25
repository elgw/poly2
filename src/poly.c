#include "poly.h"

/* Cubic spline interpolation of closed curve
 * https://mathworld.wolfram.com/CubicSpline.html
 *
 * asserts that the input domain is regularly spaced
 */


// Forward declarations for non-exported functions
static double * cbspline(double * y, int N, int f, int * np);
static double * vec_interlace(double * A, double * B, size_t N);
static double * vec_deinterlace(double * V, int stride, size_t N);
static void vec_show(double * V, int n);
static double d2(double * P0, double *P1);

static double * cbspline(double * y, int N, int f, int * np)
  /* y input vector
   * N number of elements in y
   * f upsampling factor
   * np[0] will be set to N*f
   */
{
  np[0] = N*f;

  // Solve the Toeplitz equation system. There are faster ways, see for example
  // Numerical Experience with a Superfast Real Toeplitz Solver, Gregory S. Ammar and William B. Gragg
  //

  gsl_vector * gdiag = gsl_vector_alloc(N);
  for(int kk = 0; kk<N; kk++)
    gsl_vector_set(gdiag, kk, 4);

  gsl_vector * ge = gsl_vector_alloc(N);
  for(int kk = 0; kk<N; kk++)
    gsl_vector_set(ge, kk, 1);

  gsl_vector * gb = gsl_vector_alloc(N);
  for(int kk = 0; kk<N; kk++)
  {
    gsl_vector_set(gb, kk, 3*(y[(kk+1) % N] - y[(N+kk-1) % N]));
  }

  gsl_vector * gx = gsl_vector_alloc(N);

  // int gsl_linalg_solve_symm_cyc_tridiag(const gsl_vector * diag, const gsl_vector * e, const gsl_vector * b, gsl_vector * x)
  gsl_linalg_solve_symm_cyc_tridiag(gdiag, ge, gb, gx);

  double * D = malloc(N*sizeof(double));
  for(int kk = 0; kk<N; kk++)
  {
    D[kk] = gsl_vector_get(gx, kk);
  }

  gsl_vector_free(gdiag);
  gsl_vector_free(ge);
  gsl_vector_free(gb);
  gsl_vector_free(gx);

  double * a = y;
  double * b = D;
  double * c = malloc(N*sizeof(double));
  for(int kk = 0; kk<N; kk++)
  {
    c[kk] = 3*(y[(kk+1) % N]-y[kk]) - 2*D[kk] - D[(kk+1) % N];
  }
  double * d = malloc(N*sizeof(double));
  for(int kk = 0; kk<N; kk++)
  {
    d[kk] = 2*(y[kk] - y[(kk+1) % N]) + D[kk] + D[(kk+1) % N];
  }

  double * X = malloc(N*f*sizeof(double));
  for(int kk = 0; kk<N; kk++)
  {
    for(int ll = 0; ll<f; ll++)
    {
      double xp = (double) ll/ (double) f;
      X[f*kk + ll] = a[kk] + b[kk]*xp + c[kk]*pow(xp, 2) + d[kk]*pow(xp, 3);

     // printf("x=%f, a=%f, b=%f, c=%f, d=%f -> %f\n", xp, a[kk], b[kk], c[kk], d[kk], X[f*kk+ll]);
    }
  }

  free(D);
  free(c);
  free(d);
  return X;
}




static double * vec_deinterlace(double * V, int stride, size_t N)
{
    double * D = malloc(N*sizeof(double));
    for(size_t kk = 0; kk<N; kk++)
    {
        D[kk] = V[stride*kk];
    }
    return D;
}

static double * vec_interlace(double * A, double * B, size_t N)
{
    double * V = malloc(2*N*sizeof(double));
    for(size_t kk = 0; kk<N; kk++)
    {
        V[2*kk] = A[kk];
        V[2*kk+1] = B[kk];
    }
    return V;
}

static void vec_show(double * V, int n)
{
    for(int kk = 0; kk<n; kk++)
    {
        printf("%f ", V[kk]);
    }
    printf("\n");
}


void poly_print(FILE * fid, double * P, int n)
{
    fprintf(fid, "P=[");
    for(int kk = 0; kk<n; kk++)
    {
        if(kk > 0)
        {
            fprintf(fid, ", ");
        }
        fprintf(fid, "[%f, %f]", P[2*kk], P[2*kk+1 ]);
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

static double d2(double * P0, double *P1)
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


double * poly_cbinterp(double * P, int n, int upsampling, int * N)
/* Interpolate the points in the polygon P
 * returns an interlaced array, x0, y0, x1, y1, ... with N[0] points
 * N[0] is approximately upsampling * n
 * The caller is responsible to free the returned array
 */
{
    if(n == 0){
        return NULL;
    }

    if(n < 3)
    {
        // No interpolation
        double * I = malloc(n*2*sizeof(double));
        memcpy(I, P, 2*n*sizeof(double));
        return I;
    }

    double * X = vec_deinterlace(P, 2, n);
    double * Y = vec_deinterlace(P+1, 2, n);
    assert(X[0] == P[0]);
    assert(Y[0] == P[1]);
    int np = 0;

    double * XX = cbspline(X, n, upsampling, &np);
    double * YY = cbspline(Y, n, upsampling, &np);

    N[0] = np;
    free(X);
    free(Y);
    double * R = vec_interlace(XX, YY, np);
    free(XX);
    free(YY);
    return R;

}

double * poly_bbx(double * P, int n)
{
    if(n < 1)
        return NULL;
    double * bbx = malloc(4*sizeof(double));
    double minx = P[0];
    double maxx = P[0];
    double miny = P[1];
    double maxy = P[1];
    for(int kk = 0; kk<n; kk++)
    {
        double x = P[2*kk];
        double y = P[2*kk+1];

        if(x > maxx)
        {
            maxx = x;
        }
        if(x < minx)
        {
            minx = x;
        }
        if(y>maxy)
        {
            maxy = y;
        }
        if(y<miny)
        {
            miny = y;
        }
    }
    bbx[0] = minx;
    bbx[1] = maxy;
    bbx[2] = minx;
    bbx[3] = maxy;
    return bbx;

}

double coordscale(double x, double * bbx, double w, double padding)
{
    x -= bbx[0];
    x /= (bbx[1]-bbx[0]);
    x = x*(w-2*padding) + padding;
    return x;
}

void poly_to_svg(double * P, int n, char * filename)
{

    int w = 512;
    int h = 512;
    double padding = 10;
    double * bbx = poly_bbx(P, n);

    cairo_surface_t *surface;
    cairo_t *cr;
    surface = cairo_svg_surface_create(filename, w, h);
    // surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 390, 60);
    cr = cairo_create(surface);

    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL,
                           CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 20.0);
    cairo_move_to(cr, 10.0, 20.0);
char * caption = malloc(128);
sprintf(caption, "%d", n);
    cairo_show_text(cr, caption);
free(caption);


cairo_set_source_rgb(cr, 0, 0, 0);
cairo_set_line_width(cr, 5);
double xold = P[2*(n-1)];
xold = coordscale(xold, bbx, (double) w, padding);
double yold = P[2*(n-1)+1];
yold = coordscale(yold, bbx+2, (double) h, padding);


for(int kk = 0; kk < n; kk++)
{
    printf("%f %f\n", xold, yold);
    cairo_move_to(cr, xold, yold);
    double x = P[2*kk];
    x = coordscale(x, bbx, (double) w, padding);
    double y = P[2*kk+1];
    y = coordscale(y, bbx+2, (double) h, padding);
    cairo_line_to(cr, x, y);
    xold = x;
    yold = y;
}
cairo_stroke(cr);

//cairo_surface_write_to_svg(surface, filename);

    cairo_destroy(cr);
    cairo_surface_destroy(surface);
}
