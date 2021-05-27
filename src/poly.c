#include "poly.h"

/* Cubic spline interpolation of closed curve
 * https://mathworld.wolfram.com/CubicSpline.html
 *
 * asserts that the input domain is regularly spaced
 */


// Forward declarations for non-exported functions
static double * cbspline(const double * y, int N, int f, int * np);
static double * vec_interlace(const double * A, double * B, size_t N);
static double * vec_deinterlace(const double * V, int stride, size_t N);
//static void vec_show(const double * V, int n);
static double d2(const double * P0, const double *P1);
static void eigenvector_sym_22(double a, double b, double c, double l, double * v0, double * v1);
static void eigenvalue_sym_22(double a, double b, double c, double * l0, double * l1);
double poly_orientation_with_COV(double * COV);

poly_props * poly_props_new()
{
    poly_props * props = malloc(sizeof(poly_props));
    props->Centroid = NULL;
    props->BoundingBox = NULL;
    props->measured = 0;
    props->ConvexArea = -1;
    props->Solidity = -1;
    return props;
}

void poly_props_free(poly_props ** PP)
{
    poly_props * P = PP[0];
    if(P->Centroid != NULL)
    {
        free(P->Centroid);
    }
    if(P->BoundingBox != NULL)
    {
        free(P->BoundingBox);
    }
    free(PP[0]);
    PP[0] = NULL;
}

void poly_props_print(FILE * fout, poly_props * props)
{
    fprintf(fout, "-- Polygon properties -- \n");
    if(props->measured != 1)
    {
        fprintf(fout, "  ! No measurements recorded\n");
        return;
    }
    fprintf(fout, "  Vertices: %d\n", props->nVertices);
    fprintf(fout, "  Area: %f\n", props->Area);
    if(props->Centroid != NULL)
    {
    fprintf(fout, "  Centroid: (%f, %f)\n",
            props->Centroid[0], props->Centroid[1]);
    }
    if(props->BoundingBox != NULL)
    {
    fprintf(fout, "  BoundingBox: (%f, %f, %f, %f)\n",
            props->BoundingBox[0], props->BoundingBox[1],
            props->BoundingBox[2], props->BoundingBox[3]);
    }
    assert(props->Centroid[0] >= props->BoundingBox[0]);
    assert(props->Centroid[0] <= props->BoundingBox[1]);
    assert(props->Centroid[1] >= props->BoundingBox[2]);
    assert(props->Centroid[1] <= props->BoundingBox[3]);

    fprintf(fout, "  MajorAxisLength: %f\n", props->MajorAxisLength);
    fprintf(fout, "  MinorAxisLength: %f\n", props->MinorAxisLength);
    fprintf(fout, "  Eccentricity: %f\n", props->Eccentricity);

    double ori1 = props->Orientation;
    double ori2 = ori1 + M_PI;

    fprintf(fout, "  Orientation: %f (or %f)\n", ori1, ori2);

    fprintf(fout, "  TODO: ConvexArea: %f\n", props->ConvexArea);
    fprintf(fout, "  Circularity: %f\n", props->Circularity);
    fprintf(fout, "  EquivDiameter: %f\n", props->EquivDiameter);
    fprintf(fout, "  TODO: Solidity: %f\n", props->Solidity);
    fprintf(fout, "  Perimeter: %f\n", props->Perimeter);
    return;
}

void poly_reverse(double * P, int n)
{
    // Reverse the order in P
    int mid = n/2;
    for(int from = 0; from < mid; from++)
    {
        int to = n-from-1;
        double t0 = P[2*to];
        double t1 = P[2*to + 1];
        P[2*to] = P[2*from];
        P[2*to + 1] = P[2*from + 1];
        P[2*from] = t0;
        P[2*from + 1] = t1;
    }
    return;
}

poly_props * poly_measure(const double * P, int n)
{
    poly_props * props = poly_props_new();
    props->measured = 1;
    props->nVertices = n;
    props->Area = poly_area(P, n);
    props->Perimeter = poly_circ(P, n);
    // etc ...
    props->BoundingBox = poly_bbx(P, n);


    props->Circularity = (4.0*props->Area*M_PI)/pow(props->Perimeter, 2);
    props->EquivDiameter = sqrt(4.0*props->Area/M_PI);
    //props->Solidity = props->Area/props->ConvexArea;
    double * COV = poly_cov(P, n);
    //printf("Cov = [%f, %f; %f, %f]\n", COV[0], COV[1], COV[1], COV[2]);
    props->Orientation = poly_orientation_with_COV(COV);
    props->Centroid = poly_com(P, n);

    double l0, l1;
    eigenvalue_sym_22(COV[0], COV[1], COV[2], &l0, &l1);
    props->MajorAxisLength = 4*sqrt(l0);
    props->MinorAxisLength = 4*sqrt(l1);
    props->Eccentricity = sqrt(1-l1/l0);
    free(COV);
    return props;
}

// Multiply by the value v
void poly_mult(double * P, int n, double vx, double vy)
{
    for(int kk = 0; kk<n; kk++)
    {
        P[2*kk] *= vx;
        P[2*kk+1] *= vy;
    }
}

void poly_rotate(double * P, int n , double x0, double y0, const double theta)
{
    double st = sin(theta);
    double ct = cos(theta);

    // rotation around the point (x0, y0)
    for(int kk = 0; kk<n; kk++)
    {
        double x = P[2*kk];
        double y = P[2*kk+1];
        x = x-x0;
        y = y-y0;
        double rx =  ct*x + st*y;
        double ry = -st*x + ct*y;
        rx += x0;
        ry += y0;
        P[2*kk] = rx;
        P[2*kk+1] = ry;
    }
}

void poly_translate(double * P, int n, double dx, double dy)
{
    for(int kk = 0; kk<n; kk++)
    {
        P[2*kk] += dx;
        P[2*kk+1] += dy;
    }
}

static double * cbspline(const double * y, int N, int f, int * np)
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

  const double * a = y;
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


static double * vec_deinterlace(const double * V, int stride, size_t N)
{
    double * D = malloc(N*sizeof(double));
    for(size_t kk = 0; kk<N; kk++)
    {
        D[kk] = V[stride*kk];
    }
    return D;
}

static double * vec_interlace(const double * A, double * B, size_t N)
{
    double * V = malloc(2*N*sizeof(double));
    for(size_t kk = 0; kk<N; kk++)
    {
        V[2*kk] = A[kk];
        V[2*kk+1] = B[kk];
    }
    return V;
}

/*
static void vec_show(const double * V, int n)
{
    for(int kk = 0; kk<n; kk++)
    {
        printf("%f ", V[kk]);
    }
    printf("\n");
}
*/

void poly_print(FILE * fid, const double * P, int n)
{
    fprintf(fid, "P=[");
    for(int kk = 0; kk<n; kk++)
    {
        if(kk > 0)
        {
            fprintf(fid, "; ");
        }
        fprintf(fid, "[%f, %f]", P[2*kk], P[2*kk+1 ]);
    }
    fprintf(fid, "]\n");
    //printf("Area(P) = %f\n", poly_area(P, n));
    //printf("Circ(P) = %f\n", poly_circ(P, n));
}


void com_accumulate(double * com, const double * p, const double * q)
{
    // Centre of mass, contour integral between p and q
    // Using Greens formula
    // for x:
    // dw = x dx dy
    // w = -1/2 xy dx

    // for y
    // dw = y dx dy
    // w = -1/2 xy dy

    double px = p[0];
    double py = p[1];
    double qx = q[0];
    double qy = q[1];
    double dx = qx - px;
    double dy = qy - py;

    double alpha = px*py;
    double beta = dx*py + dy*px;
    double gamma = dx*dy;

    double comx = -dx*(alpha + beta/2.0 + gamma/3.0);
    double comy = dy*(alpha + beta/2.0 + gamma/3.0);

    //printf("(%f,%f) -> (%f, %f) com: (%f, %f)\n", p[0], p[1], q[0], q[1], comx, comy);

    com[0] += comx;
    com[1] += comy;
}

double poly_accumulate_M20(const double * p, const double * q)
{
    // dP/dy = -x*x
    // P = x*x*y/2
    // 1/2 II x*x*y dx
    // x = px + t*Dx
    // y = py + t*Dy
    double px = p[0];
    double py = p[1];
    double qx = q[0];
    double qy = q[1];
    double Dx = qx - px;
    double Dy = qy - py;

    double c0 = Dx*px*px*py;
    double c1 = Dx*(2*Dx*px*py + Dy*px*px);
    double c2 = Dx*(Dx*Dx*py + 2*Dx*Dy*px);
    double c3 = Dx*Dx*Dx*Dy;
    double M20 = -(c0 + c1/2.0 + c2/3.0 + c3/4.0);
    return M20;
    }

double poly_accumulate_M11(const double * p, const double * q)
{
    // dQ/dy = x*y
    // Q = x*x*y/2
    // 1/2 II x*x*y dy
    // x = px + t*Dx
    // y = py + t*Dy
    double px = p[0];
    double py = p[1];
    double qx = q[0];
    double qy = q[1];
    double Dx = qx - px;
    double Dy = qy - py;

    double c0 = px*px*py;
    double c1 = (Dy*px*px + 2*Dx*px*py);
    double c2 = (Dx*Dx*py + 2*Dx*Dy*px);
    double c3 = Dx*Dx*Dy;
    double M11 = 0.5*Dy*(c0 + c1/2.0 + c2/3.0 + c3/4.0);
    return M11;
}

double poly_accumulate_M02(const double * p, const double * q)
{
    // dQ/dx = y*y
    // Q = x*y*y
    // II x*y*y dy
    // x = px + t*Dx
    // y = py + t*Dy
    double px = p[0];
    double py = p[1];
    double qx = q[0];
    double qy = q[1];
    double Dx = qx - px;
    double Dy = qy - py;

    double c0 = Dy*px*py*py;
    double c1 = Dy*(2*Dy*px*py + Dx*py*py);
    double c2 = Dy*(Dy*Dy*px + 2*Dy*Dx*py);
    double c3 = Dy*Dx*Dy*Dy;
    double M02 = (c0 + c1/2.0 + c2/3.0 + c3/4.0);
    return M02;
}

double * poly_com(const double * P, int n)
{
    if(n == 0)
    {
        return NULL;
    }
    double * cm = malloc(2*sizeof(double));
    cm[0] = 0; cm[1] = 0;
    if(n < 3)
    {
        for(int kk = 0; kk<n; kk++)
        {
            cm[0] += P[2*kk];
            cm[1] += P[2*kk+1];
        }
        cm[0] /= (double) n;
        cm[1] /= (double) n;
        return cm;
    }

    const double * p = P + (2*n-2);
    const double * q = P;

    com_accumulate(cm, p, q);
    for(int kk = 0; kk+1<n; kk++)
    {
        p = P + 2*kk;
        q = P + 2*(kk+1);
        com_accumulate(cm, p, q);
    }
    double a = poly_area(P, n);
    cm[0] /= a;
    cm[1] /= a;
    return cm;
}

double poly_M20(const double * P, int n)
{
    // Raw moment (2,0)
    // Corresponding to integrating x^2 over the area
    // dP/dy = -x*x, P = - x*x*y/2

    if(n < 3)
    {
        return 0;
    }

    const double * p = P + (2*n-2);
    const double * q = P;

    double M20 = poly_accumulate_M20(p, q);
    for(int kk = 0; kk+1<n; kk++)
    {
        q = P + 2*(kk+1);
        p = P + 2*kk;
        M20 += poly_accumulate_M20(p, q);
    }

    return M20;
}

double poly_M11(const double * P, int n)
{
    // Raw moment (1,1)

    if(n < 3)
    {
        return 0;
    }

    const double * p = P + (2*n-2);
    const double * q = P;

    double M11 = poly_accumulate_M11(p, q);
    for(int kk = 0; kk+1<n; kk++)
    {
        q = P + 2*(kk+1);
        p = P + 2*kk;
        M11 += poly_accumulate_M11(p, q);
    }

    return M11;
}


double poly_M02(const double * P, int n)
{
    // Raw moment (0,2)

    if(n < 3)
    {
        return 0;
    }

    const double * p = P + (2*n-2);
    const double * q = P;

    double M02 = poly_accumulate_M02(p, q);
    for(int kk = 0; kk+1<n; kk++)
    {
        q = P + 2*(kk+1);
        p = P + 2*kk;
        M02 += poly_accumulate_M02(p, q);
    }

    return M02;
}



double poly_area_rourke(const double * P, int n)
{

    // Rourke, Eq. 1.14
    // What we get by Greens formula in poly_area_green is more efficient
    // [Rourke] Joseph O'Rourke, Computational Geometry in C, Second Edition

    if(n<3)
    {
        return 0;
    }

    double a = (P[2*(n-1)] + P[0])*(P[1]-P[2*(n-1)+1]);

    for(int kk = 0; kk+1 < n; kk++)
    {
        a += (P[2*kk] + P[2*kk+2])*(P[2*kk+3]-P[2*kk+1]);
    }
    a /= 2.0;
    assert( fabs(a - poly_area(P, n)) < 1e-6);
    printf("a = %f, ag = %f\n", a, poly_area(P,n));

    return a;
}

double poly_area(const double * P, int n)
{
    // From Greens formula
    if(n<3)
    {
        return 0;
    }

    const double * p = P + (2*n-2);
    const double * q = P;
    double a = p[0]*q[1] - q[0]*p[1];

    for(int kk = 0; kk+1<n; kk++)
    {
        q = P + 2*(kk+1);
        p = P + 2*kk;
        a += p[0]*q[1] - q[0]*p[1];
    }

    a /= 2.0;

    return a;
}

double * poly_cov(const double * P, int n)
{
    /* Covariance matrix.
     * Calculated from the raw moments.
     * u20 = M20/M00 - (M10/M00)^2
     * u11 = M11/M00 - (M10*M01)/(M00)^2
     * u02 = M02/M00 - (M01/M00)^2
     * Where M00 is the same as the area
     * M10 is com_x, M01 is com_y
     */

    double M00 = poly_area(P, n);

    double * C = poly_com(P, n);
    double M10 = C[0]*M00;
    double M01 = C[1]*M00;
    free(C);
    double M20 = poly_M20(P, n);
    double M02 = poly_M02(P, n);
    double M11 = poly_M11(P, n);

    if(0){
    printf(" -- Raw moments:\n");
    printf("M00=%f\n", M00);
    printf("M10=%f, M01=%f\n", M10, M01);
    printf("M20=%f, M11=%f, M02=%f\n", M20, M11, M02);
    }

    double u11 = M11 - M10*M01/M00;
    double u20 = M20 - M10/M00*M10;
    double u02 = M02 - M01/M00*M01;
    if(0){
    printf(" -- Centered moments:\n");
    printf("u20=%f, u11=%f, u02=%f\n", u20, u11, u02);
    }
    double * COV = malloc(3*sizeof(double));
    COV[0] = u20/M00;
    COV[1] = u11/M00;
    COV[2] = u02/M00;
    return COV;
}



static void eigenvector_sym_22(double a, double b, double c, double l, double * v0, double * v1)
{
    // Eigenvector to the matrix [a b; b c] corresponding
    // to eigenvalue l
    double vx = b / sqrt(b*b + pow(l-a,2));
    double vy = b / sqrt(b*b + pow(l-c,2));
    double n = sqrt(pow(vx, 2) + pow(vy, 2));
    v0[0] = vx/n;
    v1[0] = vy/n;
    return;
}

static void eigenvalue_sym_22(double a, double b, double c, double * l0, double * l1)
{
    // Eigenvalues to the matrix [a b; b c] corresponding
    // l0 > l1 since q > 0
    double p = (a+c)/2.0;
    double q = 0.5*sqrt( pow(a-c, 2) + 4*pow(b, 2));
    l0[0] = p + q;
    l1[0] = p - q;
}



double poly_orientation_with_COV(double * COV)
{
    double a = COV[0];
    double b = COV[1];
    double c = COV[2];


    //printf("COV = [[%f, %f]; [%f, %f]]\n", a, b, b, c);
    double l0, l1;
    eigenvalue_sym_22(a, b, c, &l0, &l1);
    //printf("l0 = %f, l1 = %f\n", l0, l1);

    //printf("MajorAxisLength: %f\n", 4*sqrt(l0));
    //printf("MinorAxisLength: %f\n", 4*sqrt(l1));

    double v0, v1;
    eigenvector_sym_22(a, b, c, l0, &v0, &v1);
    double w0, w1;
    eigenvector_sym_22(a, b, c, l1, &w0, &w1);

    //printf("v0 = [%f, %f]\n", v0, v1);
    //printf("v1 = [%f, %f]\n", w0, w1);

    double orientation = atan2(v1, v0);
    //printf("Orientation: %f\n", orientation);

    return orientation;
}

double poly_orientation(const double * P, int n)
{
    double * COV = poly_cov(P, n);
    double orientation = poly_orientation_with_COV(COV);
    free(COV);
    return orientation;

}

static double d2(const double * P0, const double *P1)
{
    return sqrt( pow(P0[0] - P1[0], 2) + pow(P0[1] - P1[1], 2));
}

double poly_circ(const double * P, int n)
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


double * poly_cbinterp(const double * P, int n, int upsampling, int * N)
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

double * poly_bbx(const double * P, int n)
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
        if(y > maxy)
        {
            maxy = y;
        }
        if(y < miny)
        {
            miny = y;
        }
    }
    bbx[0] = minx;
    bbx[1] = maxx;
    bbx[2] = miny;
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
    free(bbx);
    // For clean valgrind
    FcFini();
    cairo_debug_reset_static_data();
}
