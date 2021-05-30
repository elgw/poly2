#include "poly.h"

//// Forward declarations for non-exported functions

//Cubic spline interpolation of closed curve
static double * cbspline(const double * y, int N, int f, int * np);
static double * vec_interlace(const double * A, double * B, size_t N);
static double * vec_deinterlace(const double * V, int stride, size_t N);
//static void vec_show(const double * V, int n);
static double d2(const double * P0, const double *P1);
static void eigenvector_sym_22(double a, double b, double c, double l,
                               double * v0, double * v1);
static void eigenvalue_sym_22(double a, double b, double c,
                              double * l0, double * l1);
double poly_orientation_with_COV(double * COV);
static double _is_left(const double * A, const double * B, const double * C);
static double _poly_hull_rcl(const double * A, const double * B, const double * C);

poly_props * poly_props_new()
{
    poly_props * props = malloc(sizeof(poly_props));
    props->Centroid = NULL;
    props->MajorDirection = NULL;
    props->BoundingBox = NULL;
    props->measured = 0;
    props->ConvexArea = -1;
    props->Solidity = -1;
    props->Comment = NULL;
    return props;
}

void poly_props_free(poly_props ** PP)
{
    poly_props * P = PP[0];
    if(P->MajorDirection != NULL)
    {
        free(P->MajorDirection);
    }
    if(P->Centroid != NULL)
    {
        free(P->Centroid);
    }
    if(P->BoundingBox != NULL)
    {
        free(P->BoundingBox);
    }
    if(P->Comment != NULL)
    {
        free(P->Comment);
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
    fprintf(fout, "  Simple: TODO\n");
    fprintf(fout, "  Vertices: %d\n", props->nVertices);
    if(props->VertexOrder == POLY_VERTEX_ORDER_CLOCKWISE)
    {
        fprintf(fout, "  Vertex order: clockwise\n");
    } else {
        fprintf(fout, "  Vertex order: NOT CLOCKWISE\n");
    }

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
    fprintf(fout, "  ConvexArea: %f\n", props->ConvexArea);
    fprintf(fout, "  Circularity: %f\n", props->Circularity);
    fprintf(fout, "  EquivDiameter: %f\n", props->EquivDiameter);
    fprintf(fout, "  Solidity: %f\n", props->Solidity);
    fprintf(fout, "  Perimeter: %f\n", props->Perimeter);
    if(props->Comment != NULL)
    {
        fprintf(fout, "  Comment %s\n", props->Comment);
    }
    return;
}



int poly_vertex_order(const double * P, int n)
{
    // Uses the method described in [sunday]
    // Returns < 0 if clockwise and > 0 for counter-clockwise
    // returns 0 for degenerate cases
    // note that non-simple polygons might not return 0

    // 1/ Find the vertex with smallest y,
    //    if more than one, pick the one with smallest x
    int idx = 0;
    double minx = P[0];
    double miny = P[1];

    for (int kk=1; kk<n; kk++) {
        double y = P[2*kk + 1];
        if (y > miny)
        {
            continue;
        }
        double x = P[2*kk];
        if (y == miny) {
            if (x < minx)
                continue;
        }
        idx = kk;      // a new rightmost lowest vertex
        minx = x;
        miny = y;
    }

    double val = 0;
    // 2. Return the orientation for idx-1, idx, idx+1
    if (idx == 0)
    {
        val = _is_left(P + 2*n-2, P, P+2);
    }
    else
    {
        val = _is_left(P+2*idx-2, P+2*idx, P+2*idx+2);
    }
    if(val < 0)
    {
        return POLY_VERTEX_ORDER_CLOCKWISE;
    }
    if( val > 0)
    {
        return POLY_VERTEX_ORDER_ANICLOCKWISE;
    }
    return POLY_VERTEX_ORDER_ERROR;
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
    props->VertexOrder = poly_vertex_order(P, n);
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
    props->MajorDirection = malloc(2*sizeof(double));
    eigenvector_sym_22(COV[0], COV[1], COV[2], l0,
                       props->MajorDirection, props->MajorDirection+1);
    //printf("!!?? %f, %f\n", props->MajorDirection[0], props->MajorDirection[1]);
    free(COV);

    if(1){
    int nH = 0;
    double * H = poly_hull(P, n, &nH);
    props->ConvexArea = poly_area(H, nH);
    free(H);
    props->Solidity = props->Area / props->ConvexArea;
    }
    props->measured = 1;
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
    //printf("a = %f, ag = %f\n", a, poly_area(P,n));

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
    double q = pow(b, 2) + pow(a-l, 2);
    double r = pow(c-l, 2) + pow(b, 2);
    //printf("a = %f, b = %f, c = %f\n", a, b, c);
    //printf("l = %f, q = %f, r = %f\n", l, q, r);

    if(q > r)
    {
        double n = sqrt(q);
        v0[0] = -b/n;
        v1[0] = (a-l)/n;
        return;
    }
    if(r >= q && r > 0)
    {
        double n = sqrt(r);
        v0[0] = (c-l)/n;
        v1[0] = -b/n;
        return;
    }
// Fallback for eye(2)
    v0[0] = 1;
    v1[0] = 0;
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
    //printf("l0: %f, l1: %f\n", l0[0], l1[0]);
}



double poly_orientation_with_COV(double * COV)
{
    double a = COV[0];
    double b = COV[1];
    double c = COV[2];

    //printf("COV = [[%f, %f]; [%f, %f]]\n", a, b, b, c);

    double orientation = 1.0/2.0*M_PI/2.0;
    if((a-c) != 0)
    {
        orientation = 1.0/2.0*atan(2.0*b/(a-c));
    }
    return orientation;
}
/*

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
*/

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


static double _is_left(const double * A, const double * B, const double * C)
{
    return - _poly_hull_rcl(A, B, C);
}

static double _poly_hull_rcl(const double * A, const double * B, const double * C)
{
    // One possible implementation of the (a, b, c)-function from the paper
    // used to determine if C is to the right or left of the A->B vector
    // splitting the plane
    if(C == A || C == B)
    {
        // Wouldn't work with -march=native
        // without this.
        return 0;
    }

    double qx = B[0] - A[0];
    double qy = B[1] - A[1];
    double rx = C[0] - A[0];
    double ry = C[1] - A[1];
    double val = qx*ry - qy*rx;
    return val;
}

static void _poly_hull_push(int * Q, int * t, int v)
{
    Q[++t[0]]=v;
}

static void _poly_hull_insert(int * Q, int * b, int v)
{
    Q[--b[0]]=v;
}

double * poly_hull(const double * P, int n, int * h)
{
    // Algorithm Hull from
    // On-Line Construction of the Convex Hull of a Simple Polyline
    // A. Melkman, Inf. Process. Lett., 1987 (25), pp. 11-12

    if(n<4)
    {
        // yes, ignoring also triangles.
        goto leave;
    }

    // A dequeue implemented as a redundant array.
    size_t nQ = 2*n;
    int * Q = malloc(nQ*sizeof(int));
    int b = n; // Index of bottom element
    int t = b-1; // Index of top element

/// 1. -- initialization
    if( _poly_hull_rcl(P, P+2, P+4) > 0)
    {
        _poly_hull_push(Q, &t, 0);
        _poly_hull_push(Q, &t, 1);
    } else {
        _poly_hull_push(Q, &t, 1);
        _poly_hull_push(Q, &t, 0);
    }
    _poly_hull_push(Q, &t, 2);
    _poly_hull_insert(Q, &b, 2);
    int idx = 3;


    // TODO: Consider adding a counter to abort when too many
    // iterations has passed.
label2:
    //printf("n = %d, b = %d, t = %d, idx = %d\n", n, b, t, idx);
    /// 2.

    while(! ( (_poly_hull_rcl(P+idx*2, P+Q[b]*2, P+Q[b+1]*2) < 0) ||
              (_poly_hull_rcl(P+Q[t-1]*2, P+Q[t]*2, P+idx*2) < 0)  ) )
    {
        idx++;
        // "The algorithm halts when the input is exhausted"
        if(idx >= n)
        {
            goto done;
        }
    }

    /// 3.
    while(!( _poly_hull_rcl(P+Q[t-1]*2, P+Q[t]*2, P+idx*2) > 0))
    {
        t--; // pop
        if(t - b < 2)
        {
            // Needed for non-simple or wrongly oriented polygons
            goto fail;
        }
        assert(t>b);
    }
    _poly_hull_push(Q, &t, idx);

    /// 4.
    while(!(_poly_hull_rcl(P+idx*2, P+Q[b]*2, P+Q[b+1]*2) > 0))
    {
        b++; // remove
        if(t-b < 2)
        {
            // Possibly needed for non-simple or wrongly oriented polygons
            goto fail;
        }
        assert(t>b);
    }

    _poly_hull_insert(Q, &b, idx);
    goto label2;

done: ;
    // Q[t] and Q[b] "will always refer to the same vertex"
    // hence we copy only until Q[b-1]
    assert(Q[t] == Q[b]);
    int nH = t-b;
    double * H = malloc(nH*2*sizeof(double));
    int row = 0; // output row
    for(int kk = b; kk<t; kk++)
    {
        H[row*2]     = P[Q[kk]*2];
        H[row*2 + 1] = P[Q[kk]*2 + 1];
        row++;
    }
    h[0] = nH;
    free(Q);
    return H;

fail: ;
    free(Q);

leave: ;
    h[0] = 0;
    return NULL;
}

static void draw_poly(cairo_t * cr, double * X, double * Y, int n,
                     double lineWidth, double red, double green, double blue)
{
    double alpha = 0.8;
    cairo_set_source_rgba(cr, red, green, blue, alpha);
    cairo_set_line_width(cr, lineWidth);
    cairo_move_to(cr, X[n-1], Y[n-1]);

    for(int kk = 0; kk < n; kk++)
    {
        cairo_line_to (cr, X[kk], Y[kk]);
    }
    cairo_close_path(cr);
    cairo_stroke(cr);

    for(int kk = 0; kk<n; kk++)
    {
        double pointSize = lineWidth*1.2;
        if( kk == 0)
        {
            pointSize *= 2.0;
        }
        if(kk == 1)
        {
            pointSize *= 1.5;
        }
        cairo_arc (cr, X[kk], Y[kk], pointSize, 0, 2*M_PI);
        cairo_fill (cr);
    }
    return;
}

void poly_to_svg(double * P, int n, char * filename)
{

    //// Settings
    // Width and height does not matter much since
    // it will be vector graphics.
    int w = 1024;
    int h = 512;
    double padding = 10;
    double * bbx = poly_bbx(P, n);

    double wx = bbx[1]-bbx[0];
    double wy = bbx[3]-bbx[2];
    if(wx > wy)
    {
        bbx[2] -= (wx-wy)/2;
        bbx[3] += (wx-wy)/2;
    } else {
        bbx[0] -= (wy-wx)/2;
        bbx[1] += (wy-wx)/2;
    }

    //// Initialization
    cairo_surface_t *surface;
    cairo_t *cr;
    surface = cairo_svg_surface_create(filename, w, h);
    cr = cairo_create(surface);

    //// Background
    cairo_rectangle (cr, 0, 0, w/2, h);
    cairo_set_source_rgba (cr, 1, 1, 1, 1);
    cairo_fill(cr);
    cairo_rectangle (cr, w/2, 0, w/2, h);
    cairo_set_source_rgba (cr, 0.7, 0.7, 0.7, 1);
    cairo_fill(cr);

    //// Transform polygon points
    double * X = malloc(n*sizeof(double));
    double * Y = malloc(n*sizeof(double));
    for(int kk = 0; kk<n; kk++)
    {
        X[kk] = P[2*kk];
        Y[kk] = P[2*kk+1];
        X[kk] = w/2 + coordscale(X[kk], bbx, (double) w/2, padding);
        Y[kk] = 0   + coordscale(Y[kk], bbx+2, (double) h, padding);
    }
    draw_poly(cr, X, Y, n, 4.0, 0.0, 0.0, 0.0);
    free(X);
    free(Y);

    //// Draw wire-frame for convex hull
    int nH = 0;
    double * H = poly_hull(P, n, &nH);
    if(0){
    printf("P:\n");
    poly_print(stdout, P, n);
    printf("H:\n");
    poly_print(stdout, H, nH);
    }

    if(H != NULL)
    {
    X = malloc(nH*sizeof(double));
    Y = malloc(nH*sizeof(double));
    for(int kk = 0; kk<nH; kk++)
    {
        X[kk] = H[2*kk];
        Y[kk] = H[2*kk+1];
        X[kk] = w/2 + coordscale(X[kk], bbx, (double) w/2, padding);
        Y[kk] = 0   + coordscale(Y[kk], bbx+2, (double) h, padding);
    }
    free(H);

    draw_poly(cr, X, Y, nH, 2.0, 0.0, 1.0, 0.0);

    free(X);
    free(Y);
    } else {
        printf("Could not calculate the convex hull\n");
    }

    //// Perform the feature extraction
    poly_props * props = poly_measure(P, n);
    //printf("?? %f, %f\n", props->MajorDirection[0], props->MajorDirection[1]);
    char *bp;
    size_t size;
    FILE *stream =  open_memstream (&bp, &size);
    poly_props_print(stream, props);
    fflush (stream);
    fclose(stream);


    //printf("?? %f, %f\n", props->MajorDirection[0], props->MajorDirection[1]);
    fflush(stdout);

    cairo_set_source_rgba(cr, .5, 0, .5, 1);
    cairo_set_line_width(cr, 2);
    cairo_move_to(cr, 0.75*w, 0.5*h);
    cairo_rel_line_to(cr, 200.0*props->MajorDirection[0], 200.0*props->MajorDirection[1]);
    cairo_stroke(cr);
    cairo_move_to(cr, 0.75*w, 0.5*h);
    cairo_rel_line_to(cr, -200.0*props->MajorDirection[0], -200.0*props->MajorDirection[1]);
    cairo_stroke(cr);
    poly_props_free(&props);

    //// Put some text
    PangoFontDescription* font_description=pango_font_description_new();
    pango_font_description_set_family(font_description, "Mono");
    pango_font_description_set_weight(font_description, PANGO_WEIGHT_BOLD);
    pango_font_description_set_absolute_size(font_description, 12*PANGO_SCALE);

    PangoLayout* layout=pango_cairo_create_layout(cr);
    pango_layout_set_font_description(layout, font_description);
    pango_layout_set_text(layout, bp,-1);
    cairo_set_source_rgb (cr, 0.0, 0.0, 0.0);
    cairo_move_to (cr, 10.0, 20.0);
    pango_cairo_show_layout(cr, layout);

    free(bp); // Free the feature text

    // cairo_surface_write_to_svg(surface, filename);

    //// Clean up
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
    free(bbx);
    // For clean(er) valgrind
    //FcFini();
    //cairo_debug_reset_static_data();
}
