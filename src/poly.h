/*
*
* Polygon measurements and some processing.
*
* Polygons with n vertices are stored with x and y coordinates interlaced,
* i.e. x0, y0, x1, x2, ... typically denoted P.
*
* Please see README.md for a gentle introduction.
*/

#ifndef __poly_h__
#define __poly_h__

#define _GNU_SOURCE

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* for poly_cbinterp */
#include <gsl/gsl_linalg.h>

/* for poly_to_svg */
#include <cairo.h>
#include <cairo-svg.h>
#include <fontconfig/fontconfig.h>
#include <pango/pangocairo.h>

#define POLY_VERTEX_ORDER_CLOCKWISE -1
#define POLY_VERTEX_ORDER_ERROR 0
#define POLY_VERTEX_ORDER_ANICLOCKWISE 1

typedef struct{
    int nVertices;
    int VertexOrder; // clockwise, anticlockwise or error
    int simple;
    double Area;
    double Centroid[2];
    double BoundingBox[4];
    double MajorAxisLength;
    double MinorAxisLength;
    double MajorDirection[2];
    double Eccentricity;
    double Orientation; // Using atan2 on principal axes
    //double * ConvexHull;
    double ConvexArea;
    double Circularity;
    double EquivDiameter;
    double Solidity;
    double Perimeter;
    double COV[3];
    char * Comment;
    int measured;
} poly_props;

/*
 * MEASUREMENTS
 */

/* The "High Level" interface. Measure most stuff in one go
   P is the vertices, (x0, y0), (x1, y1), ...
   n is the number of vertices
 */
poly_props * poly_measure(const double * P, int n);

void poly_props_free(poly_props**);
void poly_props_print(FILE * fout, poly_props * props);

/* Interface for specific measurements below */

/* Vertex order: we only support clockwise
 * The measurement can only be trusted for simple polygons.
 */
int poly_vertex_order(const double * P, int n);

/* Area for polygon */
double poly_area(const double * P, int n);

/* Circumference
 * Returns 0 if less than 3 points
 */
double poly_circ(const double * P, int n);

/* Bounding box of polygon, minx, maxx, miny, maxy */
double * poly_bbx(const double * P, int n);
void poly_bbx_buff(const double * P, int n, double * buff);

/* Center of mass */
double * poly_com(const double * P, int n);

/* Extracts the covariance matrix and gets orientation
 * from eigenvectors
 */
double poly_orientation(const double * P, int n);

/* Returns the convex hull of P with h points
 * Using Melkmans O(n) algorithm.
 * Returns NULL/sets h[0] to 0 when
 * less than 4 points or the algorithm fails.
 * It shouldn't fail for positively oriented simple polygons
 */
double * poly_hull(const double * P, int n, int * h);

/* Raw moments up to order 2,
 * returns M00, M10, M01, M20, M11, M02
 */
double * poly_moments_raw(const double * P, int n);

/* Individual raw moments
*/
double poly_M00(const double * P, int n);
double poly_M01(const double * P, int n);
double poly_M10(const double * P, int n);
double poly_M20(const double * P, int n);
double poly_M11(const double * P, int n);
double poly_M02(const double * P, int n);

/* covariance matrix using pre-calculated raw moments */
double * poly_cov_with_moments(const double * M);
void poly_cov_with_moments_buff(const double * M, double * buff);

/* Covariance matrix */
double * poly_cov(const double * P, int n);

/*
 * MANIPULATION
 */

/* Reverse the order of vertices in P */
void poly_reverse(double * P, int n);

/* Multiply x coordinates by vx, ... */
void poly_mult(double * P, int n, double vx, double vy);

/* Translation */
void poly_translate(double * P, int n, double dx, double dy);

/* Rotation around (x0, y0) */
void poly_rotate(double * P, int n, double x0, double y0, double theta);

/* Cubic spline interpolation
 * requires gsl
 */
double * poly_cbinterp(const double * P, int n, int upsampling, int * N);

/*
 * Check if a polygon is simple
 *
 * Returns 1 if the polygon is free from intersections.
 *  Note that this classification will always be affected by the
 *  machine precision and the results are practically undefined
 *  when edges/vertices are very close.
 *
 * This implementation uses an O(n^2) method
 * For better performance, check out:
 * - Bentley???Ottmann (to find all intersections)
 * https://en.wikipedia.org/wiki/Bentley%E2%80%93Ottmann_algorithm
*/
int poly_is_simple(const double * P, int n);
int lines_intersect(const double * a0, const double * a1,
                    const double * b0, const double * b1);
/*
 * Input / Output
 */

/* Print polygon in a MATLAB compatible way */
void poly_print(FILE *, const double * P, int n);

/* Plot the polygon and write some of the Properties
 * requires libcariro
 */
void poly_to_svg(double * P, int n, char *);

#endif
