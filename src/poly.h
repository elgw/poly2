/* Polygon processing/measurements

Polygons with n vertices are stored with x and y coordinates interlaced,
i.e. x0, y0, x1, x2, ... typically denoted P

TODO:
 - Malloc-free versions of all routines, using supplied memory buffers.
 - Principal directions.
 - Minimal bounding box (oriented according to the principal directions).
 - For convex hull, consider Graham and Yao: https://doi.org/10.1016/0196-6774(83)90013-5
 - ...
*/


#ifndef __poly_h__
#define __poly_h__

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_linalg.h>
#include <cairo.h>
#include <cairo-svg.h>
#include <fontconfig/fontconfig.h>

typedef struct{
    int nVertices;
    double Area;
    double * Centroid;
    double * BoundingBox;
    double MajorAxisLength;
    double MinorAxisLength;
    double Eccentricity;
    double Orientation; // Using atan2 on principal axes
    //double * ConvexHull;
    double ConvexArea;
    double Circularity;
    double EquivDiameter;
    double Solidity;
    double Perimeter;

    int measured;
} poly_props;

//
// MEASUREMENTS
//

// "High Level" interface, i.e. measure most stuff in one go
poly_props * poly_measure(const double * P, int n);
void poly_props_free(poly_props**);
void poly_props_print(FILE * fout, poly_props * props);

// Area for polygon
double poly_area(const double * P, int n);

// Circumference
// Returns 0 if less than 3 points
double poly_circ(const double * P, int n);

// bounding box of polygon, minx, maxx, miny, maxy
double * poly_bbx(const double * P, int n);

// Centre off mass
double * poly_com(const double * P, int n);

// Covariance matrix
// Returns comx, comy, c11, c12, c22
double * poly_cov(const double * P, int n);


// Extracts the covariance matrix and gets orientation
// from eigenvectors
double poly_orientation(const double * P, int n);

//
// MANIPULATION
//

// Multiply by v
void poly_mult(double * P, int n, double vx, double vy);

// Translation
void poly_translate(double * P, int n, double dx, double dy);
// Rotation around (x0, y0)
void poly_rotate(double * P, int n, double x0, double y0, double theta);

// Cubic spline interpolation
double * poly_cbinterp(const double * P, int n, int upsampling, int * N);

void poly_print(FILE *, const double * P, int n);

void poly_to_svg(double * P, int n, char *);

/* Cubic spline interpolation of closed curve
 * https://mathworld.wolfram.com/CubicSpline.html
 *
 * asserts that the input domain is regularly spaced
 */

// TODO
// Interpolate with max delta or to a certain number of points
// https://stats.stackexchange.com/questions/415974/how-to-find-the-covariance-matrix-of-a-polygon
// Centre of mass // -- Greens theorem
// Principal axes // -- Greens theorem
// Extent of axes (project vertices on principal axes)
// Smallest bounding box
// polygon rasterization using the even-odd rule or nonzero winding number
// Rourke Section 7.4
#endif
