/* Polygon processing/measurements

Polygons with n vertices are stored with x and y coordinates interlaced,
i.e. x0, y0, x1, x2, ... typically denoted P

TODO:
 - Malloc-free versions of all routines, using supplied memory buffers.
 - 2nd order moments.
 - Principal directions.
 - Minimal bounding box (oriented according to the principal directions).

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

//
// MEASUREMENTS
//

// Area for polygon
double poly_area(const double * P, int n);

// Circumference
// Returns 0 if less than 3 points
double poly_circ(const double * P, int n);

// bounding box of polygon, minx, maxx, miny, maxy
double * poly_bbx(const double * P, int n);

// Centre off mass
double * poly_com(const double * P, int n);

// 2nd order Moments
// Returns comx, comy, xx, xy, yy
// I.e. the moment tensor components at the end.
double * poly_moments(const double * P, int n);

//
// MANIPULATION
//

// Translation
void poly_translate(double * P, int n, double dx, double dy);

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
