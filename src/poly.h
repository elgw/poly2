/* Polygon processing

Polygons are stored with x and y coordinates interlaced,
i.e. x0, y0, x1, x2, ...

References:
[Rourke] Joseph O'Rourke, Computational Geometry in C, Second Edition

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

// Area for polygon with n vertices
double poly_area(double * P, int n);
// Circumference
// Returns 0 if less than 3 points
double poly_circ(double * P, int n);

double * poly_bbx(double * P, int n);

// Cubic spline interpolation
double * poly_cbinterp(double * P, int n, int upsampling, int * N);

void poly_print(FILE *, double * P, int n);

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
