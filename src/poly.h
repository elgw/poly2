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

// Area for polygon with n vertices
double poly_area(double * P, int n);
// Circumference
// Returns 0 if less than 3 points
double poly_circ(double * P, int n);

void poly_print(FILE *, double * P, int n);

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
