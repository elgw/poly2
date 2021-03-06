# poly2

![Screenshot](screenshot.png)

A tiny library for feature extraction of 2D polygons with focus on performance.

 - Processes > 2 M small polygons per second on a single thread (Ryzen 3800X).
 - As you would expect it is only for simple polygons,
   i.e., with no crossing edges.
 - Moments are calculated using Greens theorem, and hence the orientation
   is important.
 - The orientation can be checked with `poly_vertex_order` and changed
   with `poly_reverse` if needed.

## Dependencies:
Only tested under 64-bit Ubuntu 21.04.
For basic measurements, no special dependencies. For interpolation:
gsl, for visualizations: libcairo, and libpango.

``` shell
# sudo apt-cache search ...
sudo apt-get install libpangocairo-1.0-0
```

## Lua/Löve demo
If [Löve](https://love2d.org/) is installed, there is a small test program in `src/love/` that can be run with.
``` shell
cd src
make lpoly.so
love love/demo
```

## To do or Maybe ...
 - [ ] Tests!
 - [ ] Clean up redundant calculations.
 - [ ] Switch to CCW as the standard orientation?
 - [ ] Reduce the number of calls to malloc/free.
 - [ ] Polygon rasterization using the even-odd rule or nonzero winding number, see Rourke Section 7.4
 - [ ] Minimal bounding box (oriented according to the principal directions).
 - [ ] Smarter svg plotting with axis abstraction? Or just drop that part?

## Relevant references:
 * [Image moments](https://en.wikipedia.org/wiki/Image_moment)
 * [How to find the covariance matrix of a polygon](https://stats.stackexchange.com/questions/415974/how-to-find-the-covariance-matrix-of-a-polygon)
 * [Interactive demo of Melkman's algorithm](https://github.com/mgold/Melkmans-Algorithm-Visualized)
