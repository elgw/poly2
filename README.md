# poly2

A tiny library for feature extraction and manipulation of 2D polygons.

Only for simple polygons, i.e., with no crossing edges.

Moments are calculated using Greens theorem, and hence the orientation is important.
The orientation can be checked with `poly_vertex_order` and changed with `poly_reverse`
if needed.

If [Löve](https://love2d.org/) is installed, there is a small test program in `src/love/` that can be run with.
``` shell
love src/love/
```
Note, the library has to be build for Lua 5.1 in order for that to work.

## Dependencies:
For basic measurements, no special dependencies. For interpolation:
gsl, for visualizations: libcairo, and libpango.

``` shell
# sudo apt-cache search ...
sudo apt-get install libpangocairo-1.0-0
```

## Relevant references:
 * [https://en.wikipedia.org/wiki/Image_moment](https://en.wikipedia.org/wiki/Image_moment)
 * [https://stats.stackexchange.com/questions/415974/how-to-find-the-covariance-matrix-of-a-polygon](https://stats.stackexchange.com/questions/415974/how-to-find-the-covariance-matrix-of-a-polygon)
 * [https://maxgoldste.in/melkman/](Melkman's algorithm interactively)

## To do
See the separate [TODO](TODO.md).
