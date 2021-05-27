cflags_dbg=-Wall `pkg-config --cflags cairo gsl fontconfig` -g
ldflags_dbg=-lm `pkg-config --libs gsl cairo fontconfig`

cflags=-Wall `pkg-config --cflags cairo gsl fontconfig` -O3 -DNDEBUG
ldflags=-lm `pkg-config --libs gsl cairo fontconfig` -flto

poly_ut:
	gcc $(cflags) poly_ut.c poly.c $(ldflags) -o poly_ut

poly_bench:
	gcc $(cflags) poly_ut.c poly.c $(ldflags) -o poly_bench
