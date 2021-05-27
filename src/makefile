cflags=-Wall `pkg-config --cflags cairo gsl fontconfig` -g
ldflags=-lm `pkg-config --libs gsl cairo fontconfig`

poly_ut:
	gcc $(cflags) poly_ut.c poly.c $(ldflags) -o poly_ut
