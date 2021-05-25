cflags=-Wall `pkg-config --cflags cairo gsl`
ldflags=-lm `pkg-config --libs gsl cairo`

poly_ut:
	gcc $(cflags) poly_ut.c poly.c $(ldflags) -o poly_ut
