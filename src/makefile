cflags=-Wall
ldflags=-lm `pkg-config --libs gsl`

poly_ut:
	gcc $(cflags) poly_ut.c poly.c $(ldflags) -o poly_ut
