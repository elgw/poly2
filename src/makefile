cflags_ext=`pkg-config --cflags cairo gsl pango pangocairo fontconfig`
ldflags_ext=`pkg-config --libs gsl cairo pango pangocairo fontconfig`

cflags_dbg=-Wall -pedantic -std=c99 $(cflags_ext) -g
ldflags_dbg=-lm $(ldflags_ext)

cflags=-Wall -pedantic -std=c99 $(cflags_ext) -O3 -march=native -DNDEBUG -fno-math-errno
ldflags=-lm $(ldflags_ext) -flto

all: poly poly_ut poly_bench

poly_ut: poly_ut.c poly.c poly.h
	gcc $(cflags_dbg) poly_ut.c poly.c $(ldflags_dbg) -o poly_ut

poly_bench: poly_ut.c poly.c poly.h
	gcc $(cflags) poly_ut.c poly.c $(ldflags) -o poly_bench

tcga_test: tcga_test.c poly.c poly.h
	gcc $(cflags) tcga_test.c poly.c $(ldflags) -o tcga_test

poly: poly.c
	gcc -c $(cflags) poly.c

lpoly.so: poly.c poly_lua.c
	gcc -c $(cflags) poly.c     -fPIC -g -Wall -o fpoly.o
	gcc -c $(cflags) poly_lua.c -fPIC -g -Wall -o flpoly.o
	gcc --shared fpoly.o flpoly.o $(ldflags) -o lpoly.so
