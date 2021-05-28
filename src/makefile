cflags_ext=`pkg-config --cflags cairo gsl pango pangocairo fontconfig`
ldflags_ext=`pkg-config --libs gsl cairo pango pangocairo fontconfig`

cflags_dbg=-Wall $(cflags_ext) -g
ldflags_dbg=-lm $(ldflags_ext)

cflags=-Wall $(cflags_ext) -O3 -DNDEBUG
ldflags=-lm $(ldflags_ext) -flto

poly_ut:
	gcc $(cflags) poly_ut.c poly.c $(ldflags) -o poly_ut

poly_bench:
	gcc $(cflags) poly_ut.c poly.c $(ldflags) -o poly_bench

tcga_test:
	gcc $(cflags) tcga_test.c poly.c $(ldflags) -o tcga_test
