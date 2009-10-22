#LDFLAGS=-L/usr/local/djbfft/lib -ldjbfft -lm -Wl,-rpath=/usr/local/djbfft/lib
LDFLAGS=-lm
LDLIBS=/usr/local/djbfft/lib/djbfft.a
CFLAGS=-fomit-frame-pointer -O2 -Wextra -Wall -malign-double -I/usr/local/djbfft/include

1: