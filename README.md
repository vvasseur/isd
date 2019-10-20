# Information Set Decoding implementation

This is an implementation of Dumer's algorithm for Information Set Decoding.

It is written specifically to solve decoding challenges on
[this website](http://decodingchallenge.org/).


## Prerequisite

This will not run if your CPU does not have the AVX2 instruction set.
To check whether you have it or not you can run:
```sh
$ cat /proc/cpuinfo | grep avx2
```

## Compiling

To compile this program, first you should download and compile m4ri with OpenMP
support by doing:
```sh
$ git clone https://bitbucket.org/malb/m4ri.git
$ autoreconf --install
$ ./configure CFLAGS="-Ofast -march=native -flto" --enable-openmp
$ make
```

Executable name is `isd`.


## Parameters

Dumer parameters are chosen at compile time. They are:
- `DUMER_L` the width of the vectors used for collision in the birthday decoding
  part
- `DUMER_P` the weight of the vectors searched using birthday decoding
- `DUMER_EPS` the number of columns overlapping in the two sets
- `DUMER_BDAY` the number of attempts at choosing columns before changing the
  information set (should be 1 if DUMER_EPS is well chosen)
- `DUMER_DOOM` set to 1 to use the quasi-cyclicity of a code (all the circular
  shift of a syndrome will give the same error pattern up to blockwise circular
  shifts)
- `DUMER_LW` set to 1 to look for low-weight codeword instead of decoding

To set them the most convenient way is probably doing, for example:
```sh
$ EXTRA="-DDUMER_L=16 -DDUMER_P=4 -DDUMER_EPS=40 -DDUMER_BDAY=1 -DDUMER_DOOM=0 -DDUMER_LW=0" make -B
```

## File format

Files should respect the formats of <http://decodingchallenge.org/>.


# License

MIT (see file headers)
