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


Executable name is `isd`.


## Examples

- Syndrome decoding
```sh
$ wget http://decodingchallenge.org/Challenges/SD/SD_300_0
$ time ./isd 8 SD SD_300_0
n=300 k=150 w=38
l=16 p=4 epsilon=69 doom=0
001000010100000100000001000000000000000000000000000000001000000001010000000001000000000000000000000000000000000010000010000000000000000100000000000100100000000000011001000000000000000100000010000000000001100000010000000000000000000000001000101000000000010110001000010000010111000000000011010000000000

real	0m8.874s
user	1m8.215s
sys	0m0.156s
```
- Quasi-cyclic setting
```sh
$ wget http://decodingchallenge.org/Challenges/QC/Provider1/QC_28
$ time ./isd 8 QC QC_28
n=786 k=393 w=28
l=16 p=4 epsilon=40 doom=1
000000000000000000000000000000000000000100000000000000010010000010000000100010000000000000000000010100000100000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000000000000000001000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000010000000000010000000001000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000100000000000000000000000000000000000000010000010000000100000000000000000000000010000000000010000010000000000000000000000000000000000100000000000000000000000000000000000000100000000001010000000000000000000000000000000000000000000000000000000

real	0m3.661s
user	0m28.618s
sys	0m0.060s
```


## Parameters

Dumer parameters are chosen at compile time. They are:
- `DUMER_L` the width of the vectors used for collision in the birthday decoding
  part
- `DUMER_P` the weight of the vectors searched using birthday decoding
- `DUMER_EPS` the number of columns overlapping in the two sets
- `DUMER_DOOM` set to 1 to use the quasi-cyclicity of a code (all the circular
  shift of a syndrome will give the same error pattern up to blockwise circular
  shifts)
- `DUMER_LW` set to 1 to look for low-weight codeword instead of decoding

To set them the most convenient way is probably doing, for example:
```sh
$ EXTRA="-DDUMER_L=16 -DDUMER_P=4 -DDUMER_EPS=40 -DDUMER_DOOM=0 -DDUMER_LW=0" make -B
```

## File format

Files should respect the formats of <http://decodingchallenge.org/>.


# License

MIT (see file headers)
