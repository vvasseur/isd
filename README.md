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


## Setting parameters

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


## Choosing parameters

A python script (`optimize.py`) is provided to help choose the parameters of the algorithm.
It empirically estimates the average running time for fixed parameters.
The optimal parameters are chosen using a simple hill-climbing optimization of the running time as a function of the parameters.
Its arguments are the same as those for `isd`.
The starting values for `P`, `L`, `EPS` can be modified by appending, for example, `P=5 L=30`.
If used for low-weight codeword finding, a target weight `W` should be provided by appending it in the same fashion, for example, `W=218`.

- Example in the quasi-cyclic setting
```sh
$ wget http://decodingchallenge.org/Challenges/QC/Provider0/QC_30
$ python3 optimize.py 8 QC QC_30 P=5
    P     L   EPS Est. time
    5    12    10 0:06:02.380958
    5    12     9 0:06:06.140503
[...]
    5    19    10 0:01:52.285594
    5    19     9 0:01:53.148645
    5    19    11 0:01:59.560413
    5    20    10 0:01:58.913279
    5    20     9 0:01:55.195710
    5    20    11 0:01:53.693244
EXTRA="-DDUMER_P=5L -DDUMER_L=19L -DDUMER_EPS=10L"
```

- Example in the low-weight codeword finding setting
```sh
$ wget https://decodingchallenge.org/Challenges/LW/LW_1280_0
$ python3 optimize.py 8 LW LW_1280_0 W=230
    P     L   EPS Est. time
    4    12    10 2:43:51.936259
    4    12     9 2:39:50.257167
    4    12    11 2:37:09.806838
[...]
    4    12    12 2:37:03.774286
    4    11    12 2:43:15.992123
    4    13    12 2:48:15.027402
    4    12    13 2:37:16.885877
    4    11    13 2:47:01.252960
    4    13    13 2:51:33.162128
EXTRA="-DDUMER_P=4L -DDUMER_L=12L -DDUMER_EPS=12L"
```


## File format

Files should respect the formats of <http://decodingchallenge.org/>.


# License

MIT (see file headers)
