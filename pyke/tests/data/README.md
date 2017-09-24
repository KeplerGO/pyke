This directory contains the following test data:

* `test-tpf-all-zeros.fits`: 3x3 Target Pixel File (TPF) filled with zeros everywhere.
* `test-tpf-non-zero-center.fits`: 3x3 TPF with ones in the center pixel and zeros everywhere else.
* `test-tpf-star.fits`: 3x3 TPF as follows:
```
0 1 0
1 1 1
0 1 0
```
* `test-tpf-with-nans.fits`: 3x3 TPF as follows:
```
0   1   0
1  nan  1
0   1   0
```
* `center-mask.txt`: mask file selecting pixel [1, 1].
