# Kramers-Kronig-Relation

This script takes an arbitrary spaced (but ordered) function in the form of
x and y values (given as numpy arrays) and returns the y values of the Hilbert transformation using the same x values.

Note1: There are faster algorithms for equally spaced x-value grids.
Note2: The function uses some numpy matrices to have a fast and compact notation. If you try to transform very large arrays you might run out of memory (O(n^2))


Needed Packages: numpy (matplotlib if run as script)

Usage as script:
```
./kkr.py data.dat
```
where data.dat is a two column file readable by numpy.loadtxt

Usage as import:
```
import kkr
...
y_real = kkr(x,y_imag)
```
in case the input function should drop off like 1/x you can use:
```
y_real = kkr(x,y_imag,tail=True)
```
