# Kramers-Kronig-Relation

This script takes an arbitrary spaced (but ordered) function in the form of
x and y values (given as a numpy array) and returns the y values of the Hilbert transformation using the same x values.



Needed Packages: numpy (matplotlib if run as script)

Usage as script:
./kkr.py data.dat

Usage as import
import kkr
...
y_real = kkr(x,y_imag)

in case the input function should drop off like 1/x you can use:

y_real = kkr(x,y_imag,tail=True)