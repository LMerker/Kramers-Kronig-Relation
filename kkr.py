#!/usr/bin/env python2.7
import numpy as np


def kkr(x,y,tail=False):
    """Returns the real part for a given imaginary part of a function. Note that the reverse transformation is identical but for a  minus sign.
The given function is triangulated and then transformed using the correct limits when formally deviding by zero.
The returned values will be evaluated on the same x-values.
If the tail=True option is used, the function tries to make the outer most points into a slope that decreases with 1/x into infinity (usefull for backtransformations)"""
    # Define differences to calculate slope
    dx = x[1:] - x[:-1]
    dy = y[1:] - y[:-1]
    # Calculate slope
    a = dy/dx
    # Create a Matrix containing n rows containing  [x_0 .. x_n-2] and n rows containing [x_i,x_i...] . This way we can avoid an ugly for loop.
    xA,xy= np.meshgrid(x[:-1],x) #0.5*(x[:-1]+x[1:]))#
    # Create a Matrix containing n rows containing  [x_1 .. x_n-1] and n rows containing [x_i,x_i...] . This way we can avoid an ugly for loop.
    xB,xy= np.meshgrid(x[ 1:],x) #0.5*(x[:-1]+x[1:]))#
    # Integrate (ax + b)/(x-x_j) dx from x_i to x_i+1 = ((a*(x_j -x_i) + b)*log(| (x_i+1 - x_j)/(x_i-x_j)|)
    r = ((a*(xy-xA)) +y[:-1])*np.log(np.abs((xB-xy)/(xA-xy))) + a*(xB-xA) 
    # Special cases for x_j = x_i so int_x_i-1^x_i and int_x_i^x_i+1 make problems so you have to add them directly
    for i in xrange(len(dx)):
        
        r[  i,i] = a[i]*(xB[  i,i]-xA[  i,i])
        r[i+1,i] = a[i]*(xB[i+1,i]-xA[i+1,i])
        r[i+1,i] += y[i+1]*np.log(np.abs(dx[i+1]/dx[i])) if i+1 < len(dx) else 0.0
        
    r = np.sum(r,axis=1)/np.pi
    # there might be a A/(B-x) tail that would be neglegted. Triing to estimate it with the outer most points:
    if(tail):
        r += (x[0]-x[len(x)-1])*y[0]*y[len(y)-1]*np.log((x[0] - x)/(x[len(x)-1] -x) * y[0]/y[len(y)-1] ) /(x[0]*y[0]-x[len(x)-1]*y[len(y)-1] -x*y[0] + x*y[len(y)-1])/np.pi
    return r#0.5*(x[:-1]+x[1:]),
def kkr2(x,y):
    """ This function tries to estimate the parts of the given function as a a/(b-x) curve instead of triangulating the function.
Do not use."""
    a = (x[:-1] - x[1:])*y[:-1]*y[1:]/(y[:-1] - y[1:]) # Fitting a/(b-x) to curve
    b = (x[:-1]*y[:-1]-x[1:]*y[1:])/(y[:-1] - y[1:])
    xn = x.copy()
    xn[0] = -1e10
    xn[-1]= 1e10
    xA,xy= np.meshgrid(xn[:-1],x)#0.5*(x[:-1]+x[1:]))#
    xB,xy= np.meshgrid(xn[ 1:],x)#0.5*(x[:-1]+x[1:]))#
    r = a*np.log(np.abs( (xB-xy)*(xA-b)/((xB-b)*(xA-xy)) ))/(b-xy)
    for i in xrange(1,len(a)):
        r[i,i] = y[i]*np.log(np.abs( (x[i+1]-x[i])/(x[i]-x[i-1])*y[i+1]/y[i-1] ))
        r[i,i-1] = 0.0
    #print r,a,b,x,y
    r = np.sum(r,axis=1)/np.pi
    return r
if(__name__=='__main__'):
    import matplotlib.pyplot as plt
    import sys
    datas = []
    if (len(sys.argv) > 1):
        for filename in sys.argv[1:]:
            data = np.loadtxt(filename)
            x = data[:,0]
            y = data[:,1]
            datas.append(x)
            datas.append(kkr(x,y))
    else:
        x = np.random.rand(1000)*100 -50
        x.sort()
        y = np.zeros_like(x)
        y[np.abs(x)< 20] = 1
        datas = [x,y,'-o',x,kkr(x,y),'-x',x,-kkr(x, kkr(x,y),tail=True),'-o']
    plt.plot(*datas)
    plt.grid(True)
    plt.show()
