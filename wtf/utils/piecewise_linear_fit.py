import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

## fit piecewise curves (2 lines) and give the fitted parameters

def piecewise_linear(x, x0,y0, k1, k2):
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

def fit_piecewise_linear(f,datax,datay,plot=True):
    p , e = curve_fit(piecewise_linear, datax, datay)

    x0,y0,k1,k2=p[0],p[1],p[2],p[3]
    print 'slope1=',k1,'slope2=',k2,'y0=',y0, 'x0=',x0
    if plot==True:
        xd = np.linspace(0, 30, 100)
        plt.plot(datax, datay, "o")
        plt.plot(xd, piecewise_linear(xd, *p))
        plt.show()
    return k1,k2



if __name__=='__main__':
    ## test
    x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ,11, 12, 13, 14, 15], dtype=float)
    y = np.array([5, 7, 9, 11, 13, 15, 28.92, 42.81, 56.7, 70.59, 84.47, 98.36, 112.25, 126.14, 140.03])

    fit_piecewise_linear(piecewise_linear,x,y,plot=True)

