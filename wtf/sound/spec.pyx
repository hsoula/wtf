from libc.stdlib cimport free
from libc.string cimport memcpy

import numpy as np
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np

cdef extern from "numpy/arrayobject.h":
	ctypedef int intp
	ctypedef extern class numpy.ndarray [object PyArrayObject]:
		cdef char *data
		cdef int nd
		cdef intp *dimensions
		cdef intp *strides
		cdef int flags

cdef extern from "spectrum.h":
	double * GaussianSpectrumS(double * input, int inputsize,int increment, int winLength,int * fc, int *fft)
	void InvertAndAddS(double * spec,int frameCount,int winLength,int increment,double * data)


def pyGaussianSpectrum(ndarray input,int increment,int winLength):
	## no use of complex type everything is real and interlaced
	cdef double *inp = <double *>input.data
	cdef int in_size = input.dimensions[0]
	cdef int fc
	cdef int fft
	cdef double * res=GaussianSpectrumS(inp, in_size,increment,winLength,&fc,&fft)
	cdef np.ndarray h = np.zeros(2*fft*fc, dtype=np.double)	
	memcpy(<double*> h.data,res,2*fft*fc*sizeof(double))
	free(res)
	return h,fc,fft

def pyInvertAndAdd(ndarray spec,int frameCount,int winLength,int increment,ndarray wav):
	## simple wrapper here assuming the complex data is interlaced 
	cdef double * sp=<double*>spec.data
	cdef double * data=<double *>wav.data
	InvertAndAddS(sp,frameCount,winLength,increment,data)

