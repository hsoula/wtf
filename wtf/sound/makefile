CC=gcc -fPIC
macCFLAGS= -fPIC
files_gen=spectrum.o 

mac_include=-I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 -I/Library/Python/2.7/site-packages/numpy/core/include -I/usr/local/include -I/usr/include/python2.7 -I/anaconda2/lib/python2.7/site-packages/numpy/core/include/ -I/usr/local/include -I/anaconda2/include/python2.7

mac_libs= -lgsl -lgslcblas -lsndfile -L/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/config -ldl  
mac_libs_dir=-L./ -L/usr/local/lib

all: spec.so

spec.c: spec.pyx
	cython spec.pyx

spec.o: spec.c
	$(CC) $(macCFLAGS) $(mac_include) -c -o spec.o spec.c 

spec.so: spec.o spectrum.o
	gcc -shared spec.o spectrum.o -o spec.so $(mac_libs)  $(mac_libs_dir)
spectrum.o: spectrum.c spectrum.h
	gcc -fPIC -c spectrum.c $(mac_include) 
