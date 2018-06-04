import numpy as np
try:
    import pysndfile
except:
    print "pb with sndfile"
    pass
import sys
sys.path.append('../')
from sound import spec_wav

class cl_calls(object):
    def __init__(self,fname):
        self.sfile=pysndfile.PySndfile(fname,'r')
        if self.sfile.info.frames==0: # if number of frames in the wav file is equal to 0
            raise IOError,'%s not a correct wav file'%fname        
        self.samplerate = self.sfile.info.samplerate # sample rate = sound file sample rate
        self.frames = self.sfile.info.frames # number of frames
        self.wav = self.sfile.readf_double(self.frames)
        self.sfile.close()
    def spectrum(self,dB=0.1,f=1000.0,freq_min=200.0,freq_max=8000.0,nstd=6,fband=125):
        return spec_wav(self.wav,self.samplerate,dB=dB,f=f,freq_min=freq_min,freq_max=freq_max,nstd=nstd,fband=fband)
    

