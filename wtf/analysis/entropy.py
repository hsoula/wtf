import sys
sys.path.append('../')
import numpy as np
from scipy.signal import hilbert


def compute_entropy(spec, min_Hz, max_Hz):
    imin = np.ceil(min_Hz / spec.sample_rate * spec.fft)
    imax = np.ceil(max_Hz / spec.sample_rate * spec.fft)
    c = spec.spectrum[:, imin:imax]
    z = np.sum(c, 1)
    df = np.log2(max_Hz-min_Hz)
    v = np.array([u/z for u in c.T])
    ent = -np.sum(np.log2(v.T)*v.T, 1)/df
    return ent

def compute_hilbert(spec):
    xa = np.abs(hilbert(spec.get_envelope()))
    za = np.sum(xa)
    ## hilbert information
    # self.ana[k]=-sum(xa/za*log2(xa/za))/log2(f_fftw)