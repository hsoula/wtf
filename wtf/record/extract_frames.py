"""
    : extract_frames.py
    Created: 11/03/16
    Description:

    Extract frames should be used/extended to extract parts of a wav file and return spectrum parts
    either randomly or systematically
    This part are defined by a time window and a frequency window.

"""

import numpy as np
import sys
sys.path.append('../')
from zf import zof_chunk
#import sndfile
import pysndfile

class frames(object):
    def __init__(self, start, end, array):
        self.t = array.shape[0]
        self.f = array.shape[1]
        self.a = array[:]
        self.start = start
        self.end = end

    def __repr__(self):
        s = "%d %d %f %f" % (self.t, self.f, self.start, self.end)
        return s + " ".join(["%f" % u for u in self.a.flatten()])

class frame_record(object):
    def __init__(self, fname, offset=0):
        """
        Main init method
        :param fname: wav file name path
        :param offset: starting time
        """
        # keep only the filename
        self.fname = fname.split('/')[-1]
        self.frames = []

        # load file
        self.sfile = pysndfile.PySndfile(fname, 'r')

        self.samplerate = self.sfile.samplerate()
        f_offset = np.floor(offset * self.samplerate)
        self.sfile.seek(f_offset, 0)
        # keep trakc of chunks loaded
        self.nchunk = 0

    def __enter__(self):
        """
        Main enter method.
        Called at the beginning of a with statement.
        """
        return self

    def get_frames(self, multiplier, wdw, freq_min, freq_max, f, nstd, fband, method=0):
        f_size = multiplier * wdw * self.samplerate  # in frames
        offset = self.nchunk * multiplier * wdw  # where we are in s

        self.curr_chunk = zof_chunk(self.sfile, f_size)

        if self.curr_chunk.frames == 0:
            return 0
        self.frames = []
        self.curr_chunk.compute_spectrum(freq_Hz=f, nstd=nstd, fband=fband, method=method)
        sample_rate = self.curr_chunk.samplerate
        fft = self.curr_chunk.s.fft
        size = self.curr_chunk.spectrum.shape
        if size[0] == 0:
            return 0
        imin = int(np.ceil(float(freq_min) / sample_rate * fft))
        imax = int(np.ceil(float(freq_max) / sample_rate * fft))
        spectrum = self.curr_chunk.spectrum[:, imin:imax]
        for ix in range(multiplier):
            start = int(ix * wdw * f)  # in f * s
            end = int(start + wdw * f)  # in f * s
            array = spectrum[start:end, :]
            if array.shape[0] == 0:
                return 0
            fr = frames(start / f + offset, end / f + offset, array)
            self.frames.append(fr)
        self.nchunk += 1
        return 1

    def save_to_text(self, fname):
        """
        Save the frames in form of matrix : lines is times and columns are frequencies
        for all the frames. WARNING: can yield enormous files :)
        WARNING 2: should be called only AFTER multiple get_frames calls
        :param fname:file path for saving
        :return:
        """
        with open(fname, 'w') as iOF:
            for fr in self.frames:
                for t in fr.a:
                    s = ' '.join(['%f' % u for u in t])
                    iOF.write(s + '\n')

    def __exit__(self, exc_type, exc_value, traceback):
        """
        Main exit method.
        Called after a with statement.
        """
        del self.sfile #.close()




if __name__ == '__main__':
    import matplotlib
    import matplotlib.pyplot as plt
    fname = sys.argv[1]
    f = 1000.0
    nstd = 6
    fband = 250
    freq_max = 20000
    freq_min = 1000
    done = False
    time = 0.025  # in ms
    multiplier = 10
    f_re = frame_record(fname)
    while not done:
        ret = f_re.get_frames(multiplier, time, freq_min, freq_max, f, nstd, fband)
        if ret == 0:
            done = True

        #fn = '%f-%f.png' % (x.start, x.end)
        #plt.imsave(fn, d.T)

    f_re.save_to_text('text.txt')
