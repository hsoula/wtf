"""
    : make_mono.py
    Created: 20/04/16
    Description:

"""

import pysndfile as sndfile
import numpy as np
import sys

def make_mono(fname):
    sfile = sndfile.open(fname, 'r')
    samplerate = sfile.info.samplerate
    if sfile.info.channels == 1:
        print 'Already mono'
        exit(0)

    nchunk = 1000
    fname_r = fname[:-4] + '_R.wav'
    fname_l = fname[:-4] + '_L.wav'

    outinfo = sndfile.SF_INFO(samplerate=samplerate,
                                      channels=1,
                                      format=(sndfile.SF_FORMAT_WAV | sndfile.SF_FORMAT_PCM_16),
                                      sections=1,
                                      seekable=1)
    sfile_r = sndfile.open(fname_r, mode='w', info=outinfo)
    sfile_l = sndfile.open(fname_l, mode='w', info=outinfo)

    done = False
    while not done:
        re = sfile.readf_double(nchunk * 2)
        re_l = re[::2]
        re_r = re[1::2]
        sfile_l.writef_double(re_l)
        sfile_r.writef_double(re_r)
        if len(re) < nchunk * 2:
            done = True
    sfile.close()
    sfile_l.close()
    sfile_r.close()


if __name__ == '__main__':
    fname = sys.argv[1]
    make_mono(fname)
