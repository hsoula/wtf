"""
zof module is just a <<copy>> using sound instead of fft from numpy 
with a different interface 

"""

import numpy as np
import pysndfile
import sys

sys.path.append("../")
from sound import specs


class zof_chunk(object):
    """
    chunk class : small piece of sound file
        extracted for a large file 
    """

    def __init__(self, sfile, f_size, overlap=1024):
        """

        :param sfile: an open sndfile info on a valid file usually open by record class
        :param f_size: size in frames to read into
        :param overlap: overlap between chunk to have information before (using fft)
        :return:
        """

        #self.info = sfile.get_info()
        # information about the sound file sfile

        f_cur = sfile.seek(0, 1)
        # save the current position in the file
        # ie start of chunk
        # (2nd parameter : 1 means seek relative to the current position)

        f_rem = sfile.seek(0, 2) - f_cur
        # number of steps remaining.
        # 2nd parameter : 2 means seek relative to the file's end.
        if f_size + overlap >= f_rem:
            # if size of chunk + overlap > remaining size (ie end of file):
            f_size = f_rem - overlap # reduces f_size to remaining size
            sfile.seek(f_cur, 0)
            # return to the actual position : f_cur steps from the start of the file
            # then load data
        
            self.r = sfile.read_frames(f_size + overlap)            
        else:
            sfile.seek(f_cur, 0)  # return to the actual position
            self.r = sfile.read_frames(f_size + overlap)  # load data
            sfile.seek(-overlap, 1)  # get back to the overlap

        # print f_size
        self.f_size = f_size
        self.frames = f_size  # frame = ?
#        print(sfile.samplerate())
        self.samplerate = float(sfile.samplerate())
        self.channels = sfile.channels()
    def save_enveloppe(self, fname, offset=0):
        # save the envelope
        f = open(fname, 'w')
        for i, p in enumerate(self.env):
            f.write('%f %f\n' % (i / float(self.freq_Hz), p))
        f.close()

    def get_filtered_env(self, filter_min=0, filter_max=25000):
        ''' return the env minus the frequency below filter_min and filter_max
        cannot be called before compute_spectrum
        '''
        fft = self.fft

        min_hz = max(0, np.floor(filter_min / self.samplerate * fft))
        max_hz = min(self.fft / 2, np.ceil(filter_max / self.sampleRate * fft))
        spec = self.spectrum[:, min_hz:max_hz]
        return np.sum(spec.T, 0)

    def save_wave(self, fname):
        # file with first column between 0 and 1 : 0,1/44100,2/44100,...,44099/44100,1
        f = open(fname, 'w')
        for i, p in enumerate(self.r):
            f.write('%f %f\n' % (i / float(self.samplerate), p))
        f.close()

    def apply_stereo_method(self, method):
        ch = self.channels  # 1 = mono, 2 = stereo
        self.method = method
        # saving the method for later use
        if ch == 2:
            if method == 2:  # take the mean of r[0] and r[1], then r[2] and r[3] ...
                self.r = 0.5 * (self.r[::2] + self.r[1::2])
                # r[::2] : extract elt every two ids starting
                # by 0 (r[0],r[2],...), r[1::2] :
                # extract elt every two ids starting by 1 (r[1],r[3],...)
            else:  # take every two elt from data starting by method
                self.r = self.r[method::2]

    def compute_spectrum(self, freq_Hz=1000.0, nstd=6, fband=125, method=2):
        """
        compute the spectrum and the envelop using gaussian filtered fft
        freq_Hz: sampling frequency for the spectrum
        nstd :  number of standard deviation of the gaussian filter
        fband :  precision of frequency band
        method: if method==2 take the average of the two channels
        """

        # dealing with stereo
        self.apply_stereo_method(method)

        # houskeeping
        self.freq_Hz = freq_Hz

        # compute the spectrum
        self.s = specs(self.r,
                       self.samplerate,
                       f=freq_Hz,
                       nstd=nstd,
                       fband=fband)
        # get the power spectrum
        self.spectrum = self.s.spectrum
        # get the envelop
        self.env = np.sum(self.spectrum.T, 0)
        # we're done
        return self.env

    def extract_calls(self, threshold, filter_Hz):
        fc = self.s.fc
        spec = self.spectrum
        fft = self.s.fft
        sampleRate = float(self.samplerate)
        r_hz = np.floor(filter_Hz / sampleRate * fft)
        dB = threshold
        start = 0  # start of sound
        end = 1  # end of sound
        self.f_sounds = []

        for idx in range(fc):
            # the value from r_hz to fft/2
            z = sum(spec[idx, r_hz:fft / 2])
            # sum over y from r_hz to f_fft/2 (fftw/2 Shannon theorem)
            if z > dB and end != 0:
                start = idx + 1  # put a start for a sound
                end = 0  # now the end of the sound is set as unknown
            if z < dB and start != 0:
                end = idx + 1
                self.f_sounds.append([start, end])  # adding the new sound
                start = 0  # start a new sound
        if start != 0:  # a current sound is not finished
            self.f_sounds.append([start, -1])  # adding the sound but with -1 as end


class zof_record(object):
    def __init__(self, fname, f_size, timemax=-1):
        self.sfile = pysndfile.PySndfile(fname, 'r')
        if self.sfile.info.frames == 0:  # if number of frames in the wav file is equal to 0
            raise IOError, '%s not a correct wav file' % fname
        self.offset = 0  # current offset
        self.samplerate = self.sfile.info.samplerate  # sample rate = sound file sample rate
        self.frames = self.sfile.info.frames  # number of frames
        self.f_size = f_size
        self.last_time = -1  # flag for interrupted call
        self.nchunk = 0
        self.c_sounds = []

    def get_duration(self):
        duration = self.frames / float(self.samplerate) / 3600.
        return duration

    def _get_next_chunk(self):
        # create current chunk and increases the number of chunk if the chunk contains frames
        self.c_curr = zof_chunk(self.sfile, self.f_size)
        #                print self.c_curr.frames
        if self.c_curr.frames == 0:
            return 0
        self.nchunk += 1
        return 1

    def _get_calls(self, threshold, filter_Hz, freq_Hz=1000.0, nstd=6, fband=125, method=2):
        # fill self.c_sounds with the current chunk, dealing with overlaps
        self.c_env = self.c_curr.compute_spectrum(freq_Hz,
                                                  nstd,
                                                  fband,
                                                  method)
        self.c_curr.extract_calls(threshold, filter_Hz)

        if self.last_time != -1:
            ## call started last time : we merge the two
            if len(self.c_curr.f_sounds) > 0:
                # if the current chunk contains at least one sound (from its start to its end)
                self.c_curr.f_sounds[0][0] = self.last_time - self.offset
                # set start of the first sound of the chunk to last_time-offset
                # to add it in the c_sound list

                self.last_time = -1
                # waiting for new interrupted sound
            else:
                # empty calls the call is longer than the chunk
                print 'warning: call longer than chunk: consider increasing chunk size'
                # but nothing to do : wait for next chunk
        for sound in self.c_curr.f_sounds:
            f_start = sound[0]
            f_end = sound[1]
            if f_end != -1:
                # flag for interrupted sound (cf zf_chunk class)
                self.c_sounds.append([f_start + self.offset, f_end + self.offset])
            else:
                # we probably are at the end of the chunk:
                # self.last_time = position of the current sound
                # to finish next time (next chunk) and don't add
                # this sound to the c_sound list

                self.last_time = f_start + self.offset
                # end of the chunk of size f_size
        self.offset += self.f_size / float(self.samplerate) * freq_Hz
        # add offset to place correctly the sound in the new self.c_sounds sequence
        return self.c_sounds

    def next_chunk(self, threshold, filter_Hz, freq_Hz=1000.0, nstd=6, fband=125, method=2, verbose=0):
        # if there is another chunk, get information from it
        # (and print or save this information if verbose !=0)
        re = self._get_next_chunk()
        if re == 0:
            # we finished no more chunks
            return self.nchunk  # number of chunk

        c_sounds = self._get_calls(threshold, filter_Hz,
                                   freq_Hz=freq_Hz, nstd=nstd, fband=fband, method=method)

        # add current chunk information to c_sounds list
        if verbose > 1:
            # save spectrum file "test_spec_rec10.txt" if 10 chunks
            self.c_curr.save_spectrum('test_spec_rec%d.txt' % self.nchunk, self.offset)
        return -1

    def next_chunk_spectrum(self, freq_Hz=1000.0, nstd=6, fband=125, method=2):
        """
            This method just perform the next chunk and compute the spectrum
            IT DOES not look for the sound events!!
        :return:
        """
        re = self._get_next_chunk()
        if re == 0:
            return re

        self.c_env = self.c_curr.compute_spectrum(freq_Hz,
                                                  nstd,
                                                  fband,
                                                  method)
        return re

    def get_current_spectrum(self):
        return self.c_curr.spectrum

    def save_wave(self, fname):
        self.c_curr.save_wave(fname)


class zof_call(zof_chunk):
    def __init__(self, sfile, c_start, c_end, f_pad):
        self.sfile = sfile

        # we record sound start time
        self.c_so_start = c_start  # sound start in s
        self.c_so_end = c_end  # sound end in s

        # convert into frames
        f_start = int(np.floor(c_start * sfile.info.samplerate))
        f_end = int(np.floor(c_end * sfile.info.samplerate))

        self.f_so_start = f_start  # sound start in frames
        self.f_so_end = f_end  # sound end in frames

        if f_start - f_pad < 0:
            self.sfile.seek(0, 0)
            f_size = f_end + f_pad
        else:
            self.sfile.seek(f_start - f_pad, 0)
            f_size = f_end - f_start + 2 * f_pad

        zof_chunk.__init__(self, sfile, f_size, 0)
        # and done!
        self.f_pad = f_pad

    def extract_call(self, peak_pc):
        w = self.f_pad / float(self.samplerate) * self.freq_Hz
        m = max(self.env[w:-w])  # threshold depending on the peak height
        above = np.where(self.env > m * peak_pc)[0]
        self.nf_start = above[0]  # first position where data is higher than the threshold
        self.nf_end = above[-1]  # last position where data is higher than the threshold
        self.f_start = self.nf_start / float(self.freq_Hz) * self.samplerate + self.f_so_start
        self.f_end = self.nf_end / float(self.freq_Hz) * self.samplerate + self.f_so_start
        self.m = m

    def save_wave(self, filename):
#        outinfo = pysndfile.construct_format(samplerate=self.samplerate,
#                                  channels=1,
#                                  format=(pysndfile.SF_FORMAT_WAV | pysndfile.SF_FORMAT_PCM_16),
#                                  sections=1,
#                                  seekable=1)

        o_sfile = pysndfile.PySndfile(filename, mode="w", format='wav')
        self.sfile.seek(self.f_start, 0)
        r = self.sfile.readf_double(self.f_end - self.f_start)  # read the sound file from f_start to f_end
        r = r / max(abs(r))  # to be between 0 and 1 ?

        # deal with stereo
        ch = self.sfile.info.channels

        # we force into mono by recalling the method we used for the extraction
        if ch == 2:
            if self.method == 2:
                nr = 0.5 * (r[::2] + r[1::2])
            else:
                nr = r[self.method::2]
        else:
            nr = r
        o_sfile.write_frames(nr)
        #o_sfile.close()

    def concat(self, o_call):
        if self.f_end > o_call.f_start:
            self.f_end = o_call.f_end
        else:
            pass
