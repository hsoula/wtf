from numpy import *
import pysndfile

class zf_chunk(object):
    # chunk class : large piece of sound file
    def __init__(self, sfile, f_size, overlap=1024):  # overlap between chunk to have information before (fft calc)
        self.info = sfile.get_info()  # information about the sound file sfile
        # print "info=",self.info
        f_cur = sfile.seek(0,
                           1)  # save the current position in the file ie start of chunk (2nd parameter : 1 means seek relative to the current position)
        f_rem = sfile.seek(0,
                           2) - f_cur  # number of steps remaining. 2nd parameter : 2 means seek relative to the file's end.

        if f_size + overlap > f_rem:  # if size of chunk + overlap > remaining size (ie end of file):
            f_size = f_rem  # reduces f_size to remaining size
            sfile.seek(f_cur, 0)  # return to the actual position : f_cur steps from the start of the file
            self.r = sfile.readf_float(f_size + overlap)  ## load data

        else:
            sfile.seek(f_cur, 0)  # return to the actual position
            self.r = sfile.readf_float(f_size + overlap)  ## load data
            sfile.seek(-overlap, 1)  # get back to the overlap

        self.f_size = f_size
        self.frames = f_size  # frame = ?

    def save_spectrum(self, fname, offset=0):
        # save the envelope
        f = open(fname, 'w')
        for i, p in enumerate(self.env):
            f.write('%f %f\n' % ((i * self.f_div + self.f_fftw / 2 + offset) / float(self.info.samplerate), p))
        f.close()

    def save_wave(self, fname):
        # file with first column between 0 and 1 : 0,1/44100,2/44100,...,44099/44100,1
        f = open(fname, 'w')
        for i, p in enumerate(self.r):
            f.write('%f %f\n' % (i / float(self.info.samplerate), p))
        f.close()

    def get_env(self, f_fftw, f_div, Hz_LowPass, dB, method=2):
        ## compute the spectrum and perform the sound extraction (start and end of sounds with db higher than dB (threshold)
        ## filter > Hz_LowPass in Hz
        ## f_div the scaling divisor (in frames)
        ## f_fftw the fft window (in frames)
        ## dB power threshold 
        ## method: in case of stereo 

        ## dealing with stereo 
        ch = self.info.channels  # 1 = mono, 2 = stereo
        # saving the method for later use 
        self.method = method
        if ch == 2:
            if method == 2:  # take the mean of r[0] and r[1], then r[2] and r[3] ...
                self.r = 0.5 * (self.r[::2] + self.r[
                                              1::2])  # r[::2] : extract elt every two ids starting by 0 (r[0],r[2],...), r[1::2] : extract elt every two ids starting by 1 (r[1],r[3],...)

            else:  # take every two elt from data starting by method
                self.r = self.r[method::2]

        ## init variable 
        start = 0  # start of sound
        end = 1  # end of sound
        samplerate = float(self.info.samplerate)  # sample rate of the sound file
        self.f_sounds = []
        r_hz = int(min(f_fftw * Hz_LowPass / samplerate, f_fftw / 2 - 1))  # filter
        f_n = self.f_size / f_div  # f_n = envelope length = nb of data / division
        if f_n % f_div != 0:  # if f_n is not divisible by f_div
            f_n = f_n + 1

        self.env = zeros(f_n, dtype=np.float)
        for k in range(f_n):
            a = self.r[(k * f_div):(k * f_div) + f_fftw].copy()  # window for the fft comput.
            a = a - mean(a)  # signal over the window a  => signal - mean
            y = abs(fft.fft(a, f_fftw)) ** 2  # norm(fft(signal))^2  reel pos
            if a.shape != (f_fftw,):
                a.resize((f_fftw,))
                # print f_fftw,r_hz,f_fftw/2
            z = sum(y[r_hz:f_fftw / 2])  # sum over y from r_hz to f_fft/2 (fftw/2 Shannon theorem)
            if z > dB and end != 0:
                start = k + 1  # put a start for a sound
                end = 0  # now the end of the sound is set as unknown
            if z < dB and start != 0:  # if we have a start and z is lower than the dB threshold, we have an end for the current sound
                end = k + 1
                self.f_sounds.append([start * f_div, end * f_div])  # adding the new sound
                start = 0  # start a new sound
            self.env[k] = z
        if start != 0:  # a current sound is not finished
            self.f_sounds.append([start * f_div, -1])  # adding the sound but with -1 as end
            # we save everything (f_ftw/2)
            #	print len(self.env[self.env>dB]),dB,mean(self.env),max(self.env)
        self.f_fftw = f_fftw
        self.f_div = f_div
        return self.env


class zf_record(object):
    def __init__(self, fname, f_size, timemax=-1):
        self.sfile = pysndfile.PySndfile(fname, 'r')
        if self.sfile.info.frames == 0:  # if number of frames in the wav file is equal to 0
            raise IOError, '%s not a correct wav file' % fname
        self.offset = 0  ## current offset
        self.samplerate = self.sfile.info.samplerate  # sample rate = sound file sample rate
        self.frames = self.sfile.info.frames  # number of frames
        self.f_size = f_size
        self.last_time = -1  ## flag for interrupted call
        self.nchunk = 0
        self.c_sounds = []

    def _get_next_chunk(self):
        # create current chunk and increases the number of chunk if the chunk contains frames
        self.c_curr = zf_chunk(self.sfile, self.f_size)
        if self.c_curr.frames == 0:
            return 0
        self.nchunk += 1
        return 1

    def _get_calls(self, f_fftw, f_div, Hz_LowPass, dB, method=2):
        # fill self.c_sounds with the current chunk, dealing with overlaps
        self.c_env = self.c_curr.get_env(f_fftw, f_div, Hz_LowPass, dB, method)
        if self.last_time != -1:  ## call started last time : we merge the two
            if len(
                    self.c_curr.f_sounds) > 0:  # if the current chunk contains at least one sound (from its start to its end)
                self.c_curr.f_sounds[0][
                    0] = self.last_time - self.offset  # set start of the first sound of the chunk to last_time-offset to add it in the c_sound list
                self.last_time = -1  # waiting for new interrupted sound
            else:  ## empty calls the call is longer than the chunk
                print 'warning: call longer than chunk: consider increasing chunk size'
                ## but nothing to do : wait for next chunk
        for sound in self.c_curr.f_sounds:
            f_start = sound[0]
            f_end = sound[1]
            if f_end != -1:  # flag for interrupted sound (cf zf_chunk class)
                self.c_sounds.append([f_start + self.offset + f_fftw / 2, f_end + self.offset + f_fftw / 2])
            else:  ## we probably are at the end of the chunk: self.last_time = position of the current sound to finish next time (next chunk) and don't add this sound to the c_sound list
                self.last_time = f_start + self.offset
            # end of the chunk of size f_size
        self.offset += self.f_size  # add offset to place correctly the sound in the new self.c_sounds sequence
        return self.c_sounds

    def next_chunk(self, f_fftw, f_div, Hz_LowPass, dB, verbose=0):
        # if there is another chunk, get information from it  (and print or save this information if verbose !=0)
        re = self._get_next_chunk()
        if re == 0:
            ## we finished no more chunks
            return self.nchunk  # number of chunk
        c_sounds = self._get_calls(f_fftw, f_div, Hz_LowPass, dB)  # add current chunk information to c_sounds list
        if verbose > 1:  # save spectrum file "test_spec_rec10.txt" if 10 chunks
            self.c_curr.save_spectrum('test_spec_rec%d.txt' % self.nchunk, self.offset)
        if verbose > 0:  # print every sound (whole sequence over [0,1] : start/44100, end/44100, length of sound
            for c in c_sounds:
                print "sound : start=", c[0] / float(self.samplerate), "end=", c[1] / float(self.samplerate), "size=", (c[1] -c[0]) / float(self.samplerate)
        return -1

    def save_wave(self, fname):
        self.c_curr.save_wave(fname)


class zf_call(object):
    def __init__(self, sfile, f_start, f_end, f_pad):
        self.sfile = sfile
        ## we record sound start time
        self.f_so_start = f_start  # sound start
        self.f_so_end = f_end  # sound end
        self.ch = sfile.info.channels  # 1 = mono, 2 = stereo
        self.offset = f_pad  # offset to obtain information before sound start
        self.r_start = f_start - f_pad  # get information before start
        # we load the data
        if f_start - f_pad < 0:  ## start too early
            self.sfile.seek(0, 0)  # beginning of the file
            # print type(f_pad),type(f_start),type(self.ch)
            vpad = zeros(self.ch * (f_pad - int(f_start)),
                         dtype=np.float)  # O to fill the sequence if there is no data before (beginning of the file)
            self.r = concatenate((vpad, self.sfile.readf_float(f_end + f_pad)))  # append ? # [0,0,0,data ...]
            self.nframes = f_end - f_start + 2 * f_pad
        else:
            self.sfile.seek(f_start - f_pad, 0)  # position in f_start-f_pad
            self.nframes = f_end - f_start + 2 * f_pad  # number of frames for the whole time sequence considered
            self.r = sfile.readf_float(
                self.nframes)  # read data from f_start-f_pad to the end of the sequence (f_end+f_pad)
            # and done!
            # we can now save several stuff
        self.samplerate = self.sfile.info.samplerate

    def get_env(self, f_fftw, f_div, Hz_LowPass, method=2):
        ## compute the high pass filtered spectrum
        ## filter > Hz_LowPass in Hz
        ## f_div the scaling divisor (in frames)
        ## f_fftw the fft window (in frames)	
        ## method: in case of stereo 

        ## dealing with stereo 
        ch = self.sfile.info.channels
        ## we keep track of the method
        self.method = method
        if ch == 2:
            if method == 2:
                self.r = 0.5 * (self.r[::2] + self.r[1::2])
            else:
                self.r = self.r[method::2]
                ## init variable
        samplerate = float(self.samplerate)
        # print "samplerate",samplerate
        r_hz = int(min(f_fftw * Hz_LowPass / samplerate, f_fftw / 2 - 1))
        f_n = self.nframes / f_div
        # print f_n,self.nframes
        if f_n % f_div != 0:
            f_n = f_n + 1

        self.env = zeros(int(f_n), dtype=np.float)
        for k in range(int(f_n)):
            a = self.r[(k * f_div):(k * f_div) + f_fftw].copy()
            a = a - mean(a)
            y = abs(fft.fft(a, f_fftw)) ** 2
            if a.shape != (f_fftw,):
                a.resize((f_fftw,))
            z = sum(y[r_hz:f_fftw / 2])  # filter
            self.env[k] = z
        self.f_fftw = f_fftw
        self.f_div = f_div
        return self.env

    def save_spectrum(self, fname):
        f = open(fname, 'w')
        offset = self.r_start
        for i, p in enumerate(self.env):
            f.write('%f %f\n' % ((i * self.f_div + self.f_fftw / 2 + offset) / float(self.samplerate), p))
        f.close()

    def extract_call(self, peak_pc):

        w = self.offset / self.f_div
        # print "w=",w,self.offset
        m = max(self.env[w:-w])  # threshold depending on the peak height
        above = where(self.env > m * peak_pc)
        self.nf_start = above[0][0]  # first position where data is higher than the threshold
        self.nf_end = above[0][-1]  # last position where data is higher than the threshold
        self.f_start = self.nf_start * self.f_div + self.f_fftw / 2 + self.f_so_start - self.offset  ## ? pourquoi -offset
        self.f_end = self.nf_end * self.f_div + self.f_fftw / 2 + self.f_so_start - self.offset
        self.m = m

    def save_wave(self, filename):
        outinfo = sndfile.SF_INFO(samplerate=self.samplerate, channels=1,
                                  format=(sndfile.SF_FORMAT_WAV | sndfile.SF_FORMAT_PCM_16), sections=1, seekable=1)
        o_sfile = sndfile.open(filename, mode="w", info=outinfo)
        end = self.nf_end * self.f_div
        if end > len(self.r):  # end of file reached
            end = -1
        self.sfile.seek(self.f_start, 0)
        r = self.sfile.readf_double(self.f_end - self.f_start)  # read the sound file from f_start to f_end
        r = r / max(abs(r))  # to be between 0 and 1 ?
        # deal with stereo
        ch = self.sfile.info.channels
        ## we force into mono by recalling the method we used for the extraction		
        if ch == 2:
            if self.method == 2:
                nr = 0.5 * (r[::2] + r[1::2])
            else:
                nr = r[method::2]
        else:
            nr = r
            o_sfile.write_double(nr)
            o_sfile.close()

    def concat(self, o_call):
        if self.f_end > o_call.f_start:
            self.f_end = o_call.f_end
        else:
            pass
