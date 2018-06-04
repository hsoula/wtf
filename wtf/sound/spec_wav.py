import numpy as np
from scipy.interpolate import interp1d
from spec import pyGaussianSpectrum


class spec_wav(object):
    def __init__(self, wave_data, sampleRate=44100., dB=0.1, f=1000.0, freq_min=200.0, freq_max=8000.0, nstd=6,
                 fband=125):
        """
             Compute the spectrogram between freq_min and freq_max with nstd deviation for the gaussian kernel
             and fband spectrum.
             The temporal sampling frequency is given by f and the sample wave_data is of sampleRate.
             dB is the threshold above with the data is kept. To avoid silence.
             Should not be used for song as this!
         """
        self.wave_data = 1.0 * wave_data[:]
        twindow = nstd / (fband * 2.0 * np.pi)
        self.winLength = np.floor(twindow * sampleRate)
        self.winLength = np.round(self.winLength / 2) * 2
        self.increment = np.round(sampleRate / f)
        spec, fc, fft = pyGaussianSpectrum(wave_data, self.increment, self.winLength)
        self.whole_spec = spec[:]  ## to be used in invertandadd later
        cspec = np.reshape(spec[::2] + spec[1::2] * 1.j, (fc, fft))
        self.fc = fc
        t = np.arange(0, fc) / f
        fs = np.arange(0, fft) * sampleRate / fft
        self.spectrum = np.abs(cspec)
        a = np.sum(self.spectrum.T, 0)
        self.idx = np.abs(a) > max(np.abs(a)) * dB
        self.udx = np.where(np.abs(a) > max(np.abs(a)) * dB)[0]
        imin = np.ceil(freq_min / sampleRate * fft)
        imax = np.ceil(freq_max / sampleRate * fft)
        fs = fs[imin:imax]
        t = t[self.idx]
        self.spectrum = self.spectrum[self.idx, :]
        self.fs = fs[:]
        self.t = t[:]
        self.imin = imin
        self.imax = imax
        self.cspec = cspec
        self.a = a[self.idx]
        start = self.udx[0] * self.increment
        end = self.udx[-1] * self.increment
        self.wave_data = self.wave_data[start:end]


class specs(object):
    def __init__(self, wave_data,
                 sample_rate=44100.,
                 f=1000.0,
                 nstd=6,
                 fband=125):
        """
             Compute the spectrogram with nstd deviation for the gaussian kernel
             and fband spectrum.
             The temporal sampling frequency is given by f
             and the sample wave_data is of sampleRate.
         """
        self.wave_data = 1.0 * wave_data[:]
        self.sample_rate = sample_rate
        twindow = nstd / (fband * 2.0 * np.pi)
        self.winLength = np.floor(twindow * sample_rate)
        self.winLength = np.round(self.winLength / 2) * 2
        self.increment = np.round(sample_rate / f)
        spec, fc, fft = pyGaussianSpectrum(wave_data, self.increment, self.winLength)
        self.whole_spec = spec[:]  # to be used in invertandadd later
        cspec = np.reshape(spec[::2] + spec[1::2] * 1.j, (fc, fft))
        self.fc = fc
        self.t = np.arange(0, fc) / f
        self.fs = np.arange(0, fft) * sample_rate / fft
        self.spectrum = np.abs(cspec)
        self.fft = fft

    def get_envelope(self, freq_min=0, freq_max=22050.):
        imin = np.ceil(freq_min / self.sample_rate * self.fft)
        imax = np.ceil(freq_max / self.sample_rate * self.fft)
        c = self.spectrum[:,imin:imax]
        return np.sum(c, 1)
