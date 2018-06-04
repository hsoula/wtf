"""
    : dippers_dataset.py
    Created: 10s/03/16
    Description: Dipper dataset loading and setting up

    Needs: textfile for annotation + wave file
    fname.txt + fname.wav

    Create a dictionnary to store the data and pickle it.

"""

import numpy as np
import sys
sys.path.append('../../')
from wtf.sound import specs
from wtf.record import frames, frame_record
#import sndfile
from scipy.stats.mstats import zscore

# main parameters


dir_wav = "~/Research/birds/data/dippers/data_2014/waves/"
dir_txt = "~/Research/birds/data/dippers/data_2014/texts"

def get_frames_data(fname, step=0.025, save_file=False, sample_freq=1000.0, fband=250, time=0.025):
    txt_file = open(dir_txt + "/" + fname + '.txt')
    wav_file = dir_wav + "/" + fname + '.wav'

    all_frames = []
    all_lbl_data = []

    # f = 1000.0
    nstd = 6
    # fband = 250
    freq_max = 20000
    freq_min = 1000
    done = False
    # time = 0.025  # in ms
    multiplier = 10
    nfreq = -1

    for offset in np.arange(0, time, step):
        frames = []
        f_re = frame_record(wav_file, offset)
        while not done:
            ret = f_re.get_frames(multiplier, time, freq_min, freq_max, sample_freq, nstd, fband, nfreq)
            if ret == 0:
                done = True
            frames.extend(f_re.frames)
        done = False
        if save_file:
            xname = fname + ".spec"
            f_re.save_to_text(xname)
        ns = len(frames)
        lbl_data = np.zeros(ns)

        for line in txt_file:
            if line.startswith("\"T\"") or line.startswith("\"N\""):
                split_line = line.split()
                call_type = split_line[0][1]
                c_type = 1
                if call_type == "N":
                    c_type = 2
                start = float(split_line[1]) - offset
                end = float(split_line[2]) - offset
                nstart = int(np.floor(start / time))
                nend = int(np.floor(end / time) + 1)
                lbl_data[nstart:nend] = c_type
        all_frames.extend(frames)
    txt_file.close()
    all_lbl_data.append(lbl_data)
#    f_re.sfile.close()
    all_lbl_data = np.hstack(all_lbl_data)
    return all_frames, all_lbl_data



def predict_frames_data(fname, neural_net, save_file=False):
    txt_file = open(dir_txt + "/" + fname + '.txt')
    wav_file = dir_wav + "/" + fname + '.wav'

    f = 1000.0
    nstd = 6
    fband = 250
    freq_max = 20000
    freq_min = 1000
    done = False
    time = 0.025  # in ms
    multiplier = 10
    with  frame_record(wav_file) as f_re:
        while not done:
            ret = f_re.get_frames(multiplier, time, freq_min, freq_max, f, nstd, fband)
            if ret == 0:
                done = True
        ns = len(f_re.frames)
        lbl_data = np.zeros(ns)
        p_time = np.arange(0, ns * time, time)
        if save_file:
            xname = fname + ".spec"
            f_re.save_to_text(xname)
        f = f_re.frames[0]
        ux = (f.a / np.max(np.abs(f.a))).flatten()
        ux = ux.astype(np.float32)
        data = ux[:]
        n = len(data)
        for ix in range(1,ns):
            f = f_re.frames[ix]
            ux = (f.a / np.max(np.abs(f.a))).flatten()
            ux = ux.astype(np.float32)
            if len(ux) == n:
                data = np.vstack((data, ux))
            else:
                nux = np.concatenate((ux, np.zeros(n-len(ux))))
                data = np.vstack((data, nux.astype(np.float32)))
        preds_lbl = neural_net.predict(zscore(data))

        for line in txt_file:
            if line.startswith("\"T\"") or line.startswith("\"N\""):
                split_line = line.split()
                call_type = split_line[0][1]
                c_type = 1
                if call_type == "N":
                    c_type = 2
                start = float(split_line[1])
                end = float(split_line[2])
                nstart = np.floor(start / time)
                nend = np.floor(end / time) + 1
                lbl_data[nstart:nend] = c_type
    #f_re.sfile.close()
    return f_re.frames, lbl_data, preds_lbl, p_time


def get_all_frames(nmax=40, step=0.025, savefile=False, sample_freq=1000.0, fband=250, time=0.025):
    import os
    imx_chicks = 0
    imx_incub = 0
    lbls = []
    ix = 0
    ldir = os.listdir(dir_txt)
    txts = []
    from random import shuffle
    from random import seed
    seed(1)
    shuffle(ldir)
    for txt_file in ldir:
        if txt_file.endswith('.txt') and 'incub' in txt_file:
            print txt_file

            txts.append(txt_file)

            frames, labels = get_frames_data(txt_file[:-4],
                                             step=step,
                                             save_file=savefile,
                                             sample_freq=sample_freq,
                                             fband=fband,
                                             time=time)

            for f, l in zip(frames, labels):
                ux = (f.a / np.max(np.abs(f.a))).flatten()

                if ix == 0:
                    data = ux[:]
                    n = len(data)
                    lbls.append(l)
                    ix += 1
                else:
                    if len(ux) == n:
                        data = np.vstack((data, ux))
                        lbls.append(l)

            imx_incub += 1
            if imx_incub > nmax/2:
                break
    for txt_file in ldir:
        if txt_file.endswith('.txt') and 'chicks' in txt_file:
            print txt_file

            txts.append(txt_file)

            frames, labels = get_frames_data(txt_file[:-4], step=step, save_file=savefile)

            for f, l in zip(frames, labels):
                ux = (f.a / np.max(np.abs(f.a))).flatten()
                if len(ux) == n:
                    data = np.vstack((data, ux))
                    lbls.append(l)

            imx_chicks += 1
            if imx_chicks > nmax/2:
                break

    return data, np.array(lbls), txts


def predict_frames_conv_data(fname, neural_net, save_file=False):
    txt_file = open(dir_txt + "/" + fname + '.txt')
    wav_file = dir_wav + "/" + fname + '.wav'

    f = 1000.0
    nstd = 6
    fband = 250
    freq_max = 20000
    freq_min = 1000
    done = False
    time = 0.025  # in ms
    multiplier = 10
    f_re = frame_record(wav_file)
    nx = int(f * time)
    while not done:
        ret = f_re.get_frames(multiplier, time, freq_min, freq_max, f, nstd, fband)
        if ret == 0:
            done = True
    ns = len(f_re.frames)
    lbl_data = np.zeros(ns)
    p_time = np.arange(0, ns * time, time)
    data = []
    if save_file:
        xname = fname + ".spec"
        f_re.save_to_text(xname)
    for ix in range(ns):
        f = f_re.frames[ix]
        ux = (f.a / np.max(np.abs(f.a)))
        ux = ux.astype(np.float32)
        if ux.shape[0] < nx:
            ux = np.vstack((ux, np.zeros((nx - ux.shape[0], ux.shape[1]))))
        data.append(ux)
    data = np.reshape(data, (-1, 1, ux.shape[0], ux.shape[1])).astype(np.float32)
    preds_lbl = neural_net.predict(zscore(data))

    for line in txt_file:
        if line.startswith("\"T\"") or line.startswith("\"N\""):
            split_line = line.split()
            call_type = split_line[0][1]
            c_type = 1
            if call_type == "N":
                c_type = 2
            start = float(split_line[1])
            end = float(split_line[2])
            nstart = np.floor(start / time)
            nend = np.floor(end / time) + 1
            lbl_data[nstart:nend] = c_type
    return f_re.frames, lbl_data, preds_lbl, p_time


def predict_frames_conv_keras_data(fname, keras_model, save_file=False):
    txt_file = open(dir_txt + "/" + fname + '.txt')
    wav_file = dir_wav + "/" + fname + '.wav'

    f = 1000.0
    nstd = 6
    fband = 250
    freq_max = 20000
    freq_min = 1000
    done = False
    time = 0.025  # in ms
    multiplier = 10
    f_re = frame_record(wav_file)
    nx = int(f * time)
    while not done:
        ret = f_re.get_frames(multiplier, time, freq_min, freq_max, f, nstd, fband)
        if ret == 0:
            done = True
    ns = len(f_re.frames)
    lbl_data = np.zeros(ns)
    p_time = np.arange(0, ns * time, time)
    data = []
    if save_file:
        xname = fname + ".spec"
        f_re.save_to_text(xname)
    for ix in range(ns):
        f = f_re.frames[ix]
        ux = (f.a / np.max(np.abs(f.a)))
        ux = ux.astype(np.float32)
        if ux.shape[0] < nx:
            ux = np.vstack((ux, np.zeros((nx - ux.shape[0], ux.shape[1]))))
        data.append(ux)
    data = np.reshape(data, (-1, 1, ux.shape[0], ux.shape[1])).astype(np.float32)
    preds_lbl = keras_model.predict_classes(zscore(data))

    for line in txt_file:
        if line.startswith("\"T\"") or line.startswith("\"N\""):
            split_line = line.split()
            call_type = split_line[0][1]
            c_type = 1
            if call_type == "N":
                c_type = 2
            start = float(split_line[1])
            end = float(split_line[2])
            nstart = np.floor(start / time)
            nend = np.floor(end / time) + 1
            lbl_data[nstart:nend] = c_type
    return f_re.frames, lbl_data, preds_lbl, p_time
