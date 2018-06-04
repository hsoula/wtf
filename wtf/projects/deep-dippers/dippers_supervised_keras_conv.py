"""
    : dippers_supervised.py
    keras
    Created: 29/05/18
    Description:

"""
import numpy as np
import numpy.random as rnd
import os
import cPickle
from dippers_dataset import get_all_frames, predict_frames_data
from scipy.stats.mstats import zscore
from sklearn.metrics import confusion_matrix
import sys
from sound import specs
from record import frames, frame_record

sys.setrecursionlimit(1500)
seed = 1
rnd.seed(seed)

print "Load Data"
# get data
raw_spec_data, group, txts = get_all_frames(nmax=60, step=0.025)
n_total, n_features = raw_spec_data.shape
fraction = 0.75
spec_data = zscore(raw_spec_data.astype('float32'))

print "Size ", n_total, n_features


print "Creating training group"
idx = np.arange(n_total)
rnd.shuffle(idx)
train_idx = idx[:int(fraction * n_total)]
valid_idx = idx[int(fraction * n_total):]
train_data = spec_data[train_idx, :]
train_group = group[train_idx]
valid_data = spec_data[valid_idx, :]
valid_group = group[valid_idx]

 
from __future__ import print_function

import keras
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.optimizers import RMSprop
from keras.layers import Convolution2D, MaxPooling2D, Conv2D
from keras.utils import np_utils

num_classes = 3

# convert class vectors to binary class matrices
train_group = keras.utils.to_categorical(train_group, num_classes)
valid_group = keras.utils.to_categorical(valid_group, num_classes)


epochs = 20
batch_size = 800

ntime = 25
nfreq = 73
start = 3
train_data_f = []
for x in range(train_data.shape[0]):
    v = train_data[x, :].astype(np.float32)
    v = v.reshape((ntime, nfreq))
    train_data_f.append(v[:,start:])

valid_data_f = []
for x in range(valid_data.shape[0]):
    v = valid_data[x, :].astype(np.float32)
    v = v.reshape((ntime, nfreq)) 
    valid_data_f.append(v[:,start:])

train_data_f = np.array(train_data_f).reshape(-1,1, ntime, nfreq-start)
valid_data_f = np.array(valid_data_f).reshape(-1,1, ntime, nfreq-start)
# 7. Define model architecture
model = Sequential()
 
#model.add(Convolution2D(32, 5, 5, activation='relu', input_shape=(1,ntime,nfreq-start),dim_ordering="th"))
model.add(Conv2D(32, (5, 5), activation="relu", data_format="channels_first", input_shape=(1,ntime,nfreq-start)))

model.add(MaxPooling2D(pool_size=(2,2)))
#model.add(Convolution2D(32, 5, 5, activation='relu'))
model.add(Conv2D(32, (3, 3), activation="relu"))
model.add(MaxPooling2D(pool_size=(2,2)))
model.add(Dropout(0.5))

model.add(Flatten())
model.add(Dense(100, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(num_classes, activation='softmax'))


# 8. Compile model
model.compile(loss='categorical_crossentropy',
              optimizer='adam',
              metrics=['accuracy'])
 
# 9. Fit model on training data
model.fit(train_data_f, train_group, 
          batch_size=32, 
          epochs=20, 
          verbose=1,
          shuffle=True,
          validation_split=0.2)
 
# 10. Evaluate model on test data
score = model.evaluate(valid_data_f, valid_group, verbose=0)

print('Test loss:', score[0])
print('Test accuracy:', score[1])

from dippers_predict import prepare_data

dir_wav = "../../data/dippers/data_2014/waves/"
dir_txt = "../../data/dippers/data_2014/texts"


def prepare_data(wav_file, offset):
    f = 1000.0
    nstd = 6
    fband = 250
    freq_max = 20000
    freq_min = 1000
    done = False
    multiplier = 10
    #offset = 0
    f_re = frame_record(wav_file, offset)
    curr_time = offset
    rfr = []
    nfreq = -1
    print(wav_file)
    while not done:
        ret = f_re.get_frames(multiplier, time, freq_min, freq_max, f, nstd, fband, method=1)
        if ret == 0:
            done = True
        for ix, fr in enumerate(f_re.frames):
            nl, nf = np.shape(fr.a)
            if nl == np.floor(f * time):
                ax = (fr.a.flatten() / np.max(np.abs(fr.a))).astype(np.float32)
                rfr.append(ax)
            curr_time += time
    return rfr


time = 0.025  # in s
def predict_file(model, fname, rname):
    ret = []
    with open(rname, 'w') as iOF:
        for offset in np.arange(0, time, 0.005):
            rfr = np.array(prepare_data(fname, offset))
            rfr = rfr.reshape(-1,1, ntime,nfreq)
            rfr = rfr[:,:,:,start:]
            v = model.predict_classes(zscore(rfr).astype(np.float32))
            for i, u in enumerate(v):
                iOF.write('%f %f\n' % (i * time + offset, u))
            ret.append(v)
    iOF.close()
    return np.array(ret)

def predict_rep_all(net, rep, target_rep):
    for file in os.listdir(rep):
        if file.endswith('.wav') and not file.startswith("._"):
            print(file)
            fname = rep+'/' + file
            tname = target_rep+'/'+file
            predict_file(net, fname, tname)


predict_rep_all(model, dir_wav, 'predict_2014')

