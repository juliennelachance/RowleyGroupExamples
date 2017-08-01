#####################################
#
# Store .wav file->.mat file reader
# Author: Julienne LaChance 
#
#####################################


import wave
import struct
import numpy
import scipy.io as sio
import pylab


# Read file: 
mat_contents = sio.loadmat('normal_test.mat')
wav_dict = mat_contents['wav_dict']
val = wav_dict[0,0]
t = val['t']
samples = val['samples']

print numpy.shape(t)
print numpy.shape(samples)


pylab.plot(t[0,:], samples[0,:])
pylab.show()
