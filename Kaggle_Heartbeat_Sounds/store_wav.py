#####################################
#
# Store .wav file data as MatLab-style .mat file.
# Author: Julienne LaChance 
#
#####################################


import wave
import struct
import numpy
import scipy.io as sio

# -----------------------------------------------------
# User-defined input: 
FNAME = './normal_test.wav'	# The input .wav file.
# -----------------------------------------------------

f = wave.open(FNAME)

# frames will hold the bytestring representing all the audio frames
frames = f.readframes(-1)
samples = struct.unpack('h'*f.getnframes(), frames)
framerate = f.getframerate()
t = [float(i)/framerate for i in range(len(samples))]

wav_dict = {'t': t, 'samples': samples}
sio.savemat('normal_test.mat', {'wav_dict': wav_dict})
