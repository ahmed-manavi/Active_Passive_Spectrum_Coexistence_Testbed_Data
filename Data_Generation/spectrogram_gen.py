# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 14:48:07 2023

@author: aa2863
"""

import numpy as np
from scipy.fft import fft, fftshift

class SpectrogramGenerator:
    def __init__(self, signal, fftsize, window):
        self.signal = signal
        self.fftsize = int(fftsize)
        self.window = int(window)
        self.y, self.segment = self._calculate_spectrogram()

    def _calculate_spectrogram(self):
        y = []
        s = self.signal

        if (len(s) % self.window) != 0:
            temp = np.zeros(self.window - (self.window - (len(s) % self.window)))
            temp = temp.reshape(len(temp), 1)
            s = s.reshape(len(s), 1)
            kk = np.concatenate([s, temp], axis=0)
            s = kk.flatten()

        segment = int(len(s) / self.window)

        for i in range(segment):
            index1 = i * self.window
            index2 = (i * self.window) + self.window
            tmp_wave = s[index1:index2, ]
            temp = fft(tmp_wave, n=self.fftsize)
            temp = fftshift(temp)
            y.append(temp)

        y = np.array(y)
        y = np.transpose(y)
        return y, segment

    # def get_spectrogram(self):
    #     return self.y

    # def get_segment(self):
    #     return self.segment
