# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 15:59:44 2023

@author: aa2863
"""

import numpy as np
from scipy.fft import fft, fftshift

class PowerSpectralDensity:
    def __init__(self, signal, fftsize, window):
        self.signal = signal
        self.fftsize = int(fftsize)
        self.window = int(window)

    def compute_psd(self):
        y = []

        if len(self.signal) % self.window != 0:
            temp = np.zeros(self.window - (self.window - (len(self.signal) % self.window)))
            temp = temp.reshape(len(temp), 1)
            self.signal = self.signal.reshape(len(self.signal), 1)
            kk = np.concatenate([self.signal, temp], axis=0)
            self.signal = kk.flatten()

        segment = len(self.signal) // self.window

        for i in range(segment):
            index1 = i * self.window
            index2 = (i * self.window) + self.window
            tmp_wave = self.signal[index1:index2,]
            temp = fft(tmp_wave, n=self.fftsize)
            temp = fftshift(temp)
            y.append(temp)

        y = np.array(y)
        y = np.transpose(y)
        y = (np.abs(y)) ** 2
        y = np.sum(y, 1)
        return y,segment
