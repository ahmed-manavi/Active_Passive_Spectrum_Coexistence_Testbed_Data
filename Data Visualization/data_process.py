# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 17:32:54 2023

@author: aa2863
"""

import numpy as np
from spectrogram_gen import SpectrogramGenerator
from psd_gen import PowerSpectralDensity
from butterworth_lowpass import ButterworthLowpassFilter

class DataProcessor:
    def __init__(self, sample_data_iq, sample_fs, cutoff, fs, order):
        self.sample_data_iq = sample_data_iq
        self.sample_fs = sample_fs
        self.cutoff = cutoff
        self.fs = fs
        self.order = order
        self.fftsize1 = 1024 * 2
        self.window1 = 1024 * 2
        self.channel_data = []
        self.total_filter_data = []
        self.total_psd_data_unfilter = []
        self.spect_filter = []
        self.spect_unfilter = []

    def process_data(self):
        channel_width = int(self.sample_fs / 4)
        segment_number = int(len(self.sample_data_iq) / channel_width)

        for i in range(segment_number):
            tmp = self.sample_data_iq[i * channel_width:(channel_width + (i * channel_width)),]
            remove_data = tmp[int(1e5):int(channel_width - 1e5)]
            self.channel_data.append(remove_data)

            # PSD Calculation (unfiltered)
            psd_calculator = PowerSpectralDensity(remove_data[:], self.fftsize1, self.window1)
            psd_rad_rd, segment1 = psd_calculator.compute_psd()
            self.total_psd_data_unfilter.append(psd_rad_rd[:])

            # Spectrogram Calculation (unfiltered)
            spectrogram_generator = SpectrogramGenerator(remove_data[:], self.fftsize1, self.window1)
            sp_unfilter, sp_seg_uf = spectrogram_generator._calculate_spectrogram()
            self.spect_unfilter.append(sp_unfilter)

            # Filtering with lowpass filter
            butterworth_filter = ButterworthLowpassFilter(self.cutoff, self.fs, self.order)
            tmp_data_filter = butterworth_filter.apply_filter(remove_data)

            # PSD Calculation (filtered)
            psd_calculator = PowerSpectralDensity(tmp_data_filter[:], self.fftsize1, self.window1)
            psd_filter, self.segment1 = psd_calculator.compute_psd()
            self.total_filter_data.append(psd_filter[:])

            # Spectrogram Calculation (filtered)
            spectrogram_generator = SpectrogramGenerator(tmp_data_filter[:], self.fftsize1, self.window1)
            sp_filter, sp_seg_f = spectrogram_generator._calculate_spectrogram()
            self.spect_filter.append(sp_filter)

        return (
            self.total_psd_data_unfilter,
            self.spect_unfilter,
            self.total_filter_data,
            self.spect_filter,
            self.channel_data,
            self.segment1
        )



# Now you can access the processed data as needed.
