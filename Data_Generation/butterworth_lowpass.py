
from scipy.signal import butter, lfilter

class ButterworthLowpassFilter:
    def __init__(self, cutoff, fs, order):
        self.cutoff = cutoff
        self.fs = fs
        self.order = order
        self.b, self.a = self.butter_lowpass()

    def butter_lowpass(self):
        return butter(self.order, self.cutoff, fs=self.fs, btype='low', analog=False)

    def apply_filter(self, data):
        y = lfilter(self.b, self.a, data)
        return y
