import numpy as np

class FileReader:
    def __init__(self, file_name):
        self.file_name = file_name
        self.dt = np.dtype('<i2')

    def read_samples(self, count=-1):
        samples = np.fromfile(self.file_name, self.dt, count).astype(np.float32)
        # If you want to return complex numbers, uncomment the next line
        # samples = np.fromfile(self.file_name, self.dt, count).astype(np.float32).view(np.complex64)
        return samples
