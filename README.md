# README for 5G-Transmitted IQ Samples Dataset

## Data Storage
![Data Format](https://github.com/ahmed-manavi/Active_Passive_Spectrum_Coexistence_Testbed_Data/assets/95959590/ed6d8b6f-4fc2-4281-95be-976d9c28fca5)

## Overview
This dataset contains IQ samples received by a radiometer, derived from experiments with 5G-transmitted waveforms under various scenarios. These experiments were conducted to facilitate the research on spectrum coexistence between active wireless communication systems and passive sensing technologies.

The dataset includes **raw IQ samples**, **spectral data**, and **calibrated brightness temperature (TB) measurements**. The data is stored in the **Hierarchical Data Format (HDF5)**, which ensures efficient organization and sharing within the scientific community.

---
# Link for Datasets 
- https://bit.ly/4gYfeAV

## Folder Structure
The dataset is organized into the following folders, each representing a specific experimental scenario:
- **`inband_v1`**: Contains data from fully in-band transmission experiments where 5G waveforms were transmitted entirely within the operating bandwidth of the radiometer (1400-1427 MHz).
- **`outband_v1`**: Contains data from out-of-band experiments, where 5G transmissions were conducted outside the radiometer's operational bandwidth.
- **`transitionband_v1`**: Contains data from transition-band experiments, combining both in-band and out-of-band transmissions.

Each folder contains multiple subdirectories named systematically to reflect specific attributes of the experiments, including **resource blocks**, **modulation techniques**, **power gain levels**, and **sample numbers**.

---

## Dataset Attributes
The dataset is categorized into three levels:

1. **Level 0 (L0)**:
   - Raw IQ data split into horizontal polarization (H-pol) and vertical polarization (V-pol).
   - Includes reference sources (`ref1` and `ref2`).

2. **Level 1A (L1A)**:
   - Spectral data derived from raw antenna counts.
   - Includes:
     - Power Spectral Density (PSD).
     - Short Time Fourier Transform (STFT).
   - Data is provided for:
     - Before digital filtering.
     - After digital filtering (−10 MHz to +10 MHz).
     - Internal RFI removed (IRR).

3. **Level 1B (L1B)**:
   - Fully calibrated TB data.
   - Includes spectral time-frequency features and average TB for different polarization states.
   - Ground truth TB data (e.g., anechoic room temperature) is also recorded.

### File Naming Convention
Files are named using a systematic convention for easy identification:
`<transmission_band>_<resource_blocks>_<PowerGain>_<ModulationTechnique>_<SampleNumber>.h5`

**Example**:  
For an out-of-band transmission with 4 RBs, 0 gain, and QPSK modulation in L1B data, the filename will be:  
`outband_4RB_Gain0_QPSK_fc1_rfi_L1B_SN1.h5`

---

## Data Preprocessing
Preprocessing steps applied to the data include:
1. **Integration Time**:
   - A single sample represents 1 second of integration, divided into four 0.25-second states (H-pol, V-pol, ref1, and ref2).
2. **Digital Filtering**:
   - Butterworth low-pass filter (cutoff: 10 MHz) to remove high-frequency noise.
3. **RFI Removal**:
   - An internal RFI detection algorithm was applied to remove unwanted noise and outliers.

---

## Data Access
The dataset can be accessed from the shared OneDrive directory. The structure is as follows:
```
OneDrive/
├── inband_v1/
├── outband_v1/
└── transitionband_v1/
```

---

### Dependencies

- **Python**: Version 3.8+
- **Key Libraries**:
  - `h5py`: For HDF5 file handling.
  - `numpy`: For numerical processing.
  - `matplotlib`: For data visualization.
  - `pandas`: For data management.
  - `glob` and `os`: For file system operations.
  - `re`: For regular expressions and numerical sorting.

### Install Dependencies

Run the following command to install all the required libraries:

```bash
pip install -r requirements.txt
```

---
## Citation
If you use this dataset in your research, please cite:
> Alam, A.M., et al. (2024). *Physical Testbed and Open Dataset for Passive Sensing and Wireless Communication Spectrum Coexistence*. IEEE Access. DOI: [10.1109/ACCESS.2024.3453774](https://doi.org/10.1109/ACCESS.2024.3453774)

---

For questions or contributions, please contact the corresponding author:  
**Vuk Marojevic** - [vuk.marojevic@ece.msstate.edu](mailto:vuk.marojevic@ece.msstate.edu)





# Overall MSU Testbed
![Testbed Picture_1](https://github.com/user-attachments/assets/e8d1aff5-4294-4f7d-9338-dbc9bba128a8)


- 

## Experimental Scenario
![df701e924fi2 pdf-1](https://github.com/ahmed-manavi/Active_Passive_Spectrum_Coexistence_Testbed_Data/assets/95959590/341a1d5f-176a-4da6-819d-39bf7ff6ad9c)





# Active_Passive_Spectrum_Coexistence
![Picture1](https://github.com/ahmed-manavi/Active_Passive_Spectrum_Coexistence_Testbed_Data/assets/95959590/deb9db6d-f607-4bb9-80ba-aeb7ec957c9f)


In our rapidly expanding realm of advanced satellite and communication systems, passive radiometer sensors utilized in Earth observation, such as those for 5G, face an escalating challenge posed by radio frequency interference (RFI) stemming from human-generated signals. Effective methodologies to assess the impact of 5G on Earth observation radiometers are urgently needed. Unfortunately, the scarcity of substantial datasets in the radio frequency (RF) domain, particularly concerning active/passive coexistence, impedes progress. Our research presents a controlled test environment comprising a calibrated L-band radiometer and a 5G wireless communication system. Housed within a controlled chamber, this unique setup enables the observation and quantification of transmission effects across various frequency bands. Through the creation of a comprehensive dataset, our objective is to standardize and benchmark both wireless communication and passive sensing. Leveraging the capacity to analyze raw measurements, our test environment facilitates the detection and mitigation of RFI, promoting the harmonious coexistence of wireless communication and passive sensing technologies while establishing vital standards.

Dataset Download Link (V1): https://mstate-my.sharepoint.com/:f:/g/personal/aa2863_msstate_edu/Esvo2OU23pBDq1UmFqgFqJ8BHmAMKpDkuCrCwgNHwBEX_Q?e=2dNgtM

## Overall Testbed
![photo 2](https://github.com/ahmed-manavi/Active_Passive_Spectrum_Coexistence_Testbed_Data/assets/95959590/844fbcbc-ffb1-42e4-bf9d-3501a227e55e)

## Testbed Schematic
![photo 3](https://github.com/ahmed-manavi/Active_Passive_Spectrum_Coexistence_Testbed_Data/assets/95959590/65b5c397-ce79-4f4c-9a67-4d89d80c84bd)

## Data Storage
![Data Format](https://github.com/ahmed-manavi/Active_Passive_Spectrum_Coexistence_Testbed_Data/assets/95959590/ed6d8b6f-4fc2-4281-95be-976d9c28fca5)

## Experimental Scenario
![df701e924fi2 pdf-1](https://github.com/ahmed-manavi/Active_Passive_Spectrum_Coexistence_Testbed_Data/assets/95959590/341a1d5f-176a-4da6-819d-39bf7ff6ad9c)

## References
1. A. M. Alam, M. M. Farhad, W. Al-Qwider, A. Owfi, M. Koosha, N. Mastronarde, F. Afghah, V. Marojevic, M. Kurum, and A. C. Gurbuz, "A Physical Testbed and Open Dataset for Benchmarking of Passive Sensing and Wireless Communication Spectrum Coexistence," in IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing (Submitted)
2. A. M. Alam, M. Kurum and A. C. Gurbuz, "Radio Frequency Interference Detection for SMAP Radiometer Using Convolutional Neural Networks," in IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, vol. 15, pp. 10099-10112, 2022, doi: 10.1109/JSTARS.2022.3223198.
3. A. M. Alam, M. Kurum, M. Ogut and A. C. Gurbuz, "Microwave Radiometer Calibration Using Deep Learning With Reduced Reference Information and 2-D Spectral Features," in IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, vol. 17, pp. 748-765, 2024, doi: 10.1109/JSTARS.2023.3333268.
4. M. M. Farhad, A. M. Alam, S. Biswas, M. A. S. Rafi, A. C. Gurbuz, and M. Kurum, “Sdr-based dual polarized l-band microwave radiometer operating from small uas platforms,” IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 2023
5. W. Al-Qwider, A. M. Alam, M. Mehedi Farhad, M. Kurum, A. C. Gurbuz and V. Marojevic, "Software Radio Testbed for 5G and L-Band Radiometer Coexistence Research," IGARSS 2023 - 2023 IEEE International Geoscience and Remote Sensing Symposium, Pasadena, CA, USA, 2023, pp. 596-599, doi: 10.1109/IGARSS52108.2023.10283002.

# Contact: Ahmed Manavi Alam (aa2863@msstate.edu) for any details
