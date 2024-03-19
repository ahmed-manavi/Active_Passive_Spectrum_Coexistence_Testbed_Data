import numpy as np
import pandas as pd

class TB_Calculator:
    def __init__(self, temperature_elements, avg_power, spect_power,arguement_avg,arguement_spect):
        self.temperature_elements = temperature_elements
        self.avg_power = avg_power
        self.spect_power = spect_power
        self.arguement_avg = arguement_avg
        self.arguement_spect = arguement_spect

        # Fit a linear regression line
        self.cold_power = np.array(avg_power['cold_clean_avg'],dtype=np.float64).flatten()
        self.hot_power =  np.array(avg_power['hot_clean_avg'],dtype=np.float64).flatten()
        self.x = np.array([self.cold_power,self.hot_power]).flatten()
        self.x=np.array(self.x,dtype=np.float64)
        # print(self.x)
        
        
        #self.HS1 = np.array(temperature_elements['HS'],dtype=np.float64).flatten()
        
        self.HS1 = np.array(temperature_elements['HS'], dtype=np.float64)
        self.HS_temp = 273.00+self.HS1
        self.y = np.array([51.500000000000, self.HS_temp])
        # print(self.y)
        self.y=np.array(self.y,dtype=np.float64)
        
        
        
        self.coefficients = np.polyfit(self.x, self.y, 1)
        self.slope = self.coefficients[0]
        self.intercept = self.coefficients[1]

        self.a_2 = 1 / self.db_to_linear(0.3418)
        self.a_s = 1 / self.db_to_linear(0.22)

        self.T1 = np.array(273 + temperature_elements['Ant2']).flatten()
        self.T2 = np.array(273 + temperature_elements['Ant1']).flatten()
        self.T0 = self.y[1]
        self.TA = np.array(273 + temperature_elements['Ant1'])  # Ground Truth

        self.a_1 = 130  # empirically estimated

    def db_to_linear(self, dB_value):
        return 10 ** (dB_value / 10)

    def calculate_TB_avg(self):
        H_pol_power = np.array(self.avg_power[self.arguement_avg]).flatten()
        TA_prime = (self.slope * H_pol_power) + self.intercept

        TB_avg = (TA_prime - ((self.T1 * (1 - self.a_1) * self.a_2 * self.a_s) +
                             (self.T2 * (1 - self.a_2) * self.a_s) +
                             self.T0 * (1 - self.a_s))) / (self.a_1 * self.a_2 * self.a_s)

        return TB_avg

    def calculate_TB_spect(self):
        spect_h_power = np.array(self.spect_power[self.arguement_spect])
        spect_TA_prime = self.slope * spect_h_power + self.intercept
        spect_h_temp = (spect_TA_prime - ((self.T1 * (1 - self.a_1) * self.a_2 * self.a_s) +
                                          (self.T2 * (1 - self.a_2) * self.a_s) +
                                          self.T0 * (1 - self.a_s))) / (self.a_1 * self.a_2 * self.a_s)

        len1 = len(spect_h_temp) // 3564  # 3564 is the number of time indexes of time-frequency measurements
        TB_spect = spect_h_temp.reshape(len1, 3564)

        return TB_spect

# # Usage example:

# # Assuming you have the required DataFrames temperature_elements, avg_power, and spect_power
# # Replace 'your_file.csv' with the actual path to your CSV file

# temperature_elements = pd.read_csv("Z:/Data Files/SWIFT/Version_1/No_RFI_all/Day1/arranged_temp.csv")
# temperature_elements = temperature_elements.drop(temperature_elements.columns[0], axis=1)

# avg_power = pd.read_csv("Z:/Data Files/SWIFT/Version_1/No_RFI_all/Day1/power_files/no_rfi_1_avg_power.csv")
# spect_power = pd.read_csv('Z:/Data Files/SWIFT/Version_1/No_RFI_all/Day1/power_files/no_rfi_1_spect_power.csv')

# tb_calculator = TB_Calculator(temperature_elements, avg_power, spect_power)
# TB_avg_result = tb_calculator.calculate_TB_avg()
# TB_spect_result = tb_calculator.calculate_TB_spect()

# print(f'TB_avg: {np.mean(TB_avg_result)}')
# print(f'TB_spect - TA: {np.mean(TB_spect_result) - tb_calculator.TA}')
