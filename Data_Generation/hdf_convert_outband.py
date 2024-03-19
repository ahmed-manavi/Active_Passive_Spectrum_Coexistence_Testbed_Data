
from file_reader import FileReader
# from spectrogram_gen import SpectrogramGenerator
# from psd_gen import PowerSpectralDensity
from data_process import DataProcessor
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# from butterworth_lowpass import ButterworthLowpassFilter
import copy   
import h5py
import os
from convert_TB import TB_Calculator


loaded_array = np.load('C:/Users/Public/radiometer_calibration/Data_chamber/Saving_Data_Hdf_Format/outliers_no_rfi/outlier_mask_v3.npy')
base_path_hdf = r"Z:/Data Files/SWIFT/HDF_processed_version_1/RFI/outband5"
Base_sc16 = 'Z:/Data Files/SWIFT/version5/outband'
RB_dir = 'Four_RB'
main_folder_path = os.path.join(Base_sc16, RB_dir)

directories_all_fc = [d for d in os.listdir(main_folder_path) if os.path.isdir(os.path.join(main_folder_path, d))]
inband_directories_all_fc = [d for d in directories_all_fc if d.startswith('Four_RB')]

directories_files_dict = {}  # Dictionary to store files for each directory
temp_directories_files_dict = {}
for inband_dir in inband_directories_all_fc:
    dir_path = os.path.join(main_folder_path, inband_dir)
    sub_directories = [d for d in os.listdir(dir_path) if os.path.isdir(os.path.join(dir_path, d)) and 'outband' in d]


    files_dict = {}  # Dictionary to store files for each subdirectory
    temp_files_dict = {}
    for sub_dir in sub_directories:
        sub_dir_path = os.path.join(dir_path, sub_dir)
        files_temp = [f for f in os.listdir(sub_dir_path) if os.path.isfile(os.path.join(sub_dir_path, f)) and "arranged_temp" in f]

        files = [f for f in os.listdir(sub_dir_path) if os.path.isfile(os.path.join(sub_dir_path, f)) and f.endswith(".sc16")]
        files_dict[sub_dir] = files
        temp_files_dict[sub_dir] = files_temp

    directories_files_dict[inband_dir] = files_dict
    temp_directories_files_dict[inband_dir] = temp_files_dict

# Now directories_files_dict contains the files ending with ".sc16" for each subdirectory within each directory in inband_directories_all_fc
print(directories_files_dict)
print(temp_directories_files_dict)



empty_folders = []

# Loop through the dictionary to find any empty entries
for inband_dir, files_dict in temp_directories_files_dict.items():
    for sub_dir, files_list in files_dict.items():
        if not files_list:  # Check if the list is empty
            empty_folders.append((inband_dir, sub_dir))

if empty_folders:
    print("Empty entries found in temp_directories_files_dict:")
    for entry in empty_folders:
        print(f"Directory '{entry[0]}' and Subdirectory '{entry[1]}' has no files.")
else:
    print("No empty entries found in temp_directories_files_dict.")





# Your existing code for populating directories_files_dict

# Iterate through directories and read files
for inband_dir, files_dict in directories_files_dict.items():
    print(f"Directory: {inband_dir}")
    for sub_dir, files in files_dict.items():
        print(f"Subdirectory: {sub_dir}")
        # Extracting values from subdirectory name
        band, rb, gain, modulation, fc = sub_dir.split('_')
        file_path_L0 = os.path.join(base_path_hdf, fc, rb, gain, modulation, "L0")
        if not os.path.exists(file_path_L0):
            os.makedirs(file_path_L0)
        file_path_L1A = os.path.join(base_path_hdf, fc, rb, gain, modulation, "L1A")
        if not os.path.exists(file_path_L1A):
            os.makedirs(file_path_L1A)
        file_path_L1B = os.path.join(base_path_hdf, fc, rb, gain, modulation, "L1B")
        if not os.path.exists(file_path_L1B):
            os.makedirs(file_path_L1B)
       
        temp_files = temp_directories_files_dict[inband_dir][sub_dir]
        # temp_location = temp_files[0]
        SN_count = 1
    
    
        for file in files:
            files_path = os.path.join(Base_sc16, RB_dir, inband_dir, sub_dir, file)
            print(files_path)
            
            
            filename_L0 = sub_dir + f"_rfi_L0_SN{SN_count}.h5"
            filename_L0_path = os.path.join(file_path_L0, filename_L0)
         
            filename_L1A = sub_dir + f"_rfi_L1A_SN{SN_count}.h5"
            filename_L1A_path = os.path.join(file_path_L1A, filename_L1A)
            
            filename_L1B = sub_dir + f"_rfi_L1B_SN{SN_count}.h5"
            filename_L1B_path = os.path.join(file_path_L1B, filename_L1B)
         
         # filename_L1A = file_path_L1A + f"rfi_L1A_SN{SN_count}.h5"
         # filename_L1B = file_path_L1B + f"rfi_L1B_SN{SN_count}.h5"
        # filename_5G = path_to_files + f"in_band_8RB_10dB_5G_SN{SN_count}.h5"
         # Create an HDF5 file
            if os.path.exists(filename_L0_path):
               print(f"File '{filename_L1B_path}' already exists. Skipping...")    
               continue 
            if os.path.exists(filename_L1A_path):
               print(f"File '{filename_L1B_path}' already exists. Skipping...")    
               continue
            if os.path.exists(filename_L1B_path):
                print(f"File '{filename_L1B_path}' already exists. Skipping...")    
                continue
                
                
                
                # for file in sorted_files:
                    
                    # print(file)
                    # file_name = os.path.splitext(os.path.basename(file))[0]
                    # print(file_name)
                    
                        
            # file_path = os.path.join(directory_data_files, file)
            sample_data = FileReader(files_path)
            samples = sample_data.read_samples()
            
            
            # time = len(sample_data)// (sample_fs*2)
            # time = 3
            sample_fs = 30000000
            # sample_data_iq1=samples[0:sample_fs*2*1].view(np.complex64)
            # sample_data_iq2=samples[sample_fs*2*2:sample_fs*2*3].view(np.complex64)
            # sample_data_iq3=samples[sample_fs*2*4:sample_fs*2*5].view(np.complex64)
            # sample_data_iq4=samples[sample_fs*2*6:sample_fs*2*7].view(np.complex64)
            # sample_data_iq5=samples[sample_fs*2*8:sample_fs*2*9].view(np.complex64)
            
            
            # sample_data_iq = np.concatenate([sample_data_iq1,sample_data_iq2,sample_data_iq3,sample_data_iq4,sample_data_iq5])
            fftsize=1024*2
            window=1024*2
            sample_data_iq=[]
            num_segments=5
            for jj in range(num_segments):
                start_index = sample_fs * 2 * (2*jj)
                end_index = sample_fs * 2 * ((2*jj) + 1)
                segment_data_iq = samples[start_index:end_index].view(np.complex64)
                sample_data_iq.append(segment_data_iq)
            
            
            comp_opts=4
                
            for i1 in range(len(sample_data_iq)): 
                order = 28
                fs = 29.2*1e6      # sample rate, Hz (Although sampling rate is 30 MHz but after removing transient samples (related to switching of each state) from each channel it becomes equivalent to 29.2 MHz samples)
                cutoff = 10*1e6   
                sample_second_wise = sample_data_iq[i1]
                # Example usage:
                # Create an instance of the DataProcessor class
                data_processor = DataProcessor(sample_second_wise, sample_fs, cutoff, fs, order)
                psd_unfilter, spect_unfilter, psd_filter, spect_filter,segment_data_tr,segment_psd = data_processor.process_data()
                
                #####
                
            
                total_filter_data_psd = copy.deepcopy(psd_filter[:])
                total_filter_data_spect = copy.deepcopy(spect_filter[:])
                spect_irr = []
                outlier_quantile1 = []
                total_clean_data_with_mean1 = []
                psd_irr = []
                spect_power = []
                avg_power_psd=[]
                #loaded_array = np.load('C:/Users/Public/radiometer_calibration/Data_chamber/Saving_Data_Hdf_Format/outliers_no_rfi/outlier_mask_v2.npy')
                
                for i2 in range(len(total_filter_data_psd)):
                    
                    tmp_data_psd =copy.deepcopy(total_filter_data_psd[i2])
                    tmp_data_spect = copy.deepcopy(total_filter_data_spect[i2])
             
                    
                    spect_var = copy.deepcopy(tmp_data_spect[:,:]) 
                    
                    ccc1 =  copy.deepcopy(tmp_data_psd[:]) 
                    kkk1 = copy.deepcopy(tmp_data_psd[:])
                
                    outlier_remove_data1 = np.delete( kkk1[:] ,loaded_array[:])
                    psd_irr.append(outlier_remove_data1[:])
                    
                    no_mean_power = (np.sum((outlier_remove_data1)))/(fftsize*(fftsize-len(loaded_array))*segment_psd)
                    avg_power_psd.append(no_mean_power)
                    
                    outlier_remove_data_spect = np.delete( spect_var,loaded_array[:],axis=0)
                    spect_irr.append(outlier_remove_data_spect[:])
                    
                    
                    temp_spect_power = ((np.abs(outlier_remove_data_spect[:]))**2) / (fftsize-len(loaded_array))
                    spect_power.append(temp_spect_power.flatten())
                    
                    
                    # print(f'spect powwr : psd power {(np.mean(temp_spect_power))/no_mean_power}')
            
                    
                df_psd = pd.DataFrame({
                       'H_clean_avg': [avg_power_psd[0]],
                       'V_clean_avg': [avg_power_psd[1]],
                       'hot_clean_avg': [avg_power_psd[2]],
                       'cold_clean_avg': [avg_power_psd[3]]
                   })
                 
                
                # filename_power=f"Z:/Data Files/SWIFT/Version_1/No_RFI_all/Day1/power_files/no_rfi_{i1+1}_avg_power.csv"
                # df_psd.to_csv(filename_power, index=False)
               
                spect_data = {
                    'H_clean_spect': spect_power[0],
                    'V_clean_spect': spect_power[1],
                    'hot_clean_spect': spect_power[2],
                    'cold_clean_spect': spect_power[3]
                }
            
                df_spect = pd.DataFrame(spect_data)
                # filename_spect = f"Z:/Data Files/SWIFT/Version_1/No_RFI_all/Day1/power_files/no_rfi_{i1+1}_spect_power.csv"
                # df_spect.to_csv(filename_spect, index=False)
                
                for temp_file in temp_files:
                    temp_file_path = os.path.join(Base_sc16, RB_dir, inband_dir, sub_dir, temp_file)
                    temperature_elements = pd.read_csv(temp_file_path)  # Read the temperature file using pandas

            
              
                # temperature_elements = pd.read_csv("Z:/Data Files/SWIFT/Version_1/four_RB/temp_files/arranged_temp_4RB.csv")
                temperature_elements = temperature_elements.drop(temperature_elements.columns[0], axis=1)
                temperature_elements = temperature_elements.mean()
                TB_G = [((273 + temperature_elements['Ant1']))]
                TB_HS = [(273 + temperature_elements['HS'])]
                Tphy_ACS = [(273 + temperature_elements['ACS'])]
                TB_ACS = [51.5]
                
                
                
                avg_power = df_psd
                spect_power =  spect_data
              
                # avg_power = pd.read_csv("Z:/Data Files/SWIFT/Version_1/No_RFI_all/Day1/power_files/no_rfi_1_avg_power.csv")
                # spect_power = pd.read_csv('Z:/Data Files/SWIFT/Version_1/No_RFI_all/Day1/power_files/no_rfi_1_spect_power.csv')
                
                tb_calculator_h = TB_Calculator(temperature_elements, avg_power, spect_power,'H_clean_avg','H_clean_spect')
                TB_avg_result_h = tb_calculator_h.calculate_TB_avg()
                TB_spect_result_h = tb_calculator_h.calculate_TB_spect()
                
                
                tb_calculator_v = TB_Calculator(temperature_elements, avg_power, spect_power,'V_clean_avg','V_clean_spect')
                TB_avg_result_v = tb_calculator_v.calculate_TB_avg()
                TB_spect_result_v = tb_calculator_v.calculate_TB_spect()
                
                # print("done")
                
                
                
                # plt.figure()
                # # im1 = plt.imshow(y_res1, origin='lower', cmap='jet',  aspect='auto',vmin=298,vmax=310)
                # im1 = plt.imshow(TB_spect_result_h, origin='lower', cmap='jet',  aspect='auto',vmin=280,vmax=320)
                # plt.colorbar(label='Brightness Temperature (K)')
                # # yaxis=np.ceil(np.linspace(-15,15,11))
                # # ylabel=np.linspace(0,len1,11)
                # # xaxis=(np.linspace(0,.25,5))
                # # xlabel=np.linspace(0,len(y_res1[1]),5)
                # # plt.xlabel('Time (s)')
                # # plt.ylabel('Frequency (MHz)')
                # # plt.yticks(ylabel,yaxis)
                # # plt.xticks(xlabel,xaxis)
                # ####
                # plt.figure()
                # # im1 = plt.imshow(y_res1, origin='lower', cmap='jet',  aspect='auto',vmin=298,vmax=310)
                # im1 = plt.imshow(TB_spect_result_v, origin='lower', cmap='jet',  aspect='auto',vmin=280,vmax=320)
                # plt.colorbar(label='Brightness Temperature (K)')
            
                # plt.figure()
                # plt.plot(10*np.log10((psd_unfilter[0])))
                # # plt.plot(((total_filter_data[8])))
                # xaxis=np.ceil(np.linspace(-15,15,11))
                # # yaxis=np.ceil(np.linspace(9,15,7))
                # xlabel=np.linspace(0,len(psd_unfilter[0]),11)
                # plt.xlabel('frequency [MHz]')
                # plt.ylabel('PSD [V**2] dB')
                # plt.xticks(xlabel,xaxis)
                # # plt.yticks(yaxis)
                # plt.grid()
            
                
                
                ###############################################
                
                # path_to_files = "Z:/Data Files/SWIFT/HDF_processed_version_1/No_RFI/"
               
            
                filename_L0 = sub_dir + f"_rfi_L0_SN{SN_count}.h5"
                filename_L0_path = os.path.join(file_path_L0, filename_L0)
                
                filename_L1A = sub_dir + f"_rfi_L1A_SN{SN_count}.h5"
                filename_L1A_path = os.path.join(file_path_L1A, filename_L1A)
                
                filename_L1B = sub_dir + f"_rfi_L1B_SN{SN_count}.h5"
                filename_L1B_path = os.path.join(file_path_L1B, filename_L1B)
                
                # filename_L1A = file_path_L1A + f"rfi_L1A_SN{SN_count}.h5"
                # filename_L1B = file_path_L1B + f"rfi_L1B_SN{SN_count}.h5"
               # filename_5G = path_to_files + f"in_band_8RB_10dB_5G_SN{SN_count}.h5"
                # Create an HDF5 file
                if not os.path.exists(filename_L0_path):
                    with h5py.File(filename_L0_path, "w") as file0:
                        # Create the "L0" group
                        # l0_group = file0.create_group("L0")
                    
                        # Create the "IQ_samples" dataset inside the "L0" group
                        
                        h_pol_group_iq = file0.create_group("H-pol")
                        # Create the "V-pol" group inside "PSD"
                        v_pol_group_iq = file0.create_group("V-pol")
                        # Create the "ref1" group inside "PSD"
                        ref1_group_iq = file0.create_group("ref1")
                        # Create the "ref2" group inside "PSD"
                        ref2_group_iq = file0.create_group("ref2")
                        for p in range(0, 1):  # You mentioned psd_h_1, psd_h_2, and psd_h_3
                        # h_pol_group_iq.create_dataset("H-pol", data=sample_second_wise[0:int(sample_fs/4)],compression="gzip",compression_opts=comp_opts)
                        # v_pol_group_iq.create_dataset("V-pol", data=sample_second_wise[int(sample_fs/4):int(sample_fs/2)],compression="gzip",compression_opts=comp_opts)
                        # ref1_group_iq.create_dataset("ref1", data=sample_second_wise[int(sample_fs/2):int(sample_fs*(3/4))],compression="gzip",compression_opts=comp_opts)
                        # ref2_group_iq.create_dataset("ref2", data=sample_second_wise[int(sample_fs*(3/4)):sample_fs],compression="gzip",compression_opts=comp_opts)
                            h_pol_group_iq.create_dataset("H-pol", data=segment_data_tr[4*p],compression="gzip",compression_opts=comp_opts)
                            v_pol_group_iq.create_dataset("V-pol", data=segment_data_tr[(4*p)+1],compression="gzip",compression_opts=comp_opts)
                            ref1_group_iq.create_dataset("ref1", data=segment_data_tr[(4*p)+2],compression="gzip",compression_opts=comp_opts)
                            ref2_group_iq.create_dataset("ref2", data=segment_data_tr[(4*p)+3],compression="gzip",compression_opts=comp_opts)
                 
                else:
                    print(f"File '{filename_L0_path}' already exists. Skipping...")
                
    
                if not os.path.exists(filename_L1A_path):         
                    with h5py.File(filename_L1A_path, "w") as file1:
                        l1a_group = file1.create_group("L1A")
                    
                        # Create the "PSD" subgroup inside "L1A"
                        psd_group = l1a_group.create_group("PSD")
                        # Create the "H-pol" group inside "PSD"
                        h_pol_group = psd_group.create_group("H-pol")
                        # Create the "V-pol" group inside "PSD"
                        v_pol_group = psd_group.create_group("V-pol")
                        # Create the "ref1" group inside "PSD"
                        ref1_group = psd_group.create_group("ref1")
                        # Create the "ref2" group inside "PSD"
                        ref2_group = psd_group.create_group("ref2")
                        
                        
                        
                        # Create the "STFT" subgroup inside "L1A" (empty for now)
                        stft_group = l1a_group.create_group("STFT")
                        h_pol_group_stft = stft_group.create_group("H-pol")
                        # Create the "V-pol" group inside "PSD"
                        v_pol_group_stft = stft_group.create_group("V-pol")
                        # Create the "ref1" group inside "PSD"
                        ref1_group_stft = stft_group.create_group("ref1")
                        # Create the "ref2" group inside "PSD"
                        ref2_group_stft = stft_group.create_group("ref2")
                        
                        # print("hello")
                    
                        # Create PSD datasets
                        # for i in range(0, 1):  # You mentioned psd_h_1, psd_h_2, and psd_h_3
                        #psd_data = np.random.rand(100)  # Replace with your actual PSD data
                        i=0
                        h_pol_group.create_dataset("psd_h_unfilter", data=psd_unfilter[4*i],compression="gzip",compression_opts=comp_opts)
                        h_pol_group.create_dataset("psd_h_filter", data=psd_filter[4*i],compression="gzip",compression_opts=comp_opts)
                        h_pol_group.create_dataset("psd_h_irr", data=psd_irr[4*i],compression="gzip",compression_opts=comp_opts)
                        
                        #psd_data = np.random.rand(100)  # Replace with your actual PSD data
                        v_pol_group.create_dataset("psd_v_unfilter", data=psd_unfilter[(4*i)+1],compression="gzip",compression_opts=comp_opts)
                        v_pol_group.create_dataset("psd_v_filter", data=psd_filter[(4*i)+1],compression="gzip",compression_opts=comp_opts)
                        v_pol_group.create_dataset("psd_v_irr", data=psd_irr[(4*i)+1],compression="gzip",compression_opts=comp_opts)
                        
                        #psd_data = np.random.rand(100)  # Replace with your actual PSD data
                        ref1_group.create_dataset("psd_ref1_unfilter", data=psd_unfilter[(4*i)+2],compression="gzip",compression_opts=comp_opts)
                        ref1_group.create_dataset("psd_ref1_filter", data=psd_filter[(4*i)+2],compression="gzip",compression_opts=comp_opts)
                        ref1_group.create_dataset("psd_ref1_irr", data=psd_irr[(4*i)+2],compression="gzip",compression_opts=comp_opts)
                        
                        #psd_data = np.random.rand(100)  # Replace with your actual PSD data
                        ref2_group.create_dataset("psd_ref2_unfilter", data=psd_unfilter[(4*i)+3],compression="gzip",compression_opts=comp_opts)
                        ref2_group.create_dataset("psd_ref2_filter", data=psd_filter[(4*i)+3],compression="gzip",compression_opts=comp_opts)
                        ref2_group.create_dataset("psd_ref2_irr", data=psd_irr[(4*i)+3],compression="gzip",compression_opts=comp_opts)
                        
                        
                        
                        #psd_data = np.random.rand(100)  # Replace with your actual PSD data
                        h_pol_group_stft.create_dataset("stft_h_unfilter", data=spect_unfilter[4*i],compression="gzip",compression_opts=comp_opts)
                        h_pol_group_stft.create_dataset("stft_h_filter", data=spect_filter[4*i],compression="gzip",compression_opts=comp_opts)
                        h_pol_group_stft.create_dataset("stft_h_irr", data=spect_irr[4*i],compression="gzip",compression_opts=comp_opts)
                        
                        #psd_data = np.random.rand(100)  # Replace with your actual PSD data
                        v_pol_group_stft.create_dataset("stft_v_unfilter", data=spect_unfilter[(4*i)+1],compression="gzip",compression_opts=comp_opts)
                        v_pol_group_stft.create_dataset("stft_v_filter", data=spect_filter[(4*i)+1],compression="gzip",compression_opts=comp_opts)
                        v_pol_group_stft.create_dataset("stft_v_irr", data=spect_irr[(4*i)+1],compression="gzip",compression_opts=comp_opts)
                        
                        #psd_data = np.random.rand(100)  # Replace with your actual PSD data
                        ref1_group_stft.create_dataset("stft_ref1_unfilter", data=spect_unfilter[(4*i)+2],compression="gzip",compression_opts=comp_opts)
                        ref1_group_stft.create_dataset("stft_ref1_filter", data=spect_filter[(4*i)+2],compression="gzip",compression_opts=comp_opts)
                        ref1_group_stft.create_dataset("stft_ref1_irr", data=spect_irr[(4*i)+2],compression="gzip",compression_opts=comp_opts)
                        
                        #psd_data = np.random.rand(100)  # Replace with your actual PSD data
                        ref2_group_stft.create_dataset("stft_ref2_unfilter", data=spect_unfilter[(4*i)+3],compression="gzip",compression_opts=comp_opts)
                        ref2_group_stft.create_dataset("stft_ref2_filter", data=spect_filter[(4*i)+3],compression="gzip",compression_opts=comp_opts)
                        ref2_group_stft.create_dataset("stft_ref2_irr", data=spect_irr[(4*i)+3],compression="gzip",compression_opts=comp_opts)
                    
                else:
                    print(f"File '{filename_L1A_path}' already exists. Skipping...")
                    
                 
                if not os.path.exists(filename_L1B_path): 
                     
                    with h5py.File(filename_L1B_path, "w") as file2:
                    # Create the "L1A" group
                   
                    
                        # Create the "L1B" group (empty for now)
                        l1b_group = file2.create_group("L1B")
                        
                        # l1a_group = file1.create_group("L1A")
                    
                        # Create the "PSD" subgroup inside "L1A"
                        TB_group = l1b_group.create_group("TB")
                        # Create the "H-pol" group inside "PSD"
                        spectral = TB_group.create_group("Spectral")
                        # Create the "V-pol" group inside "PSD"
                        Avg = TB_group.create_group("Avg")
                        # Create the "ref1" group inside "PSD"
                        GT = TB_group.create_group("GT")
                        # Create the "ref2" group inside "PSD"
                        
                        TB_cal_group = l1b_group.create_group("TB_cal")
                        
                
                        spectral.create_dataset("H-pol", data=TB_spect_result_h,compression="gzip",compression_opts=comp_opts)
                        spectral.create_dataset("V-pol", data=TB_spect_result_v,compression="gzip",compression_opts=comp_opts)
                        
                        
                        Avg.create_dataset("H-pol", data=TB_avg_result_h,compression="gzip",compression_opts=comp_opts)
                        Avg.create_dataset("V-pol", data=TB_avg_result_v,compression="gzip",compression_opts=comp_opts)
                        
                        GT.create_dataset("ground_truth_temp", data=TB_G,compression="gzip",compression_opts=comp_opts)
                        
                        
                        TB_cal_group.create_dataset("TB_Ref_1_HS", data=TB_HS,compression="gzip",compression_opts=comp_opts)
                        TB_cal_group.create_dataset("Tphy_Ref_2_ACS", data=Tphy_ACS,compression="gzip",compression_opts=comp_opts)
                        TB_cal_group.create_dataset("TB_Ref_2_ACS", data=TB_ACS,compression="gzip",compression_opts=comp_opts)
                        
                                       
                else:
                    print(f"File '{filename_L1B_path}' already exists. Skipping...")     
                        
                    
                # with h5py.File(filename_5G, "w") as file3:
                # # Create the "L1A" group
               
                
                #     # Create the "L1B" group (empty for now)
                #     iq5G_wave = file3.create_group("5G")
            
                # print("HDF5 file 'your_file.h5' created successfully.")
                
                
                # import h5py
                
                # # Open the HDF5 file in read mode
                # with h5py.File("in_band_8RB_10dB.h5", "r") as file:
                #     # Recursively explore the file structure and print the keys
                #     def print_keys(name, obj):
                #         if isinstance(obj, h5py.Group):
                #             print(f"Group: {name}")
                #         elif isinstance(obj, h5py.Dataset):
                #             print(f"Dataset: {name}")
                
                #     file.visititems(print_keys)
                
                print(i1)
                SN_count = 1 + SN_count





# channel_width = int(sample_fs/4)
# segment_number = int(len(sample_data_iq)/channel_width) # Each 1s has 4 channels, so 9s data will have 36 channel info

# channel_data = []
# total_filter_data = []
# total_psd_data_unfilter=[]




# fftsize1=1024*2
# window1=1024*2
# spect_filter = []
# spect_unfilter = []

# for i in range(segment_number):
    
#     tmp = sample_data_iq[i*channel_width:(channel_width+(i*channel_width)),] # Each 1s samples has 4 channel information -- H. V, Hot, Cold
#     remove_data = tmp [int(1e5):int(channel_width-1e5)] # from each channel data we are removing first 1e5 and last 1e5 samples
#     channel_data.append(remove_data)  # storing each channel  
    
    
  
#     ##### converting channel data into psd 
#     psd_calculator = PowerSpectralDensity(remove_data[:], fftsize1, window1)
#     psd_rad_rd,segment1 = psd_calculator.compute_psd()
#     total_psd_data_unfilter.append(psd_rad_rd[:])

    
#     #####
#     spectrogram_generator = SpectrogramGenerator(remove_data[:], fftsize1, window1)
#     sp_unfilter,sp_seg_uf = spectrogram_generator._calculate_spectrogram()
#     spect_unfilter.append(sp_unfilter)

    
    
#     ######
    
#     ##### Filtering with lowpass filter -10 to 10 MHz
#     butterworth_filter = ButterworthLowpassFilter(cutoff, fs, order)
#     tmp_data_filter=butterworth_filter.apply_filter(remove_data)

    
#     psd_calculator = PowerSpectralDensity(tmp_data_filter[:], fftsize1, window1)
#     psd_filter,segment1 = psd_calculator.compute_psd()
#     total_filter_data.append(psd_filter[:])
    
#     spectrogram_generator = SpectrogramGenerator(tmp_data_filter[:], fftsize1, window1)
#     sp_filter,sp_seg_f = spectrogram_generator._calculate_spectrogram()
#     spect_filter.append(sp_filter)



# spect_1s_unfilter = np.concatenate([spect_filter[0],spect_filter[1],spect_filter[2],spect_filter[3]],axis=1)
# as_tmp = spect_unfilter[0]
# y_abs1=abs(spect_1s_unfilter)
# global_maximum = np.amax(y_abs1)
# y_res1 = 20*np.log10(y_abs1/global_maximum)
# plt.figure()
# im1 = plt.imshow(y_res1, origin='lower', cmap='jet',  aspect='auto',vmax=0,vmin=-50)
# plt.colorbar(label='dB')
# plt.xticks(xlabel,xaxis)
# plt.yticks(ylabel,yaxis)
# yaxis=np.ceil(np.linspace(-15,15,11))
# ylabel=np.linspace(0,len(psd_rad_rd),11)
# xaxis=(np.linspace(0,1,5))
# xlabel=np.linspace(0,len(spect_1s_unfilter[1]),5)
# plt.yticks(ylabel,yaxis)
# plt.xticks(xlabel,xaxis)
# plt.xlabel('Time')
# plt.ylabel('Frequency')
# plt.savefig('spect_unfilter1.png',dpi=300,bbox_inches='tight')



# spect_1s_filter = np.concatenate([spect_filter[0],spect_filter[1],spect_filter[2],spect_filter[3]],axis=1)

# # as_tmp = spect_filter[0]
# y_abs1=abs(spect_1s_filter)
# y_res1 = 20*np.log10(y_abs1/global_maximum)
# plt.figure()
# im1 = plt.imshow(y_res1, origin='lower', cmap='jet',  aspect='auto',vmax=0,vmin=-50)
# plt.colorbar(label='dB')
# yaxis=np.ceil(np.linspace(-15,15,11))
# ylabel=np.linspace(0,len(psd_rad_rd),11)
# xaxis=(np.linspace(0,1,5))
# xlabel=np.linspace(0,len(spect_1s_unfilter[1]),5)
# plt.xticks(xlabel,xaxis)
# plt.yticks(ylabel,yaxis)
# plt.xlabel('Time')
# plt.ylabel('Frequency')

###############################################






# spectrogram_generator = SpectrogramGenerator(sample_data_iq[3], fftsize, window)
# spectrogram,segment = spectrogram_generator._calculate_spectrogram()


# psd_calculator = PowerSpectralDensity(sample_data_iq[3], fftsize, window)
# psd,segment = psd_calculator.compute_psd()



# y_abs1=abs(spectrogram)

# y_res1 = 20*np.log10(y_abs1/np.amax(y_abs1))
# plt.figure()
# im1 = plt.imshow(y_res1, origin='lower', cmap='jet',  aspect='auto',vmax=0,vmin=-50)
# plt.colorbar(label='dB')
# # plt.xticks(xlabel,xaxis)
# # plt.yticks(ylabel,yaxis)
# plt.xlabel('Time')
# plt.ylabel('Frequency')

# y_res1 = 20*np.log10(y_abs1)
# plt.figure()
# im1 = plt.imshow(y_res1, origin='lower', cmap='jet',  aspect='auto')#,vmax=0,vmin=-50)
# plt.colorbar(label='dB')
# # plt.xticks(xlabel,xaxis)
# # plt.yticks(ylabel,yaxis)
# plt.xlabel('Time')
# plt.ylabel('Frequency')
            # temp_file_path = os.path.join(Base_sc16, RB_dir, inband_dir, sub_dir, temp_location)
            # temperature_elements_temp = pd.read_csv(temp_file_path) 
            
        #     # Do what you need with the file (read it, process it, etc.)
        # for temp_file in temp_files:
        #     temp_file_path = os.path.join(Base_sc16, RB_dir, inband_dir, sub_dir, temp_file)
        #     temperature_elements_temp = pd.read_csv(temp_file_path)  # Read the temperature file using pandas









################################################ Practice Do not uncomment

# -*- coding: utf-8 -*-
# """
# Created on Fri Nov  3 18:32:45 2023

# @author: aa2863
# """

# import os

# Base_sc16 = 'Z:/Data Files/SWIFT/Version_1/out_band'
# RB_dir = 'Eight_RB'
# main_folder_path = os.path.join(Base_sc16, RB_dir)

# directories_all_fc = [d for d in os.listdir(main_folder_path) if os.path.isdir(os.path.join(main_folder_path, d))]
# inband_directories_all_fc = [d for d in directories_all_fc if d.startswith('Eight_RB')]

# for inband_dir in inband_directories_all_fc:
#     dir_path = os.path.join(main_folder_path, inband_dir)
#     sub_directories = [d for d in os.listdir(dir_path) if os.path.isdir(os.path.join(dir_path, d))]
#     print(f"Subdirectories in '{inband_dir}': {sub_directories}")


# import os

# Base_sc16 = 'Z:/Data Files/SWIFT/Version_1/out_band'
# RB_dir = 'Eight_RB'
# main_folder_path = os.path.join(Base_sc16, RB_dir)

# directories_all_fc = [d for d in os.listdir(main_folder_path) if os.path.isdir(os.path.join(main_folder_path, d))]
# inband_directories_all_fc = [d for d in directories_all_fc if d.startswith('Eight_RB')]

# subdirectories_dict = {}  # Dictionary to store subdirectories for each directory

# for inband_dir in inband_directories_all_fc:
#     dir_path = os.path.join(main_folder_path, inband_dir)
#     sub_directories = [d for d in os.listdir(dir_path) if os.path.isdir(os.path.join(dir_path, d))]
#     sub_directories = [d for d in sub_directories if d.startswith('outband')]
#     subdirectories_dict[inband_dir] = sub_directories

# # Now subdirectories_dict contains the subdirectories for each directory in inband_directories_all_fc
# print(subdirectories_dict)


# import os

# Base_sc16 = 'Z:/Data Files/SWIFT/Version_1/out_band'
# RB_dir = 'Eight_RB'
# main_folder_path = os.path.join(Base_sc16, RB_dir)

# directories_all_fc = [d for d in os.listdir(main_folder_path) if os.path.isdir(os.path.join(main_folder_path, d))]
# inband_directories_all_fc = [d for d in directories_all_fc if d.startswith('Eight_RB')]

# directories_files_dict = {}  # Dictionary to store files for each directory

# for inband_dir in inband_directories_all_fc:
#     dir_path = os.path.join(main_folder_path, inband_dir)
#     sub_directories = [d for d in os.listdir(dir_path) if os.path.isdir(os.path.join(dir_path, d))]

#     files_dict = {}  # Dictionary to store files for each subdirectory

#     for sub_dir in sub_directories:
#         sub_dir_path = os.path.join(dir_path, sub_dir)
#         files = [f for f in os.listdir(sub_dir_path) if os.path.isfile(os.path.join(sub_dir_path, f))]
#         files_dict[sub_dir] = files

#     directories_files_dict[inband_dir] = files_dict

# # Now directories_files_dict contains the files for each subdirectory within each directory in inband_directories_all_fc
# print(directories_files_dict)







