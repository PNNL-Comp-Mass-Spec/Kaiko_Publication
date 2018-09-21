"""Prepare for kaiko results analysis
Finally, it generates a pickle file, which contains all output results and the number of longest matches.

USAGE:
python analysis_results.py --result_dir /Volumes/MSSHARE/Joonyong/DeepNovoRun_0-57_50/deepnovo/mgf_test/ --mgf_list_file /Volumes/MSSHARE/Joonyong/mgf_list.log --list_index_for_train_files 57
"""

import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

import argparse

###############################################################
# arguments
###############################################################
parser = argparse.ArgumentParser()

parser.add_argument(
    '--result_dir', type=str, default='/Volumes/MSSHARE/Joonyong/DeepNovoRun/deepnovo/',
    help='a directory including test result files')

parser.add_argument(
    '--mgf_list_file', type=str, default='/Volumes/MSSHARE/Joonyong/mgf_list.log',
    help='mgf list file')

parser.add_argument(
    '--last_index_for_train_files', type=int, default=234,
    help='the last index in a mgf target list for the train files from 0')
FLAGS = parser.parse_args()
###############################################################


result_dir = FLAGS.result_dir # '/Volumes/MSSHARE/Joonyong/DeepNovoRun/deepnovo/'
print("result_dir:", result_dir)
out_files = glob.glob(result_dir + '/*_out.txt')
print('num of files:', len(out_files))

mgf_tab = pd.read_csv(FLAGS.mgf_list_file, sep='\t')

# read all
###############################################################
print("----- Reading all output files -----")
all_records = []
for of in out_files:
    tab = pd.read_csv(of, sep='\t')
    tab[['file_id', "scan_number"]] = pd.DataFrame(tab.scan.str.split(':').tolist(), columns=['file_id', "scan_number"], dtype="int32")
    
    # check out scan ids
    sample_file = os.path.basename(of).rsplit('_out', 1)[0]
    true_scan_id = mgf_tab[mgf_tab.mgf_file == sample_file].id.tolist()[0]
    cur_scan_id = tab.iloc[0].file_id

    if true_scan_id > FLAGS.last_index_for_train_files:
      if cur_scan_id != true_scan_id:
          print('Scan ID is incorrect:', cur_scan_id, true_scan_id)
          tab.file_id = true_scan_id
          tab["scan"] = str(true_scan_id) + ":" + tab['scan_number'].apply(str)
      print("Read file id:", true_scan_id, sample_file)
      tab['match_score'] = tab.accuracy_AA+tab.accuracy_AA_lbyl
      all_records += tab.to_dict('records')

all_records = pd.DataFrame(all_records)
print('num of peptide seqs:', all_records.shape[0])
print('num of unique peptide seqs:', len(all_records.target_seq.drop_duplicates().tolist()))
###############################################################
###############################################################
def count_len_output_AA(row):
    seq = row['output_seq']
    try:
        return len(seq.split(','))
    except:
        return 0
print("----- Adding length of target sequences -----")
all_records['len_output_AA'] = all_records.apply(count_len_output_AA, axis=1)
all_records['len_output_AA'].astype(int)
###############################################################
###############################################################
print("----- Finding longest subsequences of matching with Novor scoring method -----")

mass_H = 1.0078
mass_H2O = 18.0106
mass_NH3 = 17.0265
mass_N_terminus = 1.0078
mass_C_terminus = 17.0027
mass_CO = 27.9949

mass_AA = {'_PAD': 0.0,
           '_GO': mass_N_terminus-mass_H,
           '_EOS': mass_C_terminus+mass_H,
           'A': 71.03711, # 0
           'R': 156.10111, # 1
           'N': 114.04293, # 2
           'Nmod': 115.02695,
           'D': 115.02694, # 3
           #~ 'C': 103.00919, # 4
           'Cmod': 160.03065, # C(+57.02)
           #~ 'Cmod': 161.01919, # C(+58.01) # orbi
           'E': 129.04259, # 5
           'Q': 128.05858, # 6
           'Qmod': 129.0426,
           'G': 57.02146, # 7
           'H': 137.05891, # 8
           'I': 113.08406, # 9
           'L': 113.08406, # 10
           'K': 128.09496, # 11
           'M': 131.04049, # 12
           'Mmod': 147.0354,
           'F': 147.06841, # 13
           'P': 97.05276, # 14
           'S': 87.03203, # 15
           'T': 101.04768, # 16
           'W': 186.07931, # 17
           'Y': 163.06333, # 18
           'V': 99.06841, # 19
          }

def get_longest_match(output_seq, target_seq, len_AA, inv=False):
    if (type(output_seq) == float) | (output_seq == ""):
        return (0, 0)
    if inv:
        output = output_seq.split(',')[::-1]
        decoder_input = target_seq.split(',')[::-1]
    else:
        output = output_seq.split(',')
        decoder_input = target_seq.split(',')
    
    decoder_input_len = len(decoder_input)
    output_len = len(output)
    decoder_input_mass = [mass_AA[x] for x in decoder_input]
    decoder_input_mass_cum = np.cumsum(decoder_input_mass)
    output_mass = [mass_AA[x] for x in output]
    output_mass_cum = np.cumsum(output_mass)
    
    length = min(decoder_input_len, output_len)
    
    num_match = 0
    i = 0
    j = 0
    match_list = [0 for i in range(decoder_input_len)]
    
    while i < decoder_input_len and j < output_len:
        if abs(decoder_input_mass_cum[i] - output_mass_cum[j]) < 0.5:
            if abs(decoder_input_mass[i] - output_mass[j]) < 0.1:
                num_match += 1
                match_list[i] = 1
            i += 1
            j += 1
        elif decoder_input_mass_cum[i] < output_mass_cum[j]:
            i += 1
        else:
            j += 1
    
    count_list = match_list
    for i in range(length-1):
        if count_list[length - i - 2] == 0:
            count_list[length - i - 2] = 0
        else:
            count_list[length - i - 2] += count_list[length - i - 1]
    
    max_idx = np.argmax(count_list)
    return (max_idx, count_list[max_idx])
    
# test
idx = 10
print(get_longest_match("V,M,Mmod,E,A,S,H,H,K,P,G,K", "V,Mmod,E,A,S,H,H,K,P,G,K", 11, True))
print(get_longest_match("V,M,Mmod,E,A,S,H,H,K,P,G,K", "V,Mmod,E,A,S,H,H,K,P,G,K", 11, False))

all_records[['longest_match_idx_with_novor', 'longest_match_length_with_novor']] = all_records.apply(
    lambda row: pd.Series(get_longest_match(row['output_seq'], row['target_seq'], row['len_AA'], inv=False), dtype='int32'), axis=1)
all_records[['longest_match_idx_inv', 'longest_match_length_inv']] = all_records.apply(
    lambda row: pd.Series(get_longest_match(row['output_seq'], row['target_seq'], row['len_AA'], inv=True), dtype='int32'), axis=1)
all_records['longest_match'] = all_records[['longest_match_length_inv', 'longest_match_length_with_novor']].max(axis=1)
###############################################################

# store
print("----- Storing pickle -----")
all_records.to_pickle(result_dir + '/all_results_df.pkl')
