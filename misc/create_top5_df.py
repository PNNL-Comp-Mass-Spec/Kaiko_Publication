import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import argparse

# ###############################################################
# # arguments
# ###############################################################
parser = argparse.ArgumentParser()

parser.add_argument(
    '--test_dir', type=str, default='/Volumes/MSSHARE/Joonyong/DeepNovoRun/deepnovo/',
    help='a directory including test result files')

parser.add_argument(
    '--mgf_list_file', type=str, default='/Volumes/MSSHARE/Joonyong/mgf_list.log',
    help='mgf list file')

parser.add_argument(
    '--last_index_for_train_files', type=int, default=100,
    help='an index in a mgf target list for the train files from 0')

FLAGS = parser.parse_args()
# ###############################################################
test_dir = FLAGS.test_dir

mgf_path = '{0}/{1}_out.txt'
mgf_list_file = FLAGS.mgf_list_file
pkl_path = test_dir+'/top5_all.pkl'
###############################################################
###############################################################
all_records = []
mgf_list = pd.read_csv(mgf_list_file, sep='\t')
valid_files = list(mgf_list[mgf_list.id>FLAGS.last_index_for_train_files].mgf_file)
for i,of in enumerate(valid_files):
    print(i, of)
    tab = pd.read_csv(mgf_path.format(test_dir, of), sep='\t')
    tab['match_score'] = tab.accuracy_AA+tab.accuracy_AA_lbyl
    tab['rank'] = tab.groupby('scan')['output_score'].rank(ascending=False, method='min', na_option='top').astype('int8')
    all_records.append(tab)
    print('num of peptide seqs:', tab.shape)
top5_all = pd.concat(all_records, ignore_index=True)
print('total num of peptide seqs:', top5_all.shape)
###############################################################
###############################################################
top5_all.len_AA = top5_all.len_AA.astype('int8')
top5_all.exact_match = top5_all.exact_match.astype('int8')
top5_all.accuracy_AA = top5_all.accuracy_AA.astype('int8')
top5_all.accuracy_AA_lbyl = top5_all.accuracy_AA_lbyl.astype('int8')
top5_all.match_score = top5_all.match_score.astype('int8')
###############################################################
###############################################################

def count_len_output_AA(row):
    seq = row['output_seq']
    try:
        return len(seq.split(','))
    except:
        return 0
print("3. Adding length of target sequences")
top5_all['len_output_AA'] = top5_all.apply(count_len_output_AA, axis=1)
top5_all['len_output_AA'].astype('int8')

def get_longest_match(output_seq, target_seq, exact_match, len_AA):
    if exact_match > 0:
        return (0, len_AA)
    if type(output_seq) == float:
        return (0, 0)
    out = output_seq.split(',')
    tar = target_seq.split(',')

    length = min(len(out), len(tar))

    match_list = [1 if out[i] == tar[i] else 0 for i in range(length)]
    count_list = match_list

    # print(match_list)
    for i in range(length-1):
        if count_list[length - i - 2] == 0:
            count_list[length - i - 2] = 0
        else:
            count_list[length - i - 2] += count_list[length - i - 1]
    # print(count_list)
    max_idx = np.argmax(count_list)
    return (max_idx, count_list[max_idx])

print("4. Finding longest subsequences of matching")

# test
print(get_longest_match(top5_all.output_seq.iloc[7], top5_all.target_seq.iloc[7], top5_all.exact_match.iloc[7], top5_all.len_AA.iloc[7]))

top5_all[['longest_match_idx', 'longest_match_length']] = top5_all.apply(
    lambda row: pd.Series(get_longest_match(row['output_seq'], row['target_seq'], row['exact_match'], row['len_AA']), dtype='int8'), axis=1)
###############################################################
###############################################################
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

# print(mass_AA)

def get_longest_match_with_novor(output_seq, target_seq, exact_match, len_AA):
    if exact_match > 0:
        return (0, len_AA)
    if type(output_seq) == float:
        return (0, 0)
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
    # print('ou tput_seq:', output)
    # print('target_seq:', decoder_input)
    # print(match_list)
    count_list = match_list
    for i in range(length-1):
        if count_list[length - i - 2] == 0:
            count_list[length - i - 2] = 0
        else:
            count_list[length - i - 2] += count_list[length - i - 1]
    # print(count_list)
    max_idx = np.argmax(count_list)
    return (max_idx, count_list[max_idx])
    return 
    
# test
print("5. Finding longest subsequences of matching with Novor scoring method")

idx = 2
print(get_longest_match_with_novor(top5_all.output_seq.iloc[idx], top5_all.target_seq.iloc[idx], top5_all.exact_match.iloc[idx], top5_all.len_AA.iloc[idx]))

top5_all[['longest_match_idx_with_novor', 'longest_match_length_with_novor']] = top5_all.apply(
    lambda row: pd.Series(get_longest_match_with_novor(row['output_seq'], row['target_seq'], row['exact_match'], row['len_AA']), dtype='int8'), axis=1)
###############################################################
print(top5_all.head())

print('to_pickle')
top5_all.to_pickle(pkl_path)



