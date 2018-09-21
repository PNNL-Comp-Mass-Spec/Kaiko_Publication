"""
python plot_figures.py --result_dir /Volumes/MSSHARE/Joonyong/DeepNovoRun_0-57/deepnovo/mgf_test --mgf_list_file /Volumes/MSSHARE/Joonyong/mgf_list.log --last_index_for_train_files 57 --training_log /Volumes/MSSHARE/Joonyong/DeepNovoRun_0-57/deepnovo.train.model_0-57/log_file_caption_2dir.tab --result_pickle /Volumes/MSSHARE/Joonyong/DeepNovoRun_0-57/deepnovo/mgf_test/all_results_df.pkl

python plot_figures.py --result_dir /Volumes/MSSHARE/Joonyong/DeepNovoRun_0-17/deepnovo --mgf_list_file /Volumes/MSSHARE/Joonyong/mgf_list.log --last_index_for_train_files 17 --training_log /Volumes/MSSHARE/Joonyong/DeepNovoRun_0-17/deepnovo.train.model_0-17/log_file_caption_2dir.tab --result_pickle /Volumes/MSSHARE/Joonyong/DeepNovoRun_0-17/deepnovo/all_results_df.pkl
"""

import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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
    '--result_pickle', type=str, default='/Volumes/MSSHARE/Joonyong/DeepNovoRun/deepnovo/all_results_df.pkl',
    help='mgf list file')

parser.add_argument(
    '--last_index_for_train_files', type=int, default=100,
    help='the last index in a mgf target list for the train files from 0')

parser.add_argument(
    '--training_log', type=str, default="/Volumes/MSSHARE/Joonyong/DeepNovoRun/log_file_caption_2dir.tab",
    help='an index in a mgf target list for the train files from 0')

FLAGS = parser.parse_args()
###############################################################
print("Loading data", FLAGS.result_pickle)
all_records = pd.read_pickle(FLAGS.result_pickle)

mgf_tab = pd.read_csv(FLAGS.mgf_list_file, sep='\t')
mgf_tab.head()

train_list = list(mgf_tab[mgf_tab.id<=FLAGS.last_index_for_train_files].mgf_file)
test_list = list(mgf_tab[mgf_tab.id>FLAGS.last_index_for_train_files].mgf_file)

#### data split into training / validation
all_records[['file_id', "scan_number"]] = pd.DataFrame(all_records.scan.str.split(':').tolist(), columns=['file_id', "scan_number"], dtype="int32")

result_dir = FLAGS.result_dir
print("result_dir:", result_dir)
out_files = glob.glob(result_dir + '/*_out.txt')
print('num of files:', len(out_files))

total_log_path = result_dir + '/log.txt'
print('log file:', total_log_path)

sumary_tab = pd.read_csv(total_log_path, sep='\t')

training_log= pd.read_csv(FLAGS.training_log, sep='\t', header=0)

plt.close('all')

###############################################################
# figure 0. training trend
###############################################################
print("Figure 0. training trend")

def training_trend(training_log):
    fig, ax = plt.subplots(1, figsize=(10,8))

    ax.plot(training_log.perplexity_train, label="Training loss")
    ax.plot(training_log.perplexity_valid_0)
    ax.plot(training_log.perplexity_valid_1)
    ax.plot(training_log.perplexity_valid_2)
    plt.legend()

    plt.savefig(result_dir + "/fig0_training_loss.pdf")

    fig, ax = plt.subplots(1, figsize=(10,8))

    ax.plot(training_log.last_accuracy_peptide_train, label="Training Acurracy")
    ax.plot(training_log.accuracy_peptide_valid_0)
    ax.plot(training_log.accuracy_peptide_valid_1)
    ax.plot(training_log.accuracy_peptide_valid_2)
    plt.legend()

    plt.savefig(result_dir + "/fig0_training_acc.pdf")

training_trend(training_log)

###############################################################
# figure 1. accuracy histogram
###############################################################
print("Figure 1. accuracy histogram")

fig, ax = plt.subplots(1, 2, figsize=(12,4))
sns.distplot(sumary_tab.recall_AA*100, ax=ax[0])
sns.distplot(sumary_tab.recall_peptide*100, ax=ax[1])
plt.savefig(result_dir + "/fig1_acc_histo.pdf")

fig, ax = plt.subplots(1, 2, figsize=(12,4))
sns.distplot(sumary_tab[sumary_tab.file.isin(test_list)].recall_AA*100, ax=ax[0])
sns.distplot(sumary_tab[sumary_tab.file.isin(test_list)].recall_peptide*100, ax=ax[1])
plt.savefig(result_dir + "/fig1_acc_histo_with_validset.pdf")

###############################################################
# figure 2. accuracy via length
###############################################################
print("Figure 2. accuracy vs lengAA")

fig, ax = plt.subplots(1, figsize=(10,8))
sns.barplot(x="len_AA", y="exact_match", data=all_records, ax=ax)
plt.savefig(result_dir + "/fig2_acc_vs_len.pdf")

fig, ax = plt.subplots(1, figsize=(10,8))
sns.barplot(x="len_AA", y="exact_match", data=all_records[all_records.file_id>FLAGS.last_index_for_train_files], ax=ax)
plt.savefig(result_dir + "/fig2_acc_vs_len_with_validset.pdf")
###############################################################
# figure 3. acc_vs_len_by_is_end_with_RK
###############################################################
print("Figure 3. accuracy vs lengAA with RK info")

plt.close('all')
sns.barplot(x="len_AA", y="exact_match", hue="is_end_with_RK", data=all_records)
plt.savefig(result_dir + "/fig3_acc_vs_len_by_is_end_with_RK.pdf")
plt.close('all')
sns.barplot(x="len_AA", y="exact_match", hue="is_start_with_RK", data=all_records)
plt.savefig(result_dir + "/fig3_acc_vs_len_by_is_start_with_RK.pdf")
plt.close('all')
sns.barplot(x="len_AA", y="exact_match", hue="is_RK_inside", data=all_records)
plt.savefig(result_dir + "/fig3_acc_vs_len_by_is_RK_inside.pdf")

plt.close('all')
sns.barplot(x="len_AA", y="exact_match", hue="is_end_with_RK", data=all_records[all_records.file_id>FLAGS.last_index_for_train_files])
plt.savefig(result_dir + "/fig3_acc_vs_len_by_is_end_with_RK_with_validset.pdf")
plt.close('all')
sns.barplot(x="len_AA", y="exact_match", hue="is_start_with_RK", data=all_records[all_records.file_id>FLAGS.last_index_for_train_files])
plt.savefig(result_dir + "/fig3_acc_vs_len_by_is_start_with_RK_with_validset.pdf")
plt.close('all')
sns.barplot(x="len_AA", y="exact_match", hue="is_RK_inside", data=all_records[all_records.file_id>FLAGS.last_index_for_train_files])
plt.savefig(result_dir + "/fig3_acc_vs_len_by_is_RK_inside_with_validset.pdf")

###############################################################
# figure 4. ROC
###############################################################
print("Figure 4. ROC curve")
from sklearn.metrics import roc_curve, auc
plt.close('all')

def draw_roc(all_records, ax, color='darkorange'):
    all_rst = all_records.dropna()[['exact_match','output_score']].copy()

    ## remove -inf
    outscore_min = np.min(all_rst[all_rst.output_score != float("-inf")].output_score)
    row_idx = all_rst[all_rst.output_score == float("-inf")].index
    all_rst.loc[row_idx, ("output_score")] = outscore_min*2
    
    fpr, tpr, _ = roc_curve(list(all_rst.exact_match), list(all_rst.output_score))
    roc_auc = auc(fpr, tpr)

    lw = 2
    ax.plot(fpr, tpr, color=color, lw=lw, label='ROC curve (area = %0.3f)' % roc_auc)
    return

fig, ax = plt.subplots(1)

draw_roc(all_records, ax, color='red')

plt.plot([0, 1], [0, 1], color='navy', lw=1, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")
plt.savefig(result_dir + "/fig4_roc_curve.pdf")

plt.close('all')
fig, ax = plt.subplots(1)

draw_roc(all_records[all_records.file_id>FLAGS.last_index_for_train_files], ax, color='red')

plt.plot([0, 1], [0, 1], color='navy', lw=1, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")
plt.savefig(result_dir + "/fig4_roc_curve_with_validset.pdf")


###############################################################
# figure 5. longest subsequences
###############################################################
print("Figure 5. logest subsequences")
plt.close('all')
# Show each observation with a scatterplot
ax = sns.factorplot(x="longest_match_idx", y="longest_match_length", col="len_AA", kind="box", col_wrap=5,
                    data=all_records[all_records.longest_match_length>0])
plt.savefig(result_dir + "/fig5_longest_match_start_idx_and_length_vs_len.pdf")

ax = sns.factorplot(x="longest_match_idx_with_novor", y="longest_match_length_with_novor", col="len_AA", kind="box", col_wrap=5,
                    data=all_records[all_records.longest_match_length>0])
plt.savefig(result_dir + "/fig5_longest_match_start_idx_and_length_with_novor_vs_len.pdf")


plt.close('all')
# Show each observation with a scatterplot
ax = sns.factorplot(x="longest_match_idx", y="longest_match_length", col="len_AA", kind="box", col_wrap=5,
                    data=all_records[(all_records.file_id>FLAGS.last_index_for_train_files) & (all_records.longest_match_length>0)])
plt.savefig(result_dir + "/fig5_longest_match_start_idx_and_length_vs_len_with_validset.pdf")

ax = sns.factorplot(x="longest_match_idx_with_novor", y="longest_match_length_with_novor", col="len_AA", kind="box", col_wrap=5,
                    data=all_records[(all_records.file_id>FLAGS.last_index_for_train_files) & (all_records.longest_match_length>0)])
plt.savefig(result_dir + "/fig5_longest_match_start_idx_and_length_with_novor_vs_len_with_validset.pdf")


