import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def training_trend(ouput_file="ftraining_acc.png"):
    print('loading all training trends ...')
    # sns.set(style="ticks")

    # output_folder = ''
    output_folder = 'mgf_list_v3/'

    train_log_0 = pd.read_csv('/Volumes/MSSHARE/Joonyong/DeepNovoRun/log_file_caption_2dir.tab', sep='\t')
    train_log_17 = pd.read_csv('/Volumes/MSSHARE/Joonyong/DeepNovoRun_0-17/deepnovo.train.model_0-17/log_file_caption_2dir.tab', sep='\t')
    train_log_57 = pd.read_csv('/Volumes/MSSHARE/Joonyong/DeepNovoRun_0-57/deepnovo.train.model_0-57/log_file_caption_2dir.tab', sep='\t')
    train_log_97 = pd.read_csv('/Volumes/MSSHARE/Joonyong/DeepNovoRun_0-97/deepnovo.train.model_0-97/log_file_caption_2dir.tab', sep='\t')
    train_log_157 = pd.read_csv('/Volumes/MSSHARE/Joonyong/DeepNovoRun_0-158/deepnovo.train.model_0-158/log_file_caption_2dir.tab', sep='\t')
    # train_log_232 = pd.read_csv('/Volumes/MSSHARE/Joonyong/DeepNovoRun_232_27/deepnovo.train.model_232_27/log_file_caption_2dir.tab', sep='\t')
    train_log_232 = pd.read_csv('/Volumes/MSSHARE/Joonyong/DeepNovoRun_235_v3_27/log_file_caption_2dir.tab', sep='\t')
    train_log_232t = pd.read_csv('/Volumes/MSSHARE/Joonyong/PnnlRun3_235_v4_0001_ep_30/log_file_caption_2dir.tab', sep='\t')
    # train_logs = [train_log_0, train_log_17, train_log_57, train_log_97, train_log_157, train_log_232, train_log_232t]
    # train_models = ["Pre-trained", "Train-300k", "Train-1m", "Train-2m", "Train-3m", "Train-4m", "Train-4m+tuned"]
    train_logs = [train_log_17, train_log_57, train_log_97, train_log_157, train_log_232, train_log_232t]
    train_models = ["Train-300k", "Train-1m", "Train-2m", "Train-3m", "Train-4m", "Train-4m+tuned"]
    # train_logs = [train_log_17, train_log_97, train_log_232, train_log_232t]
    # train_models = ["Train-300k", "Train-2m", "Train-4m", "Train-4m+tuned"]

    for i, tl in enumerate(train_logs):
        tl['model'] = train_models[i]
    train_logs_all = pd.concat(train_logs)
    threshold = 20

    df = pd.melt(train_logs_all, id_vars=['epoch', 'model'], value_vars=["perplexity_train", "perplexity_valid_1"], var_name='LossType', value_name='Loss')
    g = sns.FacetGrid(df[df.epoch < threshold], col="model", hue="LossType", col_wrap=3, size=4, palette="Paired")
    g = g.map(plt.semilogy, "epoch", "Loss")
    plt.legend(["Training","Validation for len in (10,20]"],frameon=False)
    plt.savefig(output_folder+'fig1a.'+ouput_file, dpi=600)

    # df = pd.melt(train_logs_all, id_vars=['epoch', 'model'], value_vars=["perplexity_train", "perplexity_valid_0", "perplexity_valid_1", "perplexity_valid_2"], var_name='LossType', value_name='Loss')
    g = sns.FacetGrid(df[(df.epoch < threshold)&(df.LossType=='perplexity_valid_1')], hue="model", size=4, palette="Paired")
    g = g.map(plt.semilogy, "epoch", "Loss")
    plt.legend(frameon=False)
    plt.savefig(output_folder+'fig1b1.'+ouput_file, dpi=600)

    df = pd.melt(train_logs_all, id_vars=['epoch', 'model'], value_vars=["perplexity_train", "perplexity_valid_0", "perplexity_valid_1", "perplexity_valid_2"], var_name='LossType', value_name='Loss')
    g = sns.FacetGrid(df[(df.epoch < threshold)], hue="model", col="LossType", size=4, palette="Paired")
    g = g.map(plt.semilogy, "epoch", "Loss")
    plt.legend(frameon=False)
    plt.savefig(output_folder+'fig1b2.'+ouput_file, dpi=600)

    df = pd.melt(train_logs_all, id_vars=['epoch', 'model'], value_vars=["perplexity_valid_0", "perplexity_valid_1", "perplexity_valid_2"], var_name='LossType', value_name='Loss')
    g = sns.FacetGrid(df[(df.epoch < threshold)&(df.model=='Train-4m+tuned')], hue="LossType", size=5, palette="Paired")
    g = g.map(plt.semilogy, "epoch", "Loss")
    plt.legend(["Validation for len in (0,10]","Validation for len in (10,20]","Validation for len in (20,30]"],frameon=False)
    plt.savefig(output_folder+'fig1c.'+ouput_file, dpi=600)

def plot_longest_correct_subseq(target_dfs, target_labels, output_file, target_length=14, scoring='novor'):
    plt.close('all')
    fig, axs = plt.subplots(2, len(target_dfs), sharey='row')
    
    for i, target_df in enumerate(target_dfs):
        
        df = target_df[(target_df.len_AA==target_length)]
        print(df.shape)
        ########################################
        ax = axs[1][i]
        if scoring == 'lbyl':
            ## letter by letter scoring
            sns.barplot(x="longest_match_idx", y="longest_match_idx", ax=ax,
                        data=df[(df.longest_match_length>0) & (df.longest_match_idx<5)], estimator=lambda x: len(x) / len(df) * 100)
        elif scoring == 'novor':
            ## with Novor scoring
            sns.barplot(x="longest_match_idx_with_novor", y="longest_match_idx_with_novor", ax=ax,
                        data=df[(df.longest_match_length_with_novor>0) & (df.longest_match_idx_with_novor<5)], estimator=lambda x: len(x) / len(df) * 100)
        # axs[i].set(xlim=(0,5))
        ax.set(ylabel="")
        ax.set(xlabel="")
        ########################################
        # longest match
        ax = axs[0][i]
        if scoring == 'lbyl':
            sns.factorplot(x="longest_match_idx", y="longest_match_length", kind="box", ax=ax,
                           data=df[(df.longest_match_length>0) & (df.longest_match_idx<5)])
        elif scoring == 'novor':
            sns.factorplot(x="longest_match_idx_with_novor", y="longest_match_length_with_novor", kind="box", ax=ax,
                           data=df[(df.longest_match_length_with_novor>0) & (df.longest_match_idx_with_novor<5)])
        
        ax.set(xlabel=target_labels[i])
        ax.set(ylabel="")
        ########################################
        
    axs[1][0].set(ylabel="Percent(%)")
    axs[0][0].set(ylabel="Length of correct substrings")

    # fig.suptitle("Longest correct subsequences by starting position (Length={0})".format(target_length))

    fig.tight_layout()
    fig.savefig(output_file.format(scoring, target_length), dpi=600)

def validation_by_datasize():
    #############################################################################################
    # Figure 2a. plot accuracy comparisons by species (mgf files)
    #############################################################################################
    print('Making Figure 2a ...')

    sns.set(style="whitegrid")

    print('Loading summary data files ...')

    mgf_tab = pd.read_csv('/Volumes/MSSHARE/Joonyong/mgf_list_v3.log', sep='\t')

    sumary_0 = pd.read_csv('/Volumes/MSSHARE/Joonyong/DeepNovoRun/mgf_test/log.txt', sep='\t')
    sumary_17 = pd.read_csv('/Volumes/MSSHARE/Joonyong/DeepNovoRun_0-17/mgf_test/log.txt', sep='\t')
    sumary_57 = pd.read_csv('/Volumes/MSSHARE/Joonyong/DeepNovoRun_0-57/mgf_test/log.txt', sep='\t')
    # sumary_97 = pd.read_csv('/Volumes/MSSHARE/Joonyong/DeepNovoRun_0-97/mgf_test/log.txt', sep='\t')
    # sumary_157 = pd.read_csv('/Volumes/MSSHARE/Joonyong/DeepNovoRun_0-158/mgf_test/log.txt', sep='\t')
    sumary_234 = pd.read_csv('/Volumes/MSSHARE/Joonyong/PnnlRun3_235_v4_0001_ep_60/mgf_test/log.txt', sep='\t')
    # sumary_234t = pd.read_csv('/Volumes/MSSHARE/Joonyong/PnnlRun3_235_v4_0001_ep_30_top5/mgf_test/log.txt', sep='\t')

    test_data = []

    test_data.append(mgf_tab.merge(sumary_0, left_on="mgf_file", right_on="file", how="inner"))
    test_data.append(mgf_tab.merge(sumary_17, left_on="mgf_file", right_on="file", how="inner"))
    test_data.append(mgf_tab.merge(sumary_57, left_on="mgf_file", right_on="file", how="inner"))
    # test_data.append(mgf_tab.merge(sumary_97, left_on="mgf_file", right_on="file", how="inner"))
    # test_data.append(mgf_tab.merge(sumary_157, left_on="mgf_file", right_on="file", how="inner"))
    test_data.append(mgf_tab.merge(sumary_234, left_on="mgf_file", right_on="file", how="inner"))
    # test_data.append(mgf_tab.merge(sumary_234t, left_on="mgf_file", right_on="file", how="inner"))

    data_labels = ['DeepNovo', "Kaiko-300k", "Kaiko-1M", "Kaiko-4.5M"]
    # data_labels = ['DeepNovo', "Kaiko-300k", "Kaiko-1M", "Kaiko-4.5M", "Kaiko-4.5M(Top5)"]

    last_idx_for_test_data = [234 for i in data_labels]
    print(last_idx_for_test_data)
    
    n_datasets = len(test_data)
    colors = sns.color_palette("Paired", n_datasets*2)
    fig, ax = plt.subplots(1, 2, sharey=True, figsize=(10,15))
    for i in range(n_datasets):
        j = n_datasets-i-1
        sns.barplot(x="recall_AA", y="species", data=test_data[j], linewidth=0, errwidth=0.5, ax=ax[0], color=colors[j])
        sns.barplot(x="recall_peptide", y="species", data=test_data[j], linewidth=0, errwidth=0.5, ax=ax[1], color=colors[j])

    ax[0].set(xlabel="Recall AA")
    ax[1].set(xlabel="Recall Peptide")
    ax[0].tick_params(axis="y", labelleft=False)
    ax[1].tick_params(axis="y", labelright=True)
    plt.tight_layout()
    plt.savefig("fig2a.acc_histo_by_species.png", dpi=600)
    #############################################################################################
    # Figure 2b. plot accuracy histogram by species (mgf files)
    #############################################################################################
    print('Making Figure 2b ...')

    sns.set(style="white", context="talk")

    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(10,10))
    for i, tdata in enumerate(test_data):
        sns.distplot(tdata[tdata.id>last_idx_for_test_data[i]].recall_AA*100, ax=ax[0], color=colors[i])
        sns.distplot(tdata[tdata.id>last_idx_for_test_data[i]].recall_peptide*100, ax=ax[1], label=data_labels[i], color=colors[i])

    ax[0].set(xlabel="Recall AA(%)")
    ax[1].set(xlabel="Recall Peptide(%)")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig("fig2b.acc_histo_by_datasize.png", dpi=600)

    #############################################################################################
    # Figure 3. plot accuracy by peptide length
    #############################################################################################
    print('Making Figure 3 ...')

    print('Loading pkl files for all results ...')
    
    all_records_list = []
    model_dirs = ['DeepNovoRun', 'DeepNovoRun_0-17', 'DeepNovoRun_0-57']
    
    for model_dir in model_dirs:
        all_records_list.append(pd.read_pickle('/Volumes/MSSHARE/Joonyong/{0}/mgf_test/all_results_df_v3.pkl'.format(model_dir)))

    all_records_list.append(pd.read_pickle('/Volumes/MSSHARE/Joonyong/PnnlRun3_235_v4_0001_ep_60/mgf_test/all_results_df.pkl'))

    # top5
    # top5_all = pd.read_pickle('/Volumes/MSSHARE/Joonyong/PnnlRun3_235_v4_0001_ep_30_top5/mgf_test/top5_all.pkl')
    # top5_all = top5_all[top5_all['rank']<=5].sort_values(['match_score'])
    # top5_all = top5_all.drop_duplicates(subset=['scan'], keep='last')
    # top5_all = top5_all.reset_index(drop=True)
    # top5_all[['file_id', "scan_number"]] = pd.DataFrame(top5_all.scan.str.split(':').tolist(), columns=['file_id', "scan_number"], dtype="int32")
    # all_records_list.append(top5_all)
    
    data_labels = ['DeepNovo', "Kaiko-300k", "Kaiko-1M", "Kaiko-4.5M"]
    # data_labels = ['DeepNovo', "Kaiko-300k", "Kaiko-1M", "Kaiko-4.5M", "Kaiko-4.5M(Top5)"]
    
    plt.close('all')

    sns.set(style="darkgrid")
    fig, ax = plt.subplots(figsize=(12,8))

    n_records = len(all_records_list)
    # stort test records
    test_records = []
    for i, records in enumerate(all_records_list):
        records = records.reset_index(drop=True)
        test_record = records[(records.file_id > last_idx_for_test_data[i])&(records.len_AA<26)]
        test_records.append(test_record)

    for i in range(n_records):
        j = n_records-i-1
        record = test_records[j]
        print(data_labels[j], record.shape)
        print(record.groupby('len_AA').exact_match.mean())
        sns.barplot(x="len_AA", y="exact_match", data=record, ax=ax, linewidth=0, ci=None, color=colors[j],label=data_labels[j])

    ax.set(xlabel="AA length")
    ax.set(ylabel="Exact match")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig("fig3.acc_by_len.png", dpi=600)

    #############################################################################################
    # Figure 4. plot longest correct subsring by starting position
    #############################################################################################
    # print('Making Figure 4 ...')

    # plt.close('all')

    # target_dfs = [test_records[0], test_records[-2], test_records[-1]]
    # target_labels = [data_labels[0], data_labels[-2], data_labels[-1]]

    # for leng in [9, 12, 14, 15, 17]:
    #     plot_longest_correct_subseq(target_dfs, target_labels, "fig4.longest_correct_subseq_by_pos_{0}_{1}.png", target_length=leng, scoring='novor')
    #     plot_longest_correct_subseq(target_dfs, target_labels, "fig4.longest_correct_subseq_by_pos_{0}_{1}.png", target_length=leng, scoring='lbyl')


def longest_correct_substring_by_top5(top5_df, target_length=14, scoring='novor'):
    target_dfs = []
    target_labels = []
    df = top5_df[top5_df.len_AA == target_length]
    for ranking in range(5):
        # num = (ranking+1)*2
        # if ranking == 0:
        #     num = 1
        num = (ranking+1)
        if scoring=='lbyl':
            sort_tab = df[df['rank']<=num].sort_values(['longest_match_length'])
        else:
            sort_tab = df[df['rank']<=num].sort_values(['longest_match_length_with_novor'])
        sort_tab = sort_tab.drop_duplicates(subset=['scan'], keep='last')
        # print(sort_tab.head())
        target_dfs.append(sort_tab)
        target_labels.append('Top {0}'.format(ranking+1))
    plot_longest_correct_subseq(target_dfs, target_labels, "fig5.longest_correct_subseq_by_rank_{0}_{1}.png", target_length, scoring)

def exact_match_by_top5(top5_df):
    target_dfs = []
    target_labels = []
    df = top5_df
    for ranking in range(5):
        # num = (ranking+1)*2
        # if ranking == 0:
        #     num = 1
        num = (ranking+1)

        sort_tab = df[df['rank']<=num].sort_values(['exact_match'])
        sort_tab = sort_tab.drop_duplicates(subset=['scan'], keep='last')
        print('RANK {0}: num_spectra:{1}'.format(ranking+1, sort_tab.shape[0]))
        target_dfs.append(sort_tab)
        target_labels.append('Top {0}'.format(ranking+1))
    
    sns.set(style="darkgrid")
    colors = sns.color_palette("Paired")
    fig, ax = plt.subplots(figsize=(12,8))
    ntargets = len(target_dfs)
    for i in range(ntargets):
        j = ntargets-i-1
        tarf = target_dfs[j]
        sns.barplot(x="len_AA", y="exact_match", data=tarf, ax=ax, ci=None, color=colors[j],label=target_labels[j])
    plt.legend()
    plt.tight_layout()
    plt.savefig("fig6.acc_by_len_by_rank.png", dpi=600)

if __name__ == '__main__':

    ###############################################################
    # figure 1. training trend
    ###############################################################
    # print("Making Figure 1. training trend")
    # training_trend("training_acc.png")

    ###############################################################
    # figure 2-4. test results by data size
    ###############################################################
    validation_by_datasize()

    # #############################################################################################
    # # Figure 5. 
    # #############################################################################################
    # print('Making Figure 5 ...')
    # top5_all = pd.read_pickle('/Volumes/MSSHARE/Joonyong/PnnlRun3_232_top5/mgf_test/top5_all.pkl')
    # top5_all = top5_all.reset_index(drop=True)
    # top5_all[['file_id', "scan_number"]] = pd.DataFrame(top5_all.scan.str.split(':').tolist(), columns=['file_id', "scan_number"], dtype="int32")

    # longest_correct_substring_by_top5(top5_all[top5_all.file_id > 232], 14, 'lbyl')
    # longest_correct_substring_by_top5(top5_all[top5_all.file_id > 232], 14, 'novor')
    
    # #############################################################################################
    # # Figure 6. 
    # #############################################################################################
    # print('Making Figure 6 ...')
    # exact_match_by_top5(top5_all[top5_all.file_id > 232])
