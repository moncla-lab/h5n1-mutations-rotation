import pandas as pd
from scipy.stats import fisher_exact

## specify directory, files, and gene
gene = 'PB2' # gene to analyze
tsv_dir = f'/Users/jort/coding/h5n1-mutations-rotation/build3-end-aa-analysis/build3-rerun3/{gene.lower()}/' # directory containing tsv files from get-all-aa-counts.py
host1_file = 'human_all_aa_counts.tsv' # host1 (enriched in host) tsv file
host2_file = 'avian_all_aa_counts.tsv' # host2 (background host) tsv file
output_file = 'all_aa_or_pv.tsv' # output for dataframe with odds ratios and pvalues (.tsv)


##### user input above #####


def cal_enr(host1_df_col, host2_df_col):
    '''calculate the odds ratio and pvalue for each amino acid at a given position'''

    ## get list of AAs with nonzero counts in each host
    host1_df_col_nonzero_aas = host1_df_col.loc[(host1_df_col != 0)].index.values.tolist()
    host2_df_col_nonzero_aas = host2_df_col.loc[(host2_df_col != 0)].index.values.tolist()

    ## get list of unique nonzero AAs, excluding 'other'
    nonzero_aas = [x for x in set(host1_df_col_nonzero_aas + host2_df_col_nonzero_aas) if x != 'other']

    pos_scores = {}

    ## for each nonzero AA, calculate the odds ratio and pvalue using a Fisher's exact test and assign to dict
    for aa in nonzero_aas:
        presence_host1 = host1_df_col.loc[aa]
        absence_host1 = sum(host1_df_col.loc[~host1_df_col.index.isin([aa])])
        presence_host2 = host2_df_col.loc[aa]
        absence_host2 = sum(host2_df_col.loc[~host2_df_col.index.isin([aa])])
        if presence_host2 == 0:
            presence_host2 = 1
        if absence_host1 == 0:
            absence_host1 = 1

        if presence_host1 == 0:
            presence_host1 = 1
        if absence_host2 == 0:
            absence_host2 = 1
        
        oddsr, p = fisher_exact([[presence_host1, absence_host1],[presence_host2, absence_host2]], alternative='two-sided')
        pos_scores[aa] = (oddsr, p)

    return pos_scores


## load tsv files containing AA counts
host1_df_path = tsv_dir + host1_file
host2_df_path = tsv_dir + host2_file
host1_df = pd.read_csv(host1_df_path, sep='\t', header=0, index_col=0)
host2_df = pd.read_csv(host2_df_path, sep='\t', header=0, index_col=0)


## get list of all positions based on gene length and make list of all AAs
if gene == 'HA':
    aalength = 568 # 568 AAs in HA
elif gene == 'PB2':
    aalength = 759 # 759 AAs in PB2
all_positions = range(1, aalength+1)
all_aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


## create dictionary to store data from each position
all_pos_scores = {}

## get single columns from dataframes corresponding to each position,
## then get a list of odds ratios and pvalues for each AA at that position
## and add this list to dictionary with the position as its key
for pos in all_positions:
    host1_col = host1_df.loc[:,str(pos)]
    host2_col = host2_df.loc[:,str(pos)]
    all_pos_scores[pos] = cal_enr(host1_col, host2_col)


## create lists for position, AA, odds ratio, and pvalue
pos_list = []
aa_list = []
or_list = []
pv_list = []

## append data to lists for each position/AA pair
for pos in all_pos_scores:
  for aa in all_pos_scores[pos].keys():
    pos_list.append(pos)
    aa_list.append(aa)
    or_list.append(all_pos_scores[pos][aa][0])
    pv_list.append(all_pos_scores[pos][aa][1])


## generate dataframe from data and save as a tsv file
output_df = pd.DataFrame({'position': pos_list, 'aminoacid': aa_list, 'oddsratio': or_list, 'pvalue': pv_list})
output_path = tsv_dir + output_file
output_df.to_csv(output_path, sep='\t', header=True, index=False)