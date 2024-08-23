import pandas as pd
import argparse
from filter_fun import count_mutations_per_cluster, adjust_ccf, find_founders

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str)
parser.add_argument('--min_cluster_size', type=int)
parser.add_argument('--min_founder_size', type=float)
parser.add_argument('--out', type=str)

args = parser.parse_args()

# Load input
results = pd.read_csv(args.input, sep = '\t')

# Count mutations per cluster
muts_per_cluster = count_mutations_per_cluster(results)

# Exclude clusters with less than min number of mutations
exclude = muts_per_cluster[muts_per_cluster['mut_count']<args.min_cluster_size]['cluster_id'].tolist()
results_filt = (results[~results['cluster_id'].isin(exclude)]).copy()

# Define the "founding clone" to be the cluster with the highest cellular prevalence (CCF) 
# and more than min_founder_size mutations
# If multiple samples, set the founding clone to be the one with the highest CCF across all samples

# How many total mutations in this indiviual, how many mutations should at least be in the founder?
n_mutations = results.drop_duplicates('mutation_id').shape[0] # from total unfiltered
min_mutations_founder = n_mutations * args.min_founder_size

# Find the most likely founder and exclude too small potential founders
founder, small_founders_to_exclude = find_founders(results_filt, muts_per_cluster, min_mutations_founder)
results_filt = results_filt[~results_filt['cluster_id'].isin(small_founders_to_exclude)].copy()

# Make sure all clusters in all samples have CCF < CCF(founder)
# 1. Get the CCF of the founding clone in each sample
founder_ccf_df = (results_filt[results_filt['cluster_id']==founder]
                    .drop_duplicates(['sample_id','cluster_id'])
                    [['sample_id','cellular_prevalence']])
founder_ccf_dict = dict(zip(founder_ccf_df['sample_id'],
                            founder_ccf_df['cellular_prevalence']))

# 2. Adjust CCF downwards until it fits
results_filt['cellular_prevalence'] = results_filt.apply(lambda x: 
                                                                    adjust_ccf(x['cluster_id'],
                                                                                 x['cellular_prevalence'],
                                                                                 founder,
                                                                                 founder_ccf_dict[x['sample_id']]), axis = 1)
                                                                                       
# Rename the clusters so that the founding clone is always the one called '1' (necessary for ClonEvol)
rename_df = (results_filt.drop_duplicates(['sample_id','cluster_id'])
                .groupby('cluster_id', as_index = False)
                .agg({'cellular_prevalence':'sum'})
                .sort_values('cellular_prevalence', ascending = False))

rename_dict = dict(zip(rename_df['cluster_id'].tolist(), 
                       [str(i) for i in range(1,len(rename_df['cluster_id'].tolist())+1)]))

results_filt['cluster_id'] = results_filt['cluster_id'].apply(lambda x: rename_dict[x])

# Save filtered and renamed
(results_filt.sort_values(['sample_id','cluster_id']).to_csv(args.out, sep = '\t', index = False))

# Print summary to stdout / log
print('Summary:')
print(f'Excluded clusters (original ids): {", ".join([str(cluster) for cluster in exclude])}')
print(f'Excluded founders (original ids): {", ".join([str(cluster) for cluster in small_founders_to_exclude])}')
print(f'Remaining clusters renamed as follows:')
print(f'original_id\tnew_id')
for i in rename_dict.items():
    print(f'{i[0]}\t{i[1]}')
