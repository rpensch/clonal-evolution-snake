import pandas as pd
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str)
parser.add_argument('--min_cluster_size', type=int)
parser.add_argument('--out', type=str)

args = parser.parse_args()

# Load input
results = pd.read_csv(args.input, sep = '\t')

# Count mutations per cluster
muts_per_cluster = (results
                        .drop_duplicates('mutation_id')
                        .groupby('cluster_id', as_index = False)
                        .agg({'mutation_id':'count'})
                        .rename(columns = {'mutation_id':'mut_count'}))

# Exclude clusters with less than min number of mutations
exclude = muts_per_cluster[muts_per_cluster['mut_count']<args.min_cluster_size]['cluster_id'].tolist()
results_filt = (results[~results['cluster_id'].isin(exclude)]).copy()

# Define the "founding clone" to be the cluster with the highest cellular prevalence (CCF)
# If multiple samples, set the founding clone to be the one with the highest CCF across all samples
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
print(f'Remaining clusters renamed as follows:')
print(f'original_id\tnew_id')
for i in rename_dict.items():
    print(f'{i[0]}\t{i[1]}')
