import pandas as pd

def count_mutations_per_cluster(results):
    ''' Count mutations per cluster '''
     
    muts_per_cluster = (results
                            .drop_duplicates('mutation_id')
                            .groupby('cluster_id', as_index = False)
                            .agg({'mutation_id':'count'})
                            .rename(columns = {'mutation_id':'mut_count'}))
    
    return muts_per_cluster

def find_founders(results_filt, muts_per_cluster, min_mutations_founder):
    ''' Find the most likely founder and exclude too small potential founders '''

    # Rank clusters by cellular prevalence to determine the ranking of potential founding clones
    founder_ranking = (results_filt
                        .drop_duplicates(['sample_id','cluster_id'])[['cellular_prevalence','cluster_id']]
                        .groupby('cluster_id', as_index = False).sum()
                        .sort_values('cellular_prevalence', ascending = False)['cluster_id']
                        .tolist())

    # Exclude potential founders that are smaller than min founder size
    small_founders_to_exclude = []
    for i, c in enumerate(founder_ranking):
        if muts_per_cluster[muts_per_cluster['cluster_id'] == c]['mut_count'].values[0] >= min_mutations_founder:
            founder = c
            break
        else:
            small_founders_to_exclude.append(c)
        
    return founder, small_founders_to_exclude

def adjust_ccf(cluster, ccf, founder, founder_ccf):
    ''' Make sure all clusters in all samples have CCF < CCF(founder) '''

    # No clone can have CCF > CCF(founder)
    if (ccf >= founder_ccf) & (cluster != founder):
        adjusted = ccf
        adj_perc = 0.99
        while adjusted >= founder_ccf:
            adjusted = ccf * adj_perc
            adj_perc -= 0.01
        return adjusted
    else:
        return ccf
