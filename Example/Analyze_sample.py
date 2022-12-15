import sys, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics.cluster import adjusted_mutual_info_score
from collections import Counter

#from cdiversity import identify_clonal_group
import cdiversity



def main():
    neg_file = 'Data/negation_sequences.csv'
    sample_file = 'Data/sample.csv'
    df_neg = pd.read_csv(neg_file, sep='\t')   
    df_neg.drop_duplicates(subset='junctionH',inplace=True)
    df = pd.read_csv(sample_file, sep='\t') 
    BCR_array = df['VDJ_sequence'].to_numpy()
    BCR_neg_array = df_neg['VDJ_sequence'].to_numpy()
    
    print('Grouping %s BCRs into clones...\n' % len(df))
    
    # 1 - Decompose sample into clones with different methods
    clone_true = df['true_clone'].to_numpy()
    
    print('Baseline Clone identification...')
    clone_baseline, _ = cdiversity.identify_clonal_group(df, method='junction')
    print('Done\n')
    
    print('VJ-junction Clone identification...')
    print('  Computing distance to negation..')
    threshold = cdiversity.compute_negation_threshold(df['junctionH'].to_numpy(),df_neg['junctionH'].to_numpy(), metric = 'levenshtein')
    print('  Clustering sequences..')
    clone_VJJ, _ = cdiversity.identify_clonal_group(df, method='VJJ', clone_threshold = threshold)
    print('Done\n')
    
    print('Alignement free Clone identification...')
    print('  Building tf-idf embedings..')
    vector_array, vector_neg = cdiversity.get_kmer_representations_both(BCR_array, BCR_neg_array)
    print('  Computing distance to negation..')
    threshold = cdiversity.compute_negation_threshold(vector_array, vector_neg, metric = 'cosine')
    print('  Clustering sequences..')
    clone_AF, _ = cdiversity.identify_clonal_group(df, method='AF', embeddings = vector_array, clone_threshold = threshold)
    print('Done\n')
    
    
    
    
    
    # 2 - Check for the consistency of the clonal identification methods with AMI
    print('Similarities between methods')
    c_name = ['G0', 'JO', 'VJJ', 'AF']
    clusters = [clone_true, clone_baseline, clone_VJJ, clone_AF]
    
    cluster_similarity = np.ones((len(c_name), len(c_name)))
    for ci in range(len(c_name)):
        for cj in range(ci+1, len(c_name)):
            cluster_similarity[ci,cj] = AMI_sampled(clusters[ci], clusters[cj])
            cluster_similarity[cj,ci] = cluster_similarity[ci,cj]
    plt.figure(figsize=(8,7))
    plt.imshow(cluster_similarity, cmap='BuGn', vmin = 0.5, vmax=1, origin='lower') 
    plt.xticks(np.arange(len(c_name)),c_name, rotation=90)
    plt.yticks(np.arange(len(c_name)),c_name)
    c = plt.colorbar()
    c.set_label('AMI')
    plt.show()       
    
    
    
    # 3 - Compute the diversity profiles and clone accumulation curves
    print('\n\nDiversity analysis')
    curve_true = cdiversity.cal_accumulation_curve(Counter(clone_true))
    curve_baseline = cdiversity.cal_accumulation_curve(Counter(clone_baseline))
    curve_VJJ = cdiversity.cal_accumulation_curve(Counter(clone_VJJ))
    curve_AF = cdiversity.cal_accumulation_curve(Counter(clone_AF))
    
    plt.figure(figsize=(8,4))
    plt.title('Clone accumulation curve')
    plt.plot(curve_true, lw=2, label = 'G0')
    plt.plot(curve_baseline, lw=2, label = 'JO')
    plt.plot(curve_VJJ, lw=2, label = 'VJJ')
    plt.plot(curve_AF, lw=2, label = 'AF')
    plt.legend(fontsize=16)
    plt.xlabel('Sequences')
    plt.ylabel('Clones')
    plt.xlim(xmin=1)
    plt.yscale('log')
    plt.xscale('log')
    plt.show()
    
    
    
    div_profile_true, alpha_axis = cdiversity.diversity_profile(Counter(clone_true))
    div_profile_baseline, _ = cdiversity.diversity_profile(Counter(clone_baseline))
    div_profile_VJJ, _ = cdiversity.diversity_profile(Counter(clone_VJJ))
    div_profile_AF, _ = cdiversity.diversity_profile(Counter(clone_AF))
    plt.figure(figsize=(8,4))
    plt.title('Diversity profile')
    plt.plot(div_profile_true, lw=2, label = 'G0')
    plt.plot(div_profile_baseline, lw=2, label = 'JO')
    plt.plot(div_profile_VJJ, lw=2, label = 'VJJ')
    plt.plot(div_profile_AF, lw=2, label = 'AF')
    plt.yscale('log')
    plt.xlim(0,len(alpha_axis))
    plt.legend(fontsize=16)
    plt.xticks([0,int(len(alpha_axis)/2),len(alpha_axis)], ['0','1',r'$\infty$'])
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"Hill's diversity $^\alpha D$")
    plt.show()
    
    clone_dict = Counter(clone_true)
    print('True Richness: %.0f' % cdiversity.richness(clone_dict))
    print('True Richness Chao: %.0f' % cdiversity.richness_chao(clone_dict))
    print('True Shannon entropy: %.2f' % cdiversity.Shannon_entropy(clone_dict))
    print('True Shannon entropy Chao: %.2f' % cdiversity.Shannon_entropy_Chao(clone_dict))
    
    
    
    
    
    
def AMI_sampled(cluster1, cluster2, Nmax = 500, nrepeat = 5):
    Ninstance = len(cluster1)
    
    if Ninstance > Nmax:
        AMI = 0
        for i in range(nrepeat):
            new_index = np.random.randint(0,Ninstance, size = Nmax)
            nclust1 = cluster1[new_index]
            nclust2 = cluster2[new_index]
            AMI += adjusted_mutual_info_score(nclust1, nclust2)/nrepeat
    else:
        AMI = adjusted_mutual_info_score(cluster1, cluster2)
          
    return AMI
    
    
    
if __name__ == '__main__':
    plt.rcParams.update({'font.size': 20})
    main()