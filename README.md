# Cdiversity: Quantifying B-Cell Clonal Diversity In Repertoire Data

<img align="right" src="https://raw.githubusercontent.com/Aurelien-Pelissier/cdiversity/master/Images/dprofile.png" width=350>
Advances in high-throughput sequencing technologies have enabled the high-throughput characterization of B-cell receptor sequencing data. Still, the accurate identification of clonally related BCR sequences remains a difficult challenge. Importantly, different methods may lead to different clonal definitions, which in turn can affect the quantification of clonal diversity in repertoire data [1]. This library provide different tools and metrics to (i) group B-cell repertoires into clonal groups and (ii) compute diversity indices and diversity profiles from the obtained groups [1].

&nbsp;



        
        
### Running the analysis

First, you need to install cdiversity, or alternatively you can use the `cdiversity.py` file provided in the repository:

	- pip install cdiversity
	
	
Then, you can run a repertoire analysis simulation with the toy example below. For a more complete overview, you can check out `Examples/Analyze_sample.py`.
Briefly, the analysis start by grouping Bcell into clones, and then use the obtained groups to compute various diversity metrics.

#### (I) Grouping repertoire into clones

Available methods for clonal identification are `junction`, which simply group clones together only if they have the same junction. Then, there is the commonly used `VJ-junction` methods, which group together BCR with the same V and J genes, as well as some user-specificed junction similarity (clone_threshold). Finally, the last method is `alignfree`, which compute tf-idf embedings of the BCRs to perform a fast clustering without relying on the V and J germline genes alignements.

	import pandas as pd
	import cdiversity

	df = pd.read_csv('Data/sample.csv', sep='\t') 
	clones_baseline, _ = cdiversity.identify_clonal_group(df, method='junction')
	clone_VJJ, _ = cdiversity.identify_clonal_group(df, method='VJJ', clone_threshold = 0.1)
	
For the alignement free method [2], you need to compute the optimal threshold first using the negation sequence. Something important to keep in mind is that the optimal threshold will be different for each repertoire as it depends on the learned tf-idf embeddings. Please refer to `Examples/Analyze_sample.py` to see how to proceed. Note that, as the AF method needs to compute a distance matrix for all pairwise sequences, it scales as O(N^2) and becomes slow for repertoires with more than 10k sequences.
	


#### (II) Computing diversity indices
Once the clonal groups are obtained, you can compute any diversity indices or the Hill's diversity profile with a single command.
Implemented indices are richness, richness_chao, Shannon_entropy, Shannon_entropy_chao, Simpson_index, dominance, eveness.

	from collections import Counter
	
    	clone_dict = Counter(clone_VJJ)
    	diversity = cdiversity.Shannon_entropy_Chao(clone_dict)
    	div_profile, alpha_axis = cdiversity.diversity_profile(clone_dict)
	
	
<!-- #### (III) Computing Chao diversity indices
Incomplete sample information -> Show acumulation curve, show diversity indices and Chao diversity profile
Link to new reference
Add some plots and show the chao indices, explain that it's only integers -->



## References
#### references for citing cdiversity

[1] Pelissier, A, Luo, S, et al. "Quantifying B Cell Clonal Diversity In Repertoire Data". *Submitted to Frontier in immunology* (2022) [[Preprint]](https://www.biorxiv.org/content/10.1101/2022.12.12.520133)

<!--[2] Pelissier, A et Rodriguez, M "Understanding the Germinal Center Reaction Through Clonal Diversity Indices" (in preparation)-->

#### Implemented clonal identification methods

[2] Lindenbaum, Ofir, et al. "Alignment free identification of clones in B cell receptor repertoires." *Nucleic acids research* 49.4 (2021): e21-e21.

[3] Nouri, Nima, and Steven H. Kleinstein. "Optimized threshold inference for partitioning of clones from high-throughput B cell repertoire sequencing data." *Frontiers in immunology* 9 (2018).
