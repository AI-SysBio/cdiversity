# Cdiversity: Quantifying B Cell Clonal Diversity In Repertoire Data
A Python library for computing clonal relationship and diversity in B cell repertoire data

<img align="right" src="https://raw.githubusercontent.com/Aurelien-Pelissier/cdiversity/master/Images/dprofile.png" width=400>
Advances in high-throughput sequencing technologies have enabled the high-throughput characterization of B cell receptor sequencing data, however, the accurate identification of clonally related BCR sequences remains a major challenge. In this study, we compare three different clone identification methods on both simulated and experimental data, and investigate their impact on the characterization of B cell diversity. We find that different methods may lead to different clonal definitions, which in turn can affect the quantification of clonal diversity in repertoire data. Interestingly, we find the Shannon entropy to be overall the most robust diversity index in regards to different clonal identification. [1].

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
	


#### (II) Computing diversity indices
Once the clonal groups are obtained, you can compute any diversity indices or the Hill's diversity profile with a single command.
Implemented indices are richness, richness_chao, Shannon_entropy, Shannon_entropy_chao, Simpson_index, dominance, eveness.

	from collections import Counter
	
    	clone_dict = Counter(clone_VJJ)
    	diversity = cdiversity.Shannon_entropy_Chao(clone_dict)
    	div_profile, alpha_axis = cdiversity.diversity_profile(clone_dict)



## References

[1] Pélissier, A, Luo, S, et al. "Quantifying B Cell Clonal Diversity In Repertoire Data". *Submitted to Frontier in immunology* (2022) [[Preprint]](https://www.biorxiv.org/content/10.1101/2022.12.12.520133)
