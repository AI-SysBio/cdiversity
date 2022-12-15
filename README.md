# Cdiversity: Quantifying B Cell Clonal Diversity In Repertoire Data
A Python library for computing clonal relationship and diversity in B cell repertoire data

<img align="right" src="https://raw.githubusercontent.com/Aurelien-Pelissier/cdiversity/master/Figures/dprofile.png" width=400>
Advances in high-throughput sequencing technologies have enabled the high-throughput characterization of B cell receptor sequencing data, however, the accurate identification of clonally related BCR sequences remains a major challenge. In this study, we compare three different clone identification methods on both simulated and experimental data, and investigate their impact on the characterization of B cell diversity. We find that different methods may lead to different clonal definitions, which in turn can affect the quantification of clonal diversity in repertoire data. Interestingly, we find the Shannon entropy to be overall the most robust diversity index in regards to different clonal identification. [1].

&nbsp;



        
        
### Running the analysis

First, you need to install cdiversity, or you can use the `cdiversity.py` file provided in the repository:

	- pip install cdiversity
	
	
Then, you can run a repertoire analysis simulation with the toy example below. For a more complete overview, you can check out `Examples/Analyze_sample.py`.
Briefly, the analysis start by grouping Bcell into clones, and then use the obtained groups to compute various diversity metrics.

#### (I) Group repertoire into clones

	import REGIR as gil

	#Set the simulation parameters:
	class param:
		Tend = 10		#Length of the simulation
		N_simulations = 20	#The simulation results should be averaged over many trials
		unit = 'h'		#Unit of time (used for plotting only)
		timepoints = 100	#Number of timepoints to record (used for plotting only)

	r1 = 1
	r2 = 4
	r3 = 0.03
	alpha1 = 20
	alpha2 = 5
      


#### (II) Compute diversity indices



## References

[1] Pélissier, A, Luo, S, et al. "Quantifying B Cell Clonal Diversity In Repertoire Data". *Submitted to Frontier in immunology* (2022) [[Preprint]](https://www.biorxiv.org/content/10.1101/2022.12.12.520133)
