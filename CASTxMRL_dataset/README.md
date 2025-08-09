a) MiceDataset.mat: the elaboration of ASE dataset in https://figshare.com/collections/Allele-specific_expression_reveals_genetic_drivers_of_tissue_regeneration_in_mice/6025157 that has been in the paper "A two-phase clustering procedure based on Allele Specific Expression" of Roberto Pagliarini et al.

b) ClusteringResultsASEHybridMiceTable4.mat: results of the evaluation of the two-phase procedure in clusterinf samples by wound site and/or allele. We partitioned the 40 samples in 3 classes according to their cell population
(16 immune, 12 endothelial and 12 fibroblast). Then we apply our method to one class at a time. Next, we repeat the same test, applying the procedure to all 40 samples. 

c) ClusteringResultsASEHybridMiceTable4OurElab.mat: results of the application of the two-phase procedure to the dataset obtained by by merging the MRL and CAST samples associated to the same cell population, tissue and individual. Thus, each of the 20 new samples has the details for both the two alleles for each gene. Again, we partition the samples according to their cell population (8 immune, 6 endothelial and 6 fibroblast), this time aiming to cluster them only by their wound site. 
