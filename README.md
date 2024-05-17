# Two-phase-clustering-procedure
We developed a novel two-phase unsupervised clustering procedure rooted on the spectral clustering algorithm, with the goal to partition the population into groups of similar individuals, according to their allelic expression.

function [Results] = PopulationAllelesClustering(AllelesMatrix,CountsMatrix,SampleID, Plot) 

﻿Two-phase clustering procedure based on Allele Specific Expression

This function clusters the population into groups of similar individuals, 
according to their allelic expression as desscribed in [1].

Input:

 AllelesMatrix of dimension a*m, where ﻿a is the total number of alleles observed in the RNA-seq experiment,
 while m is the number of individuals in the population. 

 CountsMatrix of dimension a*m. The ﻿columns contain read counts
 associated to individuals.

 SampleID: cell vector of dimension m, which contains the IDs of the
 individuals of the populatios.

 Plot: a string for plot. Plot the Calinski-Harabasz criterion values 
 for each number of clusters tested and for every procedure step. 

Output:

 Results, a structure composed of four fields:

   ClusterLabels: numeric column vectors of dimension n representing the
   cluster indices.

   AfterStep1: clustering results after step 1 of the procedure.

   AfterStep2: clustering results after step 2 of the procedure.

   NumSubclustersStep1: number of subclusters related to each cluster
   obtained after step1.
    
References:

[1] ﻿R. Pagliarini, F. Nascimben, and A. Policriti. "﻿A two-phase clustering procedure based on Allele Specific Expression" 

