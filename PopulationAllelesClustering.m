function [Results] = PopulationAllelesClustering(AllelesMatrix,CountsMatrix,SampleID, Plot) 
%
%﻿Two-phase clustering procedure based on Allele Specific Expression
%
%This function clusters the population into groups of similar individuals, 
%according to their allelic expression as desscribed in [1].
%
%Input:
%
% AllelesMatrix of dimension a*m, where ﻿a is the total number of alleles observed in the RNA-seq experiment,
% while m is the number of individuals in the population. 
%
% CountsMatrix of dimension a*m. The ﻿columns contain read counts
% associated to individuals.
%
% SampleID: cell vector of dimension m, which contains the IDs of the
% individuals of the populatios.
%
% Plot: a string for plot. Plot the Calinski-Harabasz criterion values 
% for each number of clusters tested and for every procedure step. 
%
%Output:
%
% Results, a structure composed of four fields:
%
%   ClusterLabels: numeric column vectors of dimension n representing the
%   cluster indices.
%
%   AfterStep1: clustering results after step 1 of the procedure.
%
%   AfterStep2: clustering results after step 2 of the procedure.
%
%   NumSubclustersStep1: number of subclusters related to each cluster
%   obtained after step1.
%    
%References:
%
%[1] ﻿R. Pagliarini, F. Nascimben, and A. Policriti. "﻿A two-phase clustering procedure based on Allele Specific Expression" 


%Closing all the open figures
close all

%Compute similarity/distance matrices   
[SimMatrix, DistMatrix] = ComputingIndividualDistances(AllelesMatrix,'Jaccard');

%Computing size of similarity matrix
[nS,mS] = size(SimMatrix);

%Store Calinski-Harabasz indexes for step 1
K = [];

%Cluster validation index
cvi = 0;

%Optimal number of clusters for step 1
optimalK = 0;

%Clutering step 1
for i = 2:ceil(nS/4)
    
    %Spectral clustering step 1
    [ClusteringResults] = PopulationAllelesSimilarity(SimMatrix,i,SampleID,CountsMatrix);
    
    labels = ClusteringResults.idxSpectral;

    %Calinski-Harabasz index computation for clustering step 1
    K(i) = chindex(labels, AllelesMatrix');

    if K(i) > cvi

        cvi = K(i);

        optimalK = i;

        ClusteringResultsStep1 = ClusteringResults;

    end

end

if nargin == 4 & strcmp(Plot,'Y')

    %Plotting Calinski-Harabasz indexes for clustering step 1 
    figure
    plot(K,'LineWidth',3)
    title('Calinski-Harabasz indexes step 1')
    xlabel('Number of clusters')
    ylabel('Calinski-Harabasz index')

else

    Plot = 'N';

end

%Setting data for clustering step 2
[n,m] = size(ClusteringResultsStep1.IndexesClustering);

ClustersIndexes = ClusteringResultsStep1.IndexesClustering;

ClusterSamples = ClusteringResultsStep1.SamplesClustering;

FinalClustering = {};

FinalClusteringStep2 = {};

NumClustersStep2 = zeros(m,1);

%Clustering step 2 in order to obtain sub-clusters
for index = 1:m

    %Clustering based on read counts
    A = PopulationCountsSimilarity(ClustersIndexes,ClusterSamples,index,CountsMatrix,Plot);

    FinalClustering{index,1} = A;

    [nU,mU] = size(FinalClusteringStep2);

    [nA,mA] = size(A);

    NumClustersStep2(index,1) = mA;

    for i = 1:nA

        for j = 1:mA
 
            FinalClusteringStep2{i,mU+j} = A{i,j};

        end

    end

end

%Generating data labels for final clustering
[Labels] = generatingDataLablesFromClustering(FinalClusteringStep2,SampleID);

%Generating output structure
Results = struct('AfterStep2',{FinalClusteringStep2},'AfterStep1',{ClusteringResultsStep1},'NumSubclustersStep1',NumClustersStep2, 'ClusterLabels',Labels);