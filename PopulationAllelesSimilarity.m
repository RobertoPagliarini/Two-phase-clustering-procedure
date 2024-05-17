function [ClusteringResults] = PopulationAllelesSimilarity(SimMatrix,k,SampleID,CountsMatrix) 
%
%This function cluster the population individuals according to the first
%step of the procedure proposed in the paper.

[n,m] = size(SimMatrix);

%Spectral clustering
idx = spectralcluster(SimMatrix,k,'Distance','precomputed');

ClustersSamples = {};

ClustersIndexes = [];

for i = 1:k  

    tmp = 0;

    for j = 1:m
        
        if idx(j) == i
            
            tmp = tmp + 1;
            
            ClustersSamples{tmp,i} = SampleID{j};
            
            ClustersIndexes(tmp,i) = j;

        end
        
    end

end

ClusteringResults = struct('SamplesClustering',{ClustersSamples},'IndexesClustering',ClustersIndexes,'idxSpectral',idx);

