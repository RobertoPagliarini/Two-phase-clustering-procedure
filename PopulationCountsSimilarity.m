function [Clusters,idx,optimalK,Data,K] = PopulationCountsSimilarity(ClustersIndexes,ClusterSamples,index,CountsMatrix, Plot) 
%
%This function cluster the population individuals according to the second
%step of the procedure proposed in the paper.

Clusters = {};

if isempty(ClustersIndexes) == 0 

    Index = ClustersIndexes(find(ClustersIndexes(:,index)),index);

    Data = CountsMatrix(:,Index);

else 

    Data = CountsMatrix;

end

[m,n] = size(Data);

%Computing correlation coefficients with correspondind p-values
[rho,pval] = corrcoef(Data);

rho(pval>=0.05) = 0;

sim = rho;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Store Calinski-Harabasz indexes for step 2
K = [];

cvi = 0;

%Optimal number of clusters for step 2
optimalK = 0;

%Clutering step 2
for i = 2:ceil(n/4)
    
    %Spectral clustering step 1

    labels = spectralcluster(sim,i,'Distance','precomputed');

    %Calinski-Harabasz index computation for clustering step 1
    K(i) = chindex(labels, Data');

    if K(i) > cvi

        cvi = K(i);

        optimalK = i;

        idx = labels;

    end

end
 
if strcmp(Plot,'Y')
    
    %Plotting Calinski-Harabasz indexes for clustering step 1 
    figure
    plot(K,'LineWidth',3)
    title('Calinski-Harabasz indexes step 2')
    xlabel('Number of clusters')
    ylabel('Calinski-Harabasz index')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:optimalK  

    tmp = 0;
    
    for j = 1:n
        
        if idx(j) == i
            
            tmp = tmp + 1;
            
            Clusters{tmp,i} = ClusterSamples{j,index};
            
        end
        
    end

end

end