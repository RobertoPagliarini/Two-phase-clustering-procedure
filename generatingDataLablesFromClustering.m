function [Labels] = generatingDataLablesFromClustering(Clusters,SampleID)
%
%This functions generates the final labels of clusters.

Labels = [];

[n,m] = size(Clusters);

for i = 1:n

    for j = 1:m

        for z = 1:length(SampleID)

            if strcmp(SampleID{z}, Clusters{i,j})

                Labels(z,1) = j;

            end

        end

    end

end



