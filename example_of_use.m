load('DataCh1Leaves.mat')

AllelesMatrix = Ch1.PopulationAlleleIndexesMatrixBinaryFiltered;

CountsMatrix = Ch1.PopulationCountsMatrixFilteredOnlyIndividuals;

SampleID= Ch1.IndividualsData;

[Results] = PopulationAllelesClustering(AllelesMatrix,CountsMatrix,SampleID,'N');


