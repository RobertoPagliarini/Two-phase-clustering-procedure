function [Sim, Dist] = ComputingIndividualDistances(AllelesMatrix,Distance)
%
%This function computes pairwise individual distances starting from alleles matrix
%
%Input:
%
% AllelesMatrix of dimension a*m, where ﻿a is the total number of alleles observed in the RNA-seq experiment,
% while m is the number of individuals in the population. 
%
% Distance: distance metric, specified as a character vector
%
%Output:
%
% Sim: similarity matrix.
%
% Dist: distance matrix.

%Computing size of input matrix
[n,m] = size(AllelesMatrix);                    

Dist = squareform(pdist(AllelesMatrix',Distance));

%Compute distances
Sim = ones(m,m) - Dist;

end