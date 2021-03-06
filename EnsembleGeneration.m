%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you find the code useful for your research, please cite the paper    %
% below:                                                                  %
%                                                                         %
% D. Huang, C.-D. Wang, H. Peng, J. Lai, & C.-K. Kwoh. "Enhanced Ensemble %
% Clustering via Fast Propagation of Cluster-wise Similarities."To appear %
% in IEEE Transactions on Systems, Man, and Cybernetics: Systems.         %
% DOI: 10.1109/TSMC.2018.2876202                                          %
%                                                                         %
% The code has been tested in Matlab R2016a and Matlab R2016b.            %
%                                                                         %
% www.researchgate.net/publication/328581758                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function baseCls = EnsembleGeneration(fea, M, lowerK, upperK)
% This function generates M base clusterings by k-means.
% For each base clustering by k-means, its cluster number is randomly selected in [lowerK, upperK].
% Dong Huang. Sep. 28, 2018.

N = size(fea,1);

clsNums = randi([lowerK, upperK],M,1);

baseCls = zeros(N,M);
parfor i = 1:M
% for i = 1: M % This is a paralleled version
    baseCls(:,i) = kmeans(fea, clsNums(i),'emptyaction','singleton');
    while length(unique(baseCls(:,i))) ~= clsNums(i);
        baseCls(:,i) = kmeans(fea, clsNums(i),'emptyaction','singleton');
    end
end
