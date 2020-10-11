function [fitness, Label] = FitnessCalculate(fea, baseCls, t, K, para)

[bcs, baseClsSegs] = getAllSegs(baseCls);
clsSim = full(simxjac(baseClsSegs));
clsSimRW = computePTS_RRW(clsSim, t);
if para(1) == 1
    ECA = getECA(bcs,clsSimRW);
else
    ECA = GenCMCo(baseCls);
end
n = size(ECA, 1);
for i = 1:n
    lastW(i) = i;
end
lastW = repmat(lastW', 1, n);
NewECA = lastW.* ECA;
if para(2) == 1
    Label = SpectralClusteringSimple(NewECA, K);
elseif para(2) == 2
    dist = pdist2(NewECA,NewECA);
    [cluster_lables] = cluster_dp(dist,K);
    Label = cluster_lables';
else
    Label = litekmeans(NewECA,K);
end
V = cleval(fea, Label, NewECA);
fitness=V(1:4);
end