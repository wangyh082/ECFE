function [cluster_lables] = cluster_dp(dist,K)

%% Input and Output
% INPUT :
% dist : A nCases*nCases matrix. each dist(i, j) represents the distance
%        between the i-th datapoint and j-th datapoint.
% para : options
%        percent : The parameter to determine dc. 1.0 to 2.0 can often yield good performance.
%        method  : alternative ways to compute density. 'gaussian' or
%                  'cut_off'. 'gaussian' often yields better performance.
% OUTPUT :
% cluster_labels : a nCases vector contains the cluster labls. Lable equals to 0 means it's in the halo region
% center_idxs    : a nCluster vector contains the center indexes.

percent = 4;
N = size(dist,1);
position = round(N*(N-1)*percent/100);

tri_u = triu(dist,1);
sda = sort(tri_u(tri_u~=0));
if size(sda,2)<position
position=round(size(sda,2)/2);
end

dc = sda(position);
clear sda; clear tri_u;

%% Compute rho(density)
%fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);
%switch para.method
    % Gaussian kernel
  %  case 'gaussian'
       rho = sum(exp(-(dist./dc).^2),2)-1;
        % "Cut off" kernel
  %  case 'cut_off'
       %  rho = sum((dist-dc)<0, 2);
%end
[~,ordrho]=sort(rho,'descend');

%% Compute delta

delta = zeros(size(rho));
nneigh = zeros(size(rho));

delta(ordrho(1)) = -1;
nneigh(ordrho(1)) = 0;
for i = 2:size(dist,1)
    range = ordrho(1:i-1);
    [delta(ordrho(i)), tmp_idx] = min(dist(ordrho(i),range));
    nneigh(ordrho(i)) = range(tmp_idx); 
end
delta(ordrho(1)) = max(delta(:));

for i=1:size(dist,1)
    ind(i)=i;
    gamma(i)=rho(i)*delta(i);
end

NCLUST = K;

%% cl 为归属标志数组，cl(i)=j 表示第 i 号数据点归属于第 j 个 cluster
%% 先统一将 cl 初始化为 -1


for i=1:size(dist,1)
    cl(i)=-1;
end

[B, Index] = sort(gamma, 'descend');
%raw_cluster_lables = cluster_lables;
icl = Index(1:K);
cl(Index(1:K)) = 1:K;
for i=1:size(dist,1)
    if (cl(ordrho(i))==-1)
        cl(ordrho(i))=cl(nneigh(ordrho(i)));
    end
end
cluster_lables=cl;
end