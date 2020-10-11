function [ss1,ss2] = minimax(C, CM0, k)
%%%s1聚类间指标 越小越好 s2 聚类内指标 越大越好

    CM = PathbasedSimi(CM0, k);
    s = 0;
    ss1 = 0;
    ss2 = 0;
    NumC = length(unique(C));
    for i = 1: NumC
        a = find(C == i);
        if length(a) <= k
            ss1 = Inf;
            ss2 = 0;
            s = Inf;
            return;
        end
        s1 = max(max(CM(a, C ~= i)));
        
        CM0_a = CM0(a, a);
        CM_a = PathbasedSimi(CM0_a, k);
        
        try
            %%%%******nothing with the clustering method to generate n*1 or 1*n vector, but with the clustering algorithm******%%%%
            [C1,~] = kmeans(CM_a,2);   
            flag = 1;
        catch e
            flag = 0;  
        end
   
        if flag == 1
            s2 = min(min(CM_a(C1==1, C1==2)));
            if (size(s2,1) == 0)
                s2 = Inf;
            end
        else
            s2 = Inf;
        end
        s = s + s1 ./ s2;
        ss1 = ss1 + s1;
        ss2 = ss2 + s2;
    end
end

function CM = PathbasedSimi(W, k)
    N = size(W, 1);
    if k <= 10
        S = sumOfNeighbors(W, k);    
        W = W .* repmat(S', 1, N) .* repmat(S, N, 1);    
    end
    W = PathbasedDist(W.^(-1)) + eye(N);
    CM = W.^(-1);
end

function M = sumOfNeighbors(W, numOfNeighbours)
    M = sort(W, 'descend');
    M= sum(M(2:numOfNeighbours+1,:));
    M = M / max(M);
end