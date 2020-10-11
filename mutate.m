function newPop = mutate(PopInd, lowerK, upperK, K)
[m,n] = size(PopInd);
pm = 1/m;
for i = 1: n
    if rand <= pm
        position = randi([1,m],1,1);
        neiborN = randi([0,ceil(m/K)],1,1);
        %                if neiborN
        newLabel = randi([lowerK,upperK],1,1);
        while newLabel == PopInd(position,i)
            newLabel = randi([lowerK,upperK],1,1);
        end
        for j = 1:neiborN
            if position + j - 1 <= m
                PopInd((position + j - 1),i) = newLabel;
            else
                PopInd(mod((position + j - 1),m),i) = newLabel;
            end
        end
    end
end
newPop = PopInd;
end