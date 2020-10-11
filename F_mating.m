function MatingPool = F_mating(Population,FrontValue,CrowdDistance,number)
%交配池选择

    [total,D] = size(Population);
    N = total/number;
    
    %二元联赛选择
    MatingPool = zeros(total,D);
    Rank = randperm(N);
    Pointer = 1;
    for i = 1 : 2 : N
        %选择父母
        k = zeros(1,2);
        for j = 1 : 2
            if Pointer >= N
                Rank = randperm(N);
                Pointer = 1;
            end
            p = Rank(Pointer);
            q = Rank(Pointer+1);
            if FrontValue(p) < FrontValue(q)
                k(j) = p;
            elseif FrontValue(p) > FrontValue(q)
                k(j) = q;
            elseif CrowdDistance(p) > CrowdDistance(q)
                k(j) = p;
            else
                k(j) = q;
            end
            Pointer = Pointer+2;
        end
        MatingPool((1+number*((i)-1)):(number*(i)),:) = Population((1+number*((k(1))-1)):(number*(k(1))),:);
        MatingPool((1+number*((i+1)-1)):(number*(i+1)),:) = Population((1+number*((k(2))-1)):(number*(k(2))),:);
    end
end

