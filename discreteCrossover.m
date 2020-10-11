function  OffspringPop = discreteCrossover(PopBest1, PopBest2)
%%%uniform crossover%%%%
        [m,n] = size(PopBest1);
        OffspringPop = ones(m,n);
        for i = 1: n
            Offspring = PopBest1(:,i);
            Best1 = PopBest1(:,i);
            Best2 = PopBest2(:,i);
            nn = size(Best1,1);
            mask = randi([0,1],nn,1);
            temp1 = find(mask == 0);
            temp2 = find(mask == 1);
            if(numel(temp1))
                Offspring(temp1) = Best1(temp1);
            end
            if(numel(temp2))
                Offspring(temp2) = Best2(temp2);
            end
            OffspringPop(:,i) = Offspring;
        end
end