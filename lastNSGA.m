clear all;
close all;
clc;

MBS = 20;   % The number of base clusterings
t = 20;   % The number of steps of the random walks
%%%%NSGA-II framework%%%%
path= 'testData\';
DD = dir(fullfile(path,'*.mat'));
for dd = 1: length(DD)
    %% Load the data.
    data1 = load(fullfile(path,DD(dd).name));
    newName = DD(dd).name;
    newName(end-3:end) = [];
    gt = data1.gnd;
    fea = data1.fea;
    clear data1
    
    [number, d] = size(fea);
    %%%Removing features with low variance%%%%
    [~, Index] = sort(var(fea,0,1), 'descend');
    if d > 300
        fea = fea(:,Index(1:300));
        d = 300;
    end
    K = numel(unique(gt));
    maxIter = 1;
    s1 = 0;
    s2 = 0;
    temp = zeros(maxIter,2);
    for iter = 1: maxIter
        lowerK = 2;
        upperK = floor(sqrt(number));
        
        M = 4;
        N = 200;
        newBaseCls = zeros(number,MBS);
        Population = zeros(N*number,MBS);
        
        for i = 1:N
            baseCls = EnsembleGeneration(fea, ceil(MBS/2), lowerK, upperK);
            randBaseCls = zeros(number,MBS-ceil(MBS/2));
            for j = 1:(MBS-ceil(MBS/2))
                randBaseClsNew = randi([1,number],1,number);
                randBaseCls(:,j) = convert_labels(randBaseClsNew)';
                while  numel(unique(randBaseCls(:,j))) == 1
                    randBaseClsNew = randi([1,number],1,number);
                    randBaseCls(:,j) = convert_labels(randBaseClsNew)';
                end
            end
            sBaseCls = [baseCls,randBaseCls];          
            for j = 1:MBS
                newBaseCls(:,j) = relabel(sBaseCls(:,j));
            end
            %%%decide how many coloums are there were duplicated%%%
            %%%unique是去除重复行后的数字%%%%
            dnumber = size(unique(newBaseCls','rows'),1);
            orinumber = size(newBaseCls,2);
            minus = orinumber - dnumber;
            if (minus >= 1)
                newClusters = randi([1,K],number,minus);
                newnewBaseCls = [unique(newBaseCls','rows')', newClusters];
                Population((1+number*(i-1)):(number*i),:) = newnewBaseCls;
            else
                Population((1+number*(i-1)):(number*i),:) = newBaseCls;
            end  
        end
        Para = [randi([1,2],1,N)' randi([1,3],1,N)'];
        FunctionValue = zeros(N,M);
        POP_label = zeros(number,N);
        for i=1:N
            [FunctionValue(i,1:4), POP_label(:,i)] = FitnessCalculate(fea, Population((1+number*(i-1)):(number*i),:), t, K, Para(i,:));
        end
        FrontValue = P_sort(FunctionValue);
        CrowdDistance = F_distance(FunctionValue,FrontValue);
        Generations = 5;
        mnumber = 0;
        for runIdx = 1 : Generations
            MatingPool = F_mating(Population, FrontValue, CrowdDistance, number);
            for i = 1 : N              
                P = 1:N;               
                k = randperm(length(P));
                %%%均匀交叉 cr = 1%%%
                if i == N
                    PopBaseCls1 = MatingPool((1+number*((1)-1)):(number*(1)),:);
                    PopBaseCls2 = MatingPool((1+number*((N)-1)):(number*(N)),:);
                else
                    PopBaseCls1 = MatingPool((1+number*((i)-1)):(number*(i)),:);
                    PopBaseCls2 = MatingPool((1+number*((i+1)-1)):(number*(i+1)),:);
                end
                
                OffspringPop = discreteCrossover(PopBaseCls1, PopBaseCls2);
                OffspringPop = mutate(OffspringPop, lowerK, upperK, K);

                NewOffspringPop = ones(number, MBS);
                
                for j = 1:MBS
                    NewOffspringPop(:,j) = relabel(OffspringPop(:,j));
                    while  numel(unique(NewOffspringPop(:,j))) == 1
                        randBaseClsNew = randi([1,number],1,number);
                        NewOffspringPop(:,j) = convert_labels(randBaseClsNew)';
                    end
                end
                
                dnumber = size(unique(NewOffspringPop','rows'),1);
                orinumber = size(NewOffspringPop,2); 
                minus = orinumber - dnumber;
                if (minus >= 1)
                    newOffspring = randi([1,K],number,minus);
                    NewOffspringPop = [unique(NewOffspringPop','rows')',newOffspring];
                end
                
                OffFunValue = zeros(1,M);
                [OffFunValue(1:4),Off_label] = FitnessCalculate(fea, NewOffspringPop, t, K, Para(i,:));
                
                Population = [Population;NewOffspringPop];
                FunctionValue = [FunctionValue;OffFunValue];
                Para = [Para;Para(i,:)];
                POP_label = [POP_label,Off_label];
            end

            [FrontValue,MaxFront] = P_sort(FunctionValue,'half');
            CrowdDistance = F_distance(FunctionValue,FrontValue);
            NewPopulation = Population;
           
            %选出非支配的个体
            Next = zeros(1,N);
            NoN = numel(FrontValue,FrontValue<MaxFront);
            Next(1:NoN) = find(FrontValue<MaxFront);
            
            %选出最后一个面的个体
            Last = find(FrontValue==MaxFront);
            [~,Rank] = sort(CrowdDistance(Last),'descend');
            Next(NoN+1:N) = Last(Rank(1:N-NoN));
            %下一代种群
            for i = 1: N
                Population((1+number*(i-1)):(number*i),:) = NewPopulation((1+number*(Next(i)-1)):(number*Next(i)),:);
            end
            % clear NewPopulation
            Population((1+number*N):(number*2*N),:) = [];
            
            POP_label = POP_label(:,Next);
            FunctionValue = FunctionValue(Next,:);     
            FrontValue = FrontValue(Next);
            Para = Para(Next,:);
            CrowdDistance = CrowdDistance(Next);
        end
        nmiScores = zeros(1,N);
        ariScores = zeros(1,N);
        for i = 1:N
            nmiScores(i) = Cal_NMI(POP_label(:,i),gt);
            ariScores(i) = RandIndex(POP_label(:,i),gt);
        end
        [mnmiScores, Index] = max(nmiScores);
        nmax = nmiScores(Index);
        rmax = ariScores(Index);
        temp(iter,:) = [nmax rmax];
        s1 = s1 + nmax;
        s2 = s2 + rmax;
    end
    finalresult(dd,1) = s1/maxIter;
    finalresult(dd,2) = s2/maxIter;
    finalresult
end
