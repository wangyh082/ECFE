clear all;
close all;
clc;

MBS = 20;   % The number of base clusterings
t = 20;   % The number of steps of the random walks
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
        p1 = [99 13  7  5  4  3  3  2  3];
        p2 = [ 0  0  0  0  1  2  2  2  2];
        p1 = p1(M-1);
        p2 = p2(M-1);
        [N,W] = F_weight(p1,p2,M);
        
        W(W==0) = 0.000001;
        T = 2;
        B = zeros(N);
        for i = 1 : N
            for j = i : N
                B(i,j) = norm(W(i,:)-W(j,:));
                B(j,i) = B(i,j);
            end
        end
        [~,B] = sort(B,2);
        B = B(:,1:T);
        
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
      
        Z = min(FunctionValue);
        Generations = 9;
        for runIdx = 1 : Generations
            for i = 1 : N
                P = 1:N;

                k = randperm(length(P));
                
                PopBaseCls1 = Population((1+number*((k(1))-1)):(number*(k(1))),:);
                PopBaseCls2 = Population((1+number*((k(2))-1)):(number*(k(2))),:);
             
                OffspringPop = discreteCrossover(PopBaseCls1, PopBaseCls2);
                OffspringPop = mutate(OffspringPop, lowerK, upperK, K);

                NewOffspringPop = ones(number,MBS);
             
                for j = 1:MBS
                    NewOffspringPop(:,j) = relabel(OffspringPop(:,j));
                    while  numel(unique(NewOffspringPop(:,j))) == 1
                        randBaseClsNew = randi([1,number],1,number);
                        NewOffspringPop(:,j) = convert_labels(randBaseClsNew)';
                    end
                end

                dnumber = size(unique(NewOffspringPop','rows'),1);
                orinumber = size(NewOffspringPop,2);

                %%%重复的个数超过1个%%%
                minus = orinumber - dnumber;
                if (minus >= 1)
                    newOffspring = randi([1,K],number,minus);
                    NewOffspringPop = [unique(NewOffspringPop','rows')',newOffspring];
                end
             
                OffFunValue = zeros(1,M);
                [OffFunValue(1:4),Off_label] = FitnessCalculate(fea, NewOffspringPop, t, K, Para(i,:));
               
                Z = min(Z,OffFunValue);
                
                for j = 1 : T
                    g_old = max(abs(FunctionValue(B(i,j),:)-Z).*W(B(i,j),:));
                    g_new = max(abs(OffFunValue-Z).*W(B(i,j),:));
                    if g_new < g_old
                        Population((1+number*(B(i,j)-1)):(number*B(i,j)),:) = OffspringPop;
                        FunctionValue(B(i,j),:) = OffFunValue;
                        POP_label(:,B(i,j))= Off_label;
                        Para(B(i,j),:) = Para(i,:);
                    else
                        Para(B(i,j),:)=[ randi([1,2],1,1)' randi([1,3],1,1)'];
                    end
                end                
            end
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