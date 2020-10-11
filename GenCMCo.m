%%%Generate Co-association Matrix
function CM = GenCMCo( PI)
    N = size(PI,1);
    CM = zeros(N);
    PiNo  = size(PI, 2);      
    
    for i = 1: PiNo
        C = PI(:,i);
        for j = 1: length(unique(PI(:,i)))
            IDX = find(C==j);
            if length(IDX) <=1 
                continue;
            end                  
            n = length(IDX);       
            Pairs = combntns(1:n,2);
            Ind = [IDX(Pairs(1:n*(n-1)/2)),IDX(Pairs(n*(n-1)/2+1:n*(n-1)))];
            CM((Ind(:,2) - 1)* N + Ind(:,1)) = CM((Ind(:,2) - 1)* N + Ind(:,1)) + 1;   
        end
    end    
%     CM = CM + CM';
    CM = (CM + CM')/ size(PI,2) + eye(N);
end