function newLabel = relabel(labels)
C = unique(labels,'stable');
k = length(C); % the number of clusters
newLabel = zeros(size(labels,1),1);
number = 1;
for i = 1:k
    ind = find(labels == C(i));
    newLabel(ind) = number;
    number = number + 1;
end