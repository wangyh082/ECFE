function V = cleval(X, CR, ECA)

    [db, dunn] = valid_DbDunn(X,CR); %compute Davies-Bouldin index and Dunn index
   
    [s1, s2] = minimax(CR, ECA, 3);

    V = [s1, -s2, db, -dunn];

end