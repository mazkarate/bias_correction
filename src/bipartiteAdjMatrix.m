function adjMatrix = bipartiteAdjMatrix(workerIDs, firmIDs)
    % Ensure the input vectors are columns
    workerIDs = workerIDs(:);
    firmIDs = firmIDs(:);
    
    % Get unique worker and firm IDs
    [uniqueWorkers, ~, workerIdx] = unique(workerIDs);
    [uniqueFirms, ~, firmIdx] = unique(firmIDs);
    
    % Number of unique workers and firms
    numWorkers = length(uniqueWorkers);
    numFirms = length(uniqueFirms);
    
    %Construct Adj. Matrix of Bipartite Graph
    B=sparse(workerIdx, firmIdx, 1, numWorkers, numFirms);
    G = [ sparse( numWorkers, numWorkers ), B;
        B.', sparse(numFirms, numFirms)];
    
    adjMatrix = G;
end


