function adjMatrix = tripartiteAdjMatrix(clusterIDs, workerIDs, firmIDs)
    
    % Get unique worker, cluster, and firm IDs
    [uniqueClusters, ~, clusterIdx] = unique(clusterIDs);
    [uniqueWorkers, ~, workerIdx] = unique(workerIDs);
    [uniqueFirms, ~, firmIdx] = unique(firmIDs);
    
    % Number of unique workers, clusters, and firms
    numWorkers = length(uniqueWorkers);
    numClusters = length(uniqueClusters);
    numFirms = length(uniqueFirms);
    
    % Construct bipartite adjacency matrices
    B1 = sparse(clusterIdx, workerIdx, 1, numClusters, numWorkers);
    B2 = sparse(clusterIdx, firmIdx, 1, numClusters, numFirms);
    
    % Construct the full adjacency matrix for the tripartite graph
    G = [ sparse(numClusters, numClusters), B1, B2;
          B1.', sparse(numWorkers, numWorkers + numFirms);
          B2.', sparse(numFirms, numWorkers + numFirms)];
    
    adjMatrix = G;
end
