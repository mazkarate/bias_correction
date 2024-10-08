function  firm_worker_cluster = check_type_cluster(id, firmid, cluster)

% Step 1: Sort by cluster
[sortedCluster, sortIdx] = sort(cluster);
sortedFirmid = firmid(sortIdx);

% Step 2: Check if each cluster is specific to one firmid
firmidChange = diff(sortedFirmid) ~= 0;
clusterChange = diff(sortedCluster) ~= 0;

% Ensure that when firmid changes, cluster also changes
isClusterFirmSpecific = all(firmidChange <= clusterChange);

% Step 3: If all clusters are unique to one firmid, we're done
if isClusterFirmSpecific
    disp('Each cluster is unique to one firmid. Use tripartite graph pruning.');
    firm_worker_cluster = true;
else
    % Step 4: Check if each cluster is specific to one id using the same sorting
    sortedId = id(sortIdx);
    idChange = diff(sortedId) ~= 0;
    
    % Ensure that when id changes, cluster also changes
    isClusterIdSpecific = all(idChange <= clusterChange);
    
    if isClusterIdSpecific
        disp('Each cluster is unique to one id. Use tripartite graph pruning.');
        firm_worker_cluster = true;
    else
        disp('Not all clusters are unique to one id or one firmid. Loop over bipartite graph.');
        firm_worker_cluster = false;
    end
end

