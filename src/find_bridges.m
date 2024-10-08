function bridge_index = find_bridges(G)
    bcedges = biconncomp(G); % Gives indices indicating assignment of each edge to a biconnected component
    edgeCounts = accumarray(bcedges', 1);
    bridge_index = find(edgeCounts(bcedges) == 1); % Returns indices of edges for which exactly one edge (the same edge itself) in the graph has the same 'bin'
%     bridge_index = edgeCounts(bcedges) == 1; % Returns index indicating edges for which exactly one edge (the same edge itself) in the graph has the same 'bin'
end


    

