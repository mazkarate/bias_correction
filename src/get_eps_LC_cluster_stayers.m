function eps_LC = get_eps_LC_cluster_stayers(id, cluster, N_C, r)

% Relabel identifiers
id      = id - min(id) + 1;
cluster = cluster - min(cluster) + 1;

% Get time in sample for each stayer
T = accumarray(id,1);
T = T(id); % Get back to right size
NT = size(T, 1);


% Form sparse matrix to store inverse of M for stayers
one_matrix = sparse(1:NT, cluster, 1);
one_matrix = (one_matrix * one_matrix');

% Use Sherman-Morrison Formula to get inverse
inv_M_aux = (1./(T- N_C)).*one_matrix; %We are missing the "+ I"
eps_LC = inv_M_aux * r + r; % the "+ r" makes up for the missing "+ I"





end