function ind_keep = pruning_worker(id,firmid)
%%%% Adapted from the replication package of the Kline, P., Saggio, R., & Sølvsten, M. (2020). Leave-out estimation of variance components. Econometrica, 88(5), 1859-1898.
%%% Source: https://www.dropbox.com/s/iaj3aap3ibfhup8/Replication%20ECTA.zip?dl=1 function pruning_unbal_v3.m


%  This function finds the leave one out largest connected set. 
%  That is, a connected set where if we were to remove the entire working
%  history of a given individual, the resulting graph still remains
%  connected. See appendix B of KSS.
%     
%-------- INPUTS: 
% id: worker identifiers. Dimensions: N* x 1
% firmid: firm identifiers. Dimensions: N* x 1

%  -------- OUTPUT:    
% ind_keep: indicator to keep observation that will filter the original data. will need
% to update the identifier codes outside

%%% Keep track of removals
% id will be filtered and renamed
% id_0 will be filtered only
% id_orig will be used to generate the filtering variable of the same size as the original data
id_0    = id;
id_orig = id_0; % to come back to the original data size

n_of_bad_workers = 1;

while n_of_bad_workers >= 1

    %find movers
    J=max(firmid);
    gcs = [NaN; id(1:end-1)];
    gcs = id~=gcs;
    lagfirmid=[NaN; firmid(1:end-1)];
    lagfirmid(gcs==1)=NaN; %%first obs for each worker
    move=(firmid~=lagfirmid);
    move(gcs==1)=0;
    move=accumarray(id,move);
    move=(move>0);
    move=move(id);
    
    %id and firmids associated with movers only
    id_movers=id(move);
    firmid_mover=firmid(move);

    %need to normalize id_movers and keep a dictionary. 
    id_movers_orig=id_movers;
    [~,~,id_movers]=unique(id_movers);
    n_movers=max(id_movers);
    
    %unique pairs    
    list=[id_movers firmid_mover];
    list=unique(list,'rows');
    
    %Construct Adj. Matrix of Bipartite Graph
    B=sparse(list(:,1),list(:,2),1,n_movers,J);
    G = [ sparse( n_movers, n_movers ), B; B.', sparse(J, J)];
    
    %Inspect the resulting Graph: find the workers that constitute an
    %articulation point
    [~,artic_points] = biconncomp(graph(G)); % Matlab built in function
    bad_workers=artic_points(artic_points<=n_movers);
    
    %Now get the right index for these workers
    sel=ismember(id_movers,bad_workers);
    bad_workers=id_movers_orig(sel);
    bad_workers=unique(bad_workers);
    n_of_bad_workers=size(bad_workers,1);

    %Drop these workers
    sel=~ismember(id,bad_workers);
    id=id(sel);
    firmid=firmid(sel);
    id_0 = id_0(sel);
    
	% Rename ids
	[~,~,n]=unique(firmid);
	firmid=n;
	[~,~,n]=unique(id);
	id=n;
 
	% Adjacency matrix
    A=build_adj(id,firmid);
    
    [sindex,sz] = conncomp(graph(A)); % Matlab built in function. Get connected sets
    idx=find(sz==max(sz)); %find largest set
    firmlst=find(sindex==idx); %firms in connected set
    sel=ismember(firmid,firmlst);
    
	% Select the connected set
    firmid=firmid(sel);
    id=id(sel);
    id_0 = id_0(sel);
    
    % Rename ids
    [~,~,n] = unique(firmid);
    firmid = n;
    [~,~,n] = unique(id);
    id=n;  
 
end

% generate indicator to keep of the length of the original data
ind_keep = ismember(id_orig,id_0); % search A in B


end

    