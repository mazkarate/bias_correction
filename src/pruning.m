    function [y,id,firmid,all_covariates,period,id_old,firmid_old] = pruning(y,id,firmid,all_covariates,n_FE,period,id_old,firmid_old);
%%%% Adjusted from the replication package of the Kline, P., Saggio, R., & Sølvsten, M. (2020). Leave-out estimation of variance components. Econometrica, 88(5), 1859-1898.
%%% Source: https://www.dropbox.com/s/iaj3aap3ibfhup8/Replication%20ECTA.zip?dl=1 function pruning_unbal_v3.m


%  This function finds the leave one out largest connected set. 
%  That is, a connected set where if we were to remove the entire working
%  history of a given individual, the resulting graph still remains
%  connected. See appendix B of KSS.
%     
%-------- INPUTS: 
% y: outcome variable. Dimensions: N* x 1; N*= # of person-year observations.
% id: worker identifiers. Dimensions: N* x 1
% firmid: firm identifiers. Dimensions: N* x 1
% all_covariates: original matrix of covariates. Dimensions: N* x K
% period: vector with periods to compute the second order moments. Dimensions: N* x 1
% id_old: original worker identifiers. Dimensions: N* x 1             
% firmid_old: original firm identifiers. Dimensions: N* x 1 

%  -------- OUTPUTS:    
% y,id, firmid, all_covariates, period, id_old and firmid_old in the leave one out sample. firmid, id and other fixed effects are all 
%          relabbeled to run estimation in this new sample. firmid_old, 
%          id_old preserve their orginal denomination. 



n_of_bad_workers=1;
while n_of_bad_workers>=1
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
    [artic_points, CC] = biconnected_components(G);
    bad_workers=artic_points(artic_points<=n_movers);
    
%Now get the right index for these workers
    sel=ismember(id_movers,bad_workers);
    bad_workers=id_movers_orig(sel);
    bad_workers=unique(bad_workers);
    n_of_bad_workers=size(bad_workers,1);

%Drop these workers
    sel=~ismember(id,bad_workers);
    y=y(sel);
    id=id(sel);
    firmid=firmid(sel);
    all_covariates = all_covariates(sel,:);
    lagfirmid=lagfirmid(sel);
    firmid_old=firmid_old(sel);
    id_old=id_old(sel);
    if ~isnan(period)
        period=period(sel); % if there is a grouping variable for the corrections update
    end
    

	% Rename ids
	[~,~,n]=unique(firmid);
	firmid=n;
	[~,~,n]=unique(id);
	id=n;

	% Relabel other FE ids
	if n_FE>2
		for j=3:n_FE
			[~,~,n]=unique(all_covariates(:,j));
			all_covariates(:,j)=n;
		end
	end
 
	% Adjacency matrix
    A=build_adj(id,firmid);
    
    [sindex, sz]=components(A); %get connected sets
    idx=find(sz==max(sz)); %find largest set
    firmlst=find(sindex==idx); %firms in connected set
    sel=ismember(firmid,firmlst);
    
	% Select
    y=y(sel);
    firmid=firmid(sel);
    id=id(sel);
    firmid_old=firmid_old(sel);
    id_old=id_old(sel);
    all_covariates = all_covariates(sel,:);
    if ~isnan(period)
        period=period(sel); % if there is a grouping variable for the corrections update
    end
    
    % Rename ids
    [~,~,n]=unique(firmid);
    firmid=n;
    [~,~,n]=unique(id);
    id=n;  
    
    % Relabel other FE
    if n_FE>2
        for j=3:n_FE
            [~,~,n]=unique(all_covariates(:,j));
            all_covariates(:,j)=n;
        end
    end
 
end
end

    