function [sel]=connected(firmid,lagfirmid,sel)
 %%% This code was taken from Card, D., Heining, J., & Kline, P. (2013). Workplace heterogeneity and the rise of West German wage inequality. The Quarterly journal of economics, 128(3), 967-1015.
 %%% Source: https://eml.berkeley.edu/~pkline/papers/code_CHK.zip

    A=sparse(lagfirmid(sel),firmid(sel),1); %adjacency matrix
    %make it square
    [m,n]=size(A);
    if m>n
        A=[A,zeros(m,m-n)];
    end
    if m<n
        A=[A;zeros(n-m,n)];
    end
    A=max(A,A'); %connections are undirected

    [sindex, sz]=components(A); %get connected sets
    idx=find(sz==max(sz)); %find largest set
    s=['# of firms: ' int2str(length(A))];
    disp(s);
    s=['# connected sets:' int2str(length(sz))];
    disp(s);
    s=['Largest connected set contains ' int2str(max(sz)) ' firms'];
    disp(s);
    fprintf('\n')
    
    clear A lagfirmid
    
    firmlst=find(sindex==idx); %firms in connected set
    sel=ismember(firmid,firmlst);
    
end