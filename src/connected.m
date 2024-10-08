function [sel]=connected(firmid, lagfirmid, sel)
 %%% This code was adapted from Card, D., Heining, J., & Kline, P. (2013). Workplace heterogeneity and the rise of West German wage inequality. The Quarterly journal of economics, 128(3), 967-1015.
 %%% Source: https://eml.berkeley.edu/~pkline/papers/code_CHK.zip

    A = sparse(double(lagfirmid(sel)), double(firmid(sel)), 1); % Adjacency matrix

    % Make it square
    [m, n] = size(A);
    if m > n
        A = [A, sparse(m, m - n)];
    end
    if m < n
        A = [A; sparse(n - m, n)];
    end
    A = max(A, A'); % Connections are undirected

    [sindex, sz] = conncomp(graph(A)); % Get connected sets
    idx = find(sz == max(sz)); % Find largest set
    s = ['# of firms: ' int2str(length(A))];
    disp(s);
    s = ['# connected sets: ' int2str(length(sz))];
    disp(s);
    s = ['Largest connected set contains ' int2str(max(sz)) ' firms'];
    disp(s);
    fprintf('\n')
    
    clear A
    
    firmlst = find(sindex==idx); % Firms in connected set
    sel=ismember(firmid,firmlst);
    
end