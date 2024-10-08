function [ch_group_cell, i_group, N_id_group_cell, N_firmid_group_cell, N_obs_group_cell] = compute_group_changes(group, id, firmid, dim_fe, N_multi_group, N_u_group)

    % Initialize cell arrays to store the outputs for each grouping
    ch_group_cell = cell(N_multi_group, 1);
    N_id_group_cell = cell(N_multi_group, 1);
    N_firmid_group_cell = cell(N_multi_group, 1);
    N_obs_group_cell = cell(N_multi_group, 1);

    [group_sorted, i_group] = sort(group); % Sort the data per group
    IDX = diff(group_sorted); % Find where the subscript changes
    clear group_sorted

    % Loop over each grouping
    for i = 1:N_multi_group
        % Find indices of nonzero elements and add start and end of data
        aux = find(IDX(:, i));
        ch_group_i = [0; aux; size(group, 1)];
        
        % Store ch_group for the current grouping
        ch_group_cell{i} = ch_group_i;
        
        % Sort id and firmid for counts
        id_sorted = id(i_group(:, i));
        firmid_sorted = firmid(i_group(:, i));

        % Initialize vectors to store unique ids, firmids, and counts
        aux_id = zeros(N_u_group(i), 1);
        aux_firmid = zeros(N_u_group(i), 1);
        aux_count = zeros(N_u_group(i), 1);

        % Loop over each group and compute unique counts
        for j = 2:(N_u_group(i) + 1)
            aux_id(j-1) = numel(unique(id_sorted(ch_group_i(j-1)+1:ch_group_i(j))));
            aux_firmid(j-1) = numel(unique(firmid_sorted(ch_group_i(j-1)+1:ch_group_i(j))));
            aux_count(j-1) = numel(firmid_sorted(ch_group_i(j-1)+1:ch_group_i(j)));
        end

        % Store the results into the cell arrays
        N_id_group_cell{i} = [dim_fe(1); aux_id];
        N_firmid_group_cell{i} = [dim_fe(2) + 1; aux_firmid];
        N_obs_group_cell{i} = [sum(aux_count); aux_count];
    end
    clear aux aux_id aux_firmid aux_count

    % Do dynamic sorting: create index to sort for first grouping from
    % original positions. Then create index to sort for second grouping
    % starting from the ordered positions of the first grouping, and so
    % on...

    if N_multi_group > 1
        for i = 2:N_multi_group
            aux_g = group(:, i);
            aux_g = aux_g(i_group(:, i - 1)); %Use previous grouping order
            [~ , i_group(:, i)] = sort(aux_g);
        end
    end



end
