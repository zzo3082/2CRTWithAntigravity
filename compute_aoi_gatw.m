function AoI = compute_aoi_gatw(Is, C_idx, user_idx, K, L)
    user_seq = Is{C_idx(user_idx)};
    
    % Use Base-3 index mapping to replace slow string Map
    % Map size is 3^K. Index = 1 + c1*3^0 + c2*3^1 + ... (where c_j in {0,1,2})
    P_old = zeros(3^K, 1);
    P_old(1) = 1;
    pow3 = 3.^(0:K-1);
    
    for i = 1:length(C_idx)
        if i == user_idx
            continue;
        end
        other_seq = Is{C_idx(i)};
        
        tau_patterns = zeros(L, K);
        for tau = 0:L-1
            shifted = mod(other_seq + tau, L);
            tau_patterns(tau+1, :) = ismember(user_seq, shifted);
        end
        
        P_new = zeros(3^K, 1);
        non_zero_idx = find(P_old > 0);
        vals = P_old(non_zero_idx);
        
        v_idx = non_zero_idx - 1;
        v_mat = repmat(v_idx, 1, K);
        pow3_mat = repmat(pow3, length(non_zero_idx), 1);
        c = rem(floor(v_mat ./ pow3_mat), 3);
        
        for t = 1:L
            pat_t = tau_patterns(t, :);
            pat_mat = repmat(pat_t, length(non_zero_idx), 1);
            
            c_new = c + pat_mat;
            c_new(c_new > 2) = 2; % Cap at 2
            
            new_idx = 1 + c_new * pow3';
            P_new = P_new + accumarray(new_idx, vals, [3^K, 1]);
        end
        P_old = P_new;
    end
    
    E_W = 0;
    E_eta = 0;
    total_combs = L^(length(C_idx) - 1);
    
    non_zero_idx = find(P_old > 0);
    
    for idx_i = 1:length(non_zero_idx)
        idx = non_zero_idx(idx_i);
        count = P_old(idx);
        prob = count / total_combs;
        
        v_val = idx - 1;
        pat = rem(floor(v_val ./ pow3), 3);
        
        W = [];
        for j = 1:K
            if pat(j) <= 1
                W(end+1) = user_seq(j);
            end
        end
        
        E_W = E_W + prob * length(W);
        if ~isempty(W)
            gaps = zeros(1, length(W));
            for j = 1:length(W)-1
                gaps(j) = W(j+1) - W(j);
            end
            gaps(end) = L - W(end) + W(1);
            
            sum_gaps_sq = sum(gaps.^2);
            E_eta = E_eta + prob * sum_gaps_sq;
        end
    end
    
    % As derived by renewing combinatorics, theoretical AoI
    AoI = E_eta / (2 * L) + 0.5;
end
