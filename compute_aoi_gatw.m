function AoI = compute_aoi_gatw(Is, C_idx, user_idx, K, L)
    user_seq = Is{C_idx(user_idx)};
    
    current_dist = containers.Map('KeyType', 'char', 'ValueType', 'double');
    init_pat = char(zeros(1, K) + '0');
    current_dist(init_pat) = 1;
    
    for i = 1:length(C_idx)
        if i == user_idx
            continue;
        end
        other_seq = Is{C_idx(i)};
        
        shift_dist = containers.Map('KeyType', 'char', 'ValueType', 'double');
        for tau = 0:L-1
            pat = char(zeros(1, K) + '0');
            shifted = mod(other_seq + tau, L);
            for j = 1:K
                if ismember(user_seq(j), shifted)
                    pat(j) = '1';
                end
            end
            if isKey(shift_dist, pat)
                shift_dist(pat) = shift_dist(pat) + 1;
            else
                shift_dist(pat) = 1;
            end
        end
        
        new_dist = containers.Map('KeyType', 'char', 'ValueType', 'double');
        k1 = keys(current_dist);
        v1 = values(current_dist);
        k2 = keys(shift_dist);
        v2 = values(shift_dist);
        for i1 = 1:length(k1)
            pat1 = k1{i1};
            cnt1 = v1{i1};
            for i2 = 1:length(k2)
                pat2 = k2{i2};
                cnt2 = v2{i2};
                
                new_pat = char(zeros(1, K) + '0');
                for j = 1:K
                    new_pat(j) = char((pat1(j) - '0') + (pat2(j) - '0') + '0');
                end
                
                if isKey(new_dist, new_pat)
                    new_dist(new_pat) = new_dist(new_pat) + cnt1 * cnt2;
                else
                    new_dist(new_pat) = cnt1 * cnt2;
                end
            end
        end
        current_dist = new_dist;
    end
    
    E_W = 0;
    E_eta = 0;
    total_combs = L^(length(C_idx) - 1);
    
    k = keys(current_dist);
    v = values(current_dist);
    
    for i = 1:length(k)
        pat = k{i};
        count = v{i};
        prob = count / total_combs;
        
        W = [];
        for j = 1:K
            if (pat(j) - '0') <= 1
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
