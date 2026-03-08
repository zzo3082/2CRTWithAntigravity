function h_val = compute_h_hat(I_tau, user_idx, C_idx, Is, L, K)
    % Compute h_hat_n for Type II sequence formula 31
    % I_tau is the shifted sequence, length K
    n_len = length(I_tau);
    h_val = 0;
    
    % Loop over all non-empty subsets A of I_tau
    % Since K is at most 6, 2^6-1 = 63 combinations
    for i = 1:(2^n_len - 1)
        % Get binary representation of i
        bin_i = dec2bin(i, n_len);
        A = [];
        for j = 1:n_len
            if bin_i(j) == '1'
                A(end+1) = I_tau(j);
            end
        end
        
        size_A = length(A);
        term = (-1)^(size_A + 1);
        
        % Product over sh in C \ {sg}
        for idx = 1:length(C_idx)
            if idx ~= user_idx
                sh = Is{C_idx(idx)};
                eval_val = compute_e(A, sh, L, K);
                term = term * (1 - (size_A * K - eval_val) / L);
            end
        end
        h_val = h_val + term;
    end
end
