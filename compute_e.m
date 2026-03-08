function e_val = compute_e(A, B, L, K)
    % A is a subset, B is a full characteristic set Is_h
    % Count the number of tau such that |A intersect (B + tau)| == 2
    e_val = 0;
    for tau = 0:L-1
        shifted = mod(B + tau, L);
        overlap = length(intersect(A, shifted));
        if overlap == 2
            e_val = e_val + 1;
        end
    end
end
