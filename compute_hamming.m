function H = compute_hamming(Isi, Isj, L)
    H = 0;
    for tau = 0:L-1
        shifted = mod(Isj + tau, L);
        overlap = length(intersect(Isi, shifted));
        H = max(H, overlap);
    end
end
