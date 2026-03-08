function inv = mod_inv(a, p)
    [g, u, ~] = gcd(a, p);
    if g ~= 1
        error('Inverse does not exist');
    end
    inv = mod(u, p);
end
