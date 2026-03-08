function p = next_prime(n)
    if isprime(n)
        p = n;
    else
        p = n + 1;
        while ~isprime(p)
            p = p + 1;
        end
    end
end
