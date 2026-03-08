function [Is, L, pK] = build_2CRT(K)
    pK = next_prime(K);
    
    if isprime(K)
        % Definition 2
        L = pK^2;
        Is = cell(1, K+1);
        for g = 0:(pK-1)
            % Is{g+1} = { mod(j * (1 + g*p), p^2) : j = 0,1,...,p-1 }
            elements = mod((0:pK-1) * (1 + g*pK), L);
            Is{g+1} = sort(unique(elements));
        end
        % Is{p+1}
        elements = mod((0:pK-1) * pK, L);
        Is{pK+1} = sort(unique(elements));
    else
        % Definition 3
        L = pK * K;
        Is = cell(1, pK+1);
        for g = 0:(pK-1)
            Is{g+1} = zeros(1, K);
            for j = 0:K-1
                target = [mod(j*g, pK), mod(j*1, K)];
                for t = 0:L-1
                    cmap = crt_map(t, pK, K);
                    if cmap(1) == target(1) && cmap(2) == target(2)
                        Is{g+1}(j+1) = t;
                        break;
                    end
                end
            end
            Is{g+1} = sort(unique(Is{g+1}));
        end
        % Is{pK+1}
        Is{pK+1} = zeros(1, K);
        for j = 0:K-1
            target = [mod(j*1, pK), mod(j*0, K)];
            for t = 0:L-1
                cmap = crt_map(t, pK, K);
                if cmap(1) == target(1) && cmap(2) == target(2)
                    Is{pK+1}(j+1) = t;
                    break;
                end
            end
        end
        Is{pK+1} = sort(unique(Is{pK+1}));
    end
    
    % Assertions
    for i = 1:length(Is)
        assert(length(Is{i}) == K, 'Each sequence must have length K');
    end
end
