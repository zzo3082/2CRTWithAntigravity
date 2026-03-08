import math

def nchoosek(n, k):
    return math.comb(n, k)

K = 6
pK = 7
L = 42

Is0 = [(j*0)%7 for j in range(K)]  # Wait, sequence depends on CRT mapping. 
# For K=6, pK=7: 
# let's write build_2CRT in python
def crt_map(t, pK, K):
    return (t % pK, t % K)

Is = []
for g in range(pK):
    seq = []
    for j in range(K):
        target = ((j*g)%pK, (j*1)%K)
        for t in range(L):
            if crt_map(t, pK, K) == target:
                seq.append(t)
                break
    Is.append(sorted(seq))

seq = []
for j in range(K):
    target = ((j*1)%pK, 0)
    for t in range(L):
        if crt_map(t, pK, K) == target:
            seq.append(t)
            break
Is.append(sorted(seq))

s0 = Is[0] # s0 is Is_0. 
# Set 1 = [s0, s2, s3, s4, s5, s7] -> these are C_idx = 1, 3, 4, 5, 6, 8 (1-indexed)
C_idx = [0, 2, 3, 4, 5, 7]
print("s0:", s0)

# s0 gap:
delta = []
for j in range(K-1):
    delta.append(s0[j+1] - s0[j])
delta.append((s0[0] - s0[-1]) % L)

print("delta:", delta)

C_delta = []
for d in range(K):
    cd = 0
    for j in range(K):
        # subscript in 0-indexed: j and (j - d)%K
        idx1 = j
        idx2 = (j - d) % K
        cd += delta[idx1] * delta[idx2]
    C_delta.append(cd)

print("C_delta:", C_delta)

h_n = []
for n in range(1, K+1):
    h = 0
    for i in range(1, n+1):
        h += ((-1)**(i+1)) * nchoosek(n, i) * ((1 - i/pK)**(K-1))
    h_n.append(h)

print("h_n:", h_n)

# Let's try combining them to get E[S]!
# E[S] = sum(C_delta[d] * something). We know Periodic Set 1 s0 Avg corresponds to E[S] + 20.5 = 25.711 => E[S] ~ 5.211.
print("Target E[S]:", 25.711 - 20.5)

