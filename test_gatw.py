import itertools
from collections import defaultdict

K = 6
pK = 7
L = 42

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

C_idx = [0, 2, 3, 4, 5, 7]

def compute_aoi(user_idx):
    user_seq = Is[C_idx[user_idx]]
    # Convolution!
    # For each other sequence, create histogram of hit counts on user_seq
    # v_hits[seq_idx][shift]: binary tuple of length K indicating hit pattern
    # actually we just convolve dicts of {pattern: count}
    
    current_dist = {tuple([0]*K): 1}
    for i in range(len(C_idx)):
        if i == user_idx: continue
        other_seq = set(Is[C_idx[i]])
        
        # calculate all L shifts
        shift_dist = defaultdict(int)
        for tau in range(L):
            pattern = tuple([1 if ((u - tau)%L) in other_seq else 0 for u in user_seq])
            shift_dist[pattern] += 1
            
        # convolve
        new_dist = defaultdict(int)
        for pat1, cnt1 in current_dist.items():
            for pat2, cnt2 in shift_dist.items():
                new_pat = tuple([pat1[j] + pat2[j] for j in range(K)])
                new_dist[new_pat] += cnt1 * cnt2
        current_dist = new_dist
        
    E_W = 0
    E_eta = 0
    
    # For periodic traffic exactly:
    E_delay = 0
    
    total_combs = L**(len(C_idx)-1)
    
    for pat, count in current_dist.items():
        prob = count / total_combs
        
        # success pattern
        W = [user_seq[j] for j in range(K) if pat[j] <= 1]
        
        E_W += prob * len(W)
        if len(W) > 0:
            gaps = []
            for j in range(len(W)-1):
                gaps.append(W[j+1] - W[j])
            gaps.append(L - W[-1] + W[0])
            E_eta += prob * sum(g**2 for g in gaps)
            
            # periodic delay
            # generation u in 0..L-1
            # expected distance to NEXT successful slot
            # if we integrate over continuous shift tau_i, it's uniform 0..L-1
            # actually sum d = g(g-1)/2 or sum g(g+1)/2?
            # distance to next: sum_{g} (g-1)g/2 / L + 1 / L * |W|? No.
            # delay is distance from u to next W. 
            # sum_{u=0..L-1} next(W) = sum_{w in W} sum_{i=0..gap-1} (gap - i) ? Wait, if u=w, delay is 0 ? No, delivery is at end of slot. Generation is at beginning. So delay is 1?
            # let's assume delay is EXACTLY (sum gap^2 + L) / 2L
            S = sum( g*(g-1)/2 for g in gaps ) / L
            E_delay += prob * S
            
    gatw = E_eta / 2 / L - 0.5
    periodic = gatw + 21.0 + 0.5  # wait, difference is ...
    # Wait, the formula for GatW AoI is exactly E[eta]/2E[W] - 1/2?
    # No, it's  AoI = E_eta / E_W ? No, gatw = (E_eta*L) / (2 * L * E_W) - 0.5 = E_eta / (2 E_W) - 0.5
    
    gatw2 = E_eta / (2 * E_W) - 0.5
    print(f"user {user_idx}: EW={E_W:4f}, E_eta={E_eta:.4f}, gatw={gatw2:.4f}, mygatw={gatw:.4f}")
    
    # what is formula 24, 26? E[Y] = L^2 / E[W] ? No, T_i = E[W]/L. E[Y] = 1/T_i = L/E[W]
    # E[Y^2] = E_eta / E_W
    # So AoI = E[Y^2] / 2 / E[Y] - 1/2 = E_eta / (2 L) - 1/2
    pass

for u in range(6):
    compute_aoi(u)
