import itertools
from collections import defaultdict

def run_gatw(K=8):
    pK = 11
    L = pK * K
    
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

    # K=8 Set 1: Just pick index 1 to 8? Let's check Set 1 indices user used: s2, s3, s4, s5, s6, s7, s8, s9
    # This corresponds to Is[2] to Is[9]
    C_idx = list(range(2, 10))

    def compute_aoi(user_idx):
        user_seq = Is[C_idx[user_idx]]
        
        current_dist = {tuple([0]*K): 1}
        for i in range(len(C_idx)):
            if i == user_idx: continue
            other_seq = set(Is[C_idx[i]])
            
            shift_dist = defaultdict(int)
            for tau in range(L):
                pattern = tuple([1 if ((u - tau)%L) in other_seq else 0 for u in user_seq])
                shift_dist[pattern] += 1
                
            new_dist = defaultdict(int)
            for pat1, cnt1 in current_dist.items():
                for pat2, cnt2 in shift_dist.items():
                    new_pat = tuple([pat1[j] + pat2[j] for j in range(K)])
                    new_dist[new_pat] += cnt1 * cnt2
            current_dist = new_dist
            
        E_W = 0
        E_eta = 0
        total_combs = L**(len(C_idx)-1)
        
        for pat, count in current_dist.items():
            prob = count / total_combs
            W = [user_seq[j] for j in range(K) if pat[j] <= 1]
            E_W += prob * len(W)
            if len(W) > 0:
                gaps = []
                for j in range(len(W)-1):
                    gaps.append(W[j+1] - W[j])
                gaps.append(L - W[-1] + W[0])
                E_eta += prob * sum(g**2 for g in gaps)
                
        gatw = E_eta / 2 / L - 0.5
        print(f"User s{C_idx[user_idx]}: gatw_base = {gatw:.4f}")
        return gatw

    for u in range(len(C_idx)):
        compute_aoi(u)

if __name__ == '__main__':
    run_gatw()
