from __future__ import annotations

def tri(x, a, b, c):
    if x <= a or x >= c:
        return 0.0
    if x == b:
        return 1.0
    return (x - a) / float(b - a) if x < b else (c - x) / float(c - b)

def trap(x, a, b, c, d):
    if x <= a or x >= d:
        return 0.0
    if b <= x <= c:
        return 1.0
    return (x - a) / float(b - a) if x < b else (d - x) / float(d - c)

def calibrate_memberships(stats, priors=None):
    D, H, N = stats.get('D', 0.0), stats.get('H', 0.0), stats.get('N', 1)
    Depth = stats.get('Depth', 10.0)

    P = priors or {
        'D': {'low': (0, 50, 150), 'mid': (120, 250, 420), 'high': (350, 650, 1000)},
        'H': {'low': (0.00, 0.06, 0.12), 'mid': (0.10, 0.22, 0.35), 'high': (0.30, 0.45, 0.60)},
        'N': {'small': (0, 120, 220), 'mid': (180, 320, 520), 'large': (450, 800, 1200)},
        'Depth': {'low': (0, 3, 8), 'mid': (5, 10, 15), 'high': (12, 20, 50)} 
    }
    
    mu = {
        'D_low':  tri(D, *P['D']['low']),
        'D_mid':  tri(D, *P['D']['mid']),
        'D_high': tri(D, *P['D']['high']),
        'H_low':  tri(H, *P['H']['low']),
        'H_mid':  tri(H, *P['H']['mid']),
        'H_high': tri(H, *P['H']['high']),
        'N_small': tri(N, *P['N']['small']),
        'N_mid':   tri(N, *P['N']['mid']),
        'N_large': tri(N, *P['N']['large']),
        'Depth_low':  tri(Depth, *P['Depth']['low']),
        'Depth_mid':  tri(Depth, *P['Depth']['mid']),
        'Depth_high': tri(Depth, *P['Depth']['high']), 
    }
    if Depth > 20:
        mu['Depth_high'] = 1.0
        
    return mu

def infer_hmm_overrides(mu, stats):
    intra_base = 0.80
    trans_base = 1.00

    intra_adj = (
        -0.05 * mu.get('H_high', 0.0) +
         0.03 * mu.get('H_low',  0.0) +
         0.08 * mu.get('Depth_high', 0.0) -
         0.08 * mu.get('Depth_low', 0.0)
    )
    
    intra_p = intra_base + intra_adj
    intra_p = min(0.98, max(0.55, intra_p)) 

    trans_adj = (
        +0.50 * mu.get('D_high', 0.0) +
        +0.30 * mu.get('N_large', 0.0) -
         0.20 * mu.get('D_low',  0.0)
    )
    trans_mult = trans_base + trans_adj
    trans_mult = min(5.0, max(0.2, trans_mult))

    return {'intra_p': float(intra_p), 'trans_mult': float(trans_mult)}

def infer_hmm_base(mu, stats):
    D = stats.get('D', 0.0)
    H = stats.get('H', 0.0)
    N = stats.get('N', 1)

    predicted_homogeneity = 0.85 - 0.2 * mu.get('H_high', 0.0)
    predicted_homogeneity = max(0.70, min(0.95, predicted_homogeneity))

    base_cross = 300
    if mu.get('D_high', 0.0) > 0:
        base_cross += 300 * mu['D_high']
    if mu.get('N_large', 0.0) > 0:
        base_cross += 200 * mu['N_large']
    if mu.get('D_low', 0.0) > 0:
        base_cross -= 150 * mu['D_low']
    predicted_cross_count = int(max(100, min(1200, base_cross)))

    return {
        'predicted_homogeneity': predicted_homogeneity,
        'predicted_cross_count': predicted_cross_count
    }
