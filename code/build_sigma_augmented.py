#!/usr/bin/env python3
# (c) 2025 S. Gobeaux — MIT
#
# Build the augmented seam automaton Σ_b:
#  Nodes: residues r mod 3^b (we interpret r as D ≡ r (mod 3^b))
#  For each node r, and each seam type κ∈{1,2}:
#     r' ≡ (3r - 1) * (2^κ)^(-1) (mod 3^b)          # seam transition
#     w_raw = log2(3) - κ                            # seam weight (raw)
#     δ_r   = log2( 1 + 1/(3*(2*D_min - 1)) )       # local bound (max) for that residue
#     w_aug = w_raw + δ_r                            # augmented weight
#
# Output CSV columns: src,dst,kappa,delta,w_raw,w_aug
#
# Rationale:
#  δ is the small “size split” correction log2(1 + 1/(3y)).
#  For a fixed residue r (on D), the worst-case δ occurs at the smallest D≥1 with D≡r (mod 3^b),
#  namely D_min = r (if r>0) or D_min = 3^b (if r=0). Then y_min = 2*D_min - 1.
#  Using this pessimistic δ_r yields a *conservative* augmented graph; if its min–mean is negative,
#  we have a uniform margin ε > 0 for blocks (⇒ an O(log y0) upper bound for flight time).
#
import argparse, csv, math
import os
from typing import Tuple

def inv_mod(a: int, m: int) -> int:
    # modular inverse (Python 3.8+: pow(a, -1, m) works; keep EEA for portability)
    try:
        return pow(a, -1, m)
    except TypeError:
        # Extended Euclidean Algorithm
        t, new_t = 0, 1
        r, new_r = m, a % m
        while new_r != 0:
            q = r // new_r
            t, new_t = new_t, t - q * new_t
            r, new_r = new_r, r - q * new_r
        if r != 1:
            raise ValueError(f"{a} not invertible mod {m}")
        if t < 0:
            t += m
        return t

def delta_residue(r: int, mod3b: int) -> float:
    # δ_r = log2( 1 + 1/(3*y_min) ) with y_min = 2*D_min - 1
    # D_min is the smallest positive integer ≡ r (mod 3^b)
    D_min = r if r != 0 else mod3b
    y_min = 2 * D_min - 1
    return math.log2(1.0 + 1.0/(3.0 * y_min))

def main():
    ap = argparse.ArgumentParser(description="Build augmented seam automaton Σ_b with δ by residue.")
    ap.add_argument("--b", type=int, required=True, help="window size b (mod 3^b)")
    ap.add_argument("--out", type=str, default=None, help="output CSV path")
    args = ap.parse_args()

    b = args.b
    mod = 3 ** b
    out_path = args.out or f"data/sigma_b{b}_aug.csv"

    inv2 = inv_mod(2, mod)
    inv4 = (inv2 * inv2) % mod

    log2_3 = math.log2(3.0)

    n_edges = 0
    
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["src","dst","kappa","delta","w_raw","w_aug"])
        for r in range(mod):
            delt = delta_residue(r, mod)
            # κ = 1
            dst1 = ((3*r - 1) * inv2) % mod
            w_raw1 = log2_3 - 1.0
            w_aug1 = w_raw1 + delt
            w.writerow([r, dst1, 1, f"{delt:.12f}", f"{w_raw1:.12f}", f"{w_aug1:.12f}"])
            n_edges += 1
            # κ = 2
            dst2 = ((3*r - 1) * inv4) % mod
            w_raw2 = log2_3 - 2.0
            w_aug2 = w_raw2 + delt
            w.writerow([r, dst2, 2, f"{delt:.12f}", f"{w_raw2:.12f}", f"{w_aug2:.12f}"])
            n_edges += 1

    print(f"[OK] Σ_b built for b={b} (mod=3^{b}={mod}) → edges={n_edges} ; file: {out_path}")

if __name__ == "__main__":
    main()
