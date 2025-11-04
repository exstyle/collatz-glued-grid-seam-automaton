#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Glued-grid lexicographic order checker (voie 2) — k-penalty version.

On vérifie sur un échantillon (ou un intervalle) que:
    R_k(L(Σ(D))) < R_k(L(D)),
avec
    R_k(L(D)) = ceil( A*log2(2D-1) + h(D mod 3^b)/2^Q - C * k(D) ),
    k(D) = nu2( 3*(2D-1) + 1 ).

La pénalisation par k(D) est conservatrice (pire cas S = k coutures).
Recommandation: A = ceil(1/(-mu))+1,  C > A + (1 - mu).
Pour mu ≈ -0.415, A=4 et C=7 convient.
"""

import argparse, json, math, random, sys

# ---------- Collatz odd-only primitives ----------

def nu2(n: int) -> int:
    if n <= 0: raise ValueError("nu2 requires n>0")
    return ((n & -n).bit_length() - 1)

def oddize(n: int) -> int:
    while (n & 1) == 0:
        n >>= 1
    return n

def Sigma(D: int) -> int:
    T = oddize(3*D - 1)
    return (T + 1) // 2

# ---------- Dual certificate & rank ----------

def load_dual_cert(cert_path: str):
    with open(cert_path, "r", encoding="utf-8") as f:
        c = json.load(f)
    b   = int(c["b"])
    Q   = int(c["scale"])
    mu  = float(c["mu_b"])
    h_d = {int(k): int(v) for k, v in c["h"].items()}
    return b, Q, mu, h_d

def k_of_D(D: int) -> int:
    y = 2*D - 1
    return nu2(3*y + 1)

def Rk_L(D: int, A: int, C: int, b: int, Q: int, hmap: dict) -> int:
    r = D % (3**b)
    y = 2*D - 1
    kval = k_of_D(D)
    h_over_Q = hmap[r] / float(1 << Q)
    return math.ceil(A * math.log2(y) + h_over_Q - C * kval)

# ---------- Check: require Rk(L(Σ(D))) < Rk(L(D)) ----------

def check_lex_order_k(cert_path: str, b: int, A: int|None, C: int|None,
                      samples: int = 0, D_range=None, seed: int = 42, failfast: bool = False) -> bool:
    b_json, Q, mu, hmap = load_dual_cert(cert_path)
    if b_json != b:
        print(f"[WARN] b du JSON = {b_json}, argument --b = {b}", file=sys.stderr)

    if A is None:
        A = math.ceil(1.0/(-mu)) + 1   # e.g. mu≈-0.415 -> A=4
    if C is None:
        # robust: C >= A + (1 - mu) + 1 for ceiling slack
        C = math.ceil(A + (1 - mu)) + 1  # e.g. A=4, mu≈-0.415 -> C=7
    print(f"[INFO] Using A={A}, C={C}, mu={mu:.12f}, Q={Q}, b={b}")

    tested = 0
    worst_margin = -10**9  # max of (Rk(LΣ)-Rk(L)); must stay < 0
    first_bad = None

    def test_D(D: int):
        nonlocal tested, worst_margin, first_bad
        RL   = Rk_L(D,  A, C, b, Q, hmap)
        Ds   = Sigma(D)
        RL_s = Rk_L(Ds, A, C, b, Q, hmap)
        margin = RL_s - RL  # want < 0
        tested += 1
        if margin >= 0 and first_bad is None:
            first_bad = (D, Ds, RL, RL_s, margin, k_of_D(D))
        if margin > worst_margin:
            worst_margin = margin
        return margin

    if D_range:
        Dmin, Dmax = D_range
        for D in range(max(1, Dmin), max(1, Dmax)+1):
            m = test_D(D)
            if failfast and m >= 0: break
    else:
        if samples <= 0:
            print("Fournir --samples N ou bien --range Dmin Dmax", file=sys.stderr)
            return False
        random.seed(seed)
        for _ in range(samples):
            D = random.randint(1, 2_000_000_000)
            m = test_D(D)
            if failfast and m >= 0: break

    if first_bad is None:
        print(f"[OK] Voie(2,k): testés={tested}, pire marge (Rk(LΣ)-Rk(L)) = {worst_margin}  (<0 attendu).")
        return True
    else:
        D, Ds, RL, RLs, margin, kval = first_bad
        print("[FAIL] contre-exemple trouvé :")
        print(f"  D={D}, Σ(D)={Ds}, k(D)={kval}")
        print(f"  Rk(L(D))={RL}, Rk(L(Σ(D)))={RLs}, marge={margin} (doit être < 0)")
        return False

# ---------- CLI ----------

def main():
    ap = argparse.ArgumentParser(description="Checker ordre lexico (voie 2, k-penalty) sur le graphe 2-pas L→L(Σ).")
    ap.add_argument("--cert", required=True, help="JSON certificat dual (h_b*_bf.json)")
    ap.add_argument("--b",    required=True, type=int, help="paramètre b (doit matcher le JSON)")
    ap.add_argument("--A",    type=int, default=None, help="par défaut: ceil(1/(-mu))+1")
    ap.add_argument("--C",    type=int, default=None, help="par défaut: ceil(A + (1-mu))+1")
    ap.add_argument("--samples", type=int, default=0, help="Nb d'échantillons aléatoires")
    ap.add_argument("--range",   nargs=2, type=int, help="Intervalle Dmin Dmax (test déterministe)")
    ap.add_argument("--seed",    type=int, default=42)
    ap.add_argument("--failfast", action="store_true")
    args = ap.parse_args()

    if args.range:
        ok = check_lex_order_k(args.cert, args.b, args.A, args.C,
                               D_range=(args.range[0], args.range[1]), failfast=args.failfast)
    else:
        ok = check_lex_order_k(args.cert, args.b, args.A, args.C,
                               samples=args.samples, seed=args.seed, failfast=args.failfast)
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
