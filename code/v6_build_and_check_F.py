# v6_build_and_check_F.py
# - Construit F avec ancrage F(1)=0 par orbites (clamp/expand)
# - Exporte F(D) en CSV
# - Vérifie les incréments LH/SEAM sur les paires assignées

import argparse
import csv

def nu2(n: int) -> int:
    c = 0
    while n % 2 == 0:
        n //= 2
        c += 1
    return c

def oddize(n: int) -> int:
    return n >> nu2(n)

def Sigma(D: int) -> int:
    return (oddize(3*D - 1) + 1) // 2

def build_F(N: int, mode: str, cap: int | None, progress: int):
    if mode not in ("clamp", "expand"):
        raise ValueError("mode must be 'clamp' or 'expand'")
    if mode == "expand" and cap is None:
        cap = 10 * N

    F: dict[int,int] = {1: 0}   # jauge unique
    assigned = 1
    explored = 0
    cycles: list[list[int]] = []
    skips = 0

    for start in range(2, N+1):
        if start in F:
            continue

        path: list[int] = []
        seen_local: dict[int,int] = {}
        u = start

        while True:
            explored += 1
            if progress and explored % progress == 0:
                print(f"[prog] explored={explored:,}  assigned={assigned:,}  F-size={len(F):,}")

            # Politique de marche
            if mode == "clamp" and u > N:
                skips += 1
                break
            if mode == "expand" and cap is not None and u > cap:
                skips += 1
                break

            if u in F:
                # raccroche: propager en arrière
                val = F[u]
                for prev in reversed(path):
                    val = val - prev         # F(prev) = F(curr) - prev
                    F[prev] = val
                    assigned += 1
                break

            if u in seen_local:
                # cycle 2-pas non trivial (ne devrait pas arriver si acyclique)
                i0 = seen_local[u]
                cyc = path[i0:] + [u]
                cycles.append(cyc)
                break

            seen_local[u] = len(path)
            path.append(u)
            u = Sigma(u)

    return F, cycles, skips

def check_increments(F: dict[int,int], N: int):
    """Retourne des compteurs/écarts sur les paires où D et Sigma(D) sont assignés."""
    n_pairs = 0
    bad_coho = 0
    bad_lh = 0
    bad_seam = 0

    for D in range(2, N+1):
        if D not in F:
            continue
        S = Sigma(D)
        if S not in F:
            continue
        n_pairs += 1
        # Cohomologie: F(S) - F(D) ?= D
        if F[S] - F[D] != D:
            bad_coho += 1
        # Psi(L)=F(D); Psi(M)=F(D)+D+2
        psi_L = F[D]
        psi_M = F[D] + D + 2
        # LH increment expected D+2
        if psi_M - psi_L != D + 2:
            bad_lh += 1
        # SEAM increment expected -2: psi(L(S)) - psi(M)
        if F[S] - psi_M != -2:
            bad_seam += 1

    return n_pairs, bad_coho, bad_lh, bad_seam

def main():
    ap = argparse.ArgumentParser(description="Construire F et vérifier LH/SEAM, exporter CSV")
    ap.add_argument("--N", type=int, default=200000)
    ap.add_argument("--mode", choices=["clamp","expand"], default="clamp")
    ap.add_argument("--cap", type=int, default=None)
    ap.add_argument("--progress", type=int, default=200000)
    ap.add_argument("--out-csv", type=str, default="data/F_prefix.csv")
    args = ap.parse_args()

    F, cycles, skips = build_F(args.N, args.mode, args.cap, args.progress)

    print(f"\n[RES] |F| = {len(F):,} (inclut D=1)")
    print(f"[RES] cycles non triviaux détectés: {len(cycles)}")
    if cycles:
        print("  Exemple de cycle en D:", cycles[0][:20], "...")
    print(f"[RES] orbites tronquées (mode={args.mode}): {skips}")

    n_pairs, bad_coho, bad_lh, bad_seam = check_increments(F, args.N)
    print(f"[CHK] paires (D, Sigma(D)) couvertes = {n_pairs:,}")
    print(f"[CHK] violations cohomologie F(S)-F(D)=D : {bad_coho}")
    print(f"[CHK] violations LH (D+2) : {bad_lh}")
    print(f"[CHK] violations SEAM (-2) : {bad_seam}")

    # Export CSV
    out_path = args.out_csv
    # s'assurer que le dossier existe (façon simple, sans pathlib)
    import os
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["D","F(D)","Sigma(D)","F(Sigma(D))"])
        for D in range(1, args.N+1):
            if D in F:
                S = Sigma(D)
                w.writerow([D, F[D], S, F.get(S,"")])

    print(f"[WRITE] F(D) → {out_path}")

if __name__ == "__main__":
    main()
