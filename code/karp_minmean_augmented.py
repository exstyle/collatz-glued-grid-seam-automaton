#!/usr/bin/env python3
# (c) 2025 S. Gobeaux — MIT
#
# Min–mean cycle (Karp) on augmented weights, with PROGRESS reporting.
#
# Usage:
#   python code/karp_minmean_augmented.py --edges data/sigma_b10_aug.csv --report --progress
# Options:
#   --progress               : print progress during DP (row by row)
#   --progress-every-k  N    : print every N DP rows (default: 1)
#   --progress-interval  S   : min seconds between prints (default: 0.5)
#   --src-col/--dst-col/--weight-col : CSV column names (defaults: src,dst,w_aug)
#   --report                : human summary with epsilon and optional bound
#   --y0 / --Ystar / --Cclose : instantiate O(log y0) bound if epsilon>0
#
import argparse, csv, math, sys, time
from typing import List, Tuple

def read_edges(path: str, src_col="src", dst_col="dst", weight_col="w_aug"):
    with open(path, "r", encoding="utf-8") as f:
        rd = csv.DictReader(f)
        rows = list(rd)
    nodes = set()
    edges = []
    for row in rows:
        try:
            u = int(row[src_col]); v = int(row[dst_col])
        except Exception as e:
            raise SystemExit(f"[ERR] bad src/dst in row: {row}") from e
        try:
            w = float(row[weight_col])
        except Exception:
            w = float(row.get("weight", row.get("w", row.get("w_aug"))))
        nodes.add(u); nodes.add(v)
        edges.append((u, v, w))
    nodes = sorted(nodes)
    mapping = {node:i for i, node in enumerate(nodes)}
    n = len(nodes)
    # incoming adjacency for Karp DP
    adj_in = [[] for _ in range(n)]
    m_edges = 0
    for (u, v, w) in edges:
        uu = mapping[u]; vv = mapping[v]
        adj_in[vv].append((uu, w))
        m_edges += 1
    return n, adj_in, mapping, m_edges

def karp_min_mean_progress(n: int,
                           adj_in: List[List[Tuple[int,float]]],
                           progress: bool = False,
                           progress_every_k: int = 1,
                           progress_interval: float = 0.5):
    """
    Karp's algorithm (min cycle mean) via DP on negated weights, with progress prints.
    Memory: O(n^2) for dp table (n+1 rows); suitable for tailles ~a few thousands.
    """
    NEG_INF = -1e300
    # dp[k][v] = best (max on neg. weights) value for walks of length k ending at v
    dp = [[NEG_INF]*n for _ in range(n+1)]
    for v in range(n):
        dp[0][v] = 0.0

    t0 = time.time()
    last_print = t0
    total_updates = 0

    if progress:
        print(f"[karp] start DP: rows=0..{n}  nodes={n}", flush=True)

    # DP filling
    for k in range(1, n+1):
        row_updates = 0
        for v in range(n):
            best = NEG_INF
            # transition on NEGATED weights => maximize
            inc = 0
            for (u, w) in adj_in[v]:
                val = dp[k-1][u] + (-w)
                if val > best:
                    best = val
                inc += 1
            dp[k][v] = best
            row_updates += inc
        total_updates += row_updates

        # progress prints
        if progress:
            now = time.time()
            if (k % progress_every_k == 0) and (now - last_print >= progress_interval):
                elapsed = now - t0
                rows_done = k
                rows_total = n
                pct = 100.0 * rows_done / rows_total
                # crude throughput estimate: transitions per second
                tps = total_updates / max(1e-9, elapsed)
                print(f"[karp] k={rows_done}/{rows_total}  ({pct:5.1f}%)  "
                      f"elapsed={elapsed:7.2f}s  updates≈{total_updates:,}  "
                      f"rate≈{tps:,.0f} trans/s",
                      flush=True)
                last_print = now

    # Compute μ_max on negated weights:
    mu_max_neg = -1e300
    dn = dp[n]
    for v in range(n):
        best_v = 1e300
        for k in range(n):
            denom = (n - k)
            if denom <= 0:  # skip k=n
                continue
            val = (dn[v] - dp[k][v]) / denom
            if val < best_v:
                best_v = val
        if best_v > mu_max_neg:
            mu_max_neg = best_v

    mu_min = -mu_max_neg  # minimum cycle mean on original weights

    if progress:
        elapsed = time.time() - t0
        print(f"[karp] done. μ*≈{mu_min:.12f}  total_elapsed={elapsed:.2f}s", flush=True)

    return mu_min

def main():
    ap = argparse.ArgumentParser(description="Karp min–mean on augmented Σ_b edges (with progress).")
    ap.add_argument("--edges", type=str, required=True, help="CSV with src,dst,w_aug (default col)")
    ap.add_argument("--src-col", type=str, default="src")
    ap.add_argument("--dst-col", type=str, default="dst")
    ap.add_argument("--weight-col", type=str, default="w_aug", help="column name (default: w_aug)")
    ap.add_argument("--report", action="store_true", help="print a human summary with epsilon and a bound template")
    ap.add_argument("--y0", type=int, default=None, help="optional y0 to instantiate the bound (needs epsilon>0)")
    ap.add_argument("--Ystar", type=int, default=1, help="threshold Y* for the bound")
    ap.add_argument("--Cclose", type=int, default=0, help="finite close-out constant C(Y*)")
    # Progress options
    ap.add_argument("--progress", action="store_true", help="print progress during DP")
    ap.add_argument("--progress-every-k", type=int, default=1, help="print every K rows (default: 1)")
    ap.add_argument("--progress-interval", type=float, default=0.5, help="min seconds between prints (default: 0.5)")
    args = ap.parse_args()

    n, adj_in, mapping, m_edges = read_edges(args.edges, args.src_col, args.dst_col, args.weight_col)

    if args.progress:
        print(f"[load] nodes={n}  edges={m_edges}  weight_col={args.weight_col}", flush=True)

    mu_aug = karp_min_mean_progress(
        n, adj_in,
        progress=args.progress,
        progress_every_k=max(1, args.progress_every_k),
        progress_interval=max(0.0, args.progress_interval)
    )

    print(f"[RES] mu_aug* ≈ {mu_aug:.12f}")

    if args.report:
        if mu_aug < 0:
            eps = -mu_aug
            print(f"[OK] epsilon = {eps:.12f}  (strictly positive margin)")
            if args.y0 is not None and args.y0 > 0:
                num = max(0.0, math.log2(args.y0/args.Ystar))
                m_bound = math.ceil(num / eps) + args.Cclose
                print(f"[BND] For y0={args.y0}, Y*={args.Ystar}, C={args.Cclose}:  m(y0) ≤ {m_bound}")
            print("[INFO] Generic bound:  m(y0) ≤ ceil( log2(y0/Y*) / epsilon ) + C(Y*)")
        else:
            print("[WARN] mu_aug* ≥ 0 → no uniform negative margin at this b/δ refinement.")
            print("       Try larger b or refine δ per residue and/or include potential h search.")

if __name__ == "__main__":
    main()
