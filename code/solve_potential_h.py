#!/usr/bin/env python3
# (c) 2025 S. Gobeaux — MIT
#
# Solve for a dual potential h on the augmented seam automaton:
#  1) Read edges (CSV) with columns: src,dst,w_aug (default names configurable).
#  2) Compute μ_aug* (minimum cycle mean) with Karp.
#  3) Build reduced weights: w_red(u->v) = w_aug(u->v) - μ_aug*.
#  4) Solve difference constraints h(v) - h(u) ≤ - w_red(u->v) via Bellman-Ford
#     (add super-source s -> every node with weight 0), returning one feasible h.
#  5) Check max slack:  max_e [ w_red(e) + h(u) - h(v) ]  (should be ≤ tiny ε_num).
#
# If μ_aug* < 0, then ε := -μ_aug* > 0 gives a UNIFORM MARGIN, hence the O(log y0) upper bound
# m(y0) ≤ ceil( log2(y0/Y*) / ε ) + C(Y*).
#
import argparse, csv, math, sys
from typing import List, Tuple, Dict

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
    # incoming adjacency for Karp
    adj_in = [[] for _ in range(n)]
    # forward edges for Bellman-Ford later
    adj_fw = [[] for _ in range(n)]
    for (u, v, w) in edges:
        uu = mapping[u]; vv = mapping[v]
        adj_in[vv].append((uu, w))
        adj_fw[uu].append((vv, w))
    return n, adj_in, adj_fw, mapping

def karp_min_mean(n: int, adj_in: List[List[Tuple[int,float]]]) -> float:
    NEG_INF = -1e300
    dp = [[NEG_INF]*n for _ in range(n+1)]
    for v in range(n):
        dp[0][v] = 0.0
    for k in range(1, n+1):
        for v in range(n):
            best = NEG_INF
            for (u, w) in adj_in[v]:
                val = dp[k-1][u] + (-w)  # maximize on negated weights
                if val > best:
                    best = val
            dp[k][v] = best
    mu_max_neg = -1e300
    for v in range(n):
        best_v = 1e300
        dn = dp[n][v]
        for k in range(n):
            denom = (n - k)
            if denom <= 0: continue
            val = (dn - dp[k][v]) / denom
            if val < best_v:
                best_v = val
        if best_v > mu_max_neg:
            mu_max_neg = best_v
    mu_min = -mu_max_neg
    return mu_min

def bellman_ford_potential(n: int, adj_fw: List[List[Tuple[int,float]]], w_red: List[List[Tuple[int,float]]], tol=1e-12):
    # Build a graph with a super-source s connected to all nodes with 0-weight edges.
    # Solve: h[v] ≤ h[u] + (-w_red[u->v])  i.e., h[v] - h[u] ≤ -w_red[u->v]
    # We'll compute shortest paths with edges cost c(u->v) = -w_red[u->v].
    INF = 1e300
    h = [0.0] * n   # init with 0 (super-source trick)
    # Relax n-1 times
    for _ in range(n-1):
        updated = False
        for u in range(n):
            hu = h[u]
            for (v, w) in w_red[u]:
                c = -w  # cost = - w_red
                if h[v] > hu + c:
                    h[v] = hu + c
                    updated = True
        if not updated:
            break
    # Check for negative cycle in reduced graph (should not happen if μ is min cycle mean)
    for u in range(n):
        hu = h[u]
        for (v, w) in w_red[u]:
            c = -w
            if h[v] > hu + c - tol:
                raise SystemExit("[ERR] Negative cycle detected in reduced graph (numerical issue?)")
    # Normalize h (optional: shift so min(h)=0)
    mn = min(h); h = [x - mn for x in h]
    return h

def max_slack(n: int, w_red: List[List[Tuple[int,float]]], h: List[float]) -> float:
    # slack(e) = w_red(u->v) + h[u] - h[v]
    mx = -1e300
    for u in range(n):
        hu = h[u]
        for (v, w) in w_red[u]:
            s = w + hu - h[v]
            if s > mx:
                mx = s
    return mx

def main():
    ap = argparse.ArgumentParser(description="Solve dual potential h for augmented min–mean; verify edge slacks.")
    ap.add_argument("--edges", type=str, required=True, help="CSV with src,dst,w_aug by default")
    ap.add_argument("--src-col", type=str, default="src")
    ap.add_argument("--dst-col", type=str, default="dst")
    ap.add_argument("--weight-col", type=str, default="w_aug")
    ap.add_argument("--out-h", type=str, default=None, help="optional CSV output for potentials h (node, h)")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--tol", type=float, default=1e-10, help="tolerance for slack check")
    args = ap.parse_args()

    n, adj_in, adj_fw, mapping = read_edges(args.edges, args.src_col, args.dst_col, args.weight_col)
    m_edges = sum(len(L) for L in adj_fw)
    mu_aug = karp_min_mean(n, adj_in)

    # Build reduced weights w_red = w_aug - mu_aug
    w_red = [[] for _ in range(n)]
    # We need to re-read edges to get the weights again in forward form
    # (could be avoided by carrying both structures from read_edges)
    with open(args.edges, "r", encoding="utf-8") as f:
        rd = csv.DictReader(f)
        for row in rd:
            u = int(row[args.src_col]); v = int(row[args.dst_col])
            w = float(row[args.weight_col])
            uu = mapping[u]; vv = mapping[v]
            w_red[uu].append((vv, w - mu_aug))

    # Solve for h with Bellman-Ford on reduced graph (cost = -w_red)
    h = bellman_ford_potential(n, adj_fw, w_red, tol=args.tol)

    # Check max slack (should be ≤ ~0)
    slack = max_slack(n, w_red, h)

    print(f"[OK] edges read: nodes={n}  edges={m_edges}  weight_col={args.weight_col}")
    print(f"[RES] mu_aug* ≈ {mu_aug:.12f}")
    if mu_aug < 0:
        eps = -mu_aug
        print(f"[OK] epsilon = {-mu_aug:.12f}  (strict uniform margin)")
        print("     ⇒ upper bound: m(y0) ≤ ceil( log2(y0/Y*) / epsilon ) + C(Y*)")
    else:
        print("[WARN] mu_aug* ≥ 0 → no strict margin (at this b / δ). Try larger b or refine δ/h.")

    print(f"[CHK] max edge slack ≈ {slack:.3e} (≤ ~1e-10 expected)")

    if args.out_h:
        # invert mapping to get original node labels if needed
        inv_map = {i:node for (node,i) in mapping.items()}
        with open(args.out_h, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["node","h"])
            for i, val in enumerate(h):
                w.writerow([inv_map[i], f"{val:.12f}"])
        print(f"[WROTE] potentials h → {args.out_h}")

if __name__ == "__main__":
    main()
