#!/usr/bin/env python3
# (c) 2025 S. Gobeaux — MIT
#
# Sweep over b for the augmented seam automaton:
#  - Build Σ_b (augmented) with δ per residue
#  - Howard min–mean to get mu_aug* and ε = -mu_aug* if negative
#  - Write summary CSV with (b,n,m,mu_aug,epsilon,delta_to_log2_4_over_3,elapsed,
#    min_red,max_red, [optional files])
#
# Usage (from repo root):
#   python code/sweep_b.py --b-min 9 --b-max 12 --out data/sweep_b9_12.csv --persist --progress
#
import argparse, csv, math, os, time
from collections import deque

LOG2_3 = math.log2(3.0)
LOG2_4_OVER_3 = math.log2(4.0/3.0)   # ≈ 0.415037499279

# ---------- helpers: arithmetic ----------
def inv_mod(a: int, m: int) -> int:
    try:
        return pow(a, -1, m)
    except TypeError:
        # EEA fallback (older Python)
        t, new_t = 0, 1
        r, new_r = m, a % m
        while new_r != 0:
            q = r // new_r
            t, new_t = new_t, t - q * new_t
            r, new_r = new_r, r - q * new_r
        if r != 1:
            raise ValueError(f"{a} not invertible mod {m}")
        if t < 0: t += m
        return t

def delta_residue(r: int, mod3b: int) -> float:
    # δ_r = log2(1 + 1/(3*y_min)), with y_min = 2*D_min - 1 and D_min ≡ r (mod 3^b), D_min>=1
    D_min = r if r != 0 else mod3b
    y_min = 2 * D_min - 1
    return math.log2(1.0 + 1.0/(3.0 * y_min))

# ---------- build Σ_b (augmented) ----------
def build_sigma_augmented(b: int):
    mod = 3 ** b
    inv2 = inv_mod(2, mod)
    inv4 = (inv2 * inv2) % mod
    n = mod
    out_edges = [[] for _ in range(n)]  # entries: (v, w_aug)
    rows = []
    for r in range(mod):
        dlt = delta_residue(r, mod)
        # κ = 1
        dst1 = ((3*r - 1) * inv2) % mod
        w1 = (LOG2_3 - 1.0) + dlt
        out_edges[r].append((dst1, w1))
        rows.append((r, dst1, 1, dlt, (LOG2_3 - 1.0), w1))
        # κ = 2
        dst2 = ((3*r - 1) * inv4) % mod
        w2 = (LOG2_3 - 2.0) + dlt
        out_edges[r].append((dst2, w2))
        rows.append((r, dst2, 2, dlt, (LOG2_3 - 2.0), w2))
    return out_edges, rows

def persist_edges_csv(path: str, rows):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["src","dst","kappa","delta","w_raw","w_aug"])
        for r, dst, kappa, dlt, wraw, waug in rows:
            w.writerow([r, dst, kappa, f"{dlt:.12f}", f"{wraw:.12f}", f"{waug:.12f}"])

# ---------- Howard min–mean (policy iteration) ----------
def initial_policy(n, out_edges):
    nxt = [-1]*n; wpi = [0.0]*n
    for u in range(n):
        if not out_edges[u]:
            raise RuntimeError(f"node {u} has no outgoing edge")
        v,w = min(out_edges[u], key=lambda t: t[1])
        nxt[u] = v; wpi[u] = w
    return nxt, wpi

def decompose_policy(nxt):
    n = len(nxt); vis = [0]*n; inx = [-1]*n; cycles = []
    for s in range(n):
        if vis[s]: continue
        u = s; stack = []; pos = {}
        while not vis[u]:
            vis[u] = 1; pos[u] = len(stack); stack.append(u); u = nxt[u]
        if vis[u] == 1:
            start = pos[u]; cyc = stack[start:]; cycles.append(cyc)
            for i,x in enumerate(cyc): inx[x] = i
        for x in stack: vis[x] = 2
    return cycles, inx

def policy_mean_of_cycle(cyc, wpi):
    return sum(wpi[u] for u in cyc)/len(cyc)

def build_predecessors(nxt):
    n = len(nxt); pred = [[] for _ in range(n)]
    for u in range(n): pred[nxt[u]].append(u)
    return pred

def evaluate_policy(nxt, wpi, tol=1e-12):
    cycles, _ = decompose_policy(nxt)
    means = [policy_mean_of_cycle(c, wpi) for c in cycles]
    mu = min(means)
    # dual potential h
    n = len(nxt); h = [0.0]*n
    pred = build_predecessors(nxt)
    crit = [False]*n
    for cyc, m in zip(cycles, means):
        if abs(m - mu) <= tol:
            base = cyc[0]; h[base] = 0.0; u = base
            while True:
                v = nxt[u]; h[v] = h[u] + (wpi[u] - mu); crit[u] = True
                u = v
                if u == base: crit[u] = True; break
    dq = deque([i for i,c in enumerate(crit) if c]); seen = crit[:]
    while dq:
        v = dq.popleft()
        for u in pred[v]:
            if not seen[u]:
                h[u] = h[v] - (wpi[u] - mu); seen[u] = True; dq.append(u)
    return mu, h

def improve_policy(out_edges, nxt, wpi, mu, h, tol=1e-12):
    improved = 0
    for u, edges in enumerate(out_edges):
        best_v, best_w = nxt[u], wpi[u]
        best_red = best_w - mu + h[nxt[u]] - h[u]
        for (v,w) in edges:
            red = w - mu + h[v] - h[u]
            if red < best_red - tol:
                best_red = red; best_v, best_w = v, w
        if best_v != nxt[u]:
            nxt[u] = best_v; wpi[u] = best_w; improved += 1
    return improved

def howard_min_mean(out_edges, progress=False, tol=1e-12, max_iter=10**6):
    n = len(out_edges)
    nxt, wpi = initial_policy(n, out_edges)
    it = 0; t0 = time.time()
    while True:
        it += 1
        mu, h = evaluate_policy(nxt, wpi, tol=tol)
        if progress:
            elapsed = time.time() - t0
            print(f"[howard] iter={it}  mu≈{mu:.12f}  elapsed={elapsed:.2f}s", flush=True)
        imp = improve_policy(out_edges, nxt, wpi, mu, h, tol=tol)
        if progress:
            print(f"[howard]   improvements={imp}", flush=True)
        if imp == 0 or it >= max_iter:
            return mu, h, it, (time.time()-t0)

def reduced_cost_stats(out_edges, mu, h):
    # red(u->v) = (w - mu) + h[v] - h[u]
    min_red = 1e300
    max_red = -1e300
    for u, edges in enumerate(out_edges):
        hu = h[u]
        for (v,w) in edges:
            red = (w - mu) + (h[v] - hu)
            if red < min_red: min_red = red
            if red > max_red: max_red = red
    return min_red, max_red

# ---------- main sweep ----------
def main():
    ap = argparse.ArgumentParser(description="Sweep b for augmented Σ_b + Howard min–mean + summary CSV.")
    ap.add_argument("--b-min", type=int, default=9)
    ap.add_argument("--b-max", type=int, default=12)
    ap.add_argument("--out", type=str, required=True, help="summary CSV path")
    ap.add_argument("--persist", action="store_true", help="write edges CSV and h CSV")
    ap.add_argument("--progress", action="store_true")
    ap.add_argument("--tol", type=float, default=1e-12)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as fsum:
        wsum = csv.writer(fsum)
        wsum.writerow([
            "b","nodes","edges",
            "mu_aug","epsilon","delta_to_log2_4_over_3",
            "elapsed_s","min_red","max_red",
            "edges_csv","h_csv"
        ])

        for b in range(args.b_min, args.b_max+1):
            n = 3**b
            if args.progress:
                print(f"\n=== b={b}  (n=3^{b}={n}, edges={2*n}) ===", flush=True)
                if n > 1_000_000:
                    print("[warn] n is large; Python may be slow / memory heavy.", flush=True)

            t0 = time.time()
            out_edges, rows = build_sigma_augmented(b)
            t_build = time.time() - t0
            if args.progress:
                print(f"[build] Σ_b augmented built in {t_build:.2f}s", flush=True)

            mu, h, it, elapsed = howard_min_mean(out_edges, progress=args.progress, tol=args.tol)
            eps = -mu if mu < 0 else 0.0
            delta = LOG2_4_OVER_3 - eps
            min_red, max_red = reduced_cost_stats(out_edges, mu, h)

            edges_csv = h_csv = ""
            if args.persist:
                edges_csv = os.path.join("data", f"sigma_b{b}_aug.csv")
                h_csv     = os.path.join("data", f"h_b{b}_howard.csv")
                persist_edges_csv(edges_csv, rows)
                os.makedirs(os.path.dirname(h_csv) or ".", exist_ok=True)
                with open(h_csv, "w", newline="", encoding="utf-8") as fh:
                    wh = csv.writer(fh); wh.writerow(["node","h"])
                    for i,val in enumerate(h):
                        wh.writerow([i, f"{val:.12f}"])
                if args.progress:
                    print(f"[write] {edges_csv} ; {h_csv}", flush=True)

            wsum.writerow([
                b, n, 2*n,
                f"{mu:.12f}", f"{eps:.12f}", f"{delta:.6e}",
                f"{elapsed:.2f}", f"{min_red:.3e}", f"{max_red:.3e}",
                edges_csv, h_csv
            ])

            if args.progress:
                print(
                    f"[res] mu≈{mu:.12f}  eps≈{eps:.12f}  Δ(ε,log2(4/3))≈{delta:.2e}  "
                    f"iters={it}  time={elapsed:.2f}s  "
                    f"min_red≈{min_red:.3e}  max_red≈{max_red:.3e}",
                    flush=True
                )

    print(f"\n[OK] Summary → {args.out}")

if __name__ == "__main__":
    main()
