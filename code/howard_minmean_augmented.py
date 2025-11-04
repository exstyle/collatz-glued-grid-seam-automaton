#!/usr/bin/env python3
# (c) 2025 S. Gobeaux — MIT
#
# Minimum mean cycle via Howard's policy iteration (augmented weights).
# - Input CSV (from build_sigma_augmented.py): src,dst,kappa,delta,w_raw,w_aug
# - Uses w_aug by default. Outputs mu_aug*, a dual h (optional), and progress.
#
# Usage:
#   python code/howard_minmean_augmented.py --edges data/sigma_b10_aug.csv --report --progress
#   python code/howard_minmean_augmented.py --edges data/sigma_b10_aug.csv --out-h data/h_b10_howard.csv
#
import argparse, csv, math, time
from collections import defaultdict, deque

def read_edges(path, src_col="src", dst_col="dst", weight_col="w_aug"):
    with open(path, "r", encoding="utf-8") as f:
        rd = csv.DictReader(f)
        rows = list(rd)
    nodes = set()
    Edges = []
    for r in rows:
        u = int(r[src_col]); v = int(r[dst_col])
        try:
            w = float(r[weight_col])
        except Exception:
            w = float(r.get("weight", r.get("w", r.get("w_aug"))))
        nodes.add(u); nodes.add(v)
        Edges.append((u, v, w))
    nodes = sorted(nodes)
    idx = {x:i for i,x in enumerate(nodes)}
    n = len(nodes)
    out_edges = [[] for _ in range(n)]
    for (u,v,w) in Edges:
        out_edges[idx[u]].append((idx[v], w))
    return n, nodes, out_edges

def initial_policy(n, out_edges):
    # choose cheapest outgoing edge per node
    nxt = [-1]*n
    wpi = [0.0]*n
    for u in range(n):
        if not out_edges[u]:
            raise SystemExit(f"[ERR] node {u} has no outgoing edge")
        v,w = min(out_edges[u], key=lambda t: t[1])
        nxt[u] = v; wpi[u] = w
    return nxt, wpi

def decompose_policy(nxt):
    # return cycles (as lists of nodes), plus per-node index in its cycle (-1 if not on a cycle)
    n = len(nxt)
    vis = [0]*n        # 0:unseen, 1:in-stack, 2:done
    inx = [-1]*n
    cycles = []
    for s in range(n):
        if vis[s]: continue
        u = s
        stack = []
        pos = {}
        while not vis[u]:
            vis[u] = 1
            pos[u] = len(stack)
            stack.append(u)
            u = nxt[u]
        if vis[u] == 1:
            # found a cycle: from pos[u] to end-1
            start = pos[u]
            cyc = stack[start:]
            cycles.append(cyc)
            for i,x in enumerate(cyc):
                inx[x] = i
        # mark path done
        for x in stack:
            vis[x] = 2
    return cycles, inx

def policy_mean_of_cycle(cyc, nxt, wpi):
    s = 0.0
    for u in cyc:
        s += wpi[u]
    return s / len(cyc)

def build_predecessors(nxt):
    n = len(nxt)
    pred = [[] for _ in range(n)]
    for u in range(n):
        pred[nxt[u]].append(u)
    return pred

def evaluate_policy(nxt, wpi, tol=1e-12):
    # 1) cycles in policy and their means
    cycles, inx = decompose_policy(nxt)
    means = [policy_mean_of_cycle(c, nxt, wpi) for c in cycles]
    mu = min(means)
    # 2) potentials h:
    #    set equality along all cycles with mean ~ mu (critical cycles),
    #    then propagate to trees feeding them: h[u] = h[v] - (wpi[u] - mu), where v = nxt[u]
    n = len(nxt)
    h = [0.0]*n
    pred = build_predecessors(nxt)
    # init queue with all nodes on critical cycles; set their h by walking the cycle
    crit_flags = [False]*n
    for cyc, m in zip(cycles, means):
        if abs(m - mu) <= tol:
            # fix h on this cycle by equality
            base = cyc[0]
            h[base] = 0.0
            u = base
            while True:
                v = nxt[u]
                h[v] = h[u] + (wpi[u] - mu)  # equality on cycle edges
                crit_flags[u] = True
                u = v
                if u == base:
                    crit_flags[u] = True
                    break
    # now BFS backwards from all critical nodes
    dq = deque([u for u in range(n) if crit_flags[u]])
    seen = [False]*n
    for u in dq:
        seen[u] = True
    while dq:
        v = dq.popleft()
        for u in pred[v]:
            if not seen[u]:
                # tree edge u -> v on the policy graph
                h[u] = h[v] - (wpi[u] - mu)
                seen[u] = True
                dq.append(u)
    return mu, h

def improve_policy(n, out_edges, nxt, wpi, mu, h, tol=1e-12):
    improved = 0
    for u in range(n):
        best_v, best_w = nxt[u], wpi[u]
        # reduced cost = w - mu + h[v] - h[u]
        best_red = best_w - mu + h[nxt[u]] - h[u]
        for (v,w) in out_edges[u]:
            red = w - mu + h[v] - h[u]
            if red < best_red - tol:
                best_red = red
                best_v, best_w = v, w
        if best_v != nxt[u]:
            nxt[u] = best_v
            wpi[u] = best_w
            improved += 1
    return improved

def howard_min_mean(n, out_edges, progress=False, tol=1e-12, max_iter=10**6):
    nxt, wpi = initial_policy(n, out_edges)
    it = 0
    t0 = time.time()
    mu = float("inf")
    while True:
        it += 1
        mu, h = evaluate_policy(nxt, wpi, tol=tol)
        if progress:
            elapsed = time.time() - t0
            print(f"[howard] iter={it}  mu≈{mu:.12f}  elapsed={elapsed:.2f}s", flush=True)
        imp = improve_policy(n, out_edges, nxt, wpi, mu, h, tol=tol)
        if progress:
            print(f"[howard]   improvements={imp}", flush=True)
        if imp == 0 or it >= max_iter:
            return mu, h, it, (time.time() - t0)

def main():
    ap = argparse.ArgumentParser(description="Howard min–mean on augmented Σ_b (w_aug).")
    ap.add_argument("--edges", required=True, help="CSV with src,dst,w_aug")
    ap.add_argument("--src-col", default="src")
    ap.add_argument("--dst-col", default="dst")
    ap.add_argument("--weight-col", default="w_aug")
    ap.add_argument("--progress", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--out-h", default=None, help="optional CSV (node,h)")
    ap.add_argument("--y0", type=int, default=None)
    ap.add_argument("--Ystar", type=int, default=1)
    ap.add_argument("--Cclose", type=int, default=0)
    ap.add_argument("--tol", type=float, default=1e-12)
    args = ap.parse_args()

    n, nodes, out_edges = read_edges(args.edges, args.src_col, args.dst_col, args.weight_col)
    if args.progress:
        print(f"[load] nodes={n}  edges={sum(len(L) for L in out_edges)}", flush=True)
    mu, h, it, elapsed = howard_min_mean(n, out_edges, progress=args.progress, tol=args.tol)
    print(f"[RES] mu_aug* ≈ {mu:.12f}  iters={it}  elapsed={elapsed:.2f}s")

    if args.report:
        if mu < 0:
            eps = -mu
            print(f"[OK] epsilon = {eps:.12f}")
            if args.y0:
                num = max(0.0, math.log2(args.y0/args.Ystar))
                m_bound = math.ceil(num / eps) + args.Cclose
                print(f"[BND] For y0={args.y0}, Y*={args.Ystar}, C={args.Cclose}:  m(y0) ≤ {m_bound}")
            print("[INFO] Bound form:  m(y0) ≤ ceil( log2(y0/Y*) / epsilon ) + C(Y*)")
        else:
            print("[WARN] mu_aug* ≥ 0 → no uniform margin (try larger b or refine δ/h).")

    if args.out_h:
        # --- quick slack check on reduced costs ---
        # reduced(u->v) = (w_aug - mu) + h[v] - h[u]  should be >= 0  (dual feasibility)
        mx = -1e300
        import csv as _csv
        with open(args.edges, "r", encoding="utf-8") as _f:
            _rd = _csv.DictReader(_f)
            for _r in _rd:
                u = int(_r[args.src_col]); v = int(_r[args.dst_col])
                w = float(_r[args.weight_col])
                uu = nodes.index(u); vv = nodes.index(v)
                red = (w - mu) + (h[vv] - h[uu])
                if red > mx: mx = red
        print(f"[CHK] max reduced cost (should be >= 0) ≈ {mx:.3e}")

        with open(args.out_h, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f); w.writerow(["node","h"])
            for i,val in enumerate(h):
                w.writerow([nodes[i], f"{val:.12f}"])
        print(f"[WROTE] potentials h → {args.out_h}")

if __name__ == "__main__":
    main()
