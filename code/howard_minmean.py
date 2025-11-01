#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import csv
import math
import sys
import time
from collections import defaultdict, deque

def read_edges(path, src_col, dst_col, w_col):
    edges = []
    with open(path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                u = row[src_col]
                v = row[dst_col]
                w = float(row[w_col])
            except KeyError as e:
                raise SystemExit(f"Missing column in CSV: {e}. Available columns: {list(row.keys())}")
            edges.append((u, v, w))
    return edges

def build_graph(edges):
    nodes = {}
    def idx(x):
        if x not in nodes:
            nodes[x] = len(nodes)
        return nodes[x]
    E = []
    for u, v, w in edges:
        E.append((idx(u), idx(v), w))
    n = len(nodes)
    adj = [[] for _ in range(n)]
    for (u,v,w) in E:
        adj[u].append((v, w))
    return n, adj, nodes

def scc_kosaraju(n, adj):
    radj = [[] for _ in range(n)]
    for u in range(n):
        for v,_ in adj[u]:
            radj[v].append(u)
    seen = [False]*n
    order = []
    sys.setrecursionlimit(max(1000000, n*2))
    def dfs1(u):
        seen[u] = True
        for v,_ in adj[u]:
            if not seen[v]:
                dfs1(v)
        order.append(u)
    for u in range(n):
        if not seen[u]:
            dfs1(u)
    comp_id = [-1]*n
    def dfs2(u, cid):
        comp_id[u] = cid
        for v in radj[u]:
            if comp_id[v] == -1:
                dfs2(v, cid)
    cid = 0
    for u in reversed(order):
        if comp_id[u] == -1:
            dfs2(u, cid)
            cid += 1
    comps = [[] for _ in range(cid)]
    for u,c in enumerate(comp_id):
        comps[c].append(u)
    return comps

def ensure_outgoing_within_comp(adj, comp):
    comp_set = set(comp)
    out = {}
    for u in comp:
        outs = [(v,w) for (v,w) in adj[u] if v in comp_set]
        if outs:
            out[u] = outs
    valid_nodes = sorted(out.keys())
    return valid_nodes, out

def initial_policy(valid_nodes, out):
    pi = {}
    for u in valid_nodes:
        best = min(out[u], key=lambda t: (t[1], t[0]))
        pi[u] = best
    return pi

def eval_policy_mu_and_bias(valid_nodes, pi, tol):
    visited = {u: -1 for u in valid_nodes}
    mu = math.inf
    best_cycles = []
    for start in valid_nodes:
        if visited[start] != -1:
            continue
        u = start
        while visited[u] == -1:
            visited[u] = start
            u = pi[u][0]
        if visited[u] == start:
            cyc = [u]
            k = pi[u][0]
            while k != u:
                cyc.append(k)
                k = pi[k][0]
            total_w = 0.0
            for a in cyc:
                b, w = pi[a]
                total_w += w
            m = total_w / len(cyc)
            if m < mu - tol:
                mu = m
                best_cycles = [cyc]
            elif abs(m - mu) <= tol:
                best_cycles.append(cyc)
    if not best_cycles:
        mu = 0.0
        best_cycles = [[valid_nodes[0]]]
    h = {u: None for u in valid_nodes}
    crit_nodes = set()
    for cyc in best_cycles:
        rep = min(cyc)
        h[rep] = 0.0
        u = rep
        while True:
            v, w = pi[u]
            if v == rep:
                break
            h[v] = h[u] + w - mu
            u = v
        crit_nodes.update(cyc)
    rev = defaultdict(list)
    for u in valid_nodes:
        v,_w = pi[u]
        rev[v].append(u)
    dq = deque([u for u in valid_nodes if h[u] is not None])
    while dq:
        v = dq.popleft()
        for u in rev[v]:
            if h[u] is None:
                w = pi[u][1]
                h[u] = h[v] - w + mu
                dq.append(u)
    for u in valid_nodes:
        if h[u] is None:
            h[u] = 0.0
    return mu, h, crit_nodes, best_cycles

def improve_policy(valid_nodes, out, pi, mu, h, tol, mode='batch', freeze_critical=False, crit_nodes=None):
    if freeze_critical and crit_nodes is None:
        crit_nodes = set()
    best_edge = {}
    candidates = []
    for u in valid_nodes:
        if freeze_critical and u in crit_nodes:
            continue
        current_v, current_w = pi[u]
        base = - h[u]
        best = (current_v, current_w)
        min_rc = current_w - mu + h[current_v] + base
        for (v, w) in out[u]:
            rc = w - mu + h[v] + base
            if rc < min_rc - tol or (abs(rc - min_rc) <= tol and (v, w) < best):
                min_rc = rc
                best = (v, w)
        adv = min_rc - (current_w - mu + h[current_v] + base)
        best_edge[u] = best
        if adv < -tol:
            candidates.append(u)
    if not candidates:
        return False, 0
    if mode == 'single':
        u = min(candidates)
        pi[u] = best_edge[u]
        return True, 1
    elif mode == 'best':
        def advantage(u):
            v, w = best_edge[u]
            current_v, current_w = pi[u]
            return (w - mu + h[v] - h[u]) - (current_w - mu + h[current_v] - h[u])
        u = min(candidates, key=lambda x: (advantage(x), x))
        pi[u] = best_edge[u]
        return True, 1
    else:
        for u in candidates:
            pi[u] = best_edge[u]
        return True, len(candidates)

def extract_one_min_cycle(pi, crit_cycles):
    if not crit_cycles:
        return []
    crit_cycles = sorted(crit_cycles, key=lambda cyc: (len(cyc), tuple(sorted(cyc))))
    return crit_cycles[0]

def run_howard_on_component(adj, comp, args, comp_idx, num_comps):
    valid_nodes, out = ensure_outgoing_within_comp(adj, comp)

    # PATCH: only skip if there is strictly no internal edge.
    if len(valid_nodes) == 0:
        if args.progress:
            print(f"[howard][comp {comp_idx}/{num_comps} | |SCC|={len(comp)}] skipped (no internal edges)")
        return None

    # If singleton, keep only if it has a self-loop (that is a genuine cycle).
    if len(valid_nodes) == 1 and len(comp) == 1:
        u = valid_nodes[0]
        has_self = any(v == u for (v, _w) in out.get(u, []))
        if not has_self:
            if args.progress:
                print(f"[howard][comp {comp_idx}/{num_comps} | |SCC|=1] skipped (no self-loop)")
            return None

    pi = initial_policy(valid_nodes, out)
    stall_hist = deque(maxlen=max(10, min(200, args.stall_iters//2 if args.stall_iters else 50)))
    improve_mode = args.improve
    t0 = time.time()
    it = 0
    while True:
        it += 1
        mu, h, crit_nodes, crit_cycles = eval_policy_mu_and_bias(valid_nodes, pi, args.tol)
        if args.progress:
            elapsed = time.time() - t0
            print(f"[howard][comp {comp_idx}/{num_comps} | |SCC|={len(comp)}] iter={it}  mu={mu}  t={elapsed:.1f}s")
        changed, switched = improve_policy(valid_nodes, out, pi, mu, h, args.tol,
                                           mode=improve_mode,
                                           freeze_critical=args.freeze_critical,
                                           crit_nodes=crit_nodes)
        if not changed:
            break
        if args.stall_iters and it % args.stall_iters == 0 and improve_mode != 'single':
            uniq = set(round(x, 12) for x in stall_hist) if stall_hist else set()
            if len(uniq) <= 3:
                improve_mode = 'single'
        stall_hist.append(mu)
        if args.max_iters_per_comp and it >= args.max_iters_per_comp:
            break
    t_total = time.time() - t0
    if args.progress:
        cyc = extract_one_min_cycle(pi, crit_cycles)
        print(f"[howard][comp {comp_idx}/{num_comps}] done: iters={it}, mu*={mu}, |cycle|={len(cyc)}, elapsed={t_total:.1f}s")
    return {
        "mu": mu,
        "iters": it,
        "cycle": extract_one_min_cycle(pi, crit_cycles),
        "policy": pi,
        "crit_nodes": crit_nodes,
        "h": h,
        "valid_nodes": valid_nodes,
        "comp": comp,
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", help="input CSV path")
    ap.add_argument("--src-col", default="src", help="source column name (default: src)")
    ap.add_argument("--dst-col", default="dst", help="dest column name (default: dst)")
    ap.add_argument("--weight-col", default="weight", help="weight column name (default: weight)")
    ap.add_argument("--progress", action="store_true", help="print per-iteration logs")
    ap.add_argument("--stall-iters", type=int, default=0, help="iterations before auto anti-cycling fallback (0=off)")
    ap.add_argument("--max-iters-per-comp", type=int, default=0, help="hard cap on iterations per SCC (0=off)")
    ap.add_argument("--tol", type=float, default=1e-12, help="numerical tolerance for tie/zero checks")
    ap.add_argument("--improve", choices=["batch","single","best"], default="batch",
                    help="policy improvement mode (batch: switch all, single: switch one, best: switch the most negative)")
    ap.add_argument("--freeze-critical", action="store_true",
                    help="never switch states that lie on cycles currently achieving mu (degeneracy fix)")
    ap.add_argument("--export-potentials", type=str, default="",
                    help="optional path to write potentials CSV (id,h)")
    args = ap.parse_args()

    edges = read_edges(args.csv, args.src_col, args.dst_col, args.weight_col)
    n, adj, nodes = build_graph(edges)
    comps = [c for c in scc_kosaraju(n, adj) if len(c) >= 1]
    comps.sort(key=lambda c: (-len(c), c[0]))
    results = []
    for i, comp in enumerate(comps, 1):
        res = run_howard_on_component(adj, comp, args, i, len(comps))
        if res is not None:
            results.append(res)
    if not results:
        print("No result.")
        return
    overall_mu = min(r["mu"] for r in results)
    print(f"global mu* (min over SCCs) = {overall_mu}")

    # Optional export of potentials (assemble per-node)
    if args.export_potentials:
        H = [0.0]*n
        for res in results:
            for u in res["valid_nodes"]:
                H[u] = res["h"][u]
        idx_to_label = {i:s for s,i in nodes.items()}
        with open(args.export_potentials, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["id","h"])
            for i in range(n):
                w.writerow([idx_to_label[i], f"{H[i]:.12f}"])
        print(f"[write] potentials → {args.export_potentials}")

if __name__ == "__main__":
    main()
