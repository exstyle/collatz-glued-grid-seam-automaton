#!/usr/bin/env python3
# dual_bellman_ford_certificate.py
# usage:
#   python dual_bellman_ford_certificate.py edges.csv --b 8 --scale 48 --out cert.json
import argparse, csv, json, math
from collections import defaultdict

def mu_b(b: int) -> float:
    # μ_b = log2( (3/4) * (1 + 1/(6*3^b - 9)) )  < 0
    den = 6*(3**b) - 9
    return math.log2((3.0/4.0) * (1.0 + 1.0/den))

def y_min_of_residue(r: int, b: int) -> int:
    # D_min(r) = r for r∈[1,3^b-1], and D_min(0) = 3^b
    if r == 0:
        return 2*(3**b) - 1
    return 2*r - 1

def ceil_fixed(x: float, scale: int) -> int:
    # borne sup int pour 2^Q * x, avec marge de sécurité 2^-Q
    # (pour garantir "upper bound" malgré l’erreur float)
    eps = 2.0**(-scale)
    return math.ceil((x + eps) * (1 << scale))

def load_edges(path):
    edges = []
    with open(path, newline='', encoding='utf-8') as f:
        rd = csv.DictReader(f)
        # auto-map colonnes
        def pick(*names):
            for n in names:
                if n in rd.fieldnames: return n
            raise ValueError(f"Colonnes manquantes. Vu: {rd.fieldnames}")
        c_r  = pick('r','src','from','u')
        c_k  = pick('kappa','k','type')
        c_rp = pick('rp','r2','dst','to','v')
        for row in rd:
            r  = int(row[c_r])
            k  = int(row[c_k])
            rp = int(row[c_rp])
            edges.append((r,k,rp))
    return edges

def build_ce(edges, b: int, Q: int):
    mu = mu_b(b)             # < 0
    log2_3 = math.log2(3.0)
    ce = []
    nodes = set()
    for (r,k,rp) in edges:
        nodes.add(r); nodes.add(rp)
        y_min = y_min_of_residue(r, b)
        delta = math.log2(1.0 + 1.0/(3.0*y_min))  # δ(r) = log2(1 + 1/(3 y_min(r)))
        w_aug = (log2_3 - k) + delta              # w_aug(r->r')
        c = ceil_fixed(w_aug - mu, Q)             # c_e = ceil( 2^Q * (w_aug - μ) )
        ce.append((r, rp, c))
    return sorted(nodes), ce, mu

def bellman_ford_difference_constraints(nodes, edges_ub):
    # edges_ub: list of (u,v,c) meaning  h[v] - h[u] <= c
    # super-source trick: init all h=0, relax |V|-1 times, detect negative cycle by any further improv.
    idx = {v:i for i,v in enumerate(nodes)}
    n = len(nodes)
    h = [0]*(n)  # all zeros is a valid starting point
    for _ in range(n-1):
        changed = False
        for (u,v,c) in edges_ub:
            iu, iv = idx[u], idx[v]
            if h[iv] > h[iu] + c:
                h[iv] = h[iu] + c
                changed = True
        if not changed:
            break
    # detect negative cycle (infeasible system)
    for (u,v,c) in edges_ub:
        iu, iv = idx[u], idx[v]
        if h[iv] > h[iu] + c:
            raise RuntimeError("Infeasible constraints: negative cycle detected with upper bounds.")
    # normalize to make min(h)=0 for readability
    mn = min(h)
    h = [x - mn for x in h]
    return {nodes[i]: h[i] for i in range(n)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", help="edges of Σ_b as CSV with columns r,kappa,rp (or src,kappa,dst)")
    ap.add_argument("--b", type=int, required=True, help="window parameter b (mod 3^b)")
    ap.add_argument("--scale", type=int, default=48, help="dyadic scale Q (default 48)")
    ap.add_argument("--out", type=str, default="certificate.json", help="output JSON")
    args = ap.parse_args()

    edges = load_edges(args.csv)
    nodes, edges_ub, mu = build_ce(edges, args.b, args.scale)
    h = bellman_ford_difference_constraints(nodes, edges_ub)

    # oscillation and L_b
    osc_int = max(h.values()) - min(h.values())
    eps = -mu                                 # ε_b = -μ_b > 0
    L_b = int(math.floor( (osc_int / (1<<args.scale)) / eps ))

    out = {
        "b": args.b,
        "scale": args.scale,
        "mu_b": mu,
        "eps_b": eps,
        "L_b": L_b,
        "osc_int": osc_int,
        "osc": osc_int / float(1<<args.scale),
        "node_count": len(nodes),
        "edge_count": len(edges),
        "h": h,                                # integers in dyadic scale
        "hash_hint": {
            "csv": args.csv
        }
    }
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print(f"[OK] Wrote {args.out}  (|V|={len(nodes)}, |E|={len(edges)}, L_b≈{L_b})")

if __name__ == "__main__":
    main()
