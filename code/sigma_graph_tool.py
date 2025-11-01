
import math
import csv
from pathlib import Path
from typing import List, Tuple, Dict
import argparse
from collections import defaultdict

LOG2_3 = math.log2(3.0)

def modinv(a: int, m: int) -> int:
    return pow(a, -1, m)

def generate_sigma_edges(b: int, ks=(1,2)):
    M = 3 ** b
    inv2_cache = {k: pow(pow(2, k, M), -1, M) for k in ks}
    edges = []
    for r in range(M):
        for k in ks:
            inv2k = inv2_cache[k]
            rp = ((3 * r - 1) * inv2k) % M
            w = LOG2_3 - k
            edges.append((r, rp, w, k))
    return edges, M

def write_states(path: Path, M: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id", "residue_mod_3^b"])
        for r in range(M):
            w.writerow([r, r])

def write_edges(path: Path, edges) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["src", "dst", "weight", "k"])
        for (u, v, wt, k) in edges:
            w.writerow([u, v, f"{wt:.12f}", k])

def write_edges_min(path: Path, edges) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["src", "dst", "weight"])
        for (u, v, wt, k) in edges:
            w.writerow([u, v, f"{wt:.12f}"])

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["gen"])
    ap.add_argument("--b", type=int, required=True)
    ap.add_argument("--out-prefix", type=str, required=True)
    args = ap.parse_args()
    if args.cmd == "gen":
        edges, M = generate_sigma_edges(args.b, ks=(1,2))
        write_states(Path(args.out_prefix + "_states.csv"), M)
        write_edges(Path(args.out_prefix + "_edges.csv"), edges)
        write_edges_min(Path(args.out_prefix + "_edges.min.csv"), edges)
        print(f"[GEN] b={args.b}  M=3^b={M}")
        print(f"[GEN] wrote {args.out_prefix}_states.csv ; {args.out_prefix}_edges.csv ; {args.out_prefix}_edges.min.csv")

if __name__ == "__main__":
    main()
