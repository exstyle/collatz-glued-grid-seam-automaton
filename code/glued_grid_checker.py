#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Glued-grid / seam-automaton checker for the odd-only Collatz map.

Ce script vérifie, pour un grand échantillon :
  1) Identités "LH + SEAM" exactes :
     - T(y) = oddize(3D - 1) avec D = (y+1)//2
     - psi(T(y)) = Sigma(D) = (oddize(3D - 1) + 1)//2
  2) Décomposition de k = nu2(3y+1) en coutures κ ∈ {1,2}
     - somme(κ) = k
     - S = total coutures sur m pas, S ≥ m et ceil(k/2) ≤ S ≤ k pas-à-pas
  3) "Size split" : Δlog2 y = (log2 3 - k) + δ(y), avec 0 < δ(y) ≤ log2(4/3)
  4) Automate Σ_b (mod 3^b), arêtes κ=1,2 avec poids w = log2(3) - κ
     - Existence boucle κ=2 en r ≡ -1 (mod 3^b)
     - Min–moyenne (Karp) = log2(3) - 2 (confirm.)
"""

from math import log
from collections import defaultdict

# ------------------------------
# Arithmétique de base
# ------------------------------

def nu2(n: int) -> int:
    """2-adic valuation: highest e s.t. 2^e divides n (n>0)."""
    if n <= 0:
        raise ValueError("nu2 requires n>0")
    e = 0
    while n % 2 == 0:
        n //= 2
        e += 1
    return e

def oddize(n: int) -> int:
    """Remove all factors of 2."""
    return n // (1 << nu2(n)) if n > 0 else 0

def T_odd(y: int) -> int:
    """Odd-only Collatz step T(y) = oddize(3y+1), y odd > 0."""
    assert y > 0 and (y % 2 == 1)
    return oddize(3*y + 1)

def psi(y: int) -> int:
    """Distance labeling: D = (y+1)//2."""
    assert y > 0 and (y % 2 == 1)
    return (y + 1) // 2

def Sigma(D: int) -> int:
    """Seam target: Σ(D) = (oddize(3D-1) + 1)//2, D>=1."""
    assert D >= 1
    return (oddize(3*D - 1) + 1) // 2

# ------------------------------
# 1) Identités LH + SEAM
# ------------------------------

def check_two_step_identity(Ymax: int = 100000) -> None:
    """Vérifie T(y) = oddize(3D-1) et psi(T(y)) = Sigma(D) sur y<=Ymax impairs."""
    for y in range(1, Ymax+1, 2):
        D = psi(y)
        lhs = T_odd(y)
        rhs = oddize(3*D - 1)
        if lhs != rhs:
            raise AssertionError(f"Fail T(y)=oddize(3D-1) at y={y} (lhs={lhs}, rhs={rhs})")
        if psi(lhs) != Sigma(D):
            raise AssertionError(f"Fail psi(T(y))=Sigma(D) at y={y}")
    print(f"[OK] Two-step identities checked for all odd y<= {Ymax}")

# ------------------------------
# 2) Décomposition en coutures κ∈{1,2}
# ------------------------------

def seam_kappas(k: int):
    """
    Décompose k = sum κi avec κi ∈ {1,2}.
    Choix canonique: le plus de 2 possible, puis 1 si k impair.
    """
    if k <= 0:
        return []
    ks = []
    while k >= 2:
        ks.append(2)
        k -= 2
    if k == 1:
        ks.append(1)
    return ks

def simulate_block_and_count_seams(y0: int, steps: int):
    """
    Simule 'steps' pas impairs, collecte:
      - la liste k_i = nu2(3y_i+1),
      - la décomposition en κ ∈ {1,2} pour chaque pas,
      - S total (nombre de coutures),
      - la somme des κ (doit valoir sum k_i),
      - vérifie S >= steps.
    """
    assert y0 > 0 and y0 % 2 == 1
    y = y0
    Ks = []
    kappas = []
    for _ in range(steps):
        k = nu2(3*y + 1)
        Ks.append(k)
        kappas.extend(seam_kappas(k))
        y = T_odd(y)
    S = len(kappas)
    assert S >= steps, f"Chaque pas impair doit avoir au moins une couture: S={S}, steps={steps}"
    assert sum(kappas) == sum(Ks), f"Somme κ ({sum(kappas)}) != somme k_i ({sum(Ks)})"
    return Ks, kappas, S

# ------------------------------
# 3) Size split: borne δ(y)
# ------------------------------

LOG2 = log(2.0)
LOG2_4_over_3 = log(4.0/3.0) / LOG2  # ~ 0.4150

def delta_of(y: int) -> float:
    """δ(y) = log2(1 + 1/(3y))."""
    return log(1.0 + 1.0/(3.0*y)) / LOG2

def check_size_split_bounds(Ymax: int = 100000) -> None:
    for y in range(1, Ymax+1, 2):
        d = delta_of(y)
        if not (0.0 < d <= LOG2_4_over_3 + 1e-15):
            raise AssertionError(f"delta bound failed at y={y}, δ={d}")
    print(f"[OK] δ(y) ∈ (0, log2(4/3)] checked for odd y<= {Ymax}")

# ------------------------------
# 4) Automate Σ_b, arêtes κ=1,2 et Karp min–moyenne
# ------------------------------

def egcd(a: int, b: int):
    """Extended GCD: return (g,x,y) with ax+by=g."""
    if b == 0:
        return (a, 1, 0)
    g, x1, y1 = egcd(b, a % b)
    return (g, y1, x1 - (a // b) * y1)

def inv_mod(a: int, m: int) -> int:
    """Inverse de a modulo m (suppose gcd(a,m)=1)."""
    g, x, _ = egcd(a, m)
    if g != 1:
        raise ZeroDivisionError("No inverse")
    return x % m

def build_sigma_graph(b: int):
    """
    Construit l’automate Σ_b (mod 3^b).
    Noeuds: r = 0..M-1 (M=3^b).
    Arêtes: pour κ=1 et κ=2, r -> r' = (3r-1) * inv_mod(2^κ, M).
    Poids: w = log2(3) - κ.
    Retourne: (M, edges) avec edges[r] = [(r', w, κ), ...]
    """
    M = pow(3, b)
    inv2 = inv_mod(2 % M, M)
    inv4 = (inv2 * inv2) % M
    edges = [[] for _ in range(M)]
    w1 = log(3.0)/LOG2 - 1.0
    w2 = log(3.0)/LOG2 - 2.0
    for r in range(M):
        r1 = ((3*r - 1) * inv2) % M
        r2 = ((3*r - 1) * inv4) % M
        edges[r].append((r1, w1, 1))
        edges[r].append((r2, w2, 2))
    return M, edges

def has_kappa2_selfloop_minus1(b: int) -> bool:
    """Vérifie la boucle κ=2 en r≡-1(mod 3^b)."""
    M = pow(3, b)
    r = (M - 1) % M  # -1 mod M
    _, edges = build_sigma_graph(b)
    for (r2, w, kappa) in edges[r]:
        if r2 == r and kappa == 2:
            return abs(w - (log(3.0)/LOG2 - 2.0)) < 1e-12
    return False

def karp_min_mean(M: int, edges):
    """
    Karp min–mean cycle algorithm.
    - Nodes: 0..M-1
    - edges[u] = list of (v, weight, kappa)
    Retourne: (mu_star, cycle_example)
    """
    INF = 1e100
    # dp[k][v] = min weight of a path of length k ending at v
    dp = [[INF]*M for _ in range(M+1)]
    # Init: zero-length paths cost 0
    for v in range(M):
        dp[0][v] = 0.0
    # Relax M times
    for k in range(1, M+1):
        for u in range(M):
            best = dp[k][u]
            # Predecessors? We'll do forward relaxation:
            # for u->v, dp[k][v] = min(dp[k][v], dp[k-1][u] + w)
            for (v, w, _) in edges[u]:
                if dp[k-1][u] + w < dp[k][v]:
                    dp[k][v] = dp[k-1][u] + w
    # Compute min mean via Karp formula
    mu_star = float('inf')
    for v in range(M):
        max_avg = -float('inf')
        for k in range(M):
            avg = (dp[M][v] - dp[k][v]) / (M - k) if dp[k][v] < INF else -float('inf')
            if avg > max_avg:
                max_avg = avg
        if max_avg < mu_star:
            mu_star = max_avg
    return mu_star

def check_sigma_b_min_mean(b_list=(3,4,5,6,7,8,9)) -> None:
    """Confirme mu* = log2(3) - 2 pour plusieurs b (3^b nœuds, 2 arêtes/noeud)."""
    target = log(3.0)/LOG2 - 2.0
    for b in b_list:
        M, edges = build_sigma_graph(b)
        ok_loop = has_kappa2_selfloop_minus1(b)
        mu = karp_min_mean(M, edges)
        print(f"[b={b}] mu*≈{mu:.12f}  (loop κ=2 @ -1: {ok_loop})")
        if abs(mu - target) > 1e-8:
            raise AssertionError(f"mu* differs from log2(3)-2 at b={b}")

# ------------------------------
# 5) Démo "bloc" (non-divergence empirique)
# ------------------------------

def block_delta_log2(y0: int, steps: int) -> float:
    """
    Calcule sum Δlog2 y_i = sum(log2((3y_i+1)/(2^{k_i} y_i))).
    Retourne la somme + vérifie δ(y_i) ≤ log2(4/3).
    """
    y = y0
    s = 0.0
    for _ in range(steps):
        k = nu2(3*y + 1)
        delta = delta_of(y)
        assert 0.0 < delta <= LOG2_4_over_3 + 1e-15
        s += (log(3.0)/LOG2 - k) + delta
        y = T_odd(y)
    return s

# ------------------------------
# 6) Main tests (tu peux ajuster les bornes)
# ------------------------------

def main():
    print("== Vérifs algébriques LH+SEAM ==")
    check_two_step_identity(Ymax=200000)

    print("\n== Size split (borne δ) ==")
    check_size_split_bounds(Ymax=200000)

    print("\n== Décomposition κ ∈ {1,2} & S≥m ==")
    Ks, kappas, S = simulate_block_and_count_seams(y0=41, steps=2000)
    print(f"Exemple: steps=2000, S={S}, sum(kappas)={sum(kappas)}, sum(Ks)={sum(Ks)}")

    print("\n== Automate Σ_b : boucle κ=2 @ r≡-1 & mu* ==")
    check_sigma_b_min_mean(b_list=(3,4,5,6,7))

    print("\n== Démo bloc Δlog2 ≤ 0 (quelques exemples) ==")
    for y0 in (1, 3, 5, 41, 97, 1001):
        val = block_delta_log2(y0=y0, steps=2000)
        print(f"y0={y0:>4} -> sum Δlog2 over 2000 steps ≈ {val:.6f}")

if __name__ == "__main__":
    main()
