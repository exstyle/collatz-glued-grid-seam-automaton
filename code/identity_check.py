import argparse

def nu2(n: int) -> int:
    c = 0
    while n % 2 == 0:
        n //= 2
        c += 1
    return c

def oddize(n: int) -> int:
    return n >> nu2(n)

def T(y: int) -> int:
    assert y % 2 == 1 and y > 0
    return oddize(3*y + 1)

def psi(y: int) -> int:
    return (y + 1) // 2

def Sigma(D: int) -> int:
    return (oddize(3*D - 1) + 1) // 2

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-odd", type=int, default=100000)
    args = ap.parse_args()

    # identities ψ∘T = Σ∘ψ
    for y in range(1, args.max_odd+1, 2):
        D = psi(y)
        if psi(T(y)) != Sigma(D):
            raise SystemExit(f"Mismatch at y={y}")
    print(f"[OK] ψ(T(y)) = Σ(ψ(y)) for all odd y <= {args.max_odd}")

    # δ(y) bound: log2((3y+1)/(2^k y)) = log2(3)-k + log2(1+1/(3y))
    # Here we just check the small correction upper bound numerically.
    import math
    worst = 0.0
    for y in range(1, args.max_odd+1, 2):
        k = nu2(3*y + 1)
        delta = math.log2(1 + 1/(3*y))
        worst = max(worst, delta)
        if not (0.0 < delta <= math.log2(4/3) + 1e-12):
            raise SystemExit(f"delta bound failed at y={y}")
        if k < 1:
            raise SystemExit(f"k>=1 violated at y={y}")
    print(f"[OK] δ(y) ∈ (0, log2(4/3)] checked up to odd y = {args.max_odd}; worst≈{worst:.9f}")

if __name__ == "__main__":
    main()
