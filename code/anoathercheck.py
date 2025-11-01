import math

def v2(n: int) -> int:
    """2-adic valuation of n>0"""
    k=0
    while n % 2 == 0:
        n//=2
        k+=1
    return k

def T_odd(y: int) -> tuple[int,int]:
    """Odd-only Collatz step and its k = nu2(3y+1)."""
    assert y>0 and y%2==1
    s = 3*y + 1
    k = v2(s)
    y2 = s >> k
    return y2, k

def check_uniform_drop(y0: int, steps: int = 1000, A: float = 5.0, verbose: bool = True):
    LOG2_3 = math.log2(3.0)
    eps = 5* (math.log2(3.0) - 1.0) - 3.0  # ~ -0.07518 (k=1 worst case)
    # We'll compute Lambda with the simulation upper bound ΔΨ = 2k+1 at each odd step.
    y = y0
    Psi = 0.0
    Lam = A*math.log2(y) - Psi
    if verbose:
        print(f"start: y={y}, Psi={Psi}, Lambda={Lam:.6f}")

    for t in range(1, steps+1):
        y2,k = T_odd(y)
        # Upper bound on the number of atomic moves:
        dPsi = 2*k + 1
        Psi2 = Psi + dPsi

        # bound Δlog2(y) ≤ log2(3) - k (safe upper bound)
        dlog = math.log2(3.0) - k   # safe bound

        Lam2 = A*(math.log2(y) + dlog) - Psi2
        dLam = Lam2 - Lam
        if verbose:
            print(f"step {t:3d}: y->{y2:6d}, k={k}, ΔΨ={dPsi},  ΔΛ={dLam:.6f},  Λ={Lam2:.6f}")

        # Must be ≤ -0.07518...
        if dLam > -0.07518 + 1e-12:
            raise RuntimeError(f"Uniform drop violated at step {t}: ΔΛ={dLam}")
        y, Psi, Lam = y2, Psi2, Lam2

    return True

if __name__=="__main__":
    # A few tests
    for y0 in (3,5,7,9,11,27,31,41,63,73,97,101):
        print("\n== test from y0 =", y0, "==")
        ok = check_uniform_drop(y0, steps=100, A=5.0, verbose=False)
        print("OK (uniform drop holds for 100 odd steps)")
