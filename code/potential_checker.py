import csv, math

def load_h(path):
    h=[]
    with open(path, newline='', encoding='utf-8') as f:
        rd=csv.reader(f); next(rd,None)
        for row in rd:
            i=int(row[0]); val=float(row[1])
            if i>=len(h): h.extend([0.0]*(i-len(h)+1))
            h[i]=val
    return h

def dual_check(h):
    N=len(h)
    b=None
    for cand in range(1,20):
        if 3**cand==N: b=cand; break
    assert b is not None, f"Size {N} is not 3^b"
    M=3**b
    LOG2_3=math.log2(3.0)
    mu=LOG2_3-2.0
    inv2={k: pow(pow(2,k,M), -1, M) for k in (1,2)}

    viol=0; min_gap=1e9; max_gap=-1e9
    for r in range(M):
        for k in (1,2):
            rp=((3*r-1)*inv2[k])%M
            w=LOG2_3-k
            # check: h[rp] <= h[r] + w - mu
            gap=(h[r] + (w-mu)) - h[rp]
            if gap < -1e-9:
                viol+=1
            if gap<min_gap: min_gap=gap
            if gap>max_gap: max_gap=gap
    return b, viol, min_gap, max_gap

h = load_h("potentials.csv")
b, viol, min_gap, max_gap = dual_check(h)
print(f"b={b}, violations={viol}, min_gap={min_gap:.3e}, max_gap={max_gap:.3e}")
