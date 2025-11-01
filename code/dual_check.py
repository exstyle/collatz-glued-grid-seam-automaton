import sys, csv, math
def dual_check(pot_file):
    h=[]
    with open(pot_file, newline="", encoding="utf-8") as f:
        rd=csv.reader(f); next(rd,None)
        for i,val in rd:
            h.append(float(val))
    N=len(h)
    b=None
    for cand in range(1,20):
        if 3**cand==N: b=cand; break
    assert b is not None, f"Size {N} != 3^b"
    M=3**b
    LOG2_3=math.log2(3.0); mu=LOG2_3-2.0
    inv2={k: pow(pow(2,k,M), -1, M) for k in (1,2)}
    viol=0; min_gap=1e9; max_gap=-1e9
    for r in range(M):
        for k in (1,2):
            rp=((3*r-1)*inv2[k])%M
            w=LOG2_3-k
            gap=(h[r] + (w - mu)) - h[rp]   # want gap >= 0
            if gap < -1e-9: viol+=1
            min_gap=min(min_gap,gap); max_gap=max(max_gap,gap)
    print(f"b={b}, violations={viol}, min_gap={min_gap:.3e}, max_gap={max_gap:.3e}")
if __name__=="__main__":
    for f in sys.argv[1:]:
        dual_check(f)
