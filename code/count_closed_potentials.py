import csv, math
h=[]; 
with open("potentials.csv", newline="", encoding="utf-8") as f:
    rd=csv.reader(f); next(rd,None)
    for i,val in rd: 
        h.append(float(val))
M=len(h); b=round(math.log(M,3))
LOG2_3=math.log2(3.0); mu=LOG2_3-2.0
inv2={k: pow(pow(2,k,3**b), -1, 3**b) for k in (1,2)}
tight0=tight1=0
for r in range(M):
    for k in (1,2):
        rp=((3*r-1)*inv2[k])%(M)
        w=LOG2_3-k
        gap=(h[r]+(w-mu))-h[rp]
        if abs(gap) < 1e-9: tight0 += 1      # arêtes saturées (surtout k=2)
        if abs(gap-1.0) < 1e-9: tight1 += 1  # arêtes avec marge 1 (k=1)
print("tight (gap=0):", tight0, "  tight (gap=1):", tight1)
