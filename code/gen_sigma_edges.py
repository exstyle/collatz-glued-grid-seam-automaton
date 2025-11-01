import sys, math, csv

def gen_sigma_edges(b:int, out_path:str):
    M = 3**b
    LOG2_3 = math.log2(3.0)
    inv2 = {k: pow(pow(2, k, M), -1, M) for k in (1,2)}  # (2^k)^(-1) mod M
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        wcsv = csv.writer(f)
        wcsv.writerow(["src","dst","weight"])
        for r in range(M):
            for k in (1,2):
                rp = ((3*r - 1) * inv2[k]) % M
                w = LOG2_3 - k
                wcsv.writerow([r, rp, f"{w:.12f}"])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python gen_sigma_edges.py <b> <out.csv>")
        sys.exit(1)
    b = int(sys.argv[1])
    out_csv = sys.argv[2]
    gen_sigma_edges(b, out_csv)
