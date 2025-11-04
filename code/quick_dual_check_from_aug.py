# quick_dual_check_from_aug.py
import json, csv, math, sys
edges_csv, cert_json = sys.argv[1], sys.argv[2]
with open(cert_json,'r',encoding='utf-8') as f:
    cert = json.load(f)
mu, Q, h = cert["mu_b"], cert["scale"], cert["h"]
ok=True
with open(edges_csv,newline='',encoding='utf-8') as f:
    rd = csv.DictReader(f)
    for row in rd:
        r  = int(row.get('src') or row.get('r') or row.get('from') or row.get('u'))
        rp = int(row.get('dst') or row.get('rp') or row.get('to') or row.get('v'))
        w_aug = float(row.get('w_aug') or row.get('waug') or row.get('w_augmented'))
        lhs = (h[str(rp)] - h[str(r)]) / (1<<Q)
        rhs = w_aug - mu
        if lhs > rhs + 1e-12:
            ok=False; print("FAIL:", r,"->",rp, "lhs",lhs, ">", "rhs",rhs); break
print("OK" if ok else "NOT OK")
