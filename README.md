# Collatz via a Glued Grid and a Seam Automaton

[![DOI](https://zenodo.org/badge/1087917470.svg)](https://doi.org/10.5281/zenodo.17504109)

**Status**: preprint-style research notes (community review invited)  
**Author**: Sgo Gobeaux  
**Repo**: code + figures + minimal scripts for identity checks, min–mean, and potential-building

---

## TL;DR — what’s proved vs. what’s conditional

- **Exact 2-move simulation (LH + SEAM)** for each odd Collatz step under the relabelling \(D=\psi(y)=(y+1)/2\):
  
  \(3y+1=2(3D-1),\quad
  T(y)=\mathrm{oddize}(3y+1)=\mathrm{oddize}(3D-1),\quad
  \psi(T(y))=\dfrac{\mathrm{oddize}(3D-1)+1}{2}=:\Sigma(D).\)

  So one odd step equals the two-move path \(L(D)\xrightarrow{\text{LH}}M(D)\xrightarrow{\text{SEAM}}L(\Sigma(D))\).

- **Seam automaton** \(\Sigma_b\) (mod \(3^b\)) with types \(\kappa\in\{1,2\}\) and weights \(w=\log_2(3)-\kappa\).  
  A \(\kappa=2\) loop at \(r\equiv -1\pmod{3^b}\) yields the **exact** min–mean
  \(\mu^\*(b)=\log_2(3)-2\in(-1,0)\) (independent of \(b\)).

- **Augmented weights** \(\tilde w=(\log_2 3-\kappa)+\delta(r)\) with
  \(\delta(r)=\log_2\!\big(1/(3y_{\min}(r))+1\big)\)
  give a **strictly negative** augmented min–mean:
  \(\mu^\*_{\text{aug}}(b)\le \log_2(5/6)<0\).
  → **Robust block inequality**:
  \(\sum_{i=0}^{m-1}\Delta\log_2 y_i \le \mu^\*_{\text{aug}}\,m + \mathrm{osc}(h)\).

- **Flight time** (odd-only): with \(\varepsilon:=-\mu^\*_{\text{aug}}>0\),
  \(m(y_0)\le \lceil \log_2(y_0/Y^\*)/\varepsilon\rceil + C(Y^\*)\).

- **Conditional acyclicity (explicit)**:  
  If there exists a node potential \(\Psi\) with
  \(\Psi(M(D))-\Psi(L(D))=D+2\) and \(\Psi(L(\Sigma(D)))-\Psi(M(D))=-2\),
  then the 2-move dynamics is acyclic ⇒ no non-trivial odd cycles.  
  Combined with the block inequality, **odd trajectories reach 1** (and with even steps, full Collatz).

> **AI assistance note.** Parts of the exposition (structuring, wording, packaging) were assisted by an LLM.  
> All mathematics, computations, and any mistakes remain the author’s responsibility.

---

## Core identities (one place)

For odd \(y=2D-1\):
  
\(3y+1=2(3D-1),\quad
T(y)=\mathrm{oddize}(3y+1)=\mathrm{oddize}(3D-1),\quad
\psi(T(y))=\dfrac{\mathrm{oddize}(3D-1)+1}{2}=:\Sigma(D).\)

Hence: \(L(D)\to M(D)\to L(\Sigma(D))\) simulates one odd step exactly.

---

## Repository layout (suggested)

- `code/` — Python scripts (builders, min–mean, potential F, sweeps)
- `data/` — Generated CSVs (edges, potentials, sweeps) — gitignored except samples
- `figures/` — Diagrams (glued grid, seam automaton)
- `doc/` — TeX/notes (optional)
- `CITATION.cff` — GitHub “Cite this repository”
- `README.md` — this file

---

## Quickstart (Python 3.10+)

**Windows tip (PowerShell):** create the `data` folder first.

    mkdir data

### 1) Build augmented Σ_b and compute min–mean (Howard)

Build Σ_b (augmented, b=10):

    python code/build_sigma_augmented.py --b 10 --out data/sigma_b10_aug.csv

Min–mean (augmented):

    python code/howard_minmean_augmented.py --edges data/sigma_b10_aug.csv --progress --report --out-h data/h_b10_howard.csv

You should see \(\mu\approx \log_2(3)-2\) (negative). The script also prints \(\varepsilon=-\mu>0\).

### 2) Sweep several b

    python code/sweep_b.py --b-min 9 --b-max 12 --out data/sweep_b9_12.csv --persist --progress

### 3) Non-augmented edges (positional args)

Generator usage is **positional**:

    python code/gen_sigma_edges.py 10 data/sigma_b10_edges.csv

Then run non-augmented min–mean:

    python code/howard_minmean.py data/sigma_b10_edges.csv --src src --dst dst --weight w --progress

### 4) Evidence for 2-step acyclicity (cohomological F)

We solve
\(F(\Sigma(D))-F(D)=D \ (D>1),\ F(1)=0\) on large prefixes with **no violations** and **no detected 2-step cycles**.

Local prefix (“clamp”):

    python code/v6_build_and_check_F.py --N 200000 --mode clamp --out-csv data/F_built_clamp.csv

Expanded coverage:

    python code/v6_build_and_check_F.py --N 1000000 --mode expand --cap 50000000000 --out-csv data/F_expand_50B.csv

Example output (typical):

- \(|F|\) in the millions; 0 non-trivial cycles detected  
- Checks passed: cohomology \(F(S)-F(D)=D\), LH \(=D+2\), SEAM \(=-2\)

> This is **evidence**, not a proof, for the Part I acyclicity condition.

---

## Flight-time upper bound (odd-only)

With \(\mu^\*_{\text{aug}}<0\) and \(\varepsilon:=-\mu^\*_{\text{aug}}>0\),

\(m(y_0)\ \le\ \left\lceil \dfrac{\log_2(y_0/Y^\*)}{\varepsilon} \right\rceil + C(Y^\*).\)

Scripts `howard_minmean_augmented.py` and `sweep_b.py` report \(\varepsilon\).

---

## What is still conditional?

Part I’s “no-cycle” uses a **node potential** \(\Psi\) with increments
\(\Delta\Psi(\mathrm{LH})=D+2\) and \(\Delta\Psi(\mathrm{SEAM})=-2\).
Once such a \(\Psi\) is **exhibited** (or replaced by a well-founded order),
Part II (strictly negative augmented min–mean + robust block inequality) completes:
**bounded + no cycle ⇒ reach 1**.

---

## Citing this work

- **Zenodo**: see DOI badge above.  
- **GitHub**: `CITATION.cff` is included so GitHub renders “Cite this repository”.

