# Collatz via a Glued Grid and a Seam Automaton

**Status**: preprint-style research notes (community review invited)  
**Author**: Sgo Gobeaux  
**Repo**: code + figures + minimal scripts for identity checks and min–mean insights

## TL;DR (claims)
- **Two-step exact simulation (LH + SEAM)** of each odd Collatz step under the re-labelling \(D=\psi(y)=(y+1)/2\).
- A simple **Lyapunov** \(\Psi\) strictly increases per simulated step ⇒ **no non-trivial odd cycles**.
- On the **seam automaton** \(\Sigma_b\) (mod \(3^b\)), each seam has weight \(w=\log_2(3)-\kappa\), \(\kappa\in\{1,2\}\).
- There is a \(\kappa=2\) **self-loop** at \(r\equiv -1 \ (\mathrm{mod}\ 3^b)\) giving the **exact** min–mean value  
  \(\mu^\*=\log_2(3)-2\) (negative and **independent of** \(b\)).
- A short **block inequality** (with \(S\ge m\) seams for \(m\) odd steps and \(\delta(y)\le\log_2(4/3)\)) yields **no divergence**.
- Together: **bounded + no cycle ⇒ reach 1** (odd-only dynamics); even steps give full Collatz.

> **AI assistance note.** Parts of the exposition (structuring, wording, packaging) were assisted by an LLM.
> All mathematical claims and any mistakes remain the responsibility of the author.

## Core identities
For odd \(y=2D-1\):
\[
3y+1 = 2\,(3D-1),\quad
T(y)=\mathrm{oddize}(3y+1)=\mathrm{oddize}(3D-1),\quad
\psi\!\big(T(y)\big)=\frac{T(y)+1}{2}=\frac{\mathrm{oddize}(3D-1)+1}{2}=:\Sigma(D).
\]
Hence the **two-move** simulation \(L(D)\xrightarrow{\mathrm{LH}}M(D)\xrightarrow{\mathrm{SEAM}}L\big(\Sigma(D)\big)\) matches one odd step.

## Repo layout (suggested)



## Quick reproduce
```bash
# Python 3.10+
python code/identity_check.py --max-odd 200000
# prints: "OK identities up to ...", sample δ(y) bounds, etc.

