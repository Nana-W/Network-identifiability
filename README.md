**Live tool:** https://nana-w.github.io/Network-identifiability/

## What it decides

For a dynamic network w = G(q)w + Br, y = Cw with known topology,
some modules known and the rest unknown, the tool decides, for your
choice of excited and measured nodes:

- **generic local identifiability** — the rank test
  rank Dh = m at a random point;
- the **identifiability degree d** — how many parameter vectors are
  compatible with a generic observation;
- the **trichotomy verdict** — not locally identifiable (d undefined,
  positive-dimensional ambiguity), locally but not globally
  identifiable (d ≥ 2), or generically globally identifiable (d = 1);
- when d ≥ 2, the **repairs** — which added measurement or
  excitation separates the candidates and the degree each leaves,
  and which module priors select the true network.

## Quick start

1. Open the live tool; the paper's 8-node example loads with d = 2.
2. Press **Run identifiability analysis**: the rank test passes, the
   two compatible parameter vectors appear with their residuals, and
   the repairs list shows that measuring node c or d, or exciting
   node b, restores d = 1. Expand *How the Jacobian was computed* and
   *How the fiber was computed* for the full numerical account.
3. Click a repair's **apply & rerun** to see the degree drop.
4. Build your own network: click the canvas (nodes cycle through
   excited → measured → both → none; drag to arrange), or paste an
   edge list in **Define by text**:


   Cyclic networks are supported; self-loops are not (G has zero
   diagonal).

## How it computes

Everything runs in the browser, no server. The identification map
h(θ) = vec C(I−G(θ))⁻¹B is evaluated by complex LU solves; the
Jacobian is exact via the resolvent identity
∂h(c,b)/∂θ(i→j) = M[c,j]·M[i,b] with M = (I−G)⁻¹; the rank comes
from a Jacobi eigenvalue sweep of JᵀJ, maximized over random
points; the fiber of h(θ) = h(θ*) is found by damped Gauss–Newton
from a few hundred complex Gaussian starts, clustered; separating
nodes compare the latent responses X = (I−G)⁻¹B and the transfers
C(I−G)⁻¹e_v across the fiber points, as in the paper's reduction
section.

## What to trust

- The rank test is exact arithmetic up to floating point and holds
  with probability one in the random draw.
- The reported d is a **numerical lower bound**: multistart Newton
  can miss fiber points. It is reliable for small networks, roughly
  m ≤ 8 unknown modules. The paper's certified route — Gröbner
  elimination, or monodromy with the trace test — is what a proof
  needs; this page is for exploration.
- Agreement of two runs is a cheap robustness check (starts are
  randomized; the probe θ* is fixed).

## Reproducing the paper's numbers

| script | reproduces |
|---|---|
| `groebner_example.py` | the lex Gröbner basis, the eliminant (3x+1)(7x−3)/21 and both fiber points of the example |
| `lifted_system.py` | the 16-variable lifted bilinear system, its two solutions, and the multihomogeneous Bézout bound 16 |
| `mixed_volume.py` | the mixed volume 2 of the eliminated form, with the normalization check |

Run with Python 3 after `pip install sympy numpy scipy`.

## Citing

If you use this tool, please cite the paper above. BibTeX will be
added once the paper is published.

## License

MIT.
