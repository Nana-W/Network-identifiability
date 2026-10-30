# Network Identifiability Lab

Interactive companion to the paper

> N. Wang, J. M. Hendrickx, D. V. Dimarogonas,
> *From Local to Global Identifiability of Dynamic Networks:
> Computing and Reducing the Identifiability Degree.*

**Live tool:** https://nana-w.github.io/Network-identifiability/

## What it decides

For a dynamic network w = G(q)w + Br, y = Cw with known topology,
some modules known and the rest unknown, the tool decides, for your
choice of excited and measured nodes:

- **generic local identifiability** — the rank test
  rank Dh = m at a random point;
- the **identifiability degree d** — how many parameter vectors are
  compatible with a generic observation;
- the **trichotomy verdict** — not locally identifiable
  (positive-dimensional ambiguity), locally but not globally
  identifiable (d ≥ 2), or generically globally identifiable (d = 1);
- when d ≥ 2, the **repairs** — which added measurement or
  excitation separates the candidates and the degree each leaves,
  and which module priors select the true network.

## Four methods for the degree

The tool implements the methods of the paper's computation section,
selectable before each run; each result carries an expandable
*How this was computed* panel with the full numerical account.

| method | guarantee | cost | notes |
|---|---|---|---|
| Multistart Newton | none — lower bound on d | number of starts | fast first look |
| Monodromy + trace test | trace test checks completeness | ∝ d paths per loop | preferred when d is small relative to the bounds |
| Bézout homotopy | every isolated solution is a path endpoint, for a generic start (Morgan–Sommese 1989) | the multihomogeneous Bézout count of the lifted bilinear system | runs **two independent start systems** and reports whether they agree; square systems (m = p), up to 400 paths |
| Gröbner (exact) | all solutions by construction, exact arithmetic | exponential in m | generates a ready-to-run sympy script for the network as drawn |

On the paper's 8-node example: d = 2, the multihomogeneous count is
6 paths, the coarser bound of the paper is 16, and the generated
exact script returns the two fiber points as exact rationals.

## Quick start

1. Open the live tool; the paper's example loads with d = 2.
2. Press **Run identifiability analysis**. The rank test passes, the
   two compatible parameter vectors appear with residuals, and the
   repairs list shows that measuring node c or d, or exciting node
   b, restores d = 1.
3. Switch the **Degree method** and rerun to compare the methods;
   expand the *How* panels to see each method's evidence.
4. Click a repair's **apply & rerun** to watch the degree drop.
5. Build your own network on the canvas (click cycles
   excited → measured → both → none; drag to arrange), or paste an
   edge list in **Define by text**:

   ```
   e1 -> a : 1      # known module with gain 1
   a -> c           # unknown module
   excite: e1, e2
   measure: m1, m2
   ```

   Cyclic networks are supported; self-loops are not (G has zero
   diagonal).

## How it computes

Everything runs in the browser, no server. The identification map
h(θ) = C(I−G(θ))⁻¹B is evaluated by complex LU solves; the Jacobian
is exact via the resolvent identity ∂h(c,b)/∂θ(i→j) = M[c,j]·M[i,b];
the rank comes from a Jacobi eigenvalue sweep maximized over random
points; the Bézout homotopy tracks a linear-product start system
respecting the 2-homogeneous structure of the lifted bilinear form
(parameters grouped by head node, latent responses by excitation
column) with the gamma trick; monodromy moves the observation around
random complex loops and applies the numerical trace test;
separating nodes compare latent responses and transfers across the
fiber, as in the paper's reduction section.

## What to trust

- The rank test holds with probability one in the random draw.
- Multistart Newton is a lower bound with no stopping rule.
- The Bézout homotopy carries a probability-one completeness
  argument for generic starts; agreement of its two independent
  runs is the reported consistency check. Monodromy's trace test
  checks completeness of the points found. Both remain floating
  point, so a proof-grade count needs the exact route.
- The generated sympy script computes d in exact rational
  arithmetic and is the certificate-grade option for small m.

## Reproducing the paper's numbers

| script | reproduces |
|---|---|
| `groebner_example.py` | the lex Gröbner basis, the eliminant (3x+1)(7x−3)/21 and both fiber points of the example |
| `lifted_system.py` | the 16-variable lifted bilinear system, its two solutions, and the multihomogeneous Bézout bound 16 |
| `mixed_volume.py` | the mixed volume 2 of the eliminated form, with the normalization check |

Run with Python 3 after `pip install sympy numpy scipy`. The tool's
Gröbner option generates the analogous exact script for any network
you draw.

## Citing

If you use this tool, please cite the paper above. BibTeX will be
added once the paper is published.

## License

MIT.
