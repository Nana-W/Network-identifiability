# Network Identifiability Lab

Interactive companion to the paper

> N. Wang, J. M. Hendrickx, D. V. Dimarogonas,
> *From Local to Global Identifiability of Dynamic Networks:
> Computing and Reducing the Identifiability Degree.*

**Live tool:** https://USERNAME.github.io/REPONAME/

Build a dynamic network in the browser (or paste it as an edge
list), mark which nodes are excited and measured, and the tool

1. runs the **rank test** for generic local identifiability at
   random points,
2. counts the parameter vectors compatible with a generic
   observation — the **identifiability degree d** — by multistart
   Newton on the fiber system `h(θ) = h(θ*)`,
3. reports the **trichotomy verdict** (not locally identifiable /
   locally but not globally, d ≥ 2 / globally identifiable, d = 1),
4. when d ≥ 2, lists the **separating repairs**: which added
   measurement or excitation reduces the degree to what, and which
   module priors select the true network.

The paper's example ships as the default preset (d = 2).

### Text format

```
e1 -> a : 1      # known module with gain 1
a -> c           # unknown module
excite: e1, e2
measure: m1, m2
```

### Caveats

The rank test holds with probability one at a random point. The
degree is a numerical count: multistart Newton can miss fiber
points, so the reported d is a lower bound, reliable for small
networks (roughly m ≤ 8 unknown modules). For larger networks use
the exact and certified methods described in the paper; the
`scripts/` folder reproduces the paper's computations.

### Scripts

| file | reproduces |
|---|---|
| `scripts/groebner_example.py` | the Gröbner basis (eq. (gb)), eliminant and fiber of the example (sympy) |
| `scripts/lifted_system.py` | the 16-variable lifted bilinear system and its two solutions |
| `scripts/mixed_volume.py` | the mixed volume 2 of the eliminated form and the Bézout bound 16 |

Each script runs with Python 3 and `pip install sympy numpy scipy`.

### License

MIT.
