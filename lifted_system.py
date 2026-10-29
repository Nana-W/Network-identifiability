"""Lifted bilinear form (I-G)X = B, CX = Y* of the paper's example.

Builds the 16-equation, 16-variable system after substituting the
measured rows, solves it exactly, and confirms the two fiber points
and the multihomogeneous Bezout bound 16.
"""
import sympy as sp
from collections import Counter

names = ['e1', 'e2', 'a', 'b', 'c', 'd', 'm1', 'm2']
I = {s: i for i, s in enumerate(names)}
n = len(names)
x, y, w, z = sp.symbols('x y w z')

G = sp.zeros(n, n)
def edge(t, h, g): G[I[h], I[t]] = g
edge('e1', 'a', 1); edge('e1', 'b', 1); edge('e2', 'a', 1)
edge('a', 'c', x);  edge('c', 'b', 1);  edge('b', 'd', y)
edge('c', 'm1', 1); edge('d', 'm1', 1); edge('c', 'm2', 2); edge('d', 'm2', 1)
edge('e2', 'm1', w); edge('e1', 'm2', z)

B = sp.zeros(n, 2); B[I['e1'], 0] = 1; B[I['e2'], 1] = 1
C = sp.zeros(2, n); C[0, I['m1']] = 1; C[1, I['m2']] = 1
measured = ('m1', 'm2')

theta_star = {x: sp.Rational(3, 7), y: sp.Rational(-5, 3),
              w: sp.Rational(2, 9), z: sp.Rational(11, 4)}
Xstar = (sp.eye(n) - G.subs(theta_star)).inv() * B
Ystar = C * Xstar

Xs = {}
for s in names:
    for r in range(2):
        if s in measured:
            Xs[(s, r)] = Ystar[measured.index(s), r]   # substituted: data
        else:
            Xs[(s, r)] = sp.Symbol(f'X_{s}_{r}')
variables = [x, y, w, z] + [Xs[(s, r)] for s in names
                            if s not in measured for r in range(2)]
eqs = []
for s in names:
    for r in range(2):
        expr = Xs[(s, r)] - B[I[s], r]
        for t in names:
            if G[I[s], I[t]] != 0:
                expr -= G[I[s], I[t]] * Xs[(t, r)]
        eqs.append(sp.expand(expr))

print(f"N_v = {len(variables)} variables, N = {len(eqs)} equations")
sols = sp.solve(eqs, variables, dict=True)
print("number of solutions d =", len(sols))
for s in sols:
    print("  theta =", [s.get(v) for v in (x, y, w, z)])

mi = Counter(hd for (tl, hd) in
             [('a','c'), ('b','d'), ('e2','m1'), ('e1','m2')])
bez = sp.prod(sp.binomial(2, k) for k in mi.values())
print("multihomogeneous Bezout bound =", bez)
