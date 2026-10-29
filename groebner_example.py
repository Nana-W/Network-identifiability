"""Eliminated form of the paper's example: Groebner basis, eliminant, fiber.

Reproduces eq. (gb): the lex basis in shape position and the eliminant
(1/21)(3x+1)(7x-3), confirming identifiability degree d = 2.
"""
import sympy as sp

x, y, w, z = sp.symbols('x y w z')
h = [x + y + x*y,          # T_{m1 e1}
     2*x + x*y,            # T_{m2 e2}
     x + x*y + w,          # T_{m1 e2}
     2*x + y + x*y + z]    # T_{m2 e1}

theta_star = {x: sp.Rational(3, 7), y: sp.Rational(-5, 3),
              w: sp.Rational(2, 9), z: sp.Rational(11, 4)}
y_star = [e.subs(theta_star) for e in h]
print("observation y* =", y_star)

G = sp.groebner([e - v for e, v in zip(h, y_star)], w, z, y, x, order='lex')
print("\nlex Groebner basis (w > z > y > x):")
for g in G.exprs:
    print("  ", sp.factor(g))

elim = [g for g in G.exprs if g.free_symbols == {x}][0]
print("\neliminant:", sp.factor(elim), "  degree d =", sp.degree(elim, x))

sols = sp.solve([e - v for e, v in zip(h, y_star)], [x, y, w, z], dict=True)
print("\nfiber points:")
for s in sols:
    print("  ", [s[v] for v in (x, y, w, z)])
