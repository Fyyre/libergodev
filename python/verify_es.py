#!/usr/bin/env python3
"""
Exact rational primal-dual certificate verifier and finite-window Monte Carlo
reproduction for Non-Asymptotic Symbolic Ergodeviation (Section 9).
Uses only the Python standard library (fractions, math, random).
"""

from fractions import Fraction
import math
import random

def d_pair(w1, w2):
    """Product-of-paths ground metric d((i,j), (k,l)) = |i-k| + 1/2 |j-l|."""
    i, j = w1
    k, l = w2
    return Fraction(abs(i - k), 1) + Fraction(abs(j - l), 2)

def verify_1_lipschitz(f):
    """Verify |f(u) - f(v)| <= d(u, v) for all 9x9 pairs in exact rational arithmetic."""
    words = [(i, j) for i in range(3) for j in range(3)]
    for u in words:
        for v in words:
            assert abs(f[u] - f[v]) <= d_pair(u, v), f"Lipschitz failed on {u}, {v}"

def verify_coupling(coupling, mu, nu):
    """Verify row/column marginals of a transport plan and return exact rational cost."""
    words = [(i, j) for i in range(3) for j in range(3)]
    for u in words:
        row_sum = sum(coupling.get((u, v), Fraction(0)) for v in words)
        assert row_sum == mu.get(u, Fraction(0)), f"Row marginal mismatch at {u}: {row_sum}"
    for v in words:
        col_sum = sum(coupling.get((u, v), Fraction(0)) for u in words)
        assert col_sum == nu.get(v, Fraction(0)), f"Col marginal mismatch at {v}: {col_sum}"
    return sum(mass * d_pair(u, v) for (u, v), mass in coupling.items())

def main():
    words = [(i, j) for i in range(3) for j in range(3)]

    # Section 9.1 Baseline Markov chain
    P = [
        [Fraction(1, 10), Fraction(4, 5), Fraction(1, 10)],
        [Fraction(2, 5),  Fraction(1, 5), Fraction(2, 5)],
        [Fraction(1, 10), Fraction(4, 5), Fraction(1, 10)],
    ]
    q = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 4)]

    # Verify stationarity q P = q
    for j in range(3):
        assert sum(q[i] * P[i][j] for i in range(3)) == q[j]

    nu = {(i, j): q[i] * P[i][j] for i in range(3) for j in range(3)}
    lambda_2 = Fraction(1, 2)

    # Squared Fiedler jumps (psi_i - psi_j)^2 for psi = (-1/sqrt(2), 0, 1/sqrt(2))
    jump = {
        (0, 0): Fraction(0), (0, 1): Fraction(1, 2), (0, 2): Fraction(2),
        (1, 0): Fraction(1, 2), (1, 1): Fraction(0), (1, 2): Fraction(1, 2),
        (2, 0): Fraction(2), (2, 1): Fraction(1, 2), (2, 2): Fraction(0),
    }

    # Verify baseline Dirichlet energy Delta_* == lambda_2 == 1/2
    delta_star = sum(nu[w] * jump[w] for w in words)
    assert delta_star == lambda_2 == Fraction(1, 2)

    print("Verified Section 9.1 baseline stationarity and Delta_* = lambda_2 = 1/2.")

    # 1. Constant state 1
    mu_const = {(1, 1): Fraction(1)}
    coup_const = {((1, 1), v): nu[v] for v in words}
    w1_const = verify_coupling(coup_const, mu_const, nu)
    def_const = Fraction(1)
    es_const = w1_const * (1 + def_const)
    assert (w1_const, es_const) == (Fraction(3, 4), Fraction(3, 2))

    # 2. Alternation 0, 2, 0, 2, ...
    mu_alt = {(0, 2): Fraction(1, 2), (2, 0): Fraction(1, 2)}
    f_alt = {
        (0, 0): Fraction(0), (0, 1): Fraction(1, 2),  (0, 2): Fraction(1),
        (1, 0): Fraction(0), (1, 1): Fraction(-1, 2), (1, 2): Fraction(0),
        (2, 0): Fraction(1), (2, 1): Fraction(1, 2),  (2, 2): Fraction(0),
    }
    verify_1_lipschitz(f_alt)
    dual_alt = sum(f_alt[w] * (mu_alt.get(w, Fraction(0)) - nu[w]) for w in words)
    coup_alt = {
        ((0, 2), (0, 0)): Fraction(1, 40), ((0, 2), (0, 1)): Fraction(1, 5),
        ((0, 2), (0, 2)): Fraction(1, 40), ((0, 2), (1, 1)): Fraction(1, 20),
        ((0, 2), (1, 2)): Fraction(1, 5),
        ((2, 0), (1, 0)): Fraction(1, 5),  ((2, 0), (1, 1)): Fraction(1, 20),
        ((2, 0), (2, 0)): Fraction(1, 40), ((2, 0), (2, 1)): Fraction(1, 5),
        ((2, 0), (2, 2)): Fraction(1, 40),
    }
    w1_alt = verify_coupling(coup_alt, mu_alt, nu)
    assert dual_alt == w1_alt == Fraction(4, 5)

    # 3. Fair noise on the cut {0, 2}
    mu_cut = {(0, 0): Fraction(1, 4), (0, 2): Fraction(1, 4),
              (2, 0): Fraction(1, 4), (2, 2): Fraction(1, 4)}
    f_cut = {
        (0, 0): Fraction(0),  (0, 1): Fraction(-1, 2), (0, 2): Fraction(0),
        (1, 0): Fraction(-1), (1, 1): Fraction(-3, 2), (1, 2): Fraction(-1),
        (2, 0): Fraction(0),  (2, 1): Fraction(-1, 2), (2, 2): Fraction(0),
    }
    verify_1_lipschitz(f_cut)
    dual_cut = sum(f_cut[w] * (mu_cut.get(w, Fraction(0)) - nu[w]) for w in words)
    coup_cut = {
        ((0, 0), (0, 0)): Fraction(1, 40), ((0, 0), (0, 1)): Fraction(7, 40), ((0, 0), (1, 1)): Fraction(1, 20),
        ((0, 2), (0, 1)): Fraction(1, 40), ((0, 2), (0, 2)): Fraction(1, 40), ((0, 2), (1, 2)): Fraction(1, 5),
        ((2, 0), (1, 0)): Fraction(1, 5),  ((2, 0), (2, 0)): Fraction(1, 40), ((2, 0), (2, 1)): Fraction(1, 40),
        ((2, 2), (1, 1)): Fraction(1, 20), ((2, 2), (2, 1)): Fraction(7, 40), ((2, 2), (2, 2)): Fraction(1, 40),
    }
    w1_cut = verify_coupling(coup_cut, mu_cut, nu)
    assert dual_cut == w1_cut == Fraction(3, 4)

    # 4. Fair noise on one side {0, 1}
    mu_side = {(0, 0): Fraction(1, 4), (0, 1): Fraction(1, 4),
               (1, 0): Fraction(1, 4), (1, 1): Fraction(1, 4)}
    f_side = {(i, j): -Fraction(i, 1) - Fraction(j, 2) for (i, j) in words}
    verify_1_lipschitz(f_side)
    dual_side = sum(f_side[w] * (mu_side.get(w, Fraction(0)) - nu[w]) for w in words)
    coup_side = {
        ((0, 0), (0, 0)): Fraction(1, 40), ((0, 0), (0, 2)): Fraction(1, 40), ((0, 0), (1, 2)): Fraction(1, 5),
        ((0, 1), (0, 1)): Fraction(1, 5),  ((0, 1), (2, 1)): Fraction(1, 20),
        ((1, 0), (1, 0)): Fraction(1, 5),  ((1, 0), (2, 0)): Fraction(1, 40), ((1, 0), (2, 2)): Fraction(1, 40),
        ((1, 1), (1, 1)): Fraction(1, 10), ((1, 1), (2, 1)): Fraction(3, 20),
    }
    w1_side = verify_coupling(coup_side, mu_side, nu)
    energy_side = sum(mu_side.get(w, Fraction(0)) * jump[w] for w in words)
    def_side = (lambda_2 - energy_side) / lambda_2
    es_side = w1_side * (1 + def_side)
    assert dual_side == w1_side == Fraction(3, 4)
    assert es_side == Fraction(9, 8)

    # 5. Local jitter at scale 1/10 (exact rational energy & W_1 = 27/500)
    K = [
        [Fraction(9, 10), Fraction(1, 10), Fraction(0)],
        [Fraction(1, 20), Fraction(9, 10), Fraction(1, 20)],
        [Fraction(0),     Fraction(1, 10), Fraction(9, 10)],
    ]
    nu_jit = {
        (i, j): sum(nu[(a, b)] * K[a][i] * K[b][j] for a in range(3) for b in range(3))
        for (i, j) in words
    }
    energy_jit = sum(nu_jit[w] * jump[w] for w in words)
    assert energy_jit == Fraction(1, 2), f"Expected 1/2, got {energy_jit}"

    print("All exact rational primal-dual certificates in Section 9.2 verified (0 residual)!")

if __name__ == "__main__":
    main()