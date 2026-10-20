#!/usr/bin/env python3
"""
Exact rational primal-dual certificate verifier and finite-window Monte Carlo
reproduction for Non-Asymptotic Symbolic Ergodeviation (Section 9).
Uses only the Python standard library (fractions, math, random).
"""

from fractions import Fraction
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
        assert row_sum == mu.get(
            u, Fraction(0)
        ), f"Row marginal mismatch at {u}: {row_sum}"
    for v in words:
        col_sum = sum(coupling.get((u, v), Fraction(0)) for u in words)
        assert col_sum == nu.get(
            v, Fraction(0)
        ), f"Col marginal mismatch at {v}: {col_sum}"
    return sum(mass * d_pair(u, v) for (u, v), mass in coupling.items())


def sinkhorn_w1(mu_counts, total_pairs, nu_float, dist_mat, reg=0.002, iters=150):
    """Fast standard-library log-domain Sinkhorn solver for finite-window Monte Carlo."""
    import math

    mu = [c / total_pairs for c in mu_counts]
    if sum(abs(mu[i] - nu_float[i]) for i in range(9)) < 1e-14:
        return 0.0
    f = [0.0] * 9
    g = [0.0] * 9
    log_mu = [math.log(m) if m > 0.0 else -1e30 for m in mu]
    log_nu = [math.log(n) for n in nu_float]

    for _ in range(iters):
        for i in range(9):
            if mu[i] <= 0.0:
                f[i] = -1e30
                continue
            vals = [(g[j] - dist_mat[i][j]) / reg for j in range(9)]
            mval = max(vals)
            f[i] = reg * (
                log_mu[i] - (mval + math.log(sum(math.exp(v - mval) for v in vals)))
            )
        for j in range(9):
            vals = [(f[i] - dist_mat[i][j]) / reg for i in range(9) if mu[i] > 0.0]
            mval = max(vals)
            g[j] = reg * (
                log_nu[j] - (mval + math.log(sum(math.exp(v - mval) for v in vals)))
            )

    w1 = 0.0
    for i in range(9):
        if mu[i] <= 0.0:
            continue
        for j in range(9):
            p_ij = math.exp((f[i] + g[j] - dist_mat[i][j]) / reg)
            w1 += p_ij * dist_mat[i][j]
    return w1


def main():
    words = [(i, j) for i in range(3) for j in range(3)]

    # Section 9.1 Baseline Markov chain
    P = [
        [Fraction(1, 10), Fraction(4, 5), Fraction(1, 10)],
        [Fraction(2, 5), Fraction(1, 5), Fraction(2, 5)],
        [Fraction(1, 10), Fraction(4, 5), Fraction(1, 10)],
    ]
    q = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 4)]

    for j in range(3):
        assert sum(q[i] * P[i][j] for i in range(3)) == q[j]

    nu = {(i, j): q[i] * P[i][j] for i in range(3) for j in range(3)}
    lambda_2 = Fraction(1, 2)

    jump = {
        (0, 0): Fraction(0),
        (0, 1): Fraction(1, 2),
        (0, 2): Fraction(2),
        (1, 0): Fraction(1, 2),
        (1, 1): Fraction(0),
        (1, 2): Fraction(1, 2),
        (2, 0): Fraction(2),
        (2, 1): Fraction(1, 2),
        (2, 2): Fraction(0),
    }

    delta_star = sum(nu[w] * jump[w] for w in words)
    assert delta_star == lambda_2 == Fraction(1, 2)
    print("1. Verified Section 9.1 baseline stationarity and Delta_* = lambda_2 = 1/2.")

    # 1. Constant state 1
    mu_const = {(1, 1): Fraction(1)}
    coup_const = {((1, 1), v): nu[v] for v in words}
    w1_const = verify_coupling(coup_const, mu_const, nu)
    es_const = w1_const * (1 + Fraction(1))
    assert (w1_const, es_const) == (Fraction(3, 4), Fraction(3, 2))

    # 2. Alternation 0, 2, 0, 2, ...
    mu_alt = {(0, 2): Fraction(1, 2), (2, 0): Fraction(1, 2)}
    f_alt = {
        (0, 0): Fraction(0),
        (0, 1): Fraction(1, 2),
        (0, 2): Fraction(1),
        (1, 0): Fraction(0),
        (1, 1): Fraction(-1, 2),
        (1, 2): Fraction(0),
        (2, 0): Fraction(1),
        (2, 1): Fraction(1, 2),
        (2, 2): Fraction(0),
    }
    verify_1_lipschitz(f_alt)
    dual_alt = sum(f_alt[w] * (mu_alt.get(w, Fraction(0)) - nu[w]) for w in words)
    coup_alt = {
        ((0, 2), (0, 0)): Fraction(1, 40),
        ((0, 2), (0, 1)): Fraction(1, 5),
        ((0, 2), (0, 2)): Fraction(1, 40),
        ((0, 2), (1, 1)): Fraction(1, 20),
        ((0, 2), (1, 2)): Fraction(1, 5),
        ((2, 0), (1, 0)): Fraction(1, 5),
        ((2, 0), (1, 1)): Fraction(1, 20),
        ((2, 0), (2, 0)): Fraction(1, 40),
        ((2, 0), (2, 1)): Fraction(1, 5),
        ((2, 0), (2, 2)): Fraction(1, 40),
    }
    w1_alt = verify_coupling(coup_alt, mu_alt, nu)
    assert dual_alt == w1_alt == Fraction(4, 5)

    # 3. Fair noise on the cut {0, 2}
    mu_cut = {
        (0, 0): Fraction(1, 4),
        (0, 2): Fraction(1, 4),
        (2, 0): Fraction(1, 4),
        (2, 2): Fraction(1, 4),
    }
    f_cut = {
        (0, 0): Fraction(0),
        (0, 1): Fraction(-1, 2),
        (0, 2): Fraction(0),
        (1, 0): Fraction(-1),
        (1, 1): Fraction(-3, 2),
        (1, 2): Fraction(-1),
        (2, 0): Fraction(0),
        (2, 1): Fraction(-1, 2),
        (2, 2): Fraction(0),
    }
    verify_1_lipschitz(f_cut)
    dual_cut = sum(f_cut[w] * (mu_cut.get(w, Fraction(0)) - nu[w]) for w in words)
    coup_cut = {
        ((0, 0), (0, 0)): Fraction(1, 40),
        ((0, 0), (0, 1)): Fraction(7, 40),
        ((0, 0), (1, 1)): Fraction(1, 20),
        ((0, 2), (0, 1)): Fraction(1, 40),
        ((0, 2), (0, 2)): Fraction(1, 40),
        ((0, 2), (1, 2)): Fraction(1, 5),
        ((2, 0), (1, 0)): Fraction(1, 5),
        ((2, 0), (2, 0)): Fraction(1, 40),
        ((2, 0), (2, 1)): Fraction(1, 40),
        ((2, 2), (1, 1)): Fraction(1, 20),
        ((2, 2), (2, 1)): Fraction(7, 40),
        ((2, 2), (2, 2)): Fraction(1, 40),
    }
    w1_cut = verify_coupling(coup_cut, mu_cut, nu)
    assert dual_cut == w1_cut == Fraction(3, 4)

    # 4. Fair noise on one side {0, 1}
    mu_side = {
        (0, 0): Fraction(1, 4),
        (0, 1): Fraction(1, 4),
        (1, 0): Fraction(1, 4),
        (1, 1): Fraction(1, 4),
    }
    f_side = {(i, j): -Fraction(i, 1) - Fraction(j, 2) for (i, j) in words}
    verify_1_lipschitz(f_side)
    dual_side = sum(f_side[w] * (mu_side.get(w, Fraction(0)) - nu[w]) for w in words)
    coup_side = {
        ((0, 0), (0, 0)): Fraction(1, 40),
        ((0, 0), (0, 2)): Fraction(1, 40),
        ((0, 0), (1, 2)): Fraction(1, 5),
        ((0, 1), (0, 1)): Fraction(1, 5),
        ((0, 1), (2, 1)): Fraction(1, 20),
        ((1, 0), (1, 0)): Fraction(1, 5),
        ((1, 0), (2, 0)): Fraction(1, 40),
        ((1, 0), (2, 2)): Fraction(1, 40),
        ((1, 1), (1, 1)): Fraction(1, 10),
        ((1, 1), (2, 1)): Fraction(3, 20),
    }
    w1_side = verify_coupling(coup_side, mu_side, nu)
    energy_side = sum(mu_side.get(w, Fraction(0)) * jump[w] for w in words)
    def_side = (lambda_2 - energy_side) / lambda_2
    es_side = w1_side * (1 + def_side)
    assert dual_side == w1_side == Fraction(3, 4)
    assert es_side == Fraction(9, 8)

    # 5. Local jitter at scale 1/10: exact rational energy = 1/2 and W_1 = 27/500
    K = [
        [Fraction(9, 10), Fraction(1, 10), Fraction(0)],
        [Fraction(1, 20), Fraction(9, 10), Fraction(1, 20)],
        [Fraction(0), Fraction(1, 10), Fraction(9, 10)],
    ]
    nu_jit = {
        (i, j): sum(nu[(a, b)] * K[a][i] * K[b][j] for a in range(3) for b in range(3))
        for (i, j) in words
    }
    energy_jit = sum(nu_jit[w] * jump[w] for w in words)
    assert energy_jit == Fraction(1, 2), f"Expected 1/2, got {energy_jit}"

    f_jit = {
        (0, 0): Fraction(1, 2),
        (0, 1): Fraction(0),
        (0, 2): Fraction(1, 2),
        (1, 0): Fraction(0),
        (1, 1): Fraction(1, 2),
        (1, 2): Fraction(0),
        (2, 0): Fraction(1, 2),
        (2, 1): Fraction(0),
        (2, 2): Fraction(1, 2),
    }
    verify_1_lipschitz(f_jit)
    dual_jit = sum(f_jit[w] * (nu_jit[w] - nu[w]) for w in words)

    coup_jit = {(w, w): min(nu_jit[w], nu[w]) for w in words}
    coup_jit[((0, 0), (0, 1))] = Fraction(27, 2000)
    coup_jit[((0, 2), (0, 1))] = Fraction(27, 2000)
    coup_jit[((2, 0), (2, 1))] = Fraction(27, 2000)
    coup_jit[((2, 2), (2, 1))] = Fraction(27, 2000)
    coup_jit[((1, 1), (1, 0))] = Fraction(27, 1000)
    coup_jit[((1, 1), (1, 2))] = Fraction(27, 1000)
    w1_jit = verify_coupling(coup_jit, nu_jit, nu)
    assert dual_jit == w1_jit == Fraction(27, 500), f"Expected 27/500, got {w1_jit}"
    print(
        "2. Verified all 5 exact rational primal-dual certificates in Section 9.2 (0 residual)."
    )

    # Section 9.4 Finite-window Monte Carlo verification (N = 400, seed = 7)
    rng = random.Random(7)
    N = 400
    N_prime = N - 1
    num_trials = 200  # Fast CI check separating baseline/jitter from structured drift
    nu_float = [float(nu[w]) for w in words]
    jump_float = [float(jump[w]) for w in words]
    dist_mat = [[float(d_pair(u, v)) for v in words] for u in words]

    def sample_baseline_window():
        u = rng.random()
        s = 0 if u < 0.25 else (1 if u < 0.75 else 2)
        seq = [s]
        for _ in range(N - 1):
            r = rng.random()
            if s == 1:
                s = 0 if r < 0.4 else (1 if r < 0.6 else 2)
            else:
                s = 0 if r < 0.1 else (1 if r < 0.9 else 2)
            seq.append(s)
        return seq

    def score_seq(seq):
        counts = [0] * 9
        for t in range(N_prime):
            counts[seq[t] * 3 + seq[t + 1]] += 1
        energy = sum(counts[w] * jump_float[w] for w in range(9)) / N_prime
        defect = max(0.0, (0.5 - energy) / 0.5)
        w1 = sinkhorn_w1(counts, N_prime, nu_float, dist_mat)
        return w1 * (1.0 + defect)

    base_scores = [score_seq(sample_baseline_window()) for _ in range(num_trials)]
    side_scores = [score_seq([rng.randint(0, 1) for _ in range(N)]) for _ in range(50)]
    alt_score = score_seq([2 if (t & 1) else 0 for t in range(N)])

    max_base = max(base_scores)
    min_side = min(side_scores)
    assert max_base < 0.35, f"Baseline max too high: {max_base}"
    assert abs(alt_score - 0.80) < 0.01, f"Alternation score mismatch: {alt_score}"
    assert min_side > 0.85, f"One-sided min too low: {min_side}"
    print(
        f"3. Verified finite-window separation (N={N}): baseline max={max_base:.4f} << alt={alt_score:.4f} < one-sided min={min_side:.4f}."
    )


if __name__ == "__main__":
    main()
