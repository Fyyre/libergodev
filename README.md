# Non-Asymptotic Symbolic Ergodeviation (`libergodev`)

## James J. (Fyyre) - fyyre @ fyyre . net

A mathematical specification and zero-allocation C99 streaming engine for **Symbolic Ergodeviation ($`\mathcal{E}_S`$)** a finite-window dynamical functional designed to separate structured state-space drift from benign high-entropy jitter in discrete event telemetry.

---

## The Problem with Classical Telemetry Invariants

When monitoring discrete event streams—such as Hardware Performance Counter (HPC) deltas, syscall traces, kernel state transitions, or network flow buckets—classical invariants fail in two predictable ways:

1. **Entropy & Lempel-Ziv are blind to geometry:** Two processes supported on disjoint regions of state space (or an attacker pacing a polling loop to match baseline entropy) can share the exact same Shannon entropy rate, Kolmogorov-Sinai entropy, Lempel-Ziv complexity, and correlation dimension while occupying completely different transition sub-manifolds (Section 3.1).

2. **Total Variation & KL divergence explode on local jitter:** Under a discrete metric, harmless scheduler or counter noise that nudges a symbol into an adjacent quantizer bin looks like a novel distribution.

3. **Naive spectral penalties fail on finite windows:** Penalizing by the algebraic connectivity of the empirical visit graph diverges on non-repeating windows before the cover time ($`\lambda_2(P_n) \le 12/(n(n+1))`$), while a two-sided ratio of Fiedler energies actively *discounts* anomalies that oscillate rapidly across rare graph cuts (shrinking a transport cost of $`4/5`$ down to $`4/17`$, Section 6.2).

---

## The Formulation

Let $`\nu_\ell`$ be the reference $`\ell`$ word marginal of a stationary exponentially mixing baseline on a finite alphabet $`A`$ equipped with ground metric $`\delta`$ and exponentially discounted word metric $`d_\theta`$. Let $`\lambda_2 > 0`$ and $`\psi`$ be the algebraic connectivity and unit Fiedler vector of the flux-weighted Rauzy graph of allowed transitions under $`\nu`$. On a finite sliding window of length $`N`$($`N' = N - \ell + 1`$words) with empirical $`\ell`$ word measure $`\hat\mu_{N,\ell}`$ and empirical Fiedler Dirichlet energy $`\hat\Delta_N`$, **Symbolic Ergodeviation** is defined for $`p \ge 1, \eta \ge 0`$ as:

```math
\mathcal{E}_S = W_p(\hat\mu_{N,\ell}, \nu_\ell) \left( 1 + \eta \frac{(\lambda_2 - \hat\Delta_N)_+}{\lambda_2} \right)
```

### Proved Properties (Sections 7–8)

* **Exact Zero Set (Prop. 7.1):** $`\mathcal{E}_S = 0 \iff \hat\mu_{N,\ell} = \nu_\ell`$. The spectral factor never shifts the zero set.
* **Local Jitter Stability (Prop. 7.3):** Perturbing symbols by at most $`\eta_{\mathrm{loc}}`$ in the ground metric $`\delta`$ increases $`\mathcal{E}_S`$ by at most $`(1+\eta)\eta_{\mathrm{loc}}\frac{1-\theta^\ell}{1-\theta}`$, regardless of how much Shannon entropy the jitter injects.
* **No-Hiding & Bottleneck Amplification (Prop. 7.4):** Any drift at distance $`W_p \ge \tau`$ scores at least $`\tau`$. If the drift is confined to one side of the Fiedler bottleneck ($`\hat\Delta_N < \lambda_2`$), the score is strictly amplified up to $`(1+\eta)\tau`$.
* **Non-Asymptotic Concentration & Spectral Stability (Thm. 8.4, Cor. 8.5, Rem. 8.3):** Explicit finite-$`N`$ Azuma-Hoeffding tail bounds under Dobrushin contraction ($`\kappa < 1`$) and Davis-Kahan $`\sin\theta`$ bounds when $`\nu`$ and $`\psi`$ are estimated from a finite baseline corpus of length $`N_0`$.

---

## Certified Separation on the Worked 3-State Chain (Section 9)

For the 3-state baseline chain with stationary law $`q = (1/4, 1/2, 1/4)`$, Dobrushin coefficient $`\kappa = 3/5`$, and Rauzy algebraic connectivity $`\lambda_2 = 1/2`$ ($`\ell = 2, \eta = 1`$), every population value below is an **exact rational number** certified by an explicit primal transport coupling and a 1-Lipschitz Kantorovich-Rubinstein dual witness:

| Process / Law | Entropy Rate (bits) | $`W_1`$ Transport | Fiedler Defect | $`\mathcal{E}_S`$ (Population) | Finite Window Mean ($`N=400`$) |
| --- | --- | --- | --- | --- | --- |
| **Baseline $`\nu`$** | $`1.222`$ | $`0`$ | $`0`$ | **$`0`$** | $`0.0511`$ (99th: $`0.1479`$) |
| **Local Jitter (scale $`1/10`$)** | $`\ge 1.222`$ ($`H_1 = 1.5`$) | $`27/500 = 0.054`$ | $`0`$ | **$`0.054`$** | $`0.0769`$ (99th: $`0.1565`$) |
| **i.i.d. on Cut $`\{0,2\}`$** | $`1.000`$ | $`3/4 = 0.750`$ | $`0`$ | **$`3/4 = 0.750`$** | $`0.7500`$ |
| **Rare-Cut Alternation $`0 \leftrightarrow 2`$** | $`0.000`$ | $`4/5 = 0.800`$ | $`0`$ | **$`4/5 = 0.800`$** | $`0.8000`$ |
| **i.i.d. One-Sided $`\{0,1\}`$** | $`1.000`$ | $`3/4 = 0.750`$ | $`1/2 = 0.500`$ | **$`9/8 = 1.125`$** | $`1.1230`$ |
| **Constant State $`1`$** | $`0.000`$ | $`3/4 = 0.750`$ | $`1`$ | **$`3/2 = 1.500`$** | $`1.5000`$ |

Notice rows 3 and 5: both have entropy rate **1 bit** and transport cost **$`3/4`$**, but the process trapped on one side of the Fiedler cut (`{0,1}`) fails to carry the baseline Dirichlet energy and is penalized from **$`3/4`$up to$`9/8`$**. Meanwhile, high-entropy local jitter (`0.054`) sits an order of magnitude below every structured deviation.

---

## Streaming C99 Engine Architecture (`include/ergodev.h`, `src/ergodev.c`)

The C99 library separates corpus fitting from online sliding-window ingestion. The hot context is caller-allocated (`ergodev_ctx_bytes`) and holds the ring, the counts, the jump table, and, in the dense case, the step LUT. The distance matrix and the Sinkhorn workspace live in a caller-owned `ergodev_scratch_t`. Push and score do not call `malloc`.

* **Per-symbol ingestion (`ergodev_push_fast`):** the dense path is integer-only. It updates the rolling word, swaps one ring slot, and updates two counters. Dirichlet energy is the count-weighted mean of the jump table, evaluated when a score is requested.
* **Corpus fit (`ergodev_fit`):** estimates the word marginal, builds the flux Laplacian, and stores algebraic connectivity, the spectral gap, the Dobrushin coefficient, and the maximum squared Fiedler jump. `q_normalize = 1` selects the stationary-mass Fiedler vector. `ergodev_init` rejects a jump table whose mean under the reference is not `lambda_2`.
* **Certified transport (`ergodev_score`):** log-domain Sinkhorn returns a cost only after the marginal residual is under the caller's tolerance. Otherwise the cost is NaN. `ergodev_screen` skips the solver when the Lemma 8.1 upper bound is already under the alarm level.
* **Thresholds:** `ergodev_bounds` returns the Azuma sufficient condition, a Freedman threshold with variance proxy `B_N^2`, and the Remark 8.3 joint product. At moderate `N` the operating cutoff is `ergodev_empirical_quantile` on held-out baseline scores.
* **Sparse words:** if `m^ell` exceeds 128 but the corpus uses at most 127 distinct words, the fit keeps that support plus one sink. That is how block length 3 fits a taxonomy larger than five symbols.
* **Union alarm (`ergodev_alarm_union`):** fires when the transport score exceeds its cutoff or the symbol window has a period at most `tau_period`. Split the false-alarm budget across the two arms.

### Placement

One `ergodev_thread_slot_t` per thread id, with two contexts. `per_event` takes one symbol per syscall, so a pure rate change does not move the score. `per_quantum` takes one symbol per fixed clock, with symbol 0 meaning idle, so volume is visible. The push helpers only ingest. Score on the caller scratch, off the observation callback.

---

## Streaming Benchmark (`make bench`)

`make bench` times 10,000,000 events, query stride 16, windows 128 through 16,384. Each query is the Lemma 8.1 upper bound, one pass over the word counts. Certified Sinkhorn is checked once on the alternation (it must converge to 0.80) and is not inside the timed loop. On this host the ring stays near 2.6 ns/event while the rescan grows with `N`:

```text
N = 128    ring 355 M/s   2.82 ns    1.5x
N = 512    ring 374 M/s   2.67 ns    5.0x
N = 2048   ring 376 M/s   2.66 ns   19.2x
N = 8192   ring 386 M/s   2.59 ns   77.3x
N = 16384  ring 384 M/s   2.61 ns  153.2x
```

## Quick Start

```bash
# 1. Build the C99 shared library, static archive, and self-test binary
make

# 2. Run both the C99 streaming self-test and Python exact primal-dual verifier
make check
```

---

## Repository Layout

```text
├── Makefile                        # Builds C99 library, runs tests
├── README.md                       # Overview and quick-start guide
├── CHANGES.MD						# Recent changes
├── docs/
│   ├── SYMBOLIC_ERGODEVIATION.md   # Mathematical specification
│   └── symbolic_ergodeviation.pdf  # Manuscript PDF
├── include/
│   └── ergodev.h                   # Caller-owned context and scratch API
├── src/
│   ├── ergodev.c                   # Ring, fit, score, bounds, alarm, self-test
│   └── bench_ergodev.c             # Ring vs rescan, Lemma 8.1 screen
└── python/
    └── verify_es.py                # Exact rational certificates and residual-checked Monte Carlo

```
