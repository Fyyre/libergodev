# Non-Asymptotic Symbolic Ergodeviation (`libergodev`)

## James J. (Fyyre) - fyyre @ fyyre . net

A mathematical specification and zero-allocation C99 streaming engine for **Symbolic Ergodeviation ($`\mathcal{E}_S`$)**—a finite-window dynamical functional designed to separate structured state-space drift from benign high-entropy jitter in discrete event telemetry.

---

## The Problem with Classical Telemetry Invariants

When monitoring discrete event streams—such as Hardware Performance Counter (HPC) deltas, syscall traces, kernel state transitions, or network flow buckets—classical invariants fail in two predictable ways:

1. **Entropy & Lempel-Ziv are blind to geometry:** Two processes supported on disjoint regions of state space (or an attacker pacing a polling loop to match baseline entropy) can share the exact same Shannon entropy rate, Kolmogorov-Sinai entropy, Lempel-Ziv complexity, and correlation dimension while occupying completely different transition sub-manifolds (Section 3.1).
2. **Total Variation & KL divergence explode on local jitter:** Under a discrete metric, harmless scheduler or counter noise that nudges a symbol into an adjacent quantizer bin looks like a novel distribution.
3. **Naive spectral penalties fail on finite windows:** Penalizing by the algebraic connectivity of the empirical visit graph diverges on non-repeating windows before the cover time ($`\lambda_2(P_n) \le 12/(n(n+1))`$), while a two-sided ratio of Fiedler energies actively *discounts* anomalies that oscillate rapidly across rare graph cuts (shrinking a transport cost of $`4/5`$down to$`4/17`$, Section 6.2).

---

## The Formulation

Let $`\nu_\ell`$be the reference$`\ell`$-word marginal of a stationary exponentially mixing baseline on a finite alphabet$`A`$equipped with ground metric$`\delta`$and exponentially discounted word metric$`d_\theta`$. Let $`\lambda_2 > 0`$and$`\psi`$be the algebraic connectivity and unit Fiedler vector of the flux-weighted Rauzy graph of allowed transitions under$`\nu`$.  On a finite sliding window of length $`N`$($`N' = N - \ell + 1`$words) with empirical$`\ell`$-word measure$`\hat\mu_{N,\ell}`$and empirical Fiedler Dirichlet energy$`\hat\Delta_N`$, **Symbolic Ergodeviation** is defined for $`p \ge 1, \eta \ge 0`$ as:

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

The C99 library separates corpus fitting from online sliding-window ingestion so that per-event updates run in constant time with **zero heap allocations**:

* **Zero `malloc` / `free`:** All state lives inside a caller-allocated `ergodev_ctx_t` struct (~34 KB default footprint), suitable for ring-0 callbacks, eBPF user-space ring consumers, or embedded telemetry daemons.
* **$O(1)$ Per-Symbol Ingestion (`ergodev_push`):** Pushing a new symbol updates the base-$`m`$ rolling word index, increments the incoming $`\ell`$-word counter, decrements the outgoing $`\ell`$-word counter from $`N`$ ticks ago, and updates the rolling Fiedler Dirichlet energy via two LUT lookups in **< 10 CPU cycles**.
* **Stack-Only $`W_1`$ Evaluation (`ergodev_score_sinkhorn`):** Evaluates entropic optimal transport over the $`m^\ell`$ word space on demand without dynamic memory allocation.

---

## Streaming Benchmark (`make bench`)

Evaluating on a 10,000,000-event stream (query stride = 16 ticks) across window lengths $`N = 128`$to$`16,384`$:

* **`ergodev_push_fast()` is 100% integer-only (zero FPU instructions):** Transition indexing uses a precomputed 1-cycle state lookup table (`step_lut`), and ring-buffer updates use pure integer increments—making per-event ingestion safe for ring-0 / kernel contexts without `XSAVE`/`FXSAVE` FPU state-save overhead and completely immune to floating-point rounding drift.
* **Constant $`O(1)`$ latency (~2.2 ns / event):** Throughput stays pinned above **400 Million events/sec** on a single core regardless of window horizon $`N`$.

```text
Window N   | Method                 | Throughput   | Latency    | Relative Speed
---------------------------------------------------------------------------------------
N = 128    | Naive Window Rescan    |   181.2 M/s  |   5.52 ns  | [######                      ]
N = 128    | libergodev O(1) Ring   |   393.8 M/s  |   2.54 ns  | [#############               ] (2.2x faster)
---------------------------------------------------------------------------------------
N = 512    | Naive Window Rescan    |    52.9 M/s  |  18.90 ns  | [#                           ]
N = 512    | libergodev O(1) Ring   |   410.1 M/s  |   2.44 ns  | [##############              ] (7.7x faster)
---------------------------------------------------------------------------------------
N = 2048   | Naive Window Rescan    |    13.7 M/s  |  73.07 ns  | [#                           ]
N = 2048   | libergodev O(1) Ring   |   416.5 M/s  |   2.40 ns  | [##############              ] (30.4x faster)
---------------------------------------------------------------------------------------
N = 8192   | Naive Window Rescan    |     3.5 M/s  | 288.45 ns  | [#                           ]
N = 8192   | libergodev O(1) Ring   |   412.0 M/s  |   2.43 ns  | [##############              ] (118.8x faster)
---------------------------------------------------------------------------------------
N = 16384  | Naive Window Rescan    |     1.7 M/s  | 574.17 ns  | [#                           ]
N = 16384  | libergodev O(1) Ring   |   458.8 M/s  |   2.18 ns  | [################            ] (263.4x faster)
---------------------------------------------------------------------------------------
```

## Quick Start

```bash
# 1. Build the C99 shared library, static archive, and self-test binary
make

# 2. Run both the C99 streaming self-test and Python exact primal-dual verifier
make check

# 3. Compile the LaTeX specification PDF (requires pdflatex)
make pdf

```

---

## Repository Layout

```text
├── Makefile                        # Builds C99 library, runs tests, compiles LaTeX PDF
├── README.md                       # Overview and quick-start guide
├── SYMBOLIC_ERGODEVIATION.md       # Full Markdown mathematical specification
├── docs/
│   └── symbolic_ergodeviation.tex  # Publication-ready AMS LaTeX manuscript
├── include/
│   └── ergodev.h                   # Zero-allocation C99 API header
├── src/
│   └── ergodev.c                   # O(1) ring-buffer engine + built-in self-test
└── python/
    └── verify_es.py                # Exact rational primal-dual verifier & Monte Carlo suite

```
