# What was fixed

These changes follow the six items in `FINDINGS.MD`, Section 9, plus a later change that put the streaming speed back. The functional is the same. The engine now refuses a transport cost that is not a coupling, and the 10-million-event benchmark no longer runs that solver on every query.

`make check` passes. `make bench` times the ring against a full window rescan.

## 1. The verifier and the manuscript

`python/verify_es.py` used to stop Sinkhorn after 150 iterations at a very small regularization. On a one-sided window that returned a score of about 0.61 against a true value near 1.13, so `make check` failed. The Monte Carlo now returns a cost only after the plan's row and column residual is under \(5\times 10^{-3}\). With that check, a length-400 run separates as baseline max 0.158, alternation 0.800, one-sided min 0.999.

The joint bound in Remark 8.3 was described as crossing half the alternation score near \(N = 16400\). The product actually written there,

\[
(B_N + t)\left(1 + \eta\frac{s}{\lambda_2}\right),
\]

crosses \(2/5\) at **\(N = 27971\)**. At \(N = 16400\) it is still about 0.565. The manuscript states that.

## 2. A score that is a coupling, and a screen that is a bound

`ergodev_score` runs log-domain Sinkhorn and stops when the marginal residual is under the caller's tolerance. If the iteration cap is hit first, \(W_1\) and \(\mathcal{E}_S\) are NaN. An unfinished plan is not reported as a transport cost. One iteration at a tiny regularization on a hard window is the regression check for that.

`ergodev_screen` skips the solver when the Lemma 8.1 upper bound is already under the alarm level:

\[
\mathcal{E}_S \le \frac{D}{2}\,\|\hat\mu-\nu\|_1\,(1 + \eta\,\mathrm{defect}).
\]

`ergodev_es_upper` computes that bound in one pass over the counts. The old bench shortcut, \(0.5\|\hat\mu-\nu\|_1(1+\eta\,\mathrm{defect})\), is gone. That expression is the discrete-metric identity, and this ground metric is not discrete. On a mixture inside the old skip region it reported 0.181 for a true score of 0.278.

## 3. Fitting a corpus

`ergodev_fit` takes a symbol stream and builds the reference the score needs: the word marginal, algebraic connectivity, the spectral gap, the Dobrushin coefficient, and the largest squared Fiedler jump. `q_normalize = 1` uses the stationary-mass vector, \(L\psi = \gamma\,\mathrm{diag}(q)\,\psi\), instead of the unnormalized unit vector. That matters when one class is rare: the unnormalized vector spends its norm on that class.

`ergodev_init` rejects a jump table whose mean under \(\nu\) is not \(\lambda_2\). A scaled table used to zero the defect with no error.

## 4. Thresholds

`ergodev_bounds` returns three sufficient conditions for the worked chain:

| Bound | What it is | Sample size where it separates |
|---|---|---|
| Azuma | Theorem 8.4 | \(N = 56977\) for the half-alternation bar \(2/5\) |
| Joint product | Remark 8.3 as written | \(N = 27971\) for the same bar |
| Freedman | variance proxy \(\mathrm{Var}(W_1) \le B_N^2\) | \(N = 9103\) for the alternation score \(4/5\) |

Below those sizes the operating cutoff is `ergodev_empirical_quantile` on held-out baseline scores. The proved bounds stay valid. They are an alarm only once they sit under the alternative.

## 5. Where the context lives

The distance matrix and the Sinkhorn workspace are no longer inside the hot context. The caller owns an `ergodev_scratch_t`. The context holds the ring, the counts, the jump table, and, for a dense alphabet, the step table. Push and score do not call `malloc`.

One `ergodev_thread_slot_t` per thread id carries two contexts:

- `per_event`: one symbol per syscall. A pure rate change does not move the score.
- `per_quantum`: one symbol per fixed clock. Symbol 0 means idle, so volume is visible.

`ergodev_thread_push_event` and `ergodev_thread_push_quantum` only ingest. Scoring runs on the scratch buffer, off the observation callback.

## 6. Cycles the bigram score cannot see, and longer words

`ergodev_alarm_union` fires if \(\mathcal{E}_S\) exceeds its cutoff or the symbol window has a period at most \(\tau_Z\). The false-alarm budget is split across the two arms. The 40-symbol Eulerian cycle of the Section 9 chain has circular pair law equal to the baseline, so the transport score on a length-400 repeat stays under 0.08 and the Lemma 8.1 screen says it cannot alarm. Its minimum period is 40, so the period arm fires. A quiet aperiodic window does not.

If \(m^\ell > 128\) but the corpus uses at most 127 distinct words, the fit keeps that support plus one sink for everything else. A 6-class nearest-neighbor chain at block length 3 (216 possible triples, 45 seen) fits and scores. The dense cap is still 128 words.

## Speed, after the first cut of the benchmark

The first benchmark that called certified Sinkhorn on every 16th event lost the flat ring latency. Windows that had not settled spent thousands of iterations, and the per-event cost rose from about 3 ns to 8–16 ns.

The timed loop is the ring push plus the Lemma 8.1 upper bound, every 16 symbols. Certified Sinkhorn is still required, once, before that loop: the alternation must converge to 0.80, and the upper bound must sit above that score. It is not inside the 10 million events.

On this host, after that split:

| Window | Ring | Against a full rescan |
|---|---|---|
| 128 | 2.82 ns (355 M/s) | 1.5× |
| 512 | 2.67 ns (374 M/s) | 5.0× |
| 2048 | 2.66 ns (376 M/s) | 19.2× |
| 8192 | 2.59 ns (386 M/s) | 77.3× |
| 16384 | 2.61 ns (384 M/s) | 153.2× |

The ring latency does not grow with the window. Hot pointers for the lookup table, the ring, and the counts are filled once at init so the push does not recompute byte offsets on every symbol.

