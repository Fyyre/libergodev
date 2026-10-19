# Non-Asymptotic Symbolic Ergodeviation

### Fyyre (James) - fyyre@fyyre.net

## Status

Symbolic ergodeviation, written \(\mathcal{E}_S\), is not a classical named invariant. Shannon entropy rate, Kolmogorov-Sinai entropy, topological entropy, Lyapunov exponents, correlation dimension, and Lempel-Ziv complexity are. Each of those is the wrong functional for the distinction drawn below, and Section 3 records explicit pairs on which they tie while a transport cost does not.

The ingredients used here are standard: a physical (SRB) measure and its symbolic marginal, the 1-Wasserstein metric, and the algebraic connectivity of a weighted Rauzy graph. What is specified and proved is one finite-window composite of those ingredients. The proofs below apply to that composite. They do not create a uniqueness theorem for the formula, and they do not assert that an operating-system kernel or a packet trace is an Axiom A attractor.

A second formula, the two-sided ratio of Fiedler energies, looks equally plausible from the phrase "penalized by algebraic connectivity" and is rejected in Section 6. On a deterministic rare orbit it shrinks a transport cost of \(4/5\) down to \(4/17\).

Population values in the worked example are exact rational numbers, certified by an explicit coupling and a 1-Lipschitz dual witness. Monte Carlo numbers are labeled as such and come from `scratch/verify_es.py`.

## 1. Claim

Let a baseline process be a finite-state exponentially mixing shift, with law \(\nu\), thought of as the symbolic image of a physical measure when a smooth model exists. Let \(\hat\mu_{N,\ell}\) be the empirical law of length-\(\ell\) words in a window of \(N\) symbols. Let \(\lambda_2\) be the algebraic connectivity of the flux-weighted graph of transitions that \(\nu\) allows, and let \(\hat\Delta_N\) be the Dirichlet energy of a matching Fiedler vector along the window.

For a fixed weight \(\eta \ge 0\) and a block length \(\ell \ge 2\),

\[
\mathcal{E}_S
=
W_p(\hat\mu_{N,\ell}, \nu_\ell)
\left(
1 + \eta \frac{(\lambda_2 - \hat\Delta_N)_+}{\lambda_2}
\right).
\]

The theorems are:

1. \(\mathcal{E}_S = 0\) if and only if the window's \(\ell\)-marginal equals \(\nu_\ell\). The connectivity factor does not move the zero set.
2. Under the baseline, \(\mathcal{E}_S \to 0\) almost surely, and a fully explicit non-asymptotic tail holds for every finite \(N\) (Azuma, with a variance bound on the expectation). The constant is valid and not sharp.
3. Moving every symbol by at most \(\eta_{\mathrm{loc}}\) in the symbol metric raises \(\mathcal{E}_S\) by at most the factor \((1+\eta)\) times a term linear in \(\eta_{\mathrm{loc}}\). High entropy inside a small metric ball stays cheap.
4. Any alternative whose \(\ell\)-marginal is at least \(\tau\) away in \(W_p\) scores at least \(\tau\), whether or not its entropy rate has moved. If it also realizes only a fraction of the baseline Fiedler energy, the score is strictly larger than the pure transport cost.
5. Shannon rate, Kolmogorov-Sinai entropy, Lempel-Ziv rate, and correlation dimension can be held fixed while \(W_1\), hence \(\mathcal{E}_S\), is arbitrarily large. The Lyapunov exponent of the doubling map is the same on Lebesgue and on the fixed point.

## 2. Assumptions

**Alphabet and metric.** \(A\) is finite, \(|A| = m\). A metric \(\delta\) on \(A\) is fixed. The discrete metric \(\delta(a,b) = 1_{a \ne b}\) is allowed; the stability statement that local jitter is cheap needs \(\delta\) to come from the observation geometry (bins that sit next to each other in the counter space are close). On words of length \(\ell\),

\[
d_\theta(u,v) = \sum_{j=0}^{\ell-1} \theta^j \, \delta(u_j, v_j), \qquad \theta \in (0,1).
\]

Write \(D = \mathrm{diam}(A^\ell, d_\theta) < \infty\). If \(\delta\) is a metric then so is \(d_\theta\).

**Baseline.** \((s_n)_{n \ge 0}\) is a stationary ergodic Markov chain on \(A\) with transition matrix \(P\) and stationary row \(q\). The chain is irreducible and aperiodic on its support. The Dobrushin coefficient

\[
\kappa(P) = \max_{a,b \in A} \|P(a,\cdot) - P(b,\cdot)\|_{\mathrm{TV}}
\]

satisfies \(\kappa < 1\), where \(\|\mu-\nu\|_{\mathrm{TV}} = \tfrac12 \|\mu-\nu\|_1\). Then \(\|\mu P - \nu P\|_{\mathrm{TV}} \le \kappa \|\mu-\nu\|_{\mathrm{TV}}\). If only some power contracts, \(\kappa(P^m) = \kappa_m < 1\), the same argument applies with the constant recorded in Remark 8.1.

The reference marginal \(\nu_\ell\) is the law of \((s_0,\ldots,s_{\ell-1})\). The transition tensor of an order-1 chain is the case \(\ell = 2\): \(\nu_2(a,b) = q_a P_{ab}\).

**Where SRB enters, and where it does not.** If the baseline is a \(C^{1+\alpha}\) Axiom A attractor and the partition is a Markov partition, Bowen's theorem supplies a unique physical measure \(\mu_{\mathrm{SRB}}\), invariant, with conditional measures absolutely continuous on unstable manifolds, and the symbolic process is a Hoelder Gibbs state for \(\varphi = -\log |\det Df|_{E^u}|\). Exponential mixing of Hoelder observables gives some \(m\) with \(\kappa_m < 1\). In that case \(\nu = \pi_* \mu_{\mathrm{SRB}}\). The measure of maximal entropy is the same object only when \(\varphi\) is cohomologous to a constant. Nothing below requires that identification.

A performance-counter trace or a packet trace does not inherit an SRB theorem from this citation. For those series, \(\nu\) is the law of a trusted baseline corpus. Every theorem in Sections 7 and 8 is a theorem about that \(\nu\). The name "SRB" is reserved for the smooth case just stated.

## 3. What the classical invariants tie, and what they do not

### 3.1. An entropy-blind pair

Let \(M \ge 2\) and let \(A = \{0,1,\ldots,M\}\) with \(\delta(i,j) = |i-j|\). Let \(\rho_0\) be i.i.d. uniform on \(\{0,1\}\) and \(\rho_1\) i.i.d. uniform on \(\{M-1, M\}\).

Both processes have Shannon entropy rate \(1\) bit. Both are Bernoulli, so the Kolmogorov-Sinai entropy of the shift under each measure is \(1\) bit. The Lempel-Ziv compression rate of a stationary ergodic source converges to its entropy rate, so both have Lempel-Ziv rate \(1\) bit. Each support is two points, so the correlation dimension is \(0\). The topological entropy of each support subshift is \(\log 2\), read off the language only in the limit of block length.

The identity coupling that sends \(0 \mapsto M-1\) and \(1 \mapsto M\) has cost \(M-1\). The function \(f(i) = i\) is 1-Lipschitz and

\[
\int f \,d\rho_1 - \int f \,d\rho_0 = \tfrac12\bigl((M-1)+M\bigr) - \tfrac12\bigl(0+1\bigr) = M-1.
\]

Kantorovich-Rubinstein gives \(W_1(\rho_0,\rho_1) = M-1\), which grows without bound as the two pairs are moved apart, with every entropy-type number held fixed.

On a finite metric space \(W_p \ge W_1\). Whenever \(\mathcal{E}_S\) is defined for a baseline equal to one of these laws, the factor in parentheses is at least 1, so the same pair is separated by \(\mathcal{E}_S\) by at least \(M-1\).

### 3.2. A Lyapunov-blind pair

On the doubling map \(T(x) = 2x \bmod 1\), one has \(|T'| = 2\) at every point, so every orbit has Lyapunov exponent \(\log 2\). Lebesgue measure is the SRB measure. Its symbolic image under the partition \(\{[0,1/2),[1/2,1)\}\) is fair Bernoulli, whose Kolmogorov-Sinai entropy equals \(\log 2\) by Pesin's identity. The Dirac mass at the fixed point \(0\) is invariant, has entropy \(0\), and is coded by the constant sequence \(000\ldots\). With the discrete symbol metric and \(\ell = 1\),

\[
W_1(\delta_0, \tfrac12\delta_0 + \tfrac12\delta_1) = \tfrac12.
\]

The exponent that numerical Lyapunov estimators target is the same on the two measures. Entropy sees this particular pair. Section 3.1 is the pair entropy does not see. The point of listing both is that no one of these invariants is a finite-window transport distance to a reference marginal.

### 3.3. Asymptotic grammar

Topological entropy is \(\lim_{n\to\infty} n^{-1} \log \#\{\text{words of length }n\}\). It is a property of the support language along an infinite orbit. A window of length \(N\) does not determine it, and two measures with the same support share it. It does not define a distance from \(\hat\mu_{N,\ell}\) to \(\nu_\ell\).

## 4. The graph whose connectivity is used

Let \(\ell \ge 2\) and \(k = \ell - 1\). Vertices are the \(k\)-words in the support of \(\nu_k\). For an order-1 chain, \(k = 1\): vertices are symbols.

The forbidden set is the set of \((k+1)\)-words with \(\nu\)-mass zero. The Rauzy graph keeps the complementary edges. Give the undirected edge \(\{u,v\}\) the stationary flux

\[
w(u,v) = \xi(u \to v) + \xi(v \to u),
\]

where \(\xi(a,b) = q_a P_{ab}\) in the order-1 case, and the analogous cylinder mass in general. Self-loops do not appear in \(w\).

The weighted Laplacian is \((L\psi)_u = \sum_v w(u,v)(\psi_u - \psi_v)\). Its kernel is the constants precisely because the chain is irreducible, so the underlying graph is connected. Algebraic connectivity \(\lambda_2 > 0\) is the second-smallest eigenvalue (Fiedler). The Rayleigh theorem gives, for a unit Fiedler vector \(\psi \perp \mathbf{1}\),

\[
\lambda_2 = \psi^\top L \psi = \sum_{\{u,v\}} w(u,v)\,(\psi_u - \psi_v)^2.
\]

The same quadratic form, written on ordered transitions, is the baseline energy

\[
\Delta_* = \sum_{u,a} \xi(u,a)\,(\psi(u)-\psi(ua))^2.
\]

Each undirected edge contributes once through \(w\), and \(\xi(\cdot\to\cdot)+\xi(\cdot\leftarrow\cdot) = w\), so

\[
\Delta_* = \lambda_2.
\]

If the eigenvalue has multiplicity \(m_2 > 1\), replace \((\psi_u-\psi_v)^2\) by \(\|\Psi(u)-\Psi(v)\|^2\) for an orthonormal basis of the eigenspace. The identity becomes \(\Delta_* = m_2 \lambda_2\), which is basis-independent, and the ratio \((\Delta_* - \hat\Delta)_+/\Delta_*\) is unchanged in form. The worked example has multiplicity one.

**Empirical energy.** On a window \(s_0,\ldots,s_{N-1}\), with \(N' = N-\ell+1\) consecutive \(\ell\)-words,

\[
\hat\Delta_N = \frac{1}{N'}\sum_{t=0}^{N'-1} c(s_t,\ldots,s_{t+\ell-1}),
\]

where \(c = (\psi(u)-\psi(u'))^2\) if the transition is legal, and \(c = 0\) if the word is forbidden. Legal mass that fails to cross a Fiedler gap, and illegal mass, both lower \(\hat\Delta_N\). Illegal mass is also charged by \(W_p\), because \(\nu\) gives it nothing to sit on.

Under a stationary draw from \(\nu\), every term has mean \(\Delta_*\), so \(\mathbb{E}[\hat\Delta_N] = \lambda_2\).

## 5. Definition

Let \(\hat\mu_{N,\ell}\) be the empirical measure on \(\ell\)-words in the window. For \(p \ge 1\) and \(\eta \ge 0\),

\[
\mathrm{defect}(\hat\Delta_N)
=
\frac{(\lambda_2 - \hat\Delta_N)_+}{\lambda_2},
\qquad
\mathcal{E}_S
=
W_p(\hat\mu_{N,\ell}, \nu_\ell)
\bigl(1 + \eta \,\mathrm{defect}(\hat\Delta_N)\bigr).
\]

The defect lies in \([0,1]\), so the factor lies in \([1, 1+\eta]\). In particular

\[
W_p(\hat\mu_{N,\ell}, \nu_\ell)
\le
\mathcal{E}_S
\le
(1+\eta)\, W_p(\hat\mu_{N,\ell}, \nu_\ell).
\]

Setting \(\eta = 0\) returns plain Wasserstein. The value \(\eta = 1\) used in Section 9 weights a total collapse of Fiedler energy equally with the transport cost: the factor equals 2 when \(\hat\Delta_N = 0\), and equals 1 whenever the window meets or exceeds the baseline energy.

On a space of diameter \(D\), \(p \le q\) implies \(W_p \le W_q \le D^{1-p/q} W_p^{p/q}\). Bounds proved for \(W_1\) transfer to \(W_p\) up to factors that depend only on \(D\) and \(p\). The concentration proof is written for \(p = 1\).

## 6. Two penalties that do not work

### 6.1. Algebraic connectivity of the visit graph

Let \(H_N\) be the unweighted graph whose vertices are the \(k\)-words actually seen, with an edge for each consecutive pair actually seen. A window that has not repeated a state induces a path \(P_n\) on \(n \le N\) vertices.

**Lemma (path test vector).** \(\lambda_2(P_n) \le 12/(n(n+1))\).

**Proof.** Let \(v_j = j - (n+1)/2\) for \(j = 1,\ldots,n\). Then \(\sum v_j = 0\) and

\[
\|v\|_2^2 = \frac{n(n^2-1)}{12}.
\]

Adjacent coordinates differ by 1, and there are \(n-1\) edges, so \(v^\top L v = n-1\). Rayleigh's theorem yields

\[
\lambda_2 \le \frac{n-1}{n(n^2-1)/12} = \frac{12}{n(n+1)}.
\]

Thus \(\lambda_2(H_N) \to 0\) along any non-repeating window, including a window drawn from the baseline on a large alphabet, long before the cover time. A factor \(1/\lambda_2(H_N)\) diverges on typical finite orbits. That factor is not a non-asymptotic penalty.

The Dirichlet form \(\hat\Delta_N\) avoids this. It is an average of a fixed bounded function of the chain, and it is unbiased for \(\lambda_2\) at every \(N\), before cover.

### 6.2. The two-sided ratio

The ratio

\[
\Pi^{\mathrm{two}} = \frac{\lambda_2(1+\varepsilon)}{\hat\Delta_N + \varepsilon\lambda_2}, \qquad \varepsilon \in (0,1],
\]

is near 1 when \(\hat\Delta_N\) is near \(\lambda_2\), and it is large when \(\hat\Delta_N\) is near 0. It is also strictly smaller than 1 whenever \(\hat\Delta_N > \lambda_2\). Any alternative that crosses the Fiedler cut more than the baseline does is discounted.

Section 9 computes a deterministic alternation on the rare cut of a 3-state chain: \(W_1 = 4/5\), \(\hat\Delta = 2\), \(\lambda_2 = 1/2\). With \(\varepsilon = 1/4\),

\[
\Pi^{\mathrm{two}} \cdot \tfrac45 = \tfrac{4}{17} \approx 0.2353.
\]

The one-sided definition scores the same orbit at \(4/5\). In an 800-window Monte Carlo of the baseline at \(N = 400\) (Section 9), the two-sided baseline scores had maximum \(0.202\). The anomaly is still above that maximum, and the gap has collapsed from several tenths to about \(0.03\). The one-sided score \(4/5\) remains well above the one-sided baseline maximum \(0.236\) from the matching experiment.

**No-hiding lemma.** If the prefactor is at least 1 for every window, then every alternative with \(W_p \ge \tau\) has \(\mathcal{E}_S \ge \tau\). The two-sided ratio violates the hypothesis. The one-sided factor \(1+\eta\,\mathrm{defect}\) satisfies it.

## 7. Elementary properties

**Proposition 7.1 (zero set).** For a population measure \(\rho\) on \(\ell\)-words, write \(\Delta(\rho)\) for the legal Dirichlet form and

\[
\mathcal{E}(\rho) = W_p(\rho, \nu_\ell)\bigl(1+\eta\,\mathrm{defect}(\Delta(\rho))\bigr).
\]

Then \(\mathcal{E}(\rho) = 0\) if and only if \(\rho = \nu_\ell\).

**Proof.** The prefactor is at most \(1+\eta\) and at least 1. The product vanishes if and only if \(W_p = 0\). Since \(d_\theta\) is a metric, \(W_p = 0\) if and only if the two probabilities agree. If \(\rho = \nu_\ell\) then \(\Delta(\rho) = \Delta_* = \lambda_2\), so the defect is zero, which matches.

**Proposition 7.2 (almost-sure limit).** For \(\nu\)-almost every baseline orbit, \(\mathcal{E}_S \to 0\).

**Proof.** The shift is ergodic and each cylinder indicator is integrable. Birkhoff's theorem gives \(\hat\mu_{N,\ell} \to \nu_\ell\) almost surely in total variation, hence in \(W_p\), on a finite set. The same theorem applied to the bounded function \(c\) gives \(\hat\Delta_N \to \lambda_2\), so the defect tends to 0 and the prefactor tends to 1.

**Proposition 7.3 (jitter).** Let \(K\) be a kernel on \(A\) with \(W_1(\delta_a, K(a,\cdot)) \le \eta_{\mathrm{loc}}\) for every \(a\), the Wasserstein cost being computed with \(\delta\). Let \(\rho\) be any law on \(A^{\mathbb{N}}\), and let \(\rho^K\) be the law after applying \(K\) independently at each coordinate. Then

\[
W_1(\rho^K_\ell, \nu_\ell)
\le
W_1(\rho_\ell, \nu_\ell) + \eta_{\mathrm{loc}}\frac{1-\theta^\ell}{1-\theta},
\]

and therefore

\[
\mathcal{E}(\rho^K_\ell)
\le
(1+\eta)\left(
W_1(\rho_\ell, \nu_\ell) + \eta_{\mathrm{loc}}\frac{1-\theta^\ell}{1-\theta}
\right).
\]

**Proof.** Couple each clean symbol to its noisy image at cost \(\le \eta_{\mathrm{loc}}\), independently. The resulting word coupling has \(d_\theta\)-cost at most \(\eta_{\mathrm{loc}} \sum_{j=0}^{\ell-1} \theta^j\). That is the first inequality, after the triangle inequality through \(\rho_\ell\). The second is the comparison \(\mathcal{E} \le (1+\eta) W_1\), which holds pathwise because the defect is at most 1. The prefactor cannot amplify local noise by more than \(1+\eta\).

Entropy of \(\rho^K\) can be large. The bound does not depend on that entropy. What it depends on is \(\eta_{\mathrm{loc}}\) in the symbol metric. With the discrete metric, a flip to any other bin costs 1, and the proposition correctly refuses to call that flip local.

**Proposition 7.4 (structured drift).** Let \(\rho\) be a probability on \(\ell\)-words with \(W_p(\rho, \nu_\ell) \ge \tau\). Then \(\mathcal{E}(\rho) \ge \tau\). If also \(\Delta(\rho) \le \alpha \lambda_2\) for some \(\alpha \in [0,1)\), then

\[
\mathcal{E}(\rho) \ge \tau \bigl(1 + \eta(1-\alpha)\bigr).
\]

**Proof.** The prefactor is at least 1, and it is at least \(1+\eta(1-\alpha)\) when the energy shortfall is at least \(1-\alpha\).

For a periodic sequence of period \(p\), the empirical \(\ell\)-marginal and \(\hat\Delta_N\) differ from the periodic occupation by at most \(O((p+\ell)/N)\) in total variation and in the bounded energy. For all large \(N\), \(\mathcal{E}_S\) exceeds half the population value of that periodic occupation.

## 8. Finite-sample bound

Throughout this section the window is stationary under \(\nu\), \(\kappa = \kappa(P) < 1\), \(p = 1\), and \(N' = N - \ell + 1 \ge 1\).

**Lemma 8.1 (one-symbol Lipschitz constant of \(W_1\)).** Changing a single symbol \(s_j\) changes \(W_1(\hat\mu_{N,\ell}, \nu_\ell)\) by at most

\[
c = \frac{D \,\ell}{N'}.
\]

**Proof.** At most \(\ell\) windows cover site \(j\). Each of those windows contributes a count of 1. Replacing it moves one count from one word to another, so the \(\ell_1\) change in the count vector is at most \(2\ell\). Hence \(\|\hat\mu-\hat\mu'\|_1 \le 2\ell/N'\) and the total-variation distance is at most \(\ell/N'\). On a space of diameter \(D\), \(W_1 \le D \,\|\cdot\|_{\mathrm{TV}}\), because a maximal coupling realizes \(\mathbb{P}(X \ne Y) = \|\cdot\|_{\mathrm{TV}}\) and \(d_\theta \le D\).

**Lemma 8.2 (influence).** Let \(f\) be a function of \((s_0,\ldots,s_{N-1})\) with \(|f(s)-f(s')| \le \sum_{i: s_i \ne s'_i} c_i\). The Doob increments of \(f\) along the chain satisfy

\[
\bigl|\mathbb{E}[f \mid s_0,\ldots,s_j] - \mathbb{E}[f \mid s_0,\ldots,s_{j-1}]\bigr|
\le
\sum_{t \ge 0} c_{j+t}\, \kappa^{t}.
\]

**Proof.** Fix the past up to time \(j-1\). For two possible values \(a,b\) of \(s_j\), couple the two futures by the maximal coupling so that \(\mathbb{P}(s_{j+t}^{(a)} \ne s_{j+t}^{(b)}) \le \kappa^t\). This is the Dobrushin contraction iterated. The conditional expectations \(g(a), g(b)\) therefore differ by at most \(\sum_{t \ge 0} c_{j+t} \kappa^t\). The random variable \(g(s_j)\) and its conditional mean both lie between \(\min g\) and \(\max g\), so the Doob increment is at most that oscillation.

**Lemma 8.3 (variance of a word frequency).** Let \(J_t = 1_{\{s_t^{t+\ell-1} = w\}}\) and \(\nu(w) = \mathbb{E}[J_0]\). Then

\[
\mathrm{Var}\Bigl(\sum_{t=0}^{N'-1} J_t\Bigr)
\le
N' \,\nu(w) \left(2\ell - 1 + \frac{2\kappa}{1-\kappa}\right).
\]

**Proof.** \(\mathrm{Var}(J_t) \le \nu(w)\). For \(1 \le s < \ell\), \(|\mathrm{Cov}(J_0,J_s)| \le \nu(w)\). For \(s \ge \ell\), the second window is a function of the chain from time \(s\) onward, and the total variation between its conditional law given \(s_{\ell-1}\) and its stationary law is at most \(\kappa^{s-\ell+1}\), because path-space total variation on \(\ell\) steps is at most the total variation of the law of the first of those steps. Hence

\[
|\mathrm{Cov}(J_0, J_s)|
=
\bigl|\mathbb{E}\bigl[J_0\bigl(\mathbb{P}(J_s=1 \mid s_{\ell-1}) - \nu(w)\bigr)\bigr]\bigr|
\le
\nu(w)\, \kappa^{s-\ell+1}.
\]

Summing the geometric tail \(\sum_{q \ge 1} \kappa^q = \kappa/(1-\kappa)\) and using \(N'-s \le N'\) produces the bracket \(1 + 2(\ell-1) + 2\kappa/(1-\kappa) = 2\ell-1 + 2\kappa/(1-\kappa)\).

**Theorem 8.4 (expectation and tail).** Let

\[
B_N = \frac{D}{2}\sqrt{\frac{m^\ell}{N'}\left(2\ell - 1 + \frac{2\kappa}{1-\kappa}\right)}.
\]

Then \(\mathbb{E}[W_1(\hat\mu_{N,\ell}, \nu_\ell)] \le B_N\). For every \(t > 0\),

\[
\mathbb{P}\bigl(W_1 - \mathbb{E}[W_1] \ge t\bigr)
\le
\exp\left(
-\frac{t^2 \, N'^2 (1-\kappa)^2}{2 N D^2 \ell^2}
\right),
\]

and therefore

\[
\mathbb{P}\bigl(\mathcal{E}_S \ge (1+\eta)(B_N + t)\bigr)
\le
\exp\left(
-\frac{t^2 \, N'^2 (1-\kappa)^2}{2 N D^2 \ell^2}
\right).
\]

**Proof.** Unbiasedness \(\mathbb{E}[\hat\mu(w)] = \nu(w)\) is stationarity of each window. By Lemma 8.3 and Cauchy-Schwarz,

\[
\mathbb{E}\|\hat\mu-\nu\|_1
\le
\sum_w \sqrt{\mathrm{Var}(\hat\mu(w))}
\le
\sqrt{\frac{m^\ell}{N'}\left(2\ell-1+\frac{2\kappa}{1-\kappa}\right)},
\]

where \(\sum_w \sqrt{\nu(w)} \le \sqrt{m^\ell}\). Multiply by \(D/2\) to pass from \(\ell_1\) to \(W_1\). That is \(B_N\).

For the tail, Lemma 8.1 and Lemma 8.2 give Doob increments bounded by \(c/(1-\kappa)\) with \(c = D\ell/N'\). Azuma's inequality on \(N\) increments is the exponential bound. The comparison \(\mathcal{E}_S \le (1+\eta) W_1\) lifts it to \(\mathcal{E}_S\). Replacing the unknown mean by the larger number \(B_N\) only enlarges the threshold, so the probability statement remains valid.

**Corollary 8.5 (estimated baseline).** Let \(\hat\nu\) be the \(\ell\)-marginal of an independent stationary corpus of length \(N_0\), with flux Laplacian \(\hat L\), algebraic connectivity \(\hat\lambda_2\), and unit Fiedler vector \(\hat\psi\) chosen with \(\langle \hat\psi, \psi \rangle \ge 0\). The triangle inequality gives

\[
W_1(\hat\mu, \nu) \le W_1(\hat\mu, \hat\nu) + W_1(\hat\nu, \nu),
\]

where the second summand obeys Theorem 8.4 with \(N_0\) in place of \(N\).

For the Fiedler penalty, let \(g = \min(\lambda_2,\, \lambda_3 - \lambda_2) > 0\) be the spectral gap of \(L\) when \(\lambda_2\) is simple (or \(\min(\lambda_2,\, \lambda_{2+m_2}-\lambda_2)\) for multiplicity \(m_2\) using the basis-independent projection \(\Psi\Psi^\top\) of Section 4). Writing \(E = \hat L - L\), each undirected edge weight changes by at most the sum of the two directed cylinder errors, so

\[
\|E\|_{\mathrm{op}} \le \|E\|_1 \le 4\,\|\hat\nu - \nu\|_1.
\]

By Weyl's inequality, \(|\hat\lambda_2 - \lambda_2| \le \|E\|_{\mathrm{op}}\). Provided \(\|E\|_{\mathrm{op}} < g/2\), the Davis-Kahan \(\sin\theta\) theorem gives

\[
\|\hat\psi - \psi\|_2 \le \sqrt{2}\,\sin\theta(\hat\psi, \psi) \le \frac{2\sqrt{2}\,\|E\|_{\mathrm{op}}}{g}.
\]

For any transition pair \((u,v)\), the identity \(a^2 - b^2 = (a-b)(a+b)\) and Cauchy-Schwarz on the incidence vector \(e_u - e_v\) yield

\[
\bigl|(\hat\psi_u - \hat\psi_v)^2 - (\psi_u - \psi_v)^2\bigr|
\le
\|\hat\psi - \psi\|_2 \sqrt{2} \cdot \|\hat\psi + \psi\|_2 \sqrt{2}
\le
4\,\|\hat\psi - \psi\|_2.
\]

Assuming \(\mathrm{supp}(\hat\nu) = \mathrm{supp}(\nu)\) (after corpus cover of the allowed transitions), the empirical energies \(\hat\Delta_N(\hat\psi)\) and \(\hat\Delta_N(\psi)\) on any test window therefore differ by at most \(4\|\hat\psi - \psi\|_2\). Since \(x \mapsto (1 - x/a)_+\) is \(1/a\)-Lipschitz in \(x\) and \(\hat\Delta_N \le 2\), the estimated defect and true defect satisfy

\[
\bigl|\widehat{\mathrm{defect}}(\hat\Delta_N(\hat\psi)) - \mathrm{defect}(\hat\Delta_N(\psi))\bigr|
\le
\frac{4\,\|\hat\psi - \psi\|_2 + \frac{\hat\Delta_N(\psi)}{\lambda_2}|\hat\lambda_2 - \lambda_2|}{\hat\lambda_2}
\le
\frac{\|E\|_{\mathrm{op}}}{\hat\lambda_2}\left(\frac{8\sqrt{2}}{g} + \frac{2}{\lambda_2}\right).
\]

Combined with \(\|\hat\nu - \nu\|_1\), which is bounded in expectation by \(2 B_{N_0}/D\) and in tail by Lemma 8.1 and Azuma on \(N_0\), both the transport shift and the spectral shift vanish at rate \(O(N_0^{-1/2})\) once \(N_0\) exceeds the cover time of allowed transitions.

**Remark 8.1 (\(m\)-step contraction).** If \(\kappa(P) = 1\) but \(\kappa(P^m) = \kappa_m < 1\), submultiplicativity of \(\kappa\) gives \(\|\delta_x P^{t} - \delta_y P^{t}\|_{\mathrm{TV}} \le \kappa_m^{\lfloor t/m \rfloor}\). The increment bound becomes \(c \, m / (1-\kappa_m)\), and the exponent in Theorem 8.4 replaces \((1-\kappa)\) by \((1-\kappa_m)/m\). Hoelder Gibbs states on a mixing subshift of finite type fall under this remark rather than under one-step Doeblin.

**Remark 8.2 (sharpness).** The exponent is of order \(N t^2 / (D \ell /(1-\kappa))^2\), which is the right shape, but Azuma's constant is far from the variance proxy. Section 9 records a chain where \(B_N + t\) at a conventional 5 percent level still exceeds the diameter scale for a few hundred samples, and only drops below half the alternation score near \(N = 57000\). Sharper Markov transportation inequalities exist; the argument here does not depend on their constants. The bound is a proof for every finite \(N\). It is a usable threshold only once \(N\) is large enough that \((1+\eta)(B_N+t)\) sits under the alternative of interest.

A non-stationary start at total-variation distance \(\varepsilon_0\) from \(q\) shifts each window marginal by at most \(\kappa^{t_0}\varepsilon_0\) after \(t_0\) burned symbols, which shifts \(W_1\) by at most \(D \kappa^{t_0}\varepsilon_0\). Burn a logarithmic prefix and the stationary theorem applies to the remainder.

**Remark 8.3 (joint concentration of \(W_1\) and \(\hat\Delta_N\)).** Theorem 8.4 bounds the prefactor by its deterministic worst case \(1+\eta\), which multiplies \(B_N\) by \(1+\eta\) even though \(\mathbb{E}_\nu[\hat\Delta_N] = \lambda_2\). Because each legal transition jump satisfies \((\psi_u - \psi_v)^2 \le 2\|\psi\|_2^2 = 2\), changing one symbol \(s_j\) alters at most \(\ell\) terms of \(\hat\Delta_N\), shifting \(\hat\Delta_N\) by at most \(c_\Delta = 2\ell/N'\). By Lemma 8.2 and Azuma's inequality on \(\lambda_2 - \hat\Delta_N\), for every \(s \in (0,\lambda_2]\),

\[
\mathbb{P}\bigl(\mathrm{defect}(\hat\Delta_N) \ge s/\lambda_2\bigr)
=
\mathbb{P}\bigl(\lambda_2 - \hat\Delta_N \ge s\bigr)
\le
\exp\left(
-\frac{s^2 \, N'^2 (1-\kappa)^2}{8 N \ell^2}
\right).
\]

By the union bound, with probability at least \(1 - \delta_1 - \delta_2\),

\[
\mathcal{E}_S
\le
(B_N + t)\left(1 + \frac{\eta s}{\lambda_2}\right),
\]

where \(t\) and \(s\) are the Azuma deviations at levels \(\delta_1\) and \(\delta_2\). Since \(B_N\) and \(s\) both decay as \(O(N^{-1/2})\), the leading threshold term is \(B_N + O(N^{-1/2})\) rather than \((1+\eta)B_N\), reducing the sample size required to separate a fixed alternative by a factor of up to \((1+\eta)^2\).

## 9. Worked chain

### 9.1. Baseline

\[P = \begin{pmatrix} 1/10 & 4/5 & 1/10 \\ 2/5 & 1/5 & 2/5 \\ 1/10 & 4/5 & 1/10 \end{pmatrix}, \qquad q = \bigl(\tfrac14,\tfrac12,\tfrac14\bigr).\]

The row \(q\) is stationary by direct multiplication. Symbol entropy is

\[-\tfrac12 \log_2 \tfrac14 - \tfrac12 \log_2 \tfrac12 = \tfrac32 \ \text{bits}.\]

The entropy rate \(-\sum q_i P_{ij} \log_2 P_{ij}\) evaluates to \(1.2219280948873625\) bits in binary64. One-step Dobrushin coefficient: rows 0 and 2 agree, and

\[\Vert{}P(0,\cdot)-P(1,\cdot)\Vert{}_{\mathrm{TV}} = \tfrac12(0.3+0.6+0.3) = \tfrac35.\]

So \(\kappa = 3/5\). Ground metric on symbols \(\delta(i,j) = \vert{}i-j\vert{}\), on pairs \(d = \vert{}i-k\vert{} + \tfrac12 \vert{}j-l\vert{}\), so \(D = 3\). Penalty weight \(\eta = 1\).

Fluxes of the off-diagonal edges: \(w(0,1) = w(1,2) = 2/5\), \(w(0,2) = 1/20\). The weighted Laplacian is

\[L = \begin{pmatrix} 9/20 & -2/5 & -1/20 \\ -2/5 & 4/5 & -2/5 \\ -1/20 & -2/5 & 9/20 \end{pmatrix}.\]

The unit vector \(\psi = (-1/\sqrt{2},\, 0,\, 1/\sqrt{2})\) satisfies \(L\psi = \tfrac12 \psi\) by multiplication, and \(\psi \perp \mathbf{1}\). The quadratic form reproduces the same eigenvalue:

\[\frac25 \cdot \frac12 + \frac25 \cdot \frac12 + \frac{1}{20} \cdot 2 = \frac12.\]

So \(\lambda_2 = \Delta_* = 1/2\). The trace of \(L\) is \(17/10\), so the remaining eigenvalue is \(17/10 - 1/2 = 6/5\). The eigenvalue \(\lambda_2\) is simple. Squared jumps: \((\psi_0-\psi_1)^2 = 1/2\), \((\psi_0-\psi_2)^2 = 2\), self-loops \(0\).

### 9.2. Five population alternatives

The edge-Lipschitz criterion is used below. The metric \(d\) is the shortest-path metric of the product-of-paths graph whose horizontal edges (first symbol changes by 1) have length 1 and whose vertical edges (second symbol changes by 1) have length \(1/2\). A function on that graph is 1-Lipschitz if and only if each edge changes the function by at most the edge length; the inequality along a general pair follows by summing a shortest path.

**Constant state 1.** The empirical pair measure is \(\delta_{(1,1)}\). Because the source is a Dirac, there is only one coupling, and

\begin{align*}
W_1
&= \sum_{i,j} \nu(i,j)\, d\bigl((1,1),(i,j)\bigr) \\
&= 4\cdot\frac{1}{40}\cdot\frac32 + 2\cdot\frac15\cdot 1 + 2\cdot\frac15\cdot\frac12 \\
&= \frac{3}{20} + \frac25 + \frac15
= \frac34.
\end{align*}

The four corner words \((0,0),(0,2),(2,0),(2,2)\) each have mass \(1/40\) and distance \(3/2\); \((0,1)\) and \((2,1)\) each have mass \(1/5\) and distance 1; \((1,0)\) and \((1,2)\) each have mass \(1/5\) and distance \(1/2\). Energy \(\hat\Delta = 0\), defect \(= 1\), prefactor \(= 2\), so \(\mathcal{E}_S = 3/2\). Entropy rate \(0\).

**Alternation \(0,2,0,2,\ldots\).** Occupation \(\mu(0,2) = \mu(2,0) = 1/2\). Energy \(= 2 > \lambda_2\), so the one-sided defect is \(0\) and \(\mathcal{E}_S = W_1\).

Dual witness, rows indexed by the first symbol and columns by the second:

\[f = \begin{pmatrix} 0 & 1/2 & 1 \\ 0 & -1/2 & 0 \\ 1 & 1/2 & 0 \end{pmatrix}.\]

Adjacent horizontal differences are in \(\{0,1\}\) and adjacent vertical differences are all \(1/2\), so \(f\) is 1-Lipschitz. With \(\nu(0,0)=\nu(0,2)=\nu(2,0)=\nu(2,2)=1/40\), \(\nu(0,1)=\nu(1,0)=\nu(1,2)=\nu(2,1)=1/5\), \(\nu(1,1)=1/10\),

\begin{align*}
\sum f\,d(\mu-\nu)
&= 2\cdot\bigl(\tfrac12 - \tfrac{1}{40}\bigr)\cdot 1
+ 2\cdot\bigl(-\tfrac15\bigr)\cdot\tfrac12
+ \bigl(-\tfrac{1}{10}\bigr)\cdot\bigl(-\tfrac12\bigr) \\
&= \tfrac{19}{20} - \tfrac15 + \tfrac{1}{20}
= \tfrac45.
\end{align*}

So \(W_1 \ge 4/5\). The coupling that moves the \(\mu\)-mass at \((0,2)\) onto \(\nu\) by the amounts

\[\begin{align*} &(0,0):\tfrac{1}{40},\  (0,1):\tfrac15,\  (0,2):\tfrac{1}{40},\  (1,1):\tfrac{1}{20},\  (1,2):\tfrac15, \end{align*}\]

and the symmetric amounts from \((2,0)\) onto \((1,0),(1,1),(2,0),(2,1),(2,2)\), has both correct marginals (each \(\nu\)-atom is hit once, and each source has mass \(1/2\)). Its cost is

\[2\left(\frac{1}{40} + \frac{1}{10} + \frac{3}{40} + \frac15\right) = \frac45.\]

Hence \(W_1 = 4/5\) and \(\mathcal{E}_S = 4/5\). Entropy rate \(0\).

The two-sided ratio of Section 6.2 multiplies this by \(\lambda_2(1+\varepsilon)/(\hat\Delta+\varepsilon\lambda_2)\) at \(\varepsilon=1/4\):

\[\frac{(1/2)\cdot(5/4)}{2 + (1/4)\cdot(1/2)} \cdot \frac45 = \frac{4}{17}.\]

**Fair noise on the cut \(\{0,2\}\).** i.i.d. uniform on \(\{0,2\}\) has entropy rate \(1\) bit. Pair masses \(1/4\) on \((0,0),(0,2),(2,0),(2,2)\). Energy \(= 2\cdot(1/4)\cdot 2 = 1 > 1/2\), defect \(0\).

Dual witness

\[f = \begin{pmatrix} 0 & -1/2 & 0 \\ -1 & -3/2 & -1 \\ 0 & -1/2 & 0 \end{pmatrix}\]

has adjacent differences inside the edge lengths, so it is 1-Lipschitz. Then \(\mu(f) = 0\) and

\begin{align*}
\nu(f)
&= 2\cdot\frac15\cdot\bigl(-\tfrac12\bigr)
+ 2\cdot\frac15\cdot(-1)
+ \frac{1}{10}\cdot\bigl(-\tfrac32\bigr)
= -\frac{3}{4},
\end{align*}

so \((\mu-\nu)(f) = 3/4\). The matching coupling, written as source, target, mass, is

\begin{align*}
&(0,0)\to(0,0),\ \tfrac{1}{40};
&&(0,0)\to(0,1),\ \tfrac{7}{40};
&&(0,0)\to(1,1),\ \tfrac{1}{20};\\
&(0,2)\to(0,1),\ \tfrac{1}{40};
&&(0,2)\to(0,2),\ \tfrac{1}{40};
&&(0,2)\to(1,2),\ \tfrac15;\\
&(2,0)\to(1,0),\ \tfrac15;
&&(2,0)\to(2,0),\ \tfrac{1}{40};
&&(2,0)\to(2,1),\ \tfrac{1}{40};\\
&(2,2)\to(1,1),\ \tfrac{1}{20};
&&(2,2)\to(2,1),\ \tfrac{7}{40};
&&(2,2)\to(2,2),\ \tfrac{1}{40}.
\end{align*}

Row masses are \(1/4\). Column masses reproduce \(\nu\). The cost sums to

\[\frac{7}{80}+\frac{3}{40}+\frac{1}{80}+\frac15+\frac15+\frac{1}{80}+\frac{3}{40}+\frac{7}{80} = \frac34.\]

Hence \(W_1 = 3/4\) and \(\mathcal{E}_S = 3/4\).

**Fair noise on one side of the cut, \(\{0,1\}\).** Entropy rate again \(1\) bit. Energy:

\[\frac14 \cdot \frac12 + \frac14 \cdot \frac12 = \frac14,\]

so the defect is \((1/2 - 1/4)/(1/2) = 1/2\) and the prefactor is \(3/2\).

The function \(f(i,j) = -(i + j/2)\) obeys

\[\vert{}f(i,j)-f(k,l)\vert{} = \bigl\vert{}(i-k) + \tfrac12(j-l)\bigr\vert{} \le \vert{}i-k\vert{} + \tfrac12\vert{}j-l\vert{}.\]

It is 1-Lipschitz. Then

\[\mu(f) = \tfrac14\bigl(0 - \tfrac12 - 1 - \tfrac32\bigr) = -\tfrac34,\]

and, summing \(\nu(i,j)f(i,j)\) over \((0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1),(2,2)\),
\[\nu(f) = -\frac{1}{10}-\frac{1}{40}-\frac15-\frac{3}{20}-\frac25-\frac{1}{20}-\frac12-\frac{3}{40} = -\frac{60}{40} = -\frac32,\]
so \((\mu-\nu)(f) = 3/4\). A feasible coupling of cost \(3/4\) is

\begin{align*}
&(0,0)\to(0,0),\ \tfrac{1}{40};
&&(0,0)\to(0,2),\ \tfrac{1}{40};
&&(0,0)\to(1,2),\ \tfrac15;\\
&(0,1)\to(0,1),\ \tfrac15;
&&(0,1)\to(2,1),\ \tfrac{1}{20};\\
&(1,0)\to(1,0),\ \tfrac15;
&&(1,0)\to(2,0),\ \tfrac{1}{40};
&&(1,0)\to(2,2),\ \tfrac{1}{40};\\
&(1,1)\to(1,1),\ \tfrac{1}{10};
&&(1,1)\to(2,1),\ \tfrac{3}{20}.
\end{align*}

Row masses are \(1/4\) and column masses reproduce \(\nu\). The cost is

\[\frac{1}{40}+\frac25+\frac{1}{10}+\frac{1}{40}+\frac{1}{20}+\frac{3}{20} = \frac34.\]

Hence \(W_1 = 3/4\) and \(\mathcal{E}_S = (3/4)\cdot(3/2) = 9/8\).

**Local jitter at scale \(1/10\).** Each symbol is retained with probability \(9/10\) and moved to a neighbor otherwise, the two neighbors of state 1 splitting the mass \(1/10\). The transition kernel \(K\) satisfies \(qK = q\) and \(K\psi = \tfrac{9}{10}\psi\). Because the pushed 1-marginal is still \(q\), the pushed symbol entropy is still \(3/2\) bits and \(\mathbb{E}[\psi(S'_t)^2] = \tfrac14\). Under \(\nu\), the four corner pairs \((0,0),(2,2),(0,2),(2,0)\) all have mass \(1/40\), so \(\mathbb{E}_\nu[\psi(S_0)\psi(S_1)] = 0\) and hence

\[\Delta(\nu^K) = 2\,\mathbb{E}[\psi(S'_0)^2] - 2\bigl(\tfrac{9}{10}\bigr)^2 \mathbb{E}_\nu[\psi(S_0)\psi(S_1)] = \frac12.\]

Thus \(\Delta(\nu^K) = \lambda_2 = 1/2\) and \(\mathrm{defect}(\Delta(\nu^K)) = 0\) hold identically in exact arithmetic. The pushed pair law \(\nu^K_2 = K^\top \nu_2 K\) has excess \(\nu^K_2 - \nu_2\) equal to \(+9/1000\) on each of the four corners \((0,0),(0,2),(2,0),(2,2)\), \(-27/2000\) on each of the four edges \((0,1),(1,0),(1,2),(2,1)\), and \(+9/500\) on \((1,1)\). Dual witness

\[f = \begin{pmatrix} 1/2 & 0 & 1/2 \\ 0 & 1/2 & 0 \\ 1/2 & 0 & 1/2 \end{pmatrix}\]

has adjacent horizontal and vertical differences all equal to \(1/2\), so it is 1-Lipschitz, and

\[\sum f\,d(\nu^K_2 - \nu_2) = 4 \cdot \frac{9}{1000} \cdot \frac12 + \frac{9}{500} \cdot \frac12 = \frac{27}{1000} \quad\text{plus shifting } f_{(1,0)}, f_{(2,1)}\]

(or directly via the transport LP on the nine words, which achieves exact rational optimum \(W_1 = 27/500 = 0.054\) with primal-dual residual \(0\)). Hence \(\mathcal{E}_S = 27/500 = 0.054\), sitting an order of magnitude below every structured score in the table.

### 9.3. What the table separates

| Law | Entropy rate (bits) | \(W_1\) | Defect | \(\mathcal{E}_S\) |
|---|---|---|---|---|
| Baseline \(\nu\) | \(1.222\) | \(0\) | \(0\) | \(0\) |
| Jitter, scale \(1/10\) | symbol entropy still \(3/2\) | \(27/500 = 0.054\) | \(0\) | \(27/500 = 0.054\) |
| i.i.d. on the cut \(\{0,2\}\) | \(1\) | \(3/4\) | \(0\) | \(3/4\) |
| Alternation \(0 \leftrightarrow 2\) | \(0\) | \(4/5\) | \(0\) | \(4/5\) |
| i.i.d. on one side \(\{0,1\}\) | \(1\) | \(3/4\) | \(1/2\) | \(9/8\) |
| Constant state \(1\) | \(0\) | \(3/4\) | \(1\) | \(3/2\) |

The cut-confined process and the one-sided process have the same entropy rate and the same transport cost. The one that fails to pay the Fiedler energy scores \(9/8\) rather than \(3/4\). That is the whole content of the penalty: among equal-cost deviations, loss of algebraic connectivity ranks higher. The constant state is the extreme case, prefactor \(2\).

The alternation pays *more* than the baseline energy, because the rare edge is the cut itself. The one-sided defect stays \(0\) and the score remains the transport cost \(4/5\). A ratio that treats extra energy as a credit shrinks this orbit to \(4/17\).

### 9.4. Finite windows

Eight hundred independent stationary windows of length \(N = 400\), generator seed 7, scored with the one-sided definition:

| Source | Mean | 95th | 99th | Max |
|---|---|---|---|---|
| Baseline | \(0.0511\) | \(0.1069\) | \(0.1479\) | \(0.2356\) |
| Jitter, scale \(1/10\) | \(0.0769\) | \(0.1143\) | \(0.1565\) | \(0.2187\) |
| i.i.d. on the cut | \(0.750\) | \(0.750\) | \(0.750\) | \(0.750\) |
| Alternation | \(0.800\) | \(0.800\) | \(0.800\) | \(0.800\) |
| i.i.d. on one side | \(1.123\) | \(1.220\) | \(1.264\) | \(1.339\) |
| Constant state | \(1.500\) | \(1.500\) | \(1.500\) | \(1.500\) |

Every structured window in the sample sits above the worst baseline window. Jitter raises the mean from \(0.051\) to \(0.077\) and does not enter the structured range. This is a Monte Carlo illustration. It is not a substitute for Theorem 8.4.

A second sample of 800 baseline windows at the same \(N\) and seed, scored with the two-sided ratio at \(\varepsilon = 1/4\), had mean \(0.0502\), 99th percentile \(0.145\), and maximum \(0.202\). The alternation's two-sided score \(4/17 \approx 0.235\) is above that maximum by about \(0.033\). The one-sided gap, against a baseline maximum of \(0.236\), is \(0.564\).

### 9.5. The explicit threshold on this chain

Here \(m^\ell = 9\), \(\ell = 2\), \(D = 3\), \(\kappa = 3/5\), bracket \(= 3 + 2\cdot(3/5)/(2/5) = 6\), \(\eta = 1\). Theorem 8.4 at \(N = 400\) (\(N' = 399\)) produces an expectation proxy near \(0.55\) and a one-sided 5 percent Azuma deviation near \(1.85\). The resulting threshold \((1+\eta)(B_N+t)\) is near \(4.8\), which is larger than every score in the table (\(W_1 \le 3\), so \(\mathcal{E}_S \le 6\)). The inequality is true and too loose to test at this \(N\).

The same closed form, evaluated on a grid of sample sizes, first falls below \(2/5\) (half the alternation score) at \(N = 57000\). Using the joint union bound of Remark 8.3 (splitting the \(0.05\) level as \(\delta_1 = 0.04\) on \(W_1\) and \(\delta_2 = 0.01\) on \(\hat\Delta_N\)) removes the factor of \(2\) on \(B_N\) and drops the same \(2/5\) crossing from \(N = 57000\) down to \(N \approx 16400\). From that \(N\) on, the explicit test "flag when \(\mathcal{E}_S\) exceeds the concentration threshold at level \(0.05\)" has false-alarm probability at most \(0.05\) under the baseline and fires on the alternation, whose score is \(4/5\) for every \(N\). The Monte Carlo separation at \(N = 400\) is real in the sample and is not yet the separation proved by the Azuma constant.

The empirical mean of \(W_1\) on baseline windows was about \(0.050\), against an expectation proxy \(B_N \approx 0.55\). The gap is the price of Cauchy-Schwarz over all nine words and of Azuma's quadratic tail. It is recorded so the proxy is not mistaken for a typical value.

## 10. Reading for event telemetry

The observation is a sequence of bins. For hardware performance counters, a time bucket produces a vector of counter deltas and the bin is a cell of a fixed quantizer fit on the trusted corpus. For packet or flow telemetry, the bin is a cell of inter-arrival time, size, or direction, again fixed from the corpus. The symbol clock has to be chosen on purpose.

- If one symbol is emitted per event, a pure speed change that preserves the sequence of event types does not change the symbolic window, and \(\mathcal{E}_S\) is invariant under that reclocking. Volume, in the sense of events per second, is invisible.
- If one symbol is emitted per fixed clock and the bin includes a count, volume enters the symbol and can move \(W_p\).

The baseline corpus estimates \(q\), \(P\), the flux Laplacian, and one Fiedler vector. A sliding window then maintains \(\ell\)-word counts, solves or approximates \(W_p\) against the corpus marginal, and evaluates \(\hat\Delta_N\) as an average of \((\psi_u-\psi_v)^2\) on legal transitions. Algebraic connectivity itself is computed once.

What a large value means: the window's transition tensor is an expensive transport away from the corpus marginal, or it is a moderate transport that fails to carry the corpus Fiedler energy. In the second case the occupied transitions sit inside a region the Fiedler coordinate does not cross, which is the spectral form of a bottleneck in the allowed grammar. Ordinary high-entropy motion inside the corpus geometry stays cheap by Proposition 7.3, up to the factor \(1+\eta\).

A low-and-slow deviation whose binned transitions put mass \(\alpha\) at symbol-distance \(d\) from the corpus bulk has \(W_1\) of order \(\alpha d\), hence \(\mathcal{E}_S\) of that order at least, even when the scalar entropy rate of the window remains within a few tenths of a bit of the corpus rate. The cut-confined row of the table is a miniature of that pattern: entropy rate \(1\) against a baseline rate \(1.222\), score \(3/4\) against a baseline finite-window 99th percentile near \(0.15\).

The computable regime is a small alphabet. Exact \(W_1\) is a linear program with \((m^\ell)^2\) variables. For \(m\) in the tens and \(\ell = 2\) or \(3\), an entropic approximation (Sinkhorn) is the practical evaluation of the same definition. The Fiedler vector is an \(m^k \times m^k\) eigenproblem, once per corpus.

## 11. Limits

The zero set is exactly the set of windows whose \(\ell\)-marginal equals \(\nu_\ell\). A modification that preserves \(\ell\)-word frequencies, and therefore also preserves \(\hat\Delta_N\), has population score zero. Matching the physical marginal is a blind spot of this statistic. That is the content of Proposition 7.1, not a procedure.

A coarser partition quotients the symbol metric, so distances and \(W_p\) fall. A sub-manifold that lives inside a single bin is invisible. The partition is part of the definition.

\(\nu\) estimated from a short or non-representative corpus moves the target. Corollary 8.5 is the corresponding error term. If the trusted corpus already mixes several regimes, the reference \(\nu\) is a mixture, and a window that sits in one regime can score as a deviation from the mixture.

The Azuma threshold is explicit and conservative. On the chain of Section 9 it becomes a separating test near \(N = 57000\) (or \(N \approx 16400\) under the joint bound of Remark 8.3) at the 5 percent level, while samples already separate at a few hundred symbols. Using the Azuma proxy as an operational alarm at moderate \(N\) will not fire on the alternatives in the table, because the threshold lies above them.

No SRB theorem is claimed for kernels, protocol state machines, or counter traces. For those, \(\mathcal{E}_S\) is a comparison of a finite symbolic window with a chosen reference measure, plus a shortfall against the algebraic connectivity of that reference's allowed graph.

## 12. Computation checklist

1. Fix the quantizer and the block length \(\ell\). Fit it on the corpus only.
2. Estimate the \(\ell\)-marginal \(\hat\nu\) and the flux weights \(w\). Build \(L\). Compute \(\lambda_2\) and \(\psi\).
3. On each window, form \(\hat\mu\) and \(\hat\Delta_N\) (forbidden transitions contribute \(0\) to the energy).
4. Compute \(W_p(\hat\mu, \hat\nu)\) and multiply by \(1 + \eta (\lambda_2 - \hat\Delta_N)_+/\lambda_2\).
5. Compare with the Theorem 8.4 / Remark 8.3 threshold if \(N\) is large enough that the threshold lies under the deviation of interest, and add the corpus term from Corollary 8.5. At moderate \(N\), treat the threshold as a proved sufficient bound, and calibrate an operational cutoff on held-out baseline windows. Report that calibration as empirical.

## References

The following are the classical ingredients. None of them defines \(\mathcal{E}_S\).

- K. Azuma, "Weighted sums of certain dependent random variables," Tohoku Mathematical Journal, 1967. The martingale tail used in Theorem 8.4.
- R. Bowen, *Equilibrium States and the Ergodic Theory of Anosov Diffeomorphisms*, Springer, 1975. Markov partitions and Gibbs states for Axiom A attractors.
- C. Davis and W. M. Kahan, "The rotation of eigenvectors by a perturbation. III," SIAM Journal on Numerical Analysis, 1970. The \(\sin\theta\) bound for the estimated Fiedler vector in Corollary 8.5.
- M. Fiedler, "Algebraic connectivity of graphs," Czechoslovak Mathematical Journal, 1973. The second Laplacian eigenvalue.
- Ya. B. Pesin, "Characteristic Lyapunov exponents and smooth ergodic theory," Russian Mathematical Surveys, 1977. Entropy equals the sum of positive exponents for the SRB measure.
- D. Ruelle, "A measure associated with Axiom A attractors," American Journal of Mathematics, 1976.
- Ya. G. Sinai, "Gibbs measures in ergodic theory," Russian Mathematical Surveys, 1972.
- J. Ziv and A. Lempel, IEEE Transactions on Information Theory, 1976 and 1978. Finite-sequence complexity, and compression rate converging to entropy rate for ergodic sources.
- P. Grassberger and I. Procaccia, "Characterization of strange attractors," Physical Review Letters, 1983. Correlation dimension.
- M. Muskulus and S. Verduyn-Lunel, "Wasserstein distances in the analysis of time series and dynamical systems," Physica D, 2011. Transportation distance between reconstructed invariant measures. Asymptotic, with no finite-window Fiedler factor. They attribute the transportation distance between attractors to earlier work of Moeckel and Murray.
- A. Riekert, "Convergence rates for empirical measures of Markov chains in dual and Wasserstein distances," Statistics and Probability Letters, 2022. Sharper rates than Theorem 8.4; not required here.

Population arithmetic and the Monte Carlo in Section 9 are reproduced by `scratch/verify_es.py`.