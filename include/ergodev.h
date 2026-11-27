#ifndef ERGODEV_H
#define ERGODEV_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Symbolic ergodeviation, streaming half.
 *
 * Supported placement (one thread, two clocks):
 *   - per_event: one symbol per syscall. The score is invariant under a
 *     pure change of rate, because the symbol string does not change.
 *   - per_quantum: one symbol per fixed clock. Symbol 0 is idle (no
 *     syscall in the bucket); the other symbols are the class alphabet.
 *     Volume enters this clock.
 * Contexts are keyed by tid. ergodev_thread_push_* only ingests. Scoring
 * uses caller-owned scratch and runs off the callback that observed the
 * event. Dense push_fast is integer-only.
 *
 * The hot context holds the ring, the counts, the jump table, and (dense)
 * the step LUT. The distance matrix and Sinkhorn workspace live in
 * ergodev_scratch_t, supplied by the caller.
 *
 * Operating threshold: an empirical quantile of held-out baseline scores
 * (ergodev_empirical_quantile). ergodev_bounds returns the Azuma sufficient
 * condition of Theorem 8.4 and a Freedman threshold whose variance proxy
 * is B_N^2. Those bounds are valid for every finite N. They are an alarm
 * only once they sit under the alternative of interest.
 */

#define ERGODEV_MAX_ALPHABET 64U
#define ERGODEV_MAX_ELL      8U
#define ERGODEV_MAX_WORDS    128U
#define ERGODEV_MAX_WINDOW   16384U
#define ERGODEV_HASH_SIZE    256U

typedef struct ergodev_ctx {
    uint32_t m;
    uint32_t ell;
    uint32_t num_words;
    uint32_t window_len;
    uint32_t num_windows;
    uint32_t high_pow;
    uint32_t sym_filled;
    uint32_t active_words;
    uint32_t word_head;
    uint32_t rolling_word;
    uint32_t sparse;
    uint32_t n_vertices;
    double   lambda_2;
    double   lambda_next;
    double   spectral_gap;
    double   eta;
    double   max_jump;
    double   diameter;
    double   kappa;
    int      q_normalized;
    uint32_t off_nu;
    uint32_t off_jump;
    uint32_t off_lut;
    uint32_t off_ring;
    uint32_t off_counts;
    uint32_t off_hash_key;
    uint32_t off_hash_val;
    uint32_t off_hash_used;
    uint32_t off_word_id;
    uint8_t  hist[ERGODEV_MAX_ELL];
    uint8_t  hist_len;
    /* Hot-path bases. Set once in init/fit so push does not recompute offsets. */
    uint16_t *lut;
    uint16_t *ring;
    uint32_t *counts;
} ergodev_ctx_t;

typedef struct {
    uint32_t num_words;
    double   reg;
    double  *dist;
    double  *work;
} ergodev_scratch_t;

typedef struct {
    double w1;
    double defect;
    double es;
    double residual;
    int    converged;
} ergodev_score_t;

typedef struct {
    double B_n;
    double azuma_t;
    double azuma_threshold;
    double freedman_t;
    double freedman_threshold;
    double sigma2;
    double joint_threshold;
    int    azuma_separates;
    int    freedman_separates;
    int    joint_separates;
} ergodev_bound_t;

typedef struct {
    int      alarm;
    int      alarm_transport;
    int      alarm_period;
    double   es;
    uint32_t period;
} ergodev_alarm_t;

typedef struct {
    int      sparse;
    int      error;
    uint32_t num_words;
    uint32_t num_windows;
    size_t   ctx_bytes;
} ergodev_plan_t;

typedef struct {
    uint32_t      tid;
    ergodev_ctx_t *per_event;
    ergodev_ctx_t *per_quantum;
} ergodev_thread_slot_t;

size_t ergodev_ctx_bytes(uint32_t num_words, uint32_t m, uint32_t num_windows, int sparse);
size_t ergodev_scratch_doubles(uint32_t num_words);

int ergodev_plan(
    uint32_t m,
    uint32_t ell,
    uint32_t window_len,
    const uint8_t *corpus,
    uint32_t corpus_len,
    ergodev_plan_t *out
);

/* Manual baseline. word_dist is read to obtain the diameter and is not stored.
 * Returns 0, or a negative code. Rejects a jump table whose nu-mean is not lambda_2. */
int ergodev_init(
    ergodev_ctx_t *ctx,
    size_t ctx_bytes,
    uint32_t m,
    uint32_t ell,
    uint32_t window_len,
    double lambda_2,
    double eta,
    const double *nu_ell,
    const double *jump_cost,
    const double *word_dist
);

/* Corpus fit. dist_out receives the word metric (num_words^2 entries).
 * symbol_metric NULL means delta(i,j) = |i-j|.
 * q_normalize selects the stationary-mass Fiedler vector (generalized eigenproblem).
 * psi_out / q_out may be NULL. On a sparse alphabet, ell = 3 is legal whenever
 * the corpus support has at most ERGODEV_MAX_WORDS - 1 distinct ell-words. */
int ergodev_fit(
    ergodev_ctx_t *ctx,
    size_t ctx_bytes,
    uint32_t m,
    uint32_t ell,
    uint32_t window_len,
    double eta,
    double theta,
    const double *symbol_metric,
    const uint8_t *corpus,
    uint32_t corpus_len,
    int q_normalize,
    double *dist_out,
    uint32_t dist_cap_words,
    double *psi_out,
    uint32_t psi_cap,
    double *q_out,
    uint32_t q_cap
);

int ergodev_scratch_bind(
    ergodev_scratch_t *scratch,
    uint32_t num_words,
    double *storage,
    size_t n_doubles,
    const double *word_dist,
    double sinkhorn_reg
);

void ergodev_reset(ergodev_ctx_t *ctx);

static inline double *ergodev_nu(ergodev_ctx_t *ctx) {
    return (double *)((unsigned char *)ctx + ctx->off_nu);
}
static inline const double *ergodev_nu_c(const ergodev_ctx_t *ctx) {
    return (const double *)((const unsigned char *)ctx + ctx->off_nu);
}
static inline double *ergodev_jump(ergodev_ctx_t *ctx) {
    return (double *)((unsigned char *)ctx + ctx->off_jump);
}
static inline const double *ergodev_jump_c(const ergodev_ctx_t *ctx) {
    return (const double *)((const unsigned char *)ctx + ctx->off_jump);
}
static inline uint16_t *ergodev_lut(ergodev_ctx_t *ctx) {
    return (uint16_t *)((unsigned char *)ctx + ctx->off_lut);
}
static inline uint16_t *ergodev_ring(ergodev_ctx_t *ctx) {
    return (uint16_t *)((unsigned char *)ctx + ctx->off_ring);
}
static inline const uint16_t *ergodev_ring_c(const ergodev_ctx_t *ctx) {
    return (const uint16_t *)((const unsigned char *)ctx + ctx->off_ring);
}
static inline uint32_t *ergodev_counts(ergodev_ctx_t *ctx) {
    return (uint32_t *)((unsigned char *)ctx + ctx->off_counts);
}
static inline const uint32_t *ergodev_counts_c(const ergodev_ctx_t *ctx) {
    return (const uint32_t *)((const unsigned char *)ctx + ctx->off_counts);
}

int ergodev_push_sparse(ergodev_ctx_t *ctx, uint8_t symbol);

/* 1 window full, 0 still warming, -1 rejected (state unchanged). */
static inline int ergodev_push_fast(ergodev_ctx_t *ctx, uint8_t symbol) {
    if (__builtin_expect(!ctx || (uint32_t)symbol >= ctx->m, 0)) {
        return -1;
    }
    if (__builtin_expect(ctx->sparse, 0)) {
        return ergodev_push_sparse(ctx, symbol);
    }
    if (__builtin_expect((uint32_t)symbol >= ERGODEV_MAX_ALPHABET, 0)) {
        return -1;
    }

    if (__builtin_expect(ctx->sym_filled < ctx->ell - 1U, 0)) {
        ctx->rolling_word = ctx->rolling_word * ctx->m + (uint32_t)symbol;
        ctx->sym_filled++;
        return 0;
    }

    uint16_t w_in;
    uint16_t *ring;
    uint32_t *counts;
    uint32_t head;
    if (__builtin_expect(ctx->sym_filled == ctx->ell - 1U, 0)) {
        w_in = (uint16_t)(ctx->rolling_word * ctx->m + (uint32_t)symbol);
    } else {
        w_in = ctx->lut[ctx->rolling_word * ctx->m + symbol];
    }
    ctx->rolling_word = w_in;
    ring = ctx->ring;
    counts = ctx->counts;

    if (__builtin_expect(ctx->sym_filled >= ctx->window_len, 1)) {
        head = ctx->word_head;
        {
            uint16_t w_out = ring[head];
            ring[head] = w_in;
            counts[w_out]--;
            counts[w_in]++;
        }
        head++;
        ctx->word_head = (head == ctx->num_windows) ? 0U : head;
        return 1;
    }

    ctx->sym_filled++;
    head = ctx->word_head;
    ring[head] = w_in;
    counts[w_in]++;
    ctx->active_words++;
    head++;
    ctx->word_head = (head == ctx->num_windows) ? 0U : head;
    return (ctx->sym_filled >= ctx->window_len) ? 1 : 0;
}

int    ergodev_push(ergodev_ctx_t *ctx, uint8_t symbol);
double ergodev_energy(const ergodev_ctx_t *ctx);
double ergodev_defect(const ergodev_ctx_t *ctx);
double ergodev_l1_drift(const ergodev_ctx_t *ctx);

/* Lemma 8.1: E_S <= (D/2) * ||mu - nu||_1 * (1 + eta * defect). */
double ergodev_es_upper(const ergodev_ctx_t *ctx);
/* 1 when the upper bound is strictly below level, so the solver can be skipped. */
int    ergodev_screen(const ergodev_ctx_t *ctx, double level);

/* 0 if the plan meets residual_tol (a coupling, cost usable as W_1).
 * -2 if max_iters is exhausted first: out->w1 and out->es are NaN.
 *    A non-coupling is not reported as a transport cost. */
int ergodev_score(
    const ergodev_ctx_t *ctx,
    ergodev_scratch_t *scratch,
    uint32_t max_iters,
    double residual_tol,
    ergodev_score_t *out
);

/* Theorem 8.4 (Azuma), Freedman with the proved proxy Var(W_1) <= B_N^2,
 * and the Remark 8.3 joint product at levels 0.8*delta and 0.2*delta.
 * max_jump replaces the manuscript's jump cap of 2. alternative <= 0 skips
 * the separates flags. */
int ergodev_bounds(
    double diameter,
    uint32_t atoms,
    uint32_t ell,
    uint32_t N,
    double kappa,
    double eta,
    double lambda_2,
    double max_jump,
    double delta,
    double alternative,
    ergodev_bound_t *out
);

/* Destructive: sorts `scores`. Quantile index is ceil(p * n) - 1. */
double ergodev_empirical_quantile(double *scores, uint32_t n, double p);

int      ergodev_window_symbols(const ergodev_ctx_t *ctx, uint8_t *out, uint32_t cap);
uint32_t ergodev_min_period(const uint8_t *symbols, uint32_t n);

/* Union rule: alarm if E_S > tau_es or the symbol window has a period <= tau_period.
 * Split the false-alarm budget across the two arms when choosing the taus. */
int ergodev_alarm_union(
    double es,
    int es_valid,
    double tau_es,
    uint32_t period,
    uint32_t n_symbols,
    uint32_t tau_period,
    ergodev_alarm_t *out
);

int ergodev_thread_push_event(ergodev_thread_slot_t *slot, uint8_t syscall_class);
int ergodev_thread_push_quantum(ergodev_thread_slot_t *slot, uint8_t occupancy_symbol);

#ifdef __cplusplus
}
#endif

#endif /* ERGODEV_H */
