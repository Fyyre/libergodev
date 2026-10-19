#ifndef ERGODEV_H
#define ERGODEV_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Compile-time upper bounds for zero-allocation static workspaces.
 * Adjust ERGODEV_MAX_WORDS (m^ell) and ERGODEV_MAX_WINDOW (N) as needed. */
#define ERGODEV_MAX_ALPHABET 16U
#define ERGODEV_MAX_WORDS    64U    /* e.g., m=3, ell=2 -> 9; m=4, ell=3 -> 64 */
#define ERGODEV_MAX_WINDOW   4096U

typedef struct {
    /* Configuration parameters */
    uint32_t m;              /* Alphabet size |A| */
    uint32_t ell;            /* Block length ell >= 2 */
    uint32_t num_words;      /* m^ell */
    uint32_t window_len;     /* Window length N (N' = N - ell + 1 words) */
    uint32_t num_windows;    /* N' = N - ell + 1 */
    uint32_t high_pow;       /* m^(ell - 1), used to strip oldest symbol */

    double   lambda_2;       /* Baseline Rauzy Laplacian algebraic connectivity */
    double   eta;            /* Fiedler defect penalty weight (eta >= 0) */

    /* Precomputed baseline tables (size: num_words) */
    double   nu[ERGODEV_MAX_WORDS];         /* Reference ell-word marginal nu_ell */
    double   jump_cost[ERGODEV_MAX_WORDS];  /* c(w) = (psi_u - psi_v)^2 if legal, 0 if forbidden */

    /* Precomputed ground metric & Gibbs kernel on ell-words (num_words x num_words) */
    double   dist_mat[ERGODEV_MAX_WORDS][ERGODEV_MAX_WORDS];
    double   gibbs_k[ERGODEV_MAX_WORDS][ERGODEV_MAX_WORDS]; /* exp(-d(u,v) / reg) */

    /* O(1) Streaming Ring-Buffer State */
    uint8_t  sym_ring[ERGODEV_MAX_WINDOW];  /* Raw symbol ring buffer of length N */
    uint16_t word_ring[ERGODEV_MAX_WINDOW]; /* Ring buffer of N' active ell-word indices */
    uint32_t counts[ERGODEV_MAX_WORDS];     /* Empirical word counts in current window */

    uint32_t sym_filled;     /* Number of symbols ingested so far (up to N) */
    uint32_t word_head;      /* Write index in word_ring [0 .. N'-1] */
    uint32_t rolling_word;   /* Current base-m rolling word state */
    double   energy_sum;     /* Rolling sum of jump_cost[w] over active N' words */
} ergodev_ctx_t;

/* Initialize context with baseline marginal nu_ell, Fiedler jump costs, and word metric.
 * Returns 0 on success, negative on invalid dimensions. */
int ergodev_init(
    ergodev_ctx_t *ctx,
    uint32_t m,
    uint32_t ell,
    uint32_t window_len,
    double lambda_2,
    double eta,
    const double *nu_ell,
    const double *jump_cost,
    const double *word_dist_matrix,
    double sinkhorn_reg
);

/* Reset streaming window state without clearing precomputed baseline tables. */
void ergodev_reset(ergodev_ctx_t *ctx);

/* O(1) per-symbol ingestion. Updates rolling ell-word histogram and Dirichlet energy.
 * Returns 1 once the window of length N is full (ready to score), 0 during warm-up. */
int ergodev_push(ergodev_ctx_t *ctx, uint8_t symbol);

/* Compute current empirical Dirichlet energy \hat\Delta_N in O(1). */
double ergodev_energy(const ergodev_ctx_t *ctx);

/* Compute current one-sided Fiedler defect ((\lambda_2 - \hat\Delta_N)_+ / \lambda_2) in O(1). */
double ergodev_defect(const ergodev_ctx_t *ctx);

/* Evaluate W_1(\hat\mu_{N,\ell}, \nu_\ell) and full Symbolic Ergodeviation \mathcal{E}_S.
 * Uses stack-only Sinkhorn iterations with zero heap allocation. */
double ergodev_score_sinkhorn(
    const ergodev_ctx_t *ctx,
    uint32_t max_iters,
    double *out_w1,
    double *out_defect
);

#ifdef __cplusplus
}
#endif

#endif /* ERGODEV_H */