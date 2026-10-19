#ifndef ERGODEV_H
#define ERGODEV_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define ERGODEV_MAX_ALPHABET 16U
#define ERGODEV_MAX_WORDS    64U
#define ERGODEV_MAX_WINDOW   16384U

typedef struct {
    uint32_t m;
    uint32_t ell;
    uint32_t num_words;
    uint32_t window_len;
    uint32_t num_windows;
    uint32_t high_pow;

    double   lambda_2;
    double   eta;

    double   nu[ERGODEV_MAX_WORDS];
    double   jump_cost[ERGODEV_MAX_WORDS];
    double   dist_mat[ERGODEV_MAX_WORDS][ERGODEV_MAX_WORDS];
    double   gibbs_k[ERGODEV_MAX_WORDS][ERGODEV_MAX_WORDS];

    /* 1-cycle state transition LUT: step_lut[w][s] = ((w % m^(ell-1)) * m) + s
     * Eliminates all integer division/modulo instructions from the hot loop! */
    uint16_t step_lut[ERGODEV_MAX_WORDS][ERGODEV_MAX_ALPHABET];

    uint16_t word_ring[ERGODEV_MAX_WINDOW];
    uint32_t counts[ERGODEV_MAX_WORDS];

    uint32_t sym_filled;
    uint32_t word_head;
    uint32_t rolling_word;
    double   energy_sum;
} ergodev_ctx_t;

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

void ergodev_reset(ergodev_ctx_t *ctx);
int  ergodev_push(ergodev_ctx_t *ctx, uint8_t symbol);

/* Hot-path inline version: 100% integer operations (zero FPU instructions).
 * Safe for kernel-mode / no-FPU contexts and eliminates FP accumulator latency. */
static inline int ergodev_push_fast(ergodev_ctx_t *ctx, uint8_t symbol) {
    uint16_t w_in = ctx->step_lut[ctx->rolling_word][symbol];
    ctx->rolling_word = w_in;

    if (__builtin_expect(ctx->sym_filled >= ctx->window_len, 1)) {
        uint32_t head = ctx->word_head;
        uint16_t w_out = ctx->word_ring[head];
        ctx->word_ring[head] = w_in;
        ctx->counts[w_out]--;
        ctx->counts[w_in]++;
        head++;
        ctx->word_head = (head == ctx->num_windows) ? 0U : head;
        return 1;
    }

    ctx->sym_filled++;
    if (ctx->sym_filled >= ctx->ell) {
        uint32_t head = ctx->word_head;
        ctx->word_ring[head] = w_in;
        ctx->counts[w_in]++;
        head++;
        ctx->word_head = (head == ctx->num_windows) ? 0U : head;
    }
    return (ctx->sym_filled >= ctx->window_len) ? 1 : 0;
}

double ergodev_energy(const ergodev_ctx_t *ctx);
double ergodev_defect(const ergodev_ctx_t *ctx);

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