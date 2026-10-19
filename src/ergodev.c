#include "ergodev.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

static uint32_t u32_pow(uint32_t base, uint32_t exp) {
    uint32_t res = 1U;
    while (exp > 0U) {
        if (exp & 1U) res *= base;
        base *= base;
        exp >>= 1U;
    }
    return res;
}

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
) {
    if (!ctx || m < 2U || m > ERGODEV_MAX_ALPHABET || ell < 2U) return -1;
    uint32_t num_words = u32_pow(m, ell);
    if (num_words > ERGODEV_MAX_WORDS) return -2;
    if (window_len < ell || window_len > ERGODEV_MAX_WINDOW) return -3;
    if (lambda_2 <= 0.0 || eta < 0.0 || sinkhorn_reg <= 0.0) return -4;

    memset(ctx, 0, sizeof(*ctx));
    ctx->m           = m;
    ctx->ell         = ell;
    ctx->num_words   = num_words;
    ctx->window_len  = window_len;
    ctx->num_windows = window_len - ell + 1U;
    ctx->high_pow    = u32_pow(m, ell - 1U);
    ctx->lambda_2    = lambda_2;
    ctx->eta         = eta;

    /* Precompute 1-cycle state transition LUT: eliminates integer div/mod in hot path */
    for (uint32_t w = 0; w < num_words; ++w) {
        uint32_t suffix = (w % ctx->high_pow) * m;
        for (uint32_t s = 0; s < m; ++s) {
            ctx->step_lut[w][s] = (uint16_t)(suffix + s);
        }
    }

    for (uint32_t i = 0; i < num_words; ++i) {
        ctx->nu[i]        = nu_ell[i];
        ctx->jump_cost[i] = jump_cost[i];
        for (uint32_t j = 0; j < num_words; ++j) {
            double d = word_dist_matrix[i * num_words + j];
            ctx->dist_mat[i][j] = d;
            ctx->gibbs_k[i][j]  = exp(-d / sinkhorn_reg);
        }
    }
    return 0;
}

void ergodev_reset(ergodev_ctx_t *ctx) {
    memset(ctx->counts, 0, sizeof(ctx->counts));
    ctx->sym_filled   = 0U;
    ctx->word_head    = 0U;
    ctx->rolling_word = 0U;
    ctx->energy_sum   = 0.0;
}

int ergodev_push(ergodev_ctx_t *ctx, uint8_t symbol) {
    if (symbol >= ctx->m) symbol %= (uint8_t)ctx->m;
    return ergodev_push_fast(ctx, symbol);
}

double ergodev_energy(const ergodev_ctx_t *ctx) {
    if (ctx->sym_filled < ctx->window_len) return 0.0;
    double sum = 0.0;
    for (uint32_t w = 0; w < ctx->num_words; ++w) {
        sum += (double)ctx->counts[w] * ctx->jump_cost[w];
    }
    return sum / (double)ctx->num_windows;
}

double ergodev_defect(const ergodev_ctx_t *ctx) {
    double delta_hat = ergodev_energy(ctx);
    double diff = ctx->lambda_2 - delta_hat;
    return (diff > 0.0) ? (diff / ctx->lambda_2) : 0.0;
}

double ergodev_score_sinkhorn(
    const ergodev_ctx_t *ctx,
    uint32_t max_iters,
    double *out_w1,
    double *out_defect
) {
    const uint32_t W = ctx->num_words;
    const double inv_n = 1.0 / (double)ctx->num_windows;

    double mu[ERGODEV_MAX_WORDS];
    double u[ERGODEV_MAX_WORDS];
    double v[ERGODEV_MAX_WORDS];

    /* Check for exact match (W_1 = 0) and build empirical marginal */
    double l1_diff = 0.0;
    for (uint32_t i = 0; i < W; ++i) {
        mu[i] = (double)ctx->counts[i] * inv_n;
        l1_diff += fabs(mu[i] - ctx->nu[i]);
        u[i] = 1.0;
        v[i] = 1.0;
    }

    double defect = ergodev_defect(ctx);
    if (l1_diff < 1e-14) {
        if (out_w1) *out_w1 = 0.0;
        if (out_defect) *out_defect = defect;
        return 0.0;
    }

    /* Stack-only Sinkhorn-Knopp matrix scaling */
    for (uint32_t iter = 0; iter < max_iters; ++iter) {
        for (uint32_t i = 0; i < W; ++i) {
            if (mu[i] <= 1e-15) {
                u[i] = 0.0;
                continue;
            }
            double kv = 0.0;
            for (uint32_t j = 0; j < W; ++j) {
                kv += ctx->gibbs_k[i][j] * v[j];
            }
            u[i] = (kv > 1e-300) ? (mu[i] / kv) : 0.0;
        }
        for (uint32_t j = 0; j < W; ++j) {
            if (ctx->nu[j] <= 1e-15) {
                v[j] = 0.0;
                continue;
            }
            double ktu = 0.0;
            for (uint32_t i = 0; i < W; ++i) {
                ktu += ctx->gibbs_k[i][j] * u[i];
            }
            v[j] = (ktu > 1e-300) ? (ctx->nu[j] / ktu) : 0.0;
        }
    }

    /* Evaluate primal transport cost <P, D> */
    double w1 = 0.0;
    for (uint32_t i = 0; i < W; ++i) {
        if (u[i] == 0.0) continue;
        for (uint32_t j = 0; j < W; ++j) {
            double p_ij = u[i] * ctx->gibbs_k[i][j] * v[j];
            w1 += p_ij * ctx->dist_mat[i][j];
        }
    }

    double es = w1 * (1.0 + ctx->eta * defect);
    if (out_w1) *out_w1 = w1;
    if (out_defect) *out_defect = defect;
    return es;
}

#ifndef ERGODEV_NO_MAIN
/* =========================================================================
 * Built-in Self-Test: Reproduces Section 9 Worked Chain in O(1) Streaming
 * ========================================================================= */
int main(void) {
    /* Section 9.1: 3-state baseline chain, ell = 2 (9 words)
     * Word index w = 3*i + j for i,j in {0,1,2}. */
    const double nu_2[9] = {
        1.0/40.0, 1.0/5.0,  1.0/40.0, /* (0,0), (0,1), (0,2) */
        1.0/5.0,  1.0/10.0, 1.0/5.0,  /* (1,0), (1,1), (1,2) */
        1.0/40.0, 1.0/5.0,  1.0/40.0  /* (2,0), (2,1), (2,2) */
    };

    /* Fiedler vector psi = (-1/sqrt(2), 0, 1/sqrt(2)) -> jumps (psi_i - psi_j)^2 */
    const double jump_cost[9] = {
        0.0, 0.5, 2.0,
        0.5, 0.0, 0.5,
        2.0, 0.5, 0.0
    };

    /* Product metric d((i,j), (k,l)) = |i - k| + 0.5 * |j - l| */
    double dist_mat[81];
    for (int w1 = 0; w1 < 9; ++w1) {
        int i = w1 / 3, j = w1 % 3;
        for (int w2 = 0; w2 < 9; ++w2) {
            int k = w2 / 3, l = w2 % 3;
            dist_mat[w1 * 9 + w2] = fabs((double)(i - k)) + 0.5 * fabs((double)(j - l));
        }
    }

    ergodev_ctx_t ctx;
    /* Window N = 401 -> N' = 400 pair transitions; Sinkhorn reg = 0.005 */
    ergodev_init(&ctx, 3, 2, 401, 0.5, 1.0, nu_2, jump_cost, dist_mat, 0.005);

    printf("%-28s | %-10s | %-10s | %-10s\n", "Stream Source (N'=400)", "W_1", "Defect", "E_S Score");
    printf("-------------------------------------------------------------------\n");

    /* 1. Alternation 0, 2, 0, 2, ... (Expected: W_1 = 0.800, Defect = 0.000, E_S = 0.800) */
    ergodev_reset(&ctx);
    for (uint32_t t = 0; t < 401; ++t) {
        ergodev_push(&ctx, (t & 1U) ? 2U : 0U);
    }
    double w1, def, es;
    es = ergodev_score_sinkhorn(&ctx, 200, &w1, &def);
    printf("%-28s | %-10.4f | %-10.4f | %-10.4f\n", "Alternation 0 <-> 2", w1, def, es);

    /* 2. Constant state 1, 1, 1, ... (Expected: W_1 = 0.750, Defect = 1.000, E_S = 1.500) */
    ergodev_reset(&ctx);
    for (uint32_t t = 0; t < 401; ++t) {
        ergodev_push(&ctx, 1U);
    }
    es = ergodev_score_sinkhorn(&ctx, 200, &w1, &def);
    printf("%-28s | %-10.4f | %-10.4f | %-10.4f\n", "Constant State 1", w1, def, es);

    /* 3. Exact 4-cycle on {0,2}: 0,0,2,2,0,0,2,2... (Expected: W_1 = 0.750, Defect = 0, E_S = 0.750) */
    ergodev_reset(&ctx);
    const uint8_t cut_cycle[4] = {0, 0, 2, 2};
    for (uint32_t t = 0; t < 401; ++t) {
        ergodev_push(&ctx, cut_cycle[t & 3U]);
    }
    es = ergodev_score_sinkhorn(&ctx, 200, &w1, &def);
    printf("%-28s | %-10.4f | %-10.4f | %-10.4f\n", "Balanced Cut {0,2}", w1, def, es);

    /* 4. Exact 4-cycle on one side {0,1}: 0,0,1,1... (Expected: W_1 = 0.750, Defect = 0.5, E_S = 1.125) */
    ergodev_reset(&ctx);
    const uint8_t side_cycle[4] = {0, 0, 1, 1};
    for (uint32_t t = 0; t < 401; ++t) {
        ergodev_push(&ctx, side_cycle[t & 3U]);
    }
    es = ergodev_score_sinkhorn(&ctx, 200, &w1, &def);
    printf("%-28s | %-10.4f | %-10.4f | %-10.4f\n", "One-Sided Submanifold {0,1}", w1, def, es);

    return 0;
}
#endif /* ERGODEV_NO_MAIN */