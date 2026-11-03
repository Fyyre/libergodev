#include "ergodev.h"
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/* Overflow-safe power computation: returns 0 on success, -1 on overflow > limit */
static int u32_pow_checked(uint32_t base, uint32_t exp, uint32_t limit, uint32_t *out) {
    uint64_t res = 1ULL;
    for (uint32_t i = 0; i < exp; ++i) {
        res *= (uint64_t)base;
        if (res > (uint64_t)limit) {
            return -1;
        }
    }
    *out = (uint32_t)res;
    return 0;
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
    if (!ctx || !nu_ell || !jump_cost || !word_dist_matrix) return -1;
    if (m < 2U || m > ERGODEV_MAX_ALPHABET || ell < 2U) return -2;

    uint32_t num_words = 0U;
    if (u32_pow_checked(m, ell, ERGODEV_MAX_WORDS, &num_words) != 0 || num_words == 0U) {
        return -2;
    }

    uint32_t high_pow = 0U;
    if (u32_pow_checked(m, ell - 1U, ERGODEV_MAX_WORDS, &high_pow) != 0 || high_pow == 0U) {
        return -2;
    }

    if (window_len < ell || window_len > ERGODEV_MAX_WINDOW) return -3;
    if (!isfinite(lambda_2) || lambda_2 <= 0.0) return -4;
    if (!isfinite(eta) || eta < 0.0) return -4;
    if (!isfinite(sinkhorn_reg) || sinkhorn_reg <= 0.0) return -4;

    double nu_sum = 0.0;
    for (uint32_t i = 0; i < num_words; ++i) {
        if (!isfinite(nu_ell[i]) || nu_ell[i] < 0.0) return -5;
        if (!isfinite(jump_cost[i]) || jump_cost[i] < 0.0) return -5;
        nu_sum += nu_ell[i];
        for (uint32_t j = 0; j < num_words; ++j) {
            double d = word_dist_matrix[i * num_words + j];
            if (!isfinite(d) || d < 0.0) return -5;
        }
    }
    if (fabs(nu_sum - 1.0) > 1e-6) return -5;

    memset(ctx, 0, sizeof(*ctx));
    ctx->m            = m;
    ctx->ell          = ell;
    ctx->num_words    = num_words;
    ctx->window_len   = window_len;
    ctx->num_windows  = window_len - ell + 1U;
    ctx->high_pow     = high_pow;
    ctx->lambda_2     = lambda_2;
    ctx->eta          = eta;
    ctx->sinkhorn_reg = sinkhorn_reg;

    for (uint32_t w = 0; w < num_words; ++w) {
        uint32_t suffix = (w % high_pow) * m;
        for (uint32_t s = 0; s < m; ++s) {
            ctx->step_lut[w][s] = (uint16_t)(suffix + s);
        }
    }

    int underflow_risk = 0;
    for (uint32_t i = 0; i < num_words; ++i) {
        ctx->nu[i]        = nu_ell[i];
        ctx->jump_cost[i] = jump_cost[i];
        for (uint32_t j = 0; j < num_words; ++j) {
            double d = word_dist_matrix[i * num_words + j];
            ctx->dist_mat[i][j] = d;
            double ratio = -d / sinkhorn_reg;
            if (ratio < -100.0) {
                underflow_risk = 1;
            }
            ctx->gibbs_k[i][j] = exp(ratio);
        }
    }
    ctx->use_log_sinkhorn = underflow_risk;
    return 0;
}

void ergodev_reset(ergodev_ctx_t *ctx) {
    if (!ctx) return;
    memset(ctx->counts, 0, sizeof(ctx->counts));
    ctx->sym_filled   = 0U;
    ctx->active_words = 0U;
    ctx->word_head    = 0U;
    ctx->rolling_word = 0U;
}

int ergodev_push(ergodev_ctx_t *ctx, uint8_t symbol) {
    return ergodev_push_fast(ctx, symbol);
}

double ergodev_energy(const ergodev_ctx_t *ctx) {
    if (!ctx || ctx->active_words == 0U) return 0.0;
    double sum = 0.0;
    for (uint32_t w = 0; w < ctx->num_words; ++w) {
        sum += (double)ctx->counts[w] * ctx->jump_cost[w];
    }
    return sum / (double)ctx->active_words;
}

double ergodev_defect(const ergodev_ctx_t *ctx) {
    if (!ctx || ctx->active_words == 0U) return 0.0;
    double delta_hat = ergodev_energy(ctx);
    double diff = ctx->lambda_2 - delta_hat;
    return (diff > 0.0) ? (diff / ctx->lambda_2) : 0.0;
}

double ergodev_l1_drift(const ergodev_ctx_t *ctx) {
    if (!ctx || ctx->active_words == 0U) return 0.0;
    const double inv_n = 1.0 / (double)ctx->active_words;
    double l1 = 0.0;
    for (uint32_t i = 0; i < ctx->num_words; ++i) {
        l1 += fabs((double)ctx->counts[i] * inv_n - ctx->nu[i]);
    }
    return l1;
}

/* Log-domain stabilized Sinkhorn for small epsilon (prevents exp(-d/eps) == 0 underflow) */
static double sinkhorn_log_domain(
    const ergodev_ctx_t *ctx,
    const double *mu,
    uint32_t max_iters
) {
    const uint32_t W = ctx->num_words;
    const double eps = ctx->sinkhorn_reg;
    double f[ERGODEV_MAX_WORDS] = {0};
    double g[ERGODEV_MAX_WORDS] = {0};
    double log_mu[ERGODEV_MAX_WORDS];
    double log_nu[ERGODEV_MAX_WORDS];

    for (uint32_t i = 0; i < W; ++i) {
        log_mu[i] = (mu[i] > 0.0) ? log(mu[i]) : -1e30;
        log_nu[i] = (ctx->nu[i] > 0.0) ? log(ctx->nu[i]) : -1e30;
    }

    for (uint32_t iter = 0; iter < max_iters; ++iter) {
        /* Update f_i = eps * log(mu_i) - eps * logsumexp_j((g_j - d_ij) / eps) */
        for (uint32_t i = 0; i < W; ++i) {
            if (mu[i] <= 0.0) { f[i] = -1e30; continue; }
            double max_val = -1e300;
            for (uint32_t j = 0; j < W; ++j) {
                if (ctx->nu[j] <= 0.0) continue;
                double val = (g[j] - ctx->dist_mat[i][j]) / eps;
                if (val > max_val) max_val = val;
            }
            double sum_exp = 0.0;
            for (uint32_t j = 0; j < W; ++j) {
                if (ctx->nu[j] <= 0.0) continue;
                sum_exp += exp(((g[j] - ctx->dist_mat[i][j]) / eps) - max_val);
            }
            f[i] = eps * (log_mu[i] - (max_val + log(sum_exp)));
        }

        /* Update g_j = eps * log(nu_j) - eps * logsumexp_i((f_i - d_ij) / eps) */
        for (uint32_t j = 0; j < W; ++j) {
            if (ctx->nu[j] <= 0.0) { g[j] = -1e30; continue; }
            double max_val = -1e300;
            for (uint32_t i = 0; i < W; ++i) {
                if (mu[i] <= 0.0) continue;
                double val = (f[i] - ctx->dist_mat[i][j]) / eps;
                if (val > max_val) max_val = val;
            }
            double sum_exp = 0.0;
            for (uint32_t i = 0; i < W; ++i) {
                if (mu[i] <= 0.0) continue;
                sum_exp += exp(((f[i] - ctx->dist_mat[i][j]) / eps) - max_val);
            }
            g[j] = eps * (log_nu[j] - (max_val + log(sum_exp)));
        }
    }

    double w1 = 0.0;
    for (uint32_t i = 0; i < W; ++i) {
        if (mu[i] <= 0.0) continue;
        for (uint32_t j = 0; j < W; ++j) {
            if (ctx->nu[j] <= 0.0) continue;
            double p_ij = exp((f[i] + g[j] - ctx->dist_mat[i][j]) / eps);
            w1 += p_ij * ctx->dist_mat[i][j];
        }
    }
    return w1;
}

double ergodev_score_sinkhorn(
    const ergodev_ctx_t *ctx,
    uint32_t max_iters,
    double *out_w1,
    double *out_defect
) {
    if (!ctx || ctx->active_words == 0U) {
        if (out_w1) *out_w1 = 0.0;
        if (out_defect) *out_defect = 0.0;
        return 0.0;
    }

    const uint32_t W = ctx->num_words;
    /* Normalize by active_words so partial warm-up windows still form a valid probability measure */
    const double inv_n = 1.0 / (double)ctx->active_words;

    double mu[ERGODEV_MAX_WORDS];
    double u[ERGODEV_MAX_WORDS];
    double v[ERGODEV_MAX_WORDS];

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

    double w1 = 0.0;
    if (ctx->use_log_sinkhorn) {
        w1 = sinkhorn_log_domain(ctx, mu, max_iters);
    } else {
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

        for (uint32_t i = 0; i < W; ++i) {
            if (u[i] == 0.0) continue;
            for (uint32_t j = 0; j < W; ++j) {
                double p_ij = u[i] * ctx->gibbs_k[i][j] * v[j];
                w1 += p_ij * ctx->dist_mat[i][j];
            }
        }
    }

    double prefactor = (defect > 0.0) ? (1.0 + ctx->eta * defect) : 1.0;
    double es = w1 * prefactor;
    if (out_w1) *out_w1 = w1;
    if (out_defect) *out_defect = defect;
    return es;
}

#ifndef ERGODEV_NO_MAIN
static void require_close(const char *label, double actual, double expected, double tol) {
    if (!isfinite(actual) || fabs(actual - expected) > tol) {
        fprintf(stderr, "FAIL [%s]: expected %.6f, got %.6f (tol %.6f)\n",
                label, expected, actual, tol);
        exit(EXIT_FAILURE);
    }
}

static void require_true(const char *label, int cond) {
    if (!cond) {
        fprintf(stderr, "FAIL [%s]: assertion failed\n", label);
        exit(EXIT_FAILURE);
    }
}

int main(void) {
    const double nu_2[9] = {
        1.0/40.0, 1.0/5.0,  1.0/40.0,
        1.0/5.0,  1.0/10.0, 1.0/5.0,
        1.0/40.0, 1.0/5.0,  1.0/40.0
    };
    const double jump_cost[9] = {
        0.0, 0.5, 2.0,
        0.5, 0.0, 0.5,
        2.0, 0.5, 0.0
    };
    double dist_mat[81];
    for (int w1 = 0; w1 < 9; ++w1) {
        int i = w1 / 3, j = w1 % 3;
        for (int w2 = 0; w2 < 9; ++w2) {
            int k = w2 / 3, l = w2 % 3;
            dist_mat[w1 * 9 + w2] = fabs((double)(i - k)) + 0.5 * fabs((double)(j - l));
        }
    }

    ergodev_ctx_t ctx;
    require_true("init valid", ergodev_init(&ctx, 3, 2, 401, 0.5, 1.0, nu_2, jump_cost, dist_mat, 0.005) == 0);

    printf("%-28s | %-10s | %-10s | %-10s\n", "Stream Source (N'=400)", "W_1", "Defect", "E_S Score");
    printf("-------------------------------------------------------------------\n");

    double w1, def, es;

    /* 1. Alternation 0, 2, 0, 2, ... */
    ergodev_reset(&ctx);
    for (uint32_t t = 0; t < 401; ++t) {
        require_true("push alt", ergodev_push(&ctx, (t & 1U) ? 2U : 0U) >= 0);
    }
    es = ergodev_score_sinkhorn(&ctx, 200, &w1, &def);
    printf("%-28s | %-10.4f | %-10.4f | %-10.4f\n", "Alternation 0 <-> 2", w1, def, es);
    require_close("Alt W_1", w1, 0.80, 1e-3);
    require_close("Alt Defect", def, 0.00, 1e-6);
    require_close("Alt E_S", es, 0.80, 1e-3);

    /* 2. Constant state 1, 1, 1, ... */
    ergodev_reset(&ctx);
    for (uint32_t t = 0; t < 401; ++t) {
        require_true("push const", ergodev_push(&ctx, 1U) >= 0);
    }
    es = ergodev_score_sinkhorn(&ctx, 200, &w1, &def);
    printf("%-28s | %-10.4f | %-10.4f | %-10.4f\n", "Constant State 1", w1, def, es);
    require_close("Const W_1", w1, 0.75, 1e-3);
    require_close("Const Defect", def, 1.00, 1e-6);
    require_close("Const E_S", es, 1.50, 1e-3);

    /* 3. Balanced Cut {0,2} */
    ergodev_reset(&ctx);
    const uint8_t cut_cycle[4] = {0, 0, 2, 2};
    for (uint32_t t = 0; t < 401; ++t) {
        ergodev_push(&ctx, cut_cycle[t & 3U]);
    }
    es = ergodev_score_sinkhorn(&ctx, 200, &w1, &def);
    printf("%-28s | %-10.4f | %-10.4f | %-10.4f\n", "Balanced Cut {0,2}", w1, def, es);
    require_close("Cut W_1", w1, 0.75, 1e-3);
    require_close("Cut Defect", def, 0.00, 1e-6);
    require_close("Cut E_S", es, 0.75, 1e-3);

    /* 4. One-Sided Submanifold {0,1} */
    ergodev_reset(&ctx);
    const uint8_t side_cycle[4] = {0, 0, 1, 1};
    for (uint32_t t = 0; t < 401; ++t) {
        ergodev_push(&ctx, side_cycle[t & 3U]);
    }
    es = ergodev_score_sinkhorn(&ctx, 200, &w1, &def);
    printf("%-28s | %-10.4f | %-10.4f | %-10.4f\n", "One-Sided Submanifold {0,1}", w1, def, es);
    require_close("Side W_1", w1, 0.75, 1e-3);
    require_close("Side Defect", def, 0.50, 1e-6);
    require_close("Side E_S", es, 1.125, 1e-3);

    /* =====================================================================
     * Regression Suite for Findings #1 - #5
     * ===================================================================== */
    ergodev_ctx_t bad_ctx;
    /* #1: m^ell overflow wrap (2^32 == 0 mod 2^32) must be rejected */
    require_true("reject m^ell wrap", ergodev_init(&bad_ctx, 2, 32, 401, 0.5, 1.0, nu_2, jump_cost, dist_mat, 0.005) < 0);

    /* #2: Partial window (first 15 symbols of Constant State 1 & 10 symbols of Alternation) */
    ergodev_reset(&ctx);
    for (uint32_t t = 0; t < 15; ++t) ergodev_push(&ctx, 1U);
    require_close("Short Const 15", ergodev_score_sinkhorn(&ctx, 200, &w1, &def), 1.50, 1e-3);

    ergodev_reset(&ctx);
    for (uint32_t t = 0; t < 11; ++t) ergodev_push(&ctx, (t & 1U) ? 2U : 0U); /* 10 pairs -> 5 (0,2), 5 (2,0) */
    require_close("Short Alt 100 iters", ergodev_score_sinkhorn(&ctx, 100, NULL, NULL), 0.80, 1e-3);
    require_close("Short Alt 200 iters", ergodev_score_sinkhorn(&ctx, 200, NULL, NULL), 0.80, 1e-3);

    /* #3: Tiny epsilon = 0.0005 log-domain stabilization */
    require_true("init eps=0.0005", ergodev_init(&bad_ctx, 3, 2, 401, 0.5, 1.0, nu_2, jump_cost, dist_mat, 0.0005) == 0);
    for (uint32_t t = 0; t < 401; ++t) ergodev_push(&bad_ctx, 1U);
    require_close("Const eps=0.0005", ergodev_score_sinkhorn(&bad_ctx, 200, NULL, NULL), 1.50, 1e-3);

    /* #4: Reject NaN and +Inf in init */
    require_true("reject NaN lambda_2", ergodev_init(&bad_ctx, 3, 2, 401, NAN, 1.0, nu_2, jump_cost, dist_mat, 0.005) < 0);
    require_true("reject NaN eps", ergodev_init(&bad_ctx, 3, 2, 401, 0.5, 1.0, nu_2, jump_cost, dist_mat, NAN) < 0);
    require_true("reject Inf eta", ergodev_init(&bad_ctx, 3, 2, 401, 0.5, INFINITY, nu_2, jump_cost, dist_mat, 0.005) < 0);

    /* #5: Reject out-of-range symbols (16, 255) in both push and push_fast */
    require_true("reject sym 255", ergodev_push(&ctx, 255U) == -1);
    require_true("reject sym 16 fast", ergodev_push_fast(&ctx, 16U) == -1);

    printf("-------------------------------------------------------------------\n");
    printf("All population targets and edge-case regression assertions passed!\n");
    return 0;
}
#endif /* ERGODEV_NO_MAIN */