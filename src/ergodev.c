#include "ergodev.h"

#include <math.h>
#include <float.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

static size_t align8(size_t x) {
    return (x + 7u) & ~(size_t)7u;
}

static int u32_pow_checked(uint32_t base, uint32_t exp, uint32_t limit, uint32_t *out) {
    uint64_t res = 1ULL;
    uint32_t i;
    for (i = 0; i < exp; ++i) {
        res *= (uint64_t)base;
        if (res > (uint64_t)limit) return -1;
    }
    *out = (uint32_t)res;
    return 0;
}

typedef struct {
    uint32_t off_nu, off_jump, off_lut, off_ring, off_counts;
    uint32_t off_hash_key, off_hash_val, off_hash_used, off_word_id;
    size_t total;
} layout_t;

static void compute_layout(uint32_t num_words, uint32_t m, uint32_t num_windows, int sparse, layout_t *L) {
    size_t off = align8(sizeof(ergodev_ctx_t));
    memset(L, 0, sizeof(*L));
    L->off_nu = (uint32_t)off;
    off = align8(off + (size_t)num_words * sizeof(double));
    L->off_jump = (uint32_t)off;
    off = align8(off + (size_t)num_words * sizeof(double));
    if (!sparse) {
        L->off_lut = (uint32_t)off;
        off = align8(off + (size_t)num_words * (size_t)m * sizeof(uint16_t));
    } else {
        L->off_hash_key = (uint32_t)off;
        off = align8(off + (size_t)ERGODEV_HASH_SIZE * sizeof(uint32_t));
        L->off_hash_val = (uint32_t)off;
        off = align8(off + (size_t)ERGODEV_HASH_SIZE * sizeof(uint16_t));
        L->off_hash_used = (uint32_t)off;
        off = align8(off + (size_t)ERGODEV_HASH_SIZE * sizeof(uint8_t));
        L->off_word_id = (uint32_t)off;
        off = align8(off + (size_t)num_words * sizeof(uint32_t));
    }
    L->off_ring = (uint32_t)off;
    off = align8(off + (size_t)num_windows * sizeof(uint16_t));
    L->off_counts = (uint32_t)off;
    off = align8(off + (size_t)num_words * sizeof(uint32_t));
    L->total = off;
}

static void bind_hot(ergodev_ctx_t *ctx, const layout_t *L) {
    unsigned char *base = (unsigned char *)ctx;
    ctx->ring = (uint16_t *)(base + L->off_ring);
    ctx->counts = (uint32_t *)(base + L->off_counts);
    ctx->lut = (L->off_lut != 0U) ? (uint16_t *)(base + L->off_lut) : NULL;
}

size_t ergodev_ctx_bytes(uint32_t num_words, uint32_t m, uint32_t num_windows, int sparse) {
    layout_t L;
    if (num_words == 0U || num_words > ERGODEV_MAX_WORDS || m == 0U || num_windows == 0U) return 0U;
    compute_layout(num_words, m, num_windows, sparse, &L);
    return L.total;
}

size_t ergodev_scratch_doubles(uint32_t num_words) {
    if (num_words == 0U || num_words > ERGODEV_MAX_WORDS) return 0U;
    return (size_t)num_words * (size_t)num_words + 6U * (size_t)num_words;
}

static int shape_ok(uint32_t m, uint32_t ell, uint32_t window_len, uint32_t *num_words, uint32_t *high_pow, uint32_t *num_windows) {
    if (m < 2U || m > ERGODEV_MAX_ALPHABET || ell < 2U || ell > ERGODEV_MAX_ELL) return -2;
    if (u32_pow_checked(m, ell, ERGODEV_MAX_WORDS, num_words) != 0 || *num_words == 0U) return -2;
    if (u32_pow_checked(m, ell - 1U, ERGODEV_MAX_WORDS, high_pow) != 0 || *high_pow == 0U) return -2;
    if (window_len < ell || window_len > ERGODEV_MAX_WINDOW) return -3;
    *num_windows = window_len - ell + 1U;
    return 0;
}

static void apply_layout(ergodev_ctx_t *ctx, const layout_t *L) {
    ctx->off_nu = L->off_nu;
    ctx->off_jump = L->off_jump;
    ctx->off_lut = L->off_lut;
    ctx->off_ring = L->off_ring;
    ctx->off_counts = L->off_counts;
    ctx->off_hash_key = L->off_hash_key;
    ctx->off_hash_val = L->off_hash_val;
    ctx->off_hash_used = L->off_hash_used;
    ctx->off_word_id = L->off_word_id;
    bind_hot(ctx, L);
}

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
) {
    uint32_t num_words = 0U, high_pow = 0U, num_windows = 0U;
    layout_t L;
    double nu_sum = 0.0, delta = 0.0, diameter = 0.0, max_jump = 0.0;
    uint32_t i, j;
    int shape;

    if (!ctx || !nu_ell || !jump_cost || !word_dist) return -1;
    shape = shape_ok(m, ell, window_len, &num_words, &high_pow, &num_windows);
    if (shape != 0) return shape;
    compute_layout(num_words, m, num_windows, 0, &L);
    if (ctx_bytes < L.total) return -6;
    if (!isfinite(lambda_2) || lambda_2 <= 0.0) return -4;
    if (!isfinite(eta) || eta < 0.0) return -4;

    for (i = 0; i < num_words; ++i) {
        if (!isfinite(nu_ell[i]) || nu_ell[i] < 0.0) return -5;
        if (!isfinite(jump_cost[i]) || jump_cost[i] < 0.0) return -5;
        nu_sum += nu_ell[i];
        delta += nu_ell[i] * jump_cost[i];
        if (jump_cost[i] > max_jump) max_jump = jump_cost[i];
        for (j = 0; j < num_words; ++j) {
            double d = word_dist[(size_t)i * num_words + j];
            if (!isfinite(d) || d < 0.0) return -5;
            if (d > diameter) diameter = d;
        }
    }
    if (fabs(nu_sum - 1.0) > 1e-6) return -5;
    if (fabs(delta - lambda_2) > 1e-5 * (1.0 + fabs(lambda_2))) return -5;

    memset(ctx, 0, ctx_bytes);
    apply_layout(ctx, &L);
    ctx->m = m;
    ctx->ell = ell;
    ctx->num_words = num_words;
    ctx->window_len = window_len;
    ctx->num_windows = num_windows;
    ctx->high_pow = high_pow;
    ctx->lambda_2 = lambda_2;
    ctx->eta = eta;
    ctx->max_jump = max_jump;
    ctx->diameter = diameter;
    ctx->sparse = 0U;
    ctx->n_vertices = high_pow;

    memcpy(ergodev_nu(ctx), nu_ell, (size_t)num_words * sizeof(double));
    memcpy(ergodev_jump(ctx), jump_cost, (size_t)num_words * sizeof(double));
    for (i = 0; i < num_words; ++i) {
        uint32_t suffix = (i % high_pow) * m;
        for (j = 0; j < m; ++j) {
            ergodev_lut(ctx)[(size_t)i * m + j] = (uint16_t)(suffix + j);
        }
    }
    return 0;
}

void ergodev_reset(ergodev_ctx_t *ctx) {
    if (!ctx) return;
    memset(ergodev_counts(ctx), 0, (size_t)ctx->num_words * sizeof(uint32_t));
    ctx->sym_filled = 0U;
    ctx->active_words = 0U;
    ctx->word_head = 0U;
    ctx->rolling_word = 0U;
    ctx->hist_len = 0U;
    memset(ctx->hist, 0, sizeof(ctx->hist));
}

int ergodev_push(ergodev_ctx_t *ctx, uint8_t symbol) {
    return ergodev_push_fast(ctx, symbol);
}

static uint32_t *hash_keys(ergodev_ctx_t *ctx) {
    return (uint32_t *)((unsigned char *)ctx + ctx->off_hash_key);
}
static uint16_t *hash_vals(ergodev_ctx_t *ctx) {
    return (uint16_t *)((unsigned char *)ctx + ctx->off_hash_val);
}
static uint8_t *hash_used(ergodev_ctx_t *ctx) {
    return (unsigned char *)ctx + ctx->off_hash_used;
}
static uint32_t *word_ids(ergodev_ctx_t *ctx) {
    return (uint32_t *)((unsigned char *)ctx + ctx->off_word_id);
}
static const uint32_t *word_ids_c(const ergodev_ctx_t *ctx) {
    return (const uint32_t *)((const unsigned char *)ctx + ctx->off_word_id);
}

static int hash_lookup(const ergodev_ctx_t *ctx, uint32_t key, uint32_t *slot) {
    const uint32_t *keys = (const uint32_t *)((const unsigned char *)ctx + ctx->off_hash_key);
    const uint16_t *vals = (const uint16_t *)((const unsigned char *)ctx + ctx->off_hash_val);
    const uint8_t *used = (const unsigned char *)ctx + ctx->off_hash_used;
    uint32_t mask = ERGODEV_HASH_SIZE - 1U;
    uint32_t i = (key * 0x9E3779B9u) & mask;
    uint32_t n;
    for (n = 0; n < ERGODEV_HASH_SIZE; ++n) {
        if (!used[i]) return 0;
        if (keys[i] == key) {
            *slot = vals[i];
            return 1;
        }
        i = (i + 1U) & mask;
    }
    return 0;
}

static int hash_insert(ergodev_ctx_t *ctx, uint32_t key, uint16_t slot) {
    uint32_t *keys = hash_keys(ctx);
    uint16_t *vals = hash_vals(ctx);
    uint8_t *used = hash_used(ctx);
    uint32_t mask = ERGODEV_HASH_SIZE - 1U;
    uint32_t i = (key * 0x9E3779B9u) & mask;
    uint32_t n;
    for (n = 0; n < ERGODEV_HASH_SIZE; ++n) {
        if (!used[i] || keys[i] == key) {
            used[i] = 1U;
            keys[i] = key;
            vals[i] = slot;
            return 0;
        }
        i = (i + 1U) & mask;
    }
    return -1;
}

static int ingest_slot(ergodev_ctx_t *ctx, uint16_t w_in) {
    if (ctx->sym_filled >= ctx->window_len) {
        uint32_t head = ctx->word_head;
        uint16_t *ring = ergodev_ring(ctx);
        uint32_t *counts = ergodev_counts(ctx);
        uint16_t w_out = ring[head];
        ring[head] = w_in;
        counts[w_out]--;
        counts[w_in]++;
        head++;
        ctx->word_head = (head == ctx->num_windows) ? 0U : head;
        return 1;
    }
    ctx->sym_filled++;
    {
        uint32_t head = ctx->word_head;
        ergodev_ring(ctx)[head] = w_in;
        ergodev_counts(ctx)[w_in]++;
        ctx->active_words++;
        head++;
        ctx->word_head = (head == ctx->num_windows) ? 0U : head;
        return (ctx->sym_filled >= ctx->window_len) ? 1 : 0;
    }
}

int ergodev_push_sparse(ergodev_ctx_t *ctx, uint8_t symbol) {
    uint32_t id = 0U;
    uint32_t slot = 0U;
    uint32_t k;
    if (!ctx || !ctx->sparse || (uint32_t)symbol >= ctx->m) return -1;

    if (ctx->hist_len < ctx->ell) {
        ctx->hist[ctx->hist_len++] = symbol;
        if (ctx->hist_len < ctx->ell) {
            ctx->sym_filled++;
            return 0;
        }
    } else {
        memmove(ctx->hist, ctx->hist + 1, (size_t)ctx->ell - 1U);
        ctx->hist[ctx->ell - 1U] = symbol;
    }
    for (k = 0; k < ctx->ell; ++k) id = id * ctx->m + ctx->hist[k];
    if (!hash_lookup(ctx, id, &slot)) slot = (uint32_t)ctx->num_words - 1U;
    if (ctx->hist_len == ctx->ell && ctx->sym_filled < ctx->ell - 1U) {
        /* hist_len reaches ell together with the first word; sym_filled counts warmup. */
    }
    /* Warmup symbols were counted above only while hist was short. The first
     * complete word is emitted with sym_filled still equal to ell-1. */
    if (ctx->sym_filled < ctx->ell - 1U) ctx->sym_filled = ctx->ell - 1U;
    return ingest_slot(ctx, (uint16_t)slot);
}

double ergodev_energy(const ergodev_ctx_t *ctx) {
    const uint32_t *counts;
    const double *jump;
    double sum = 0.0;
    uint32_t w;
    if (!ctx || ctx->active_words == 0U) return 0.0;
    counts = ergodev_counts_c(ctx);
    jump = ergodev_jump_c(ctx);
    for (w = 0; w < ctx->num_words; ++w) sum += (double)counts[w] * jump[w];
    return sum / (double)ctx->active_words;
}

double ergodev_defect(const ergodev_ctx_t *ctx) {
    double diff;
    if (!ctx || ctx->active_words == 0U || ctx->lambda_2 <= 0.0) return 0.0;
    diff = ctx->lambda_2 - ergodev_energy(ctx);
    return (diff > 0.0) ? (diff / ctx->lambda_2) : 0.0;
}

double ergodev_l1_drift(const ergodev_ctx_t *ctx) {
    const uint32_t *counts;
    const double *nu;
    double inv_n, l1 = 0.0;
    uint32_t i;
    if (!ctx || ctx->active_words == 0U) return 0.0;
    counts = ergodev_counts_c(ctx);
    nu = ergodev_nu_c(ctx);
    inv_n = 1.0 / (double)ctx->active_words;
    for (i = 0; i < ctx->num_words; ++i) l1 += fabs((double)counts[i] * inv_n - nu[i]);
    return l1;
}

double ergodev_es_upper(const ergodev_ctx_t *ctx) {
    const uint32_t *counts;
    const double *nu;
    const double *jump;
    double inv_n, l1, sum, energy, defect;
    uint32_t i, n, active;
    if (!ctx || ctx->active_words == 0U) return 0.0;
    active = ctx->active_words;
    n = ctx->num_words;
    counts = ergodev_counts_c(ctx);
    nu = ergodev_nu_c(ctx);
    jump = ergodev_jump_c(ctx);
    inv_n = 1.0 / (double)active;
    l1 = 0.0;
    sum = 0.0;
    for (i = 0; i < n; ++i) {
        double c = (double)counts[i];
        sum += c * jump[i];
        l1 += fabs(c * inv_n - nu[i]);
    }
    energy = sum * inv_n;
    defect = (ctx->lambda_2 - energy > 0.0) ? ((ctx->lambda_2 - energy) / ctx->lambda_2) : 0.0;
    return 0.5 * ctx->diameter * l1 * (1.0 + ctx->eta * defect);
}

int ergodev_screen(const ergodev_ctx_t *ctx, double level) {
    if (!ctx || !isfinite(level)) return 0;
    return ergodev_es_upper(ctx) < level;
}

int ergodev_score(
    const ergodev_ctx_t *ctx,
    ergodev_scratch_t *scratch,
    uint32_t max_iters,
    double residual_tol,
    ergodev_score_t *out
) {
    const uint32_t *counts;
    const double *nu;
    double *mu, *f, *g, *log_mu, *log_nu;
    double inv_n, l1_diff = 0.0, eps, defect;
    uint32_t n, i, j, iter;
    if (!out) return -1;
    out->w1 = 0.0;
    out->defect = 0.0;
    out->es = 0.0;
    out->residual = 0.0;
    out->converged = 0;
    if (!ctx || !scratch || !scratch->dist || !scratch->work || ctx->active_words == 0U) return -1;
    if (scratch->num_words != ctx->num_words || !(scratch->reg > 0.0) || !(residual_tol > 0.0)) return -1;

    n = ctx->num_words;
    counts = ergodev_counts_c(ctx);
    nu = ergodev_nu_c(ctx);
    inv_n = 1.0 / (double)ctx->active_words;
    mu = scratch->work;
    f = scratch->work + n;
    g = scratch->work + 2U * n;
    log_mu = scratch->work + 3U * n;
    log_nu = scratch->work + 4U * n;
    eps = scratch->reg;

    for (i = 0; i < n; ++i) {
        mu[i] = (double)counts[i] * inv_n;
        l1_diff += fabs(mu[i] - nu[i]);
        f[i] = 0.0;
        g[i] = 0.0;
        log_mu[i] = (mu[i] > 0.0) ? log(mu[i]) : -1e30;
        log_nu[i] = (nu[i] > 0.0) ? log(nu[i]) : -1e30;
    }
    defect = ergodev_defect(ctx);
    out->defect = defect;
    if (l1_diff < 1e-14) {
        out->converged = 1;
        out->residual = l1_diff;
        return 0;
    }

    for (iter = 0; iter < max_iters; ++iter) {
        for (i = 0; i < n; ++i) {
            double max_val, sum_exp;
            if (mu[i] <= 0.0) { f[i] = -1e30; continue; }
            max_val = -1e300;
            for (j = 0; j < n; ++j) {
                double val;
                if (nu[j] <= 0.0) continue;
                val = (g[j] - scratch->dist[(size_t)i * n + j]) / eps;
                if (val > max_val) max_val = val;
            }
            sum_exp = 0.0;
            for (j = 0; j < n; ++j) {
                if (nu[j] <= 0.0) continue;
                sum_exp += exp(((g[j] - scratch->dist[(size_t)i * n + j]) / eps) - max_val);
            }
            f[i] = eps * (log_mu[i] - (max_val + log(sum_exp)));
        }
        for (j = 0; j < n; ++j) {
            double max_val, sum_exp;
            if (nu[j] <= 0.0) { g[j] = -1e30; continue; }
            max_val = -1e300;
            for (i = 0; i < n; ++i) {
                double val;
                if (mu[i] <= 0.0) continue;
                val = (f[i] - scratch->dist[(size_t)i * n + j]) / eps;
                if (val > max_val) max_val = val;
            }
            sum_exp = 0.0;
            for (i = 0; i < n; ++i) {
                if (mu[i] <= 0.0) continue;
                sum_exp += exp(((f[i] - scratch->dist[(size_t)i * n + j]) / eps) - max_val);
            }
            g[j] = eps * (log_nu[j] - (max_val + log(sum_exp)));
        }
        if (((iter + 1U) % 25U) == 0U || iter + 1U == max_iters) {
            double residual_now = 0.0;
            for (i = 0; i < n; ++i) {
                double row = 0.0, max_val = -1e300, sum_exp = 0.0;
                if (mu[i] <= 0.0) continue;
                for (j = 0; j < n; ++j) {
                    double val;
                    if (nu[j] <= 0.0) continue;
                    val = (f[i] + g[j] - scratch->dist[(size_t)i * n + j]) / eps;
                    if (val > max_val) max_val = val;
                }
                if (max_val > -1e299) {
                    for (j = 0; j < n; ++j) {
                        if (nu[j] <= 0.0) continue;
                        sum_exp += exp(((f[i] + g[j] - scratch->dist[(size_t)i * n + j]) / eps) - max_val);
                    }
                    row = exp(max_val + log(sum_exp));
                }
                residual_now += fabs(row - mu[i]);
            }
            for (j = 0; j < n; ++j) {
                double col = 0.0, max_val = -1e300, sum_exp = 0.0;
                if (nu[j] <= 0.0) continue;
                for (i = 0; i < n; ++i) {
                    double val;
                    if (mu[i] <= 0.0) continue;
                    val = (f[i] + g[j] - scratch->dist[(size_t)i * n + j]) / eps;
                    if (val > max_val) max_val = val;
                }
                if (max_val > -1e299) {
                    for (i = 0; i < n; ++i) {
                        if (mu[i] <= 0.0) continue;
                        sum_exp += exp(((f[i] + g[j] - scratch->dist[(size_t)i * n + j]) / eps) - max_val);
                    }
                    col = exp(max_val + log(sum_exp));
                }
                residual_now += fabs(col - nu[j]);
            }
            if (residual_now <= residual_tol) break;
        }
    }

    {
        double residual = 0.0;
        for (i = 0; i < n; ++i) {
            double row = 0.0;
            if (mu[i] <= 0.0) continue;
            row = 0.0;
            {
                double max_val = -1e300, sum_exp = 0.0;
                for (j = 0; j < n; ++j) {
                    double val;
                    if (nu[j] <= 0.0) continue;
                    val = (f[i] + g[j] - scratch->dist[(size_t)i * n + j]) / eps;
                    if (val > max_val) max_val = val;
                }
                if (max_val > -1e299) {
                    for (j = 0; j < n; ++j) {
                        if (nu[j] <= 0.0) continue;
                        sum_exp += exp(((f[i] + g[j] - scratch->dist[(size_t)i * n + j]) / eps) - max_val);
                    }
                    row = exp(max_val + log(sum_exp));
                }
            }
            residual += fabs(row - mu[i]);
        }
        for (j = 0; j < n; ++j) {
            double col = 0.0, max_val = -1e300, sum_exp = 0.0;
            if (nu[j] <= 0.0) continue;
            for (i = 0; i < n; ++i) {
                double val;
                if (mu[i] <= 0.0) continue;
                val = (f[i] + g[j] - scratch->dist[(size_t)i * n + j]) / eps;
                if (val > max_val) max_val = val;
            }
            if (max_val > -1e299) {
                for (i = 0; i < n; ++i) {
                    if (mu[i] <= 0.0) continue;
                    sum_exp += exp(((f[i] + g[j] - scratch->dist[(size_t)i * n + j]) / eps) - max_val);
                }
                col = exp(max_val + log(sum_exp));
            }
            residual += fabs(col - nu[j]);
        }
        out->residual = residual;
        if (!(residual <= residual_tol)) {
            out->w1 = nan("");
            out->es = nan("");
            out->converged = 0;
            return -2;
        }
        {
            double w1 = 0.0;
            for (i = 0; i < n; ++i) {
                if (mu[i] <= 0.0) continue;
                for (j = 0; j < n; ++j) {
                    double arg, p;
                    if (nu[j] <= 0.0) continue;
                    arg = (f[i] + g[j] - scratch->dist[(size_t)i * n + j]) / eps;
                    if (arg > 700.0) {
                        out->w1 = nan("");
                        out->es = nan("");
                        out->converged = 0;
                        return -2;
                    }
                    p = (arg < -700.0) ? 0.0 : exp(arg);
                    w1 += p * scratch->dist[(size_t)i * n + j];
                }
            }
            out->w1 = w1;
            out->es = w1 * (1.0 + ctx->eta * defect);
            out->converged = 1;
            return 0;
        }
    }
}

int ergodev_scratch_bind(
    ergodev_scratch_t *scratch,
    uint32_t num_words,
    double *storage,
    size_t n_doubles,
    const double *word_dist,
    double sinkhorn_reg
) {
    size_t need, n2, i;
    if (!scratch || !storage || !word_dist) return -1;
    if (num_words == 0U || num_words > ERGODEV_MAX_WORDS) return -2;
    if (!isfinite(sinkhorn_reg) || sinkhorn_reg <= 0.0) return -4;
    need = ergodev_scratch_doubles(num_words);
    if (n_doubles < need) return -6;
    n2 = (size_t)num_words * (size_t)num_words;
    for (i = 0; i < n2; ++i) {
        if (!isfinite(word_dist[i]) || word_dist[i] < 0.0) return -5;
    }
    memcpy(storage, word_dist, n2 * sizeof(double));
    scratch->num_words = num_words;
    scratch->reg = sinkhorn_reg;
    scratch->dist = storage;
    scratch->work = storage + n2;
    return 0;
}

static double azuma_dev(uint32_t N, uint32_t Np, double kappa, double lip, double delta) {
    double coeff = ((double)Np * (double)Np) * (1.0 - kappa) * (1.0 - kappa) /
                   (2.0 * (double)N * lip * lip);
    return sqrt(-log(delta) / coeff);
}

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
) {
    uint32_t Np;
    double bracket, B, c, M, L, b, disc, tF, tA, sE, joint, dw, de;
    if (!out) return -1;
    memset(out, 0, sizeof(*out));
    if (!(diameter > 0.0) || atoms < 1U || ell < 2U || N <= ell) return -3;
    if (!(kappa >= 0.0) || !(kappa < 1.0) || !(eta >= 0.0)) return -4;
    if (!(lambda_2 > 0.0) || !(max_jump > 0.0) || !(delta > 0.0) || !(delta < 1.0)) return -4;
    Np = N - ell + 1U;
    bracket = (double)(2U * ell - 1U) + 2.0 * kappa / (1.0 - kappa);
    B = 0.5 * diameter * sqrt(((double)atoms / (double)Np) * bracket);
    c = diameter * (double)ell / (double)Np;
    M = c / (1.0 - kappa);
    tA = azuma_dev(N, Np, kappa, diameter * (double)ell, delta);
    /* sigma^2 = B^2 >= Var(W_1), because W_1 <= (D/2) sqrt(atoms) ||mu-nu||_2
     * and E||mu-nu||_2^2 = sum Var(mu(w)) <= bracket / N' by Lemma 8.3. */
    L = -log(delta);
    b = 2.0 * L * M / 3.0;
    disc = b * b + 8.0 * L * B * B;
    tF = 0.5 * (b + sqrt(disc));
    dw = 0.8 * delta;
    de = 0.2 * delta;
    sE = azuma_dev(N, Np, kappa, max_jump * (double)ell, de);
    joint = (B + azuma_dev(N, Np, kappa, diameter * (double)ell, dw)) *
            (1.0 + eta * fmin(1.0, sE / lambda_2));
    out->B_n = B;
    out->azuma_t = tA;
    out->azuma_threshold = (1.0 + eta) * (B + tA);
    out->freedman_t = tF;
    out->freedman_threshold = (1.0 + eta) * (B + tF);
    out->sigma2 = B * B;
    out->joint_threshold = joint;
    if (alternative > 0.0) {
        out->azuma_separates = out->azuma_threshold < alternative;
        out->freedman_separates = out->freedman_threshold < alternative;
        out->joint_separates = out->joint_threshold < alternative;
    }
    return 0;
}

static int cmp_double(const void *a, const void *b) {
    double x = *(const double *)a;
    double y = *(const double *)b;
    if (x < y) return -1;
    if (x > y) return 1;
    return 0;
}

double ergodev_empirical_quantile(double *scores, uint32_t n, double p) {
    double kreal;
    uint32_t k;
    if (!scores || n == 0U || !(p > 0.0) || p > 1.0) return nan("");
    qsort(scores, n, sizeof(double), cmp_double);
    kreal = ceil(p * (double)n);
    if (kreal < 1.0) kreal = 1.0;
    k = (uint32_t)kreal - 1U;
    if (k >= n) k = n - 1U;
    return scores[k];
}

int ergodev_window_symbols(const ergodev_ctx_t *ctx, uint8_t *out, uint32_t cap) {
    uint32_t nwords, start, i, nsym, id, k;
    const uint16_t *ring;
    if (!ctx || !out || ctx->active_words == 0U) return -1;
    nwords = ctx->active_words;
    nsym = nwords + ctx->ell - 1U;
    if (cap < nsym) return -6;
    ring = ergodev_ring_c(ctx);
    start = (ctx->sym_filled >= ctx->window_len) ? ctx->word_head : 0U;
    for (i = 0; i < nwords; ++i) {
        uint16_t slot = ring[(start + i) % ctx->num_windows];
        if (ctx->sparse) {
            id = word_ids_c(ctx)[slot];
            if (id == UINT32_MAX) return -5;
        } else {
            id = slot;
        }
        out[i] = (uint8_t)(id / ctx->high_pow);
        if (i + 1U == nwords) {
            uint32_t tmp = id;
            uint8_t syms[ERGODEV_MAX_ELL];
            for (k = ctx->ell; k-- > 0U; ) {
                syms[k] = (uint8_t)(tmp % ctx->m);
                tmp /= ctx->m;
            }
            for (k = 0; k < ctx->ell; ++k) out[(nwords - 1U) + k] = syms[k];
        }
    }
    return (int)nsym;
}

uint32_t ergodev_min_period(const uint8_t *symbols, uint32_t n) {
    uint32_t p, i;
    if (!symbols || n == 0U) return 0U;
    for (p = 1U; p < n; ++p) {
        for (i = 0; i + p < n; ++i) {
            if (symbols[i] != symbols[i + p]) break;
        }
        if (i + p >= n) return p;
    }
    return n;
}

int ergodev_alarm_union(
    double es,
    int es_valid,
    double tau_es,
    uint32_t period,
    uint32_t n_symbols,
    uint32_t tau_period,
    ergodev_alarm_t *out
) {
    if (!out) return -1;
    memset(out, 0, sizeof(*out));
    out->es = es;
    out->period = period;
    out->alarm_transport = (es_valid && isfinite(es) && es > tau_es) ? 1 : 0;
    out->alarm_period = (n_symbols > 0U && period > 0U && period < n_symbols && period <= tau_period) ? 1 : 0;
    out->alarm = (out->alarm_transport || out->alarm_period) ? 1 : 0;
    return 0;
}

int ergodev_thread_push_event(ergodev_thread_slot_t *slot, uint8_t syscall_class) {
    if (!slot || !slot->per_event) return -1;
    return ergodev_push_fast(slot->per_event, syscall_class);
}

int ergodev_thread_push_quantum(ergodev_thread_slot_t *slot, uint8_t occupancy_symbol) {
    if (!slot || !slot->per_quantum) return -1;
    return ergodev_push_fast(slot->per_quantum, occupancy_symbol);
}

/* ---- spectrum and corpus fit ---- */

static int jacobi_sym(double *A, double *V, double *evals, uint32_t n) {
    uint32_t sweep, p, q, r, i;
    double scale = 1.0;
    if (n == 0U) return -2;
    for (i = 0; i < n; ++i) scale += fabs(A[(size_t)i * n + i]);
    for (i = 0; i < n * n; ++i) V[i] = 0.0;
    for (i = 0; i < n; ++i) V[(size_t)i * n + i] = 1.0;
    for (sweep = 0; sweep < 80U; ++sweep) {
        double off = 0.0;
        for (p = 0; p < n; ++p) {
            for (q = p + 1U; q < n; ++q) off += A[(size_t)p * n + q] * A[(size_t)p * n + q];
        }
        if (off <= 1e-28 * scale * scale) break;
        for (p = 0; p < n; ++p) {
            for (q = p + 1U; q < n; ++q) {
                double apq = A[(size_t)p * n + q];
                double app, aqq, tau, t, c, s;
                if (fabs(apq) <= 1e-18 * scale) continue;
                app = A[(size_t)p * n + p];
                aqq = A[(size_t)q * n + q];
                tau = (aqq - app) / (2.0 * apq);
                if (tau >= 0.0) t = 1.0 / (tau + sqrt(1.0 + tau * tau));
                else t = -1.0 / (-tau + sqrt(1.0 + tau * tau));
                c = 1.0 / sqrt(1.0 + t * t);
                s = t * c;
                A[(size_t)p * n + p] = app - t * apq;
                A[(size_t)q * n + q] = aqq + t * apq;
                A[(size_t)p * n + q] = 0.0;
                A[(size_t)q * n + p] = 0.0;
                for (r = 0; r < n; ++r) {
                    double arp, arq;
                    if (r == p || r == q) continue;
                    arp = A[(size_t)r * n + p];
                    arq = A[(size_t)r * n + q];
                    A[(size_t)r * n + p] = A[(size_t)p * n + r] = c * arp - s * arq;
                    A[(size_t)r * n + q] = A[(size_t)q * n + r] = s * arp + c * arq;
                }
                for (r = 0; r < n; ++r) {
                    double vrp = V[(size_t)r * n + p];
                    double vrq = V[(size_t)r * n + q];
                    V[(size_t)r * n + p] = c * vrp - s * vrq;
                    V[(size_t)r * n + q] = s * vrp + c * vrq;
                }
            }
        }
    }
    for (i = 0; i < n; ++i) evals[i] = A[(size_t)i * n + i];
    /* Order eigenvalues, carrying eigenvector columns with them. */
    for (p = 0; p < n; ++p) {
        uint32_t best = p;
        for (q = p + 1U; q < n; ++q) if (evals[q] < evals[best]) best = q;
        if (best != p) {
            double tmp = evals[p];
            evals[p] = evals[best];
            evals[best] = tmp;
            for (r = 0; r < n; ++r) {
                double v = V[(size_t)r * n + p];
                V[(size_t)r * n + p] = V[(size_t)r * n + best];
                V[(size_t)r * n + best] = v;
            }
        }
    }
    return 0;
}

static double symbol_delta(const double *metric, uint32_t m, uint32_t a, uint32_t b) {
    if (!metric) return fabs((double)a - (double)b);
    return metric[(size_t)a * m + b];
}

static void decode_word(uint32_t id, uint32_t m, uint32_t ell, uint8_t *syms) {
    uint32_t k;
    for (k = ell; k-- > 0U; ) {
        syms[k] = (uint8_t)(id % m);
        id /= m;
    }
}

static double word_distance(const double *metric, uint32_t m, uint32_t ell, double theta,
                            uint32_t ida, uint32_t idb) {
    uint8_t sa[ERGODEV_MAX_ELL], sb[ERGODEV_MAX_ELL];
    memset(sa, 0, sizeof(sa));
    memset(sb, 0, sizeof(sb));
    double d = 0.0, th = 1.0;
    uint32_t k;
    decode_word(ida, m, ell, sa);
    decode_word(idb, m, ell, sb);
    for (k = 0; k < ell; ++k) {
        d += th * symbol_delta(metric, m, sa[k], sb[k]);
        th *= theta;
    }
    return d;
}

static int spectrum_from_flux(
    double *xi, double *q, double *P, double *Lmat, double *V, double *evals,
    double *psi, uint32_t Vn, int q_normalize, double *lambda, double *lambda_next, double *gap, double *kappa
) {
    uint32_t u, v, a, b;
    double kap = 0.0;
    for (u = 0; u < Vn; ++u) {
        q[u] = 0.0;
        for (v = 0; v < Vn; ++v) q[u] += xi[(size_t)u * Vn + v];
    }
    {
        double qs = 0.0;
        for (u = 0; u < Vn; ++u) qs += q[u];
        if (qs <= 0.0) return -5;
        for (u = 0; u < Vn; ++u) q[u] /= qs;
    }
    for (u = 0; u < Vn; ++u) {
        double out = 0.0;
        for (v = 0; v < Vn; ++v) out += xi[(size_t)u * Vn + v];
        if (out <= 0.0) {
            for (v = 0; v < Vn; ++v) P[(size_t)u * Vn + v] = (u == v) ? 1.0 : 0.0;
        } else {
            for (v = 0; v < Vn; ++v) P[(size_t)u * Vn + v] = xi[(size_t)u * Vn + v] / out;
        }
    }
    for (a = 0; a < Vn; ++a) {
        for (b = 0; b < Vn; ++b) {
            double tv = 0.0;
            for (u = 0; u < Vn; ++u) tv += fabs(P[(size_t)a * Vn + u] - P[(size_t)b * Vn + u]);
            tv *= 0.5;
            if (tv > kap) kap = tv;
        }
    }
    memset(Lmat, 0, (size_t)Vn * Vn * sizeof(double));
    for (u = 0; u < Vn; ++u) {
        for (v = u + 1U; v < Vn; ++v) {
            double w = xi[(size_t)u * Vn + v] + xi[(size_t)v * Vn + u];
            Lmat[(size_t)u * Vn + u] += w;
            Lmat[(size_t)v * Vn + v] += w;
            Lmat[(size_t)u * Vn + v] -= w;
            Lmat[(size_t)v * Vn + u] -= w;
        }
    }
    if (q_normalize) {
        for (u = 0; u < Vn; ++u) {
            for (v = 0; v < Vn; ++v) {
                double dq = sqrt(fmax(q[u], 1e-300) * fmax(q[v], 1e-300));
                Lmat[(size_t)u * Vn + v] /= dq;
            }
        }
    }
    if (jacobi_sym(Lmat, V, evals, Vn) != 0) return -5;
    if (Vn < 2U || !(evals[1] > 1e-10)) return -5;
    for (u = 0; u < Vn; ++u) psi[u] = V[(size_t)u * Vn + 1U];
    if (q_normalize) {
        double mean = 0.0, var = 0.0;
        for (u = 0; u < Vn; ++u) psi[u] /= sqrt(fmax(q[u], 1e-300));
        for (u = 0; u < Vn; ++u) mean += q[u] * psi[u];
        for (u = 0; u < Vn; ++u) psi[u] -= mean;
        for (u = 0; u < Vn; ++u) var += q[u] * psi[u] * psi[u];
        if (!(var > 0.0)) return -5;
        var = sqrt(var);
        for (u = 0; u < Vn; ++u) psi[u] /= var;
    } else {
        double mean = 0.0, norm = 0.0;
        for (u = 0; u < Vn; ++u) mean += psi[u];
        mean /= (double)Vn;
        for (u = 0; u < Vn; ++u) psi[u] -= mean;
        for (u = 0; u < Vn; ++u) norm += psi[u] * psi[u];
        if (!(norm > 0.0)) return -5;
        norm = sqrt(norm);
        for (u = 0; u < Vn; ++u) psi[u] /= norm;
    }
    if (psi[0] < 0.0) {
        for (u = 0; u < Vn; ++u) psi[u] = -psi[u];
    }
    *lambda = evals[1];
    *lambda_next = (Vn >= 3U) ? evals[2] : evals[1];
    *gap = *lambda;
    if (Vn >= 3U) {
        double sep = *lambda_next - *lambda;
        if (sep < *gap) *gap = sep;
    }
    if (*gap < 0.0) *gap = 0.0;
    *kappa = kap;
    return 0;
}

static int collect_dense_counts(uint32_t m, uint32_t ell, uint32_t high_pow, const uint8_t *corpus,
                                uint32_t corpus_len, double *counts, uint32_t *ntrans) {
    uint32_t i, filled = 0U, id = 0U, n = 0U;
    for (i = 0; i < corpus_len; ++i) {
        if ((uint32_t)corpus[i] >= m) return -5;
        if (filled < ell - 1U) {
            id = id * m + corpus[i];
            filled++;
            continue;
        }
        if (filled == ell - 1U) id = id * m + corpus[i];
        else id = (id % high_pow) * m + corpus[i];
        filled = ell;
        counts[id] += 1.0;
        n++;
    }
    *ntrans = n;
    return (n == 0U) ? -3 : 0;
}

int ergodev_plan(
    uint32_t m,
    uint32_t ell,
    uint32_t window_len,
    const uint8_t *corpus,
    uint32_t corpus_len,
    ergodev_plan_t *out
) {
    uint32_t num_words = 0U, high_pow = 0U, num_windows = 0U;
    int shape;
    if (!out) return -1;
    memset(out, 0, sizeof(*out));
    shape = shape_ok(m, ell, window_len, &num_words, &high_pow, &num_windows);
    if (shape == 0) {
        out->sparse = 0;
        out->num_words = num_words;
        out->num_windows = num_windows;
        out->ctx_bytes = ergodev_ctx_bytes(num_words, m, num_windows, 0);
        out->error = 0;
        (void)corpus;
        (void)corpus_len;
        return 0;
    }
    if (shape == -3) {
        out->error = -3;
        return -3;
    }
    /* Dense power does not fit: sparse support, counted from the corpus. */
    if (m < 2U || m > ERGODEV_MAX_ALPHABET || ell < 2U || ell > ERGODEV_MAX_ELL) {
        out->error = -2;
        return -2;
    }
    if (window_len < ell || window_len > ERGODEV_MAX_WINDOW) {
        out->error = -3;
        return -3;
    }
    if (u32_pow_checked(m, ell - 1U, UINT32_MAX, &high_pow) != 0) {
        out->error = -2;
        return -2;
    }
    if (!corpus || corpus_len < ell) {
        out->error = -3;
        return -3;
    }
    {
        uint32_t keys[ERGODEV_HASH_SIZE];
        uint8_t used[ERGODEV_HASH_SIZE];
        uint32_t nunique = 0U, i, filled = 0U, id = 0U;
        memset(used, 0, sizeof(used));
        for (i = 0; i < corpus_len; ++i) {
            uint32_t mask, h, n;
            int found = 0;
            if ((uint32_t)corpus[i] >= m) {
                out->error = -5;
                return -5;
            }
            if (filled < ell - 1U) {
                id = id * m + corpus[i];
                filled++;
                continue;
            }
            if (filled == ell - 1U) id = id * m + corpus[i];
            else id = (id % high_pow) * m + corpus[i];
            filled = ell;
            mask = ERGODEV_HASH_SIZE - 1U;
            h = (id * 0x9E3779B9u) & mask;
            for (n = 0; n < ERGODEV_HASH_SIZE; ++n) {
                if (!used[h]) break;
                if (keys[h] == id) { found = 1; break; }
                h = (h + 1U) & mask;
            }
            if (!found) {
                if (nunique + 1U >= ERGODEV_MAX_WORDS) {
                    out->error = -7;
                    return -7;
                }
                used[h] = 1U;
                keys[h] = id;
                nunique++;
            }
        }
        if (nunique == 0U) {
            out->error = -3;
            return -3;
        }
        out->sparse = 1;
        out->num_words = nunique + 1U; /* sink */
        out->num_windows = window_len - ell + 1U;
        out->ctx_bytes = ergodev_ctx_bytes(out->num_words, m, out->num_windows, 1);
        out->error = 0;
        return 0;
    }
}

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
) {
    ergodev_plan_t plan;
    int rc;
    if (!ctx || !corpus || !dist_out) return -1;
    if (!isfinite(eta) || eta < 0.0) return -4;
    if (!(theta > 0.0) || theta > 1.0) return -4;
    rc = ergodev_plan(m, ell, window_len, corpus, corpus_len, &plan);
    if (rc != 0) return rc;
    if (ctx_bytes < plan.ctx_bytes) return -6;
    if (dist_cap_words < plan.num_words) return -6;

    if (!plan.sparse) {
        uint32_t num_words = plan.num_words;
        uint32_t high_pow = 0U, num_windows = plan.num_windows, ntrans = 0U;
        double *counts = NULL, *xi = NULL, *q = NULL, *P = NULL, *Lmat = NULL, *V = NULL, *evals = NULL, *psi = NULL;
        double *nu_tmp = NULL, *jump_tmp = NULL;
        uint32_t Vn, w, pass;
        double lambda = 0.0, lambda_next = 0.0, gap = 0.0, kap = 0.0, delta = 0.0, max_jump = 0.0, diameter = 0.0;
        layout_t L;
        int ok = -5;
        u32_pow_checked(m, ell - 1U, ERGODEV_MAX_WORDS, &high_pow);
        Vn = high_pow;
        counts = (double *)calloc(num_words, sizeof(double));
        nu_tmp = (double *)calloc(num_words, sizeof(double));
        jump_tmp = (double *)calloc(num_words, sizeof(double));
        xi = (double *)calloc((size_t)Vn * Vn, sizeof(double));
        q = (double *)calloc(Vn, sizeof(double));
        P = (double *)calloc((size_t)Vn * Vn, sizeof(double));
        Lmat = (double *)calloc((size_t)Vn * Vn, sizeof(double));
        V = (double *)calloc((size_t)Vn * Vn, sizeof(double));
        evals = (double *)calloc(Vn, sizeof(double));
        psi = (double *)calloc(Vn, sizeof(double));
        if (!counts || !nu_tmp || !jump_tmp || !xi || !q || !P || !Lmat || !V || !evals || !psi) {
            ok = -6;
            goto dense_done;
        }
        if (collect_dense_counts(m, ell, high_pow, corpus, corpus_len, counts, &ntrans) != 0) {
            ok = -5;
            goto dense_done;
        }
        for (pass = 0; pass < 2U; ++pass) {
            double nrm;
            if (pass == 1U) {
                for (w = 0; w < num_words; ++w) counts[w] += 1.0;
                ntrans += num_words;
            }
            nrm = (double)ntrans;
            memset(xi, 0, (size_t)Vn * Vn * sizeof(double));
            for (w = 0; w < num_words; ++w) {
                uint32_t prefix = w / m;
                uint32_t suffix = w % high_pow;
                nu_tmp[w] = counts[w] / nrm;
                xi[(size_t)prefix * Vn + suffix] += nu_tmp[w];
            }
            if (spectrum_from_flux(xi, q, P, Lmat, V, evals, psi, Vn, q_normalize,
                                   &lambda, &lambda_next, &gap, &kap) == 0) {
                ok = 0;
                break;
            }
            ok = -5;
        }
        if (ok != 0) goto dense_done;
        max_jump = 0.0;
        delta = 0.0;
        for (w = 0; w < num_words; ++w) {
            uint32_t prefix = w / m;
            uint32_t suffix = w % high_pow;
            double jmp = (psi[prefix] - psi[suffix]) * (psi[prefix] - psi[suffix]);
            jump_tmp[w] = jmp;
            if (jmp > max_jump) max_jump = jmp;
            delta += nu_tmp[w] * jmp;
        }
        /* Identity Delta_* = lambda by using the realized Dirichlet form. */
        if (!(delta > 0.0)) { ok = -5; goto dense_done; }
        if ((psi_out && psi_cap < Vn) || (q_out && q_cap < Vn)) { ok = -6; goto dense_done; }
        lambda = delta;
        diameter = 0.0;
        for (w = 0; w < num_words; ++w) {
            uint32_t z;
            for (z = 0; z < num_words; ++z) {
                double d = word_distance(symbol_metric, m, ell, theta, w, z);
                dist_out[(size_t)w * num_words + z] = d;
                if (d > diameter) diameter = d;
            }
        }
        compute_layout(num_words, m, num_windows, 0, &L);
        memset(ctx, 0, ctx_bytes);
        apply_layout(ctx, &L);
        ctx->m = m;
        ctx->ell = ell;
        ctx->num_words = num_words;
        ctx->window_len = window_len;
        ctx->num_windows = num_windows;
        ctx->high_pow = high_pow;
        ctx->lambda_2 = lambda;
        ctx->lambda_next = lambda_next;
        ctx->spectral_gap = gap;
        ctx->eta = eta;
        ctx->max_jump = (max_jump > 0.0) ? max_jump : 1e-15;
        ctx->diameter = (diameter > 0.0) ? diameter : 1e-15;
        ctx->kappa = kap;
        ctx->q_normalized = q_normalize ? 1 : 0;
        ctx->sparse = 0U;
        ctx->n_vertices = Vn;
        memcpy(ergodev_nu(ctx), nu_tmp, (size_t)num_words * sizeof(double));
        memcpy(ergodev_jump(ctx), jump_tmp, (size_t)num_words * sizeof(double));
        for (w = 0; w < num_words; ++w) {
            uint32_t suffix = (w % high_pow) * m;
            uint32_t s;
            for (s = 0; s < m; ++s) ergodev_lut(ctx)[(size_t)w * m + s] = (uint16_t)(suffix + s);
        }
        if (psi_out) memcpy(psi_out, psi, (size_t)Vn * sizeof(double));
        if (q_out) memcpy(q_out, q, (size_t)Vn * sizeof(double));
        ok = 0;
    dense_done:
        free(counts); free(nu_tmp); free(jump_tmp); free(xi); free(q); free(P);
        free(Lmat); free(V); free(evals); free(psi);
        return ok;
    }

    /* Sparse support: distinct ell-words of the corpus, plus one sink. */
    {
        typedef struct { uint32_t id; uint32_t count; uint32_t prefix; uint32_t suffix; } item_t;
        item_t *items = NULL;
        uint32_t *verts = NULL;
        uint32_t nitems = 0U, Vn = 0U, i, filled = 0U, id = 0U, ntrans = 0U, high_pow = 0U;
        double *xi = NULL, *q = NULL, *P = NULL, *Lmat = NULL, *Vmat = NULL, *evals = NULL, *psi = NULL;
        double lambda = 0.0, lambda_next = 0.0, gap = 0.0, kap = 0.0, delta = 0.0, max_jump = 0.0, diameter = 0.0;
        layout_t L;
        int ok = -5;
        uint32_t num_words, num_windows;
        if (u32_pow_checked(m, ell - 1U, UINT32_MAX, &high_pow) != 0) return -2;
        items = (item_t *)calloc(ERGODEV_MAX_WORDS, sizeof(item_t));
        verts = (uint32_t *)calloc(ERGODEV_MAX_WORDS, sizeof(uint32_t));
        if (!items || !verts) { ok = -6; goto sparse_done; }
        for (i = 0; i < corpus_len; ++i) {
            uint32_t t;
            int found = 0;
            if ((uint32_t)corpus[i] >= m) { ok = -5; goto sparse_done; }
            if (filled < ell - 1U) {
                id = id * m + corpus[i];
                filled++;
                continue;
            }
            if (filled == ell - 1U) id = id * m + corpus[i];
            else id = (id % high_pow) * m + corpus[i];
            filled = ell;
            for (t = 0; t < nitems; ++t) if (items[t].id == id) { items[t].count++; found = 1; break; }
            if (!found) {
                if (nitems + 1U >= ERGODEV_MAX_WORDS) { ok = -7; goto sparse_done; }
                items[nitems].id = id;
                items[nitems].count = 1U;
                items[nitems].prefix = id / m;
                items[nitems].suffix = id % high_pow;
                nitems++;
            }
            ntrans++;
        }
        if (nitems == 0U || ntrans == 0U) { ok = -3; goto sparse_done; }
        for (i = 0; i < nitems; ++i) {
            uint32_t pfx = items[i].prefix, sfx = items[i].suffix, t;
            int have = 0;
            for (t = 0; t < Vn; ++t) if (verts[t] == pfx) { have = 1; break; }
            if (!have) {
                if (Vn >= ERGODEV_MAX_WORDS) { ok = -7; goto sparse_done; }
                verts[Vn++] = pfx;
            }
            have = 0;
            for (t = 0; t < Vn; ++t) if (verts[t] == sfx) { have = 1; break; }
            if (!have) {
                if (Vn >= ERGODEV_MAX_WORDS) { ok = -7; goto sparse_done; }
                verts[Vn++] = sfx;
            }
        }
        xi = (double *)calloc((size_t)Vn * Vn, sizeof(double));
        q = (double *)calloc(Vn, sizeof(double));
        P = (double *)calloc((size_t)Vn * Vn, sizeof(double));
        Lmat = (double *)calloc((size_t)Vn * Vn, sizeof(double));
        Vmat = (double *)calloc((size_t)Vn * Vn, sizeof(double));
        evals = (double *)calloc(Vn, sizeof(double));
        psi = (double *)calloc(Vn, sizeof(double));
        if (!xi || !q || !P || !Lmat || !Vmat || !evals || !psi) { ok = -6; goto sparse_done; }
        for (i = 0; i < nitems; ++i) {
            uint32_t vp = 0U, vs = 0U, t;
            for (t = 0; t < Vn; ++t) if (verts[t] == items[i].prefix) vp = t;
            for (t = 0; t < Vn; ++t) if (verts[t] == items[i].suffix) vs = t;
            xi[(size_t)vp * Vn + vs] += (double)items[i].count / (double)ntrans;
        }
        if (spectrum_from_flux(xi, q, P, Lmat, Vmat, evals, psi, Vn, q_normalize,
                               &lambda, &lambda_next, &gap, &kap) != 0) {
            ok = -5;
            goto sparse_done;
        }
        num_words = nitems + 1U;
        num_windows = window_len - ell + 1U;
        compute_layout(num_words, m, num_windows, 1, &L);
        if (ctx_bytes < L.total) { ok = -6; goto sparse_done; }
        memset(ctx, 0, ctx_bytes);
        apply_layout(ctx, &L);
        ctx->m = m;
        ctx->ell = ell;
        ctx->num_words = num_words;
        ctx->window_len = window_len;
        ctx->num_windows = num_windows;
        ctx->high_pow = high_pow;
        ctx->eta = eta;
        ctx->sparse = 1U;
        ctx->n_vertices = Vn;
        ctx->q_normalized = q_normalize ? 1 : 0;
        ctx->kappa = kap;
        ctx->lambda_next = lambda_next;
        ctx->spectral_gap = gap;
        for (i = 0; i < nitems; ++i) {
            uint32_t vp = 0U, vs = 0U, t;
            double jmp;
            for (t = 0; t < Vn; ++t) if (verts[t] == items[i].prefix) vp = t;
            for (t = 0; t < Vn; ++t) if (verts[t] == items[i].suffix) vs = t;
            ergodev_nu(ctx)[i] = (double)items[i].count / (double)ntrans;
            jmp = (psi[vp] - psi[vs]) * (psi[vp] - psi[vs]);
            ergodev_jump(ctx)[i] = jmp;
            if (jmp > max_jump) max_jump = jmp;
            delta += ergodev_nu(ctx)[i] * jmp;
            word_ids(ctx)[i] = items[i].id;
            if (hash_insert(ctx, items[i].id, (uint16_t)i) != 0) { ok = -7; goto sparse_done; }
        }
        ergodev_nu(ctx)[nitems] = 0.0;
        ergodev_jump(ctx)[nitems] = 0.0;
        word_ids(ctx)[nitems] = UINT32_MAX;
        if (!(delta > 0.0)) { ok = -5; goto sparse_done; }
        if ((psi_out && psi_cap < Vn) || (q_out && q_cap < Vn)) { ok = -6; goto sparse_done; }
        ctx->lambda_2 = delta;
        ctx->max_jump = (max_jump > 0.0) ? max_jump : 1e-15;
        diameter = 0.0;
        for (i = 0; i < num_words; ++i) {
            uint32_t z;
            for (z = 0; z < num_words; ++z) {
                double d;
                if (i == nitems || z == nitems) {
                    double step = 0.0;
                    uint32_t a, b;
                    for (a = 0; a < m; ++a) for (b = 0; b < m; ++b) {
                        double sd = symbol_delta(symbol_metric, m, a, b);
                        if (sd > step) step = sd;
                    }
                    if (fabs(theta - 1.0) < 1e-15) d = step * (double)ell;
                    else d = step * (1.0 - pow(theta, (double)ell)) / (1.0 - theta);
                    if (i == z) d = 0.0;
                } else {
                    d = word_distance(symbol_metric, m, ell, theta, items[i].id, items[z].id);
                }
                dist_out[(size_t)i * num_words + z] = d;
                if (d > diameter) diameter = d;
            }
        }
        ctx->diameter = (diameter > 0.0) ? diameter : 1e-15;
        if (psi_out) memcpy(psi_out, psi, (size_t)Vn * sizeof(double));
        if (q_out) memcpy(q_out, q, (size_t)Vn * sizeof(double));
        ok = 0;
    sparse_done:
        free(items); free(verts); free(xi); free(q); free(P); free(Lmat); free(Vmat); free(evals); free(psi);
        return ok;
    }
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
        fprintf(stderr, "FAIL [%s]\n", label);
        exit(EXIT_FAILURE);
    }
}

static ergodev_ctx_t *new_worked(const double *nu, const double *jump, const double *dist, uint32_t window, double lambda) {
    uint32_t nwin = window - 2U + 1U;
    size_t bytes = ergodev_ctx_bytes(9U, 3U, nwin, 0);
    ergodev_ctx_t *ctx = (ergodev_ctx_t *)malloc(bytes);
    int rc;
    require_true("malloc ctx", ctx != NULL);
    rc = ergodev_init(ctx, bytes, 3U, 2U, window, lambda, 1.0, nu, jump, dist);
    require_true("init worked", rc == 0);
    return ctx;
}

static int score_window(ergodev_ctx_t *ctx, const double *dist, double eps, uint32_t iters, double tol, ergodev_score_t *sc) {
    size_t nd = ergodev_scratch_doubles(ctx->num_words);
    double *store = (double *)malloc(nd * sizeof(double));
    ergodev_scratch_t scratch;
    int rc;
    require_true("malloc scratch", store != NULL);
    rc = ergodev_scratch_bind(&scratch, ctx->num_words, store, nd, dist, eps);
    require_true("scratch bind", rc == 0);
    rc = ergodev_score(ctx, &scratch, iters, tol, sc);
    free(store);
    return rc;
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
    double bad_jump[9];
    ergodev_ctx_t *ctx;
    ergodev_score_t sc;
    ergodev_bound_t bound;
    uint32_t w, t;
    int i, j;

    for (i = 0; i < 9; ++i) {
        int a = i / 3, b = i % 3;
        bad_jump[i] = 1.0;
        for (j = 0; j < 9; ++j) {
            int c = j / 3, d = j % 3;
            dist_mat[i * 9 + j] = fabs((double)(a - c)) + 0.5 * fabs((double)(b - d));
        }
    }

    ctx = new_worked(nu_2, jump_cost, dist_mat, 401U, 0.5);
    printf("%-28s | %-10s | %-10s | %-10s\n", "Stream Source (N'=400)", "W_1", "Defect", "E_S Score");
    printf("-------------------------------------------------------------------\n");

    for (t = 0; t < 401U; ++t) require_true("push alt", ergodev_push(ctx, (t & 1U) ? 2U : 0U) >= 0);
    require_true("alt conv", score_window(ctx, dist_mat, 0.02, 400U, 1e-5, &sc) == 0 && sc.converged);
    printf("%-28s | %-10.4f | %-10.4f | %-10.4f\n", "Alternation 0 <-> 2", sc.w1, sc.defect, sc.es);
    require_close("Alt W_1", sc.w1, 0.80, 1e-3);
    require_close("Alt Defect", sc.defect, 0.00, 1e-6);
    require_close("Alt E_S", sc.es, 0.80, 1e-3);
    require_true("alt not screened at 0.40", ergodev_screen(ctx, 0.40) == 0);
    require_true("alt upper exceeds score", ergodev_es_upper(ctx) + 1e-9 >= sc.es);

    ergodev_reset(ctx);
    for (t = 0; t < 401U; ++t) require_true("push const", ergodev_push(ctx, 1U) >= 0);
    require_true("const conv", score_window(ctx, dist_mat, 0.02, 400U, 1e-5, &sc) == 0);
    printf("%-28s | %-10.4f | %-10.4f | %-10.4f\n", "Constant State 1", sc.w1, sc.defect, sc.es);
    require_close("Const W_1", sc.w1, 0.75, 1e-3);
    require_close("Const Defect", sc.defect, 1.00, 1e-6);
    require_close("Const E_S", sc.es, 1.50, 1e-3);

    ergodev_reset(ctx);
    {
        const uint8_t cut_cycle[4] = {0, 0, 2, 2};
        for (t = 0; t < 401U; ++t) ergodev_push(ctx, cut_cycle[t & 3U]);
    }
    require_true("cut conv", score_window(ctx, dist_mat, 0.02, 400U, 1e-5, &sc) == 0);
    printf("%-28s | %-10.4f | %-10.4f | %-10.4f\n", "Balanced Cut {0,2}", sc.w1, sc.defect, sc.es);
    require_close("Cut E_S", sc.es, 0.75, 1e-3);
    require_close("Cut Defect", sc.defect, 0.00, 1e-6);

    ergodev_reset(ctx);
    {
        const uint8_t side_cycle[4] = {0, 0, 1, 1};
        for (t = 0; t < 401U; ++t) ergodev_push(ctx, side_cycle[t & 3U]);
    }
    require_true("side conv", score_window(ctx, dist_mat, 0.02, 400U, 1e-5, &sc) == 0);
    printf("%-28s | %-10.4f | %-10.4f | %-10.4f\n", "One-Sided Submanifold {0,1}", sc.w1, sc.defect, sc.es);
    require_close("Side W_1", sc.w1, 0.75, 1e-3);
    require_close("Side Defect", sc.defect, 0.50, 1e-6);
    require_close("Side E_S", sc.es, 1.125, 1e-3);

    /* Non-coupling is not a cost. */
    {
        size_t nd = ergodev_scratch_doubles(9U);
        double *store = (double *)malloc(nd * sizeof(double));
        ergodev_scratch_t scratch;
        require_true("scratch tiny", ergodev_scratch_bind(&scratch, 9U, store, nd, dist_mat, 1e-4) == 0);
        require_true("early stop", ergodev_score(ctx, &scratch, 1U, 1e-8, &sc) == -2);
        require_true("early not converged", sc.converged == 0 && !isfinite(sc.w1) && !isfinite(sc.es));
        free(store);
    }

    /* Short windows still match the population laws. */
    ergodev_reset(ctx);
    for (t = 0; t < 15U; ++t) ergodev_push(ctx, 1U);
    require_true("short const", score_window(ctx, dist_mat, 0.02, 400U, 1e-5, &sc) == 0);
    require_close("Short Const 15", sc.es, 1.50, 1e-3);
    ergodev_reset(ctx);
    for (t = 0; t < 11U; ++t) ergodev_push(ctx, (t & 1U) ? 2U : 0U);
    require_true("short alt", score_window(ctx, dist_mat, 0.02, 200U, 1e-5, &sc) == 0);
    require_close("Short Alt", sc.es, 0.80, 1e-3);

    /* Reject an inconsistent jump table, non-finite scalars, overflow, bad symbols. */
    {
        size_t bytes = ergodev_ctx_bytes(9U, 3U, 400U, 0);
        ergodev_ctx_t *bad = (ergodev_ctx_t *)malloc(bytes);
        require_true("reject jump mean", ergodev_init(bad, bytes, 3U, 2U, 401U, 0.5, 1.0, nu_2, bad_jump, dist_mat) < 0);
        require_true("reject m^ell wrap", ergodev_init(bad, bytes, 2U, 32U, 401U, 0.5, 1.0, nu_2, jump_cost, dist_mat) < 0);
        require_true("reject NaN lambda", ergodev_init(bad, bytes, 3U, 2U, 401U, NAN, 1.0, nu_2, jump_cost, dist_mat) < 0);
        require_true("reject Inf eta", ergodev_init(bad, bytes, 3U, 2U, 401U, 0.5, INFINITY, nu_2, jump_cost, dist_mat) < 0);
        require_true("reject tiny buffer", ergodev_init(bad, 8U, 3U, 2U, 401U, 0.5, 1.0, nu_2, jump_cost, dist_mat) < 0);
        free(bad);
    }
    require_true("reject sym 255", ergodev_push(ctx, 255U) == -1);
    require_true("reject sym 16 fast", ergodev_push_fast(ctx, 16U) == -1);

    /* ell = 3 dense ring matches an independent histogram. */
    {
        enum { M = 3, ELL = 3, NW = 27, WIN = 64, STREAM = 800 };
        double nu[NW], jump[NW], dist[NW * NW];
        uint32_t nwin = WIN - ELL + 1U;
        size_t bytes = ergodev_ctx_bytes(NW, M, nwin, 0);
        ergodev_ctx_t *c3 = (ergodev_ctx_t *)malloc(bytes);
        uint8_t stream[STREAM];
        uint32_t rng = 1U, brute[NW], sum_b = 0U, mismatches = 0U;
        const uint32_t *counts;
        require_true("ell3 malloc", c3 != NULL);
        for (w = 0; w < NW; ++w) {
            nu[w] = 1.0 / (double)NW;
            jump[w] = 0.1;
        }
        for (i = 0; i < NW * NW; ++i) dist[i] = (i % (NW + 1) == 0) ? 0.0 : 1.0;
        require_true("ell3 init", ergodev_init(c3, bytes, M, ELL, WIN, 0.1, 1.0, nu, jump, dist) == 0);
        for (i = 0; i < STREAM; ++i) {
            rng = rng * 1664525U + 1013904223U;
            stream[i] = (uint8_t)((rng >> 16) % M);
            ergodev_push_fast(c3, stream[i]);
        }
        memset(brute, 0, sizeof(brute));
        for (i = STREAM - WIN; i <= STREAM - ELL; ++i) {
            uint32_t id = 0U, k;
            for (k = 0; k < ELL; ++k) id = id * M + stream[i + k];
            brute[id]++;
        }
        counts = ergodev_counts_c(c3);
        for (w = 0; w < NW; ++w) {
            sum_b += brute[w];
            if (counts[w] != brute[w]) mismatches++;
        }
        require_true("ell3 histogram", mismatches == 0U && sum_b == c3->active_words);
        free(c3);
    }

    /* Theorem 8.4 / Remark 8.3 / Freedman on the worked chain. */
    require_true("bounds 400", ergodev_bounds(3.0, 9U, 2U, 400U, 0.6, 1.0, 0.5, 2.0, 0.05, 0.4, &bound) == 0);
    require_close("B_400", bound.B_n, 0.5518254055364693, 1e-9);
    require_close("azuma 400", bound.azuma_threshold, 4.784473112848603, 1e-9);
    require_close("freedman 400", bound.freedman_threshold, 3.8812327488306093, 1e-9);
    require_true("azuma does not separate at 400", bound.azuma_separates == 0);
    require_true("bounds 56977", ergodev_bounds(3.0, 9U, 2U, 56977U, 0.6, 1.0, 0.5, 2.0, 0.05, 0.4, &bound) == 0);
    require_true("azuma separates at 56977", bound.azuma_separates == 1);
    require_true("bounds 16400", ergodev_bounds(3.0, 9U, 2U, 16400U, 0.6, 1.0, 0.5, 2.0, 0.05, 0.4, &bound) == 0);
    require_close("joint 16400", bound.joint_threshold, 0.5649598148315592, 1e-9);
    require_true("joint does not separate at 16400", bound.joint_separates == 0);
    require_true("bounds 27971", ergodev_bounds(3.0, 9U, 2U, 27971U, 0.6, 1.0, 0.5, 2.0, 0.05, 0.4, &bound) == 0);
    require_true("joint separates at 27971", bound.joint_separates == 1);
    require_true("bounds freedman alt", ergodev_bounds(3.0, 9U, 2U, 9103U, 0.6, 1.0, 0.5, 2.0, 0.05, 0.8, &bound) == 0);
    require_true("freedman separates alternation at 9103", bound.freedman_separates == 1);
    {
        double sample[5] = {5.0, 1.0, 4.0, 2.0, 3.0};
        require_close("quantile 0.8", ergodev_empirical_quantile(sample, 5U, 0.8), 4.0, 1e-15);
    }

    /* Corpus fit, unnormalized and q-normalized, on the Section 9 chain. */
    {
        enum { CLEN = 40000 };
        uint8_t *corpus = (uint8_t *)malloc(CLEN);
        uint32_t rng = 7U;
        uint8_t state = 1U;
        double psi[8], qv[8], mean = 0.0, var = 0.0;
        ergodev_plan_t plan;
        ergodev_ctx_t *fit;
        double *dist;
        uint32_t k;
        require_true("corpus", corpus != NULL);
        for (k = 0; k < CLEN; ++k) {
            uint32_t r;
            rng = rng * 1664525U + 1013904223U;
            r = (rng >> 16) % 10U;
            if (state == 1U) state = (r < 4U) ? 0U : ((r < 6U) ? 1U : 2U);
            else state = (r < 1U) ? 0U : ((r < 9U) ? 1U : 2U);
            corpus[k] = state;
        }
        require_true("plan dense", ergodev_plan(3U, 2U, 401U, corpus, CLEN, &plan) == 0 && plan.sparse == 0);
        fit = (ergodev_ctx_t *)malloc(plan.ctx_bytes);
        dist = (double *)malloc((size_t)plan.num_words * plan.num_words * sizeof(double));
        require_true("fit unnorm", ergodev_fit(fit, plan.ctx_bytes, 3U, 2U, 401U, 1.0, 0.5, NULL,
                                                corpus, CLEN, 0, dist, plan.num_words, psi, 8U, qv, 8U) == 0);
        require_close("fitted lambda", fit->lambda_2, 0.5, 0.03);
        require_true("gap stored", fit->spectral_gap > 0.0 && fit->kappa > 0.0 && fit->kappa < 1.0);
        require_true("fit qnorm", ergodev_fit(fit, plan.ctx_bytes, 3U, 2U, 401U, 1.0, 0.5, NULL,
                                               corpus, CLEN, 1, dist, plan.num_words, psi, 8U, qv, 8U) == 0);
        require_true("qnorm flag", fit->q_normalized == 1);
        for (k = 0; k < fit->n_vertices; ++k) {
            mean += qv[k] * psi[k];
            var += qv[k] * psi[k] * psi[k];
        }
        require_close("q-mean", mean, 0.0, 1e-8);
        require_close("q-var", var, 1.0, 1e-8);
        require_true("max jump is the Azuma cap", fit->max_jump > 0.0);

        /* Pair-preserving cycle: period arm fires, transport arm stays quiet. */
        {
            const char *period = "0221212121212121212011111010101010101010";
            uint8_t syms[400];
            ergodev_alarm_t alarm;
            uint32_t p;
            int nsym;
            ergodev_reset(fit);
            for (k = 0; k < 400U; ++k) {
                syms[k] = (uint8_t)(period[k % 40U] - '0');
                ergodev_push(fit, syms[k]);
            }
            nsym = ergodev_window_symbols(fit, syms, 400U);
            require_true("window symbols", nsym == 400);
            p = ergodev_min_period(syms, (uint32_t)nsym);
            require_true("cycle period 40", p == 40U);
            require_true("cycle screened by Lemma 8.1", ergodev_screen(fit, 0.05) == 1);
            require_true("cycle score", score_window(fit, dist, 0.05, 2500U, 1e-3, &sc) == 0);
            require_true("cycle transport under 0.08", sc.es < 0.08);
            require_true("union", ergodev_alarm_union(sc.es, 1, 0.20, p, (uint32_t)nsym, 64U, &alarm) == 0);
            require_true("period arm", alarm.alarm == 1 && alarm.alarm_period == 1 && alarm.alarm_transport == 0);
            require_true("no alarm on aperiodic quiet window",
                         ergodev_alarm_union(0.01, 1, 0.20, 400U, 400U, 64U, &alarm) == 0 && alarm.alarm == 0);
            require_true("transport arm",
                         ergodev_alarm_union(1.0, 1, 0.20, 400U, 400U, 64U, &alarm) == 0 &&
                         alarm.alarm_transport == 1 && alarm.alarm == 1);
        }

        /* Per-thread, two clocks. Scoring is not invoked by the push. */
        {
            ergodev_thread_slot_t slot;
            ergodev_ctx_t *quantum;
            ergodev_plan_t qp;
            double qnu[4] = {0.7, 0.1, 0.1, 0.1};
            double qjump[4] = {0.0, 0.2, 0.2, 0.2};
            double qdist[16];
            size_t qbytes;
            int a, b;
            for (a = 0; a < 4; ++a) for (b = 0; b < 4; ++b) qdist[a * 4 + b] = fabs((double)(a - b));
            /* Jump mean is 0.1*0.2*3 = 0.06. */
            qbytes = ergodev_ctx_bytes(4U, 2U, 64U - 2U + 1U, 0);
            quantum = (ergodev_ctx_t *)malloc(qbytes);
            require_true("quantum init", ergodev_init(quantum, qbytes, 2U, 2U, 64U, 0.06, 1.0, qnu, qjump, qdist) == 0);
            slot.tid = 42U;
            slot.per_event = fit;
            slot.per_quantum = quantum;
            ergodev_reset(fit);
            for (k = 0; k < 32U; ++k) {
                require_true("thread event", ergodev_thread_push_event(&slot, (uint8_t)(k % 3U)) >= 0);
                require_true("thread quantum", ergodev_thread_push_quantum(&slot, (k % 5U == 0U) ? 0U : 1U) >= 0);
            }
            require_true("event window advanced", fit->active_words > 0U);
            require_true("quantum window advanced", quantum->active_words > 0U);
            free(quantum);
            (void)qp;
        }
        free(dist);
        free(fit);
        free(corpus);
    }

    /* Sparse ell = 3 on a 6-class nearest-neighbor taxonomy. */
    {
        enum { M = 6, ELL = 3, CLEN = 8000, WIN = 128 };
        uint8_t *corpus = (uint8_t *)malloc(CLEN);
        uint32_t rng = 3U, k;
        uint8_t state = 2U;
        ergodev_plan_t plan;
        ergodev_ctx_t *fit;
        double *dist;
        require_true("sparse corpus", corpus != NULL);
        for (k = 0; k < CLEN; ++k) {
            int delta;
            rng = rng * 1664525U + 1013904223U;
            delta = (int)((rng >> 16) % 3U) - 1;
            if ((int)state + delta < 0) delta = 0;
            if ((int)state + delta >= M) delta = 0;
            state = (uint8_t)((int)state + delta);
            corpus[k] = state;
        }
        require_true("plan sparse", ergodev_plan(M, ELL, WIN, corpus, CLEN, &plan) == 0);
        require_true("plan is sparse", plan.sparse == 1 && plan.num_words < ERGODEV_MAX_WORDS);
        require_true("dense would not fit", M * M * M > ERGODEV_MAX_WORDS);
        fit = (ergodev_ctx_t *)malloc(plan.ctx_bytes);
        dist = (double *)calloc((size_t)plan.num_words * plan.num_words, sizeof(double));
        require_true("sparse fit", ergodev_fit(fit, plan.ctx_bytes, M, ELL, WIN, 1.0, 0.5, NULL,
                                                corpus, CLEN, 1, dist, plan.num_words, NULL, 0U, NULL, 0U) == 0);
        require_true("sparse flag", fit->sparse == 1U && fit->lambda_2 > 0.0 && fit->spectral_gap >= 0.0);
        for (k = 0; k < WIN; ++k) require_true("sparse push", ergodev_push(fit, corpus[k]) >= 0);
        require_true("sparse score", score_window(fit, dist, 0.05, 4000U, 1e-2, &sc) == 0 && isfinite(sc.es));
        free(dist);
        free(fit);
        free(corpus);
    }

    free(ctx);
    printf("-------------------------------------------------------------------\n");
    printf("All population targets, bounds, fit, sparse, and alarm checks passed.\n");
    return 0;
}
#endif

