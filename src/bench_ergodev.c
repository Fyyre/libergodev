#define _POSIX_C_SOURCE 199309L
#include "ergodev.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

static uint64_t now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

/* Naive stateless approach: Rescans the N-symbol window from scratch at each query
 * to build the word histogram before evaluating the pre-filter and Sinkhorn W_1. */
static double naive_rescan_and_score(
    ergodev_ctx_t *scratch_ctx,
    const uint8_t *stream,
    size_t pos,
    uint32_t win_len
) {
    memset(scratch_ctx->counts, 0, sizeof(scratch_ctx->counts));
    uint32_t num_pairs = win_len - 1U;
    const uint8_t *win = stream + pos - win_len + 1U;

    for (uint32_t i = 0; i < num_pairs; ++i) {
        uint32_t w = (uint32_t)win[i] * 3U + (uint32_t)win[i + 1U];
        scratch_ctx->counts[w]++;
    }
    scratch_ctx->sym_filled   = win_len;
    scratch_ctx->active_words = num_pairs;

    /* By Lemma 8.1, W_1 <= (D/2) * ||\hat\mu - \nu||_1.
     * Checking both Fiedler defect AND L1 drift guarantees high-energy anomalies
     * like Alternation (defect = 0, L1 = 1.90, W_1 = 0.80) always trigger Sinkhorn. */
    double def = ergodev_defect(scratch_ctx);
    double l1  = ergodev_l1_drift(scratch_ctx);
    if (def > 0.25 || l1 > 0.35) {
        return ergodev_score_sinkhorn(scratch_ctx, 15, NULL, NULL);
    }
    return 0.5 * l1 * (1.0 + scratch_ctx->eta * def);
}

static void print_bar(double mops, double max_mops) {
    int width = 28;
    int filled = (int)((mops / max_mops) * width);
    if (filled < 1 && mops > 0.0) filled = 1;
    if (filled > width) filled = width;
    printf("[");
    for (int i = 0; i < width; ++i) {
        putchar(i < filled ? '#' : ' ');
    }
    printf("]");
}

int main(void) {
    const size_t NUM_EVENTS = 10000000ULL;
    const uint32_t EVAL_STRIDE = 16U;

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
        for (int w2 = 0; w2 < 9; ++w2) {
            dist_mat[w1 * 9 + w2] = fabs((double)(w1/3 - w2/3)) + 0.5 * fabs((double)(w1%3 - w2%3));
        }
    }

    ergodev_ctx_t ctx;
    ergodev_ctx_t scratch_ctx;

    /* Sanity check: Verify the pre-filter does not mask zero-defect Alternation (0.80) */
    uint8_t alt_buf[129];
    ergodev_init(&ctx, 3, 2, 129, 0.5, 1.0, nu_2, jump_cost, dist_mat, 0.005);
    ergodev_init(&scratch_ctx, 3, 2, 129, 0.5, 1.0, nu_2, jump_cost, dist_mat, 0.005);
    for (uint32_t t = 0; t < 129; ++t) {
        alt_buf[t] = (t & 1U) ? 2U : 0U;
        ergodev_push_fast(&ctx, alt_buf[t]);
    }
    double alt_naive = naive_rescan_and_score(&scratch_ctx, alt_buf, 128, 129);
    double alt_def   = ergodev_defect(&ctx);
    double alt_l1    = ergodev_l1_drift(&ctx);
    double alt_ergo  = (alt_def > 0.25 || alt_l1 > 0.35)
                     ? ergodev_score_sinkhorn(&ctx, 200, NULL, NULL)
                     : 0.5 * alt_l1 * (1.0 + ctx.eta * alt_def);
    if (fabs(alt_naive - 0.80) > 0.02 || fabs(alt_ergo - 0.80) > 0.01) {
        fprintf(stderr, "FAIL: Alternation pre-filter check failed (naive=%.4f, ergo=%.4f)\n",
                alt_naive, alt_ergo);
        return EXIT_FAILURE;
    }

    uint8_t *stream = (uint8_t *)malloc(NUM_EVENTS);
    if (!stream) return EXIT_FAILURE;

    /* Generate stationary baseline stream + periodic structured anomaly bursts */
    uint32_t rng = 0xC0FFEEU;
    uint8_t state = 1U;
    for (size_t i = 0; i < NUM_EVENTS; ++i) {
        rng ^= rng << 13;
        rng ^= rng >> 17;
        rng ^= rng << 5;
        uint32_t phase = (uint32_t)(i % 1000000ULL);
        if (phase > 997500ULL) {
            /* Anomaly burst A: Rare-cut alternation 0 <-> 2 (defect = 0, W_1 = 0.80) */
            stream[i] = (uint8_t)((i & 1ULL) ? 2U : 0U);
        } else if (phase > 995000ULL) {
            /* Anomaly burst B: One-sided submanifold {0, 1} (defect = 0.5, E_S = 1.125) */
            stream[i] = (uint8_t)(rng & 1U);
        } else {
            /* Baseline Markov chain transitions */
            uint32_t r = rng % 10U;
            if (state == 1U) {
                state = (r < 4U) ? 0U : ((r < 8U) ? 2U : 1U);
            } else {
                state = (r < 8U) ? 1U : ((r == 8U) ? 0U : 2U);
            }
            stream[i] = state;
        }
    }

    const uint32_t windows[] = { 128U, 512U, 2048U, 8192U, 16384U };
    const size_t num_w = sizeof(windows) / sizeof(windows[0]);

    printf("=======================================================================================\n");
    printf("  LIBERGODEV STREAMING BENCHMARK: 10,000,000 Events (Query Stride = %u ticks)\n", EVAL_STRIDE);
    printf("=======================================================================================\n\n");
    printf("%-10s | %-22s | %-12s | %-10s | %s\n",
           "Window N", "Method", "Throughput", "Latency", "Relative Speed");
    printf("---------------------------------------------------------------------------------------\n");

    volatile double sink = 0.0;

    for (size_t w_idx = 0; w_idx < num_w; ++w_idx) {
        uint32_t N = windows[w_idx];
        if (N > ERGODEV_MAX_WINDOW) continue;

        /* 1. Naive O(N) Stateless Window Rescan */
        ergodev_init(&scratch_ctx, 3, 2, N, 0.5, 1.0, nu_2, jump_cost, dist_mat, 0.05);
        uint64_t t0 = now_ns();
        double acc_naive = 0.0;
        for (size_t i = N - 1U; i < NUM_EVENTS; i += EVAL_STRIDE) {
            acc_naive += naive_rescan_and_score(&scratch_ctx, stream, i, N);
        }
        uint64_t dt_naive = now_ns() - t0;
        sink += acc_naive;

        double mops_naive = (double)NUM_EVENTS / ((double)dt_naive * 1e-3);
        double ns_naive   = (double)dt_naive / (double)NUM_EVENTS;

        /* 2. libergodev O(1) Ring-Buffer Streaming */
        ergodev_init(&ctx, 3, 2, N, 0.5, 1.0, nu_2, jump_cost, dist_mat, 0.05);
        uint64_t t1 = now_ns();
        double acc_ergo = 0.0;
        for (size_t i = 0; i < NUM_EVENTS; ++i) {
            int ready = ergodev_push_fast(&ctx, stream[i]);
            if (ready && (((i - (N - 1U)) & (EVAL_STRIDE - 1U)) == 0U)) {
                double def = ergodev_defect(&ctx);
                double l1  = ergodev_l1_drift(&ctx);
                if (def > 0.25 || l1 > 0.35) {
                    acc_ergo += ergodev_score_sinkhorn(&ctx, 15, NULL, NULL);
                } else {
                    acc_ergo += 0.5 * l1 * (1.0 + ctx.eta * def);
                }
            }
        }
        uint64_t dt_ergo = now_ns() - t1;
        sink += acc_ergo;

        /* Verify both methods produced identical accumulated scores across all 10M events */
        if (fabs(acc_naive - acc_ergo) > 1e-6) {
            fprintf(stderr, "FAIL [N=%u]: Score mismatch between naive (%.6f) and streaming (%.6f)\n",
                    N, acc_naive, acc_ergo);
            free(stream);
            return EXIT_FAILURE;
        }

        double mops_ergo = (double)NUM_EVENTS / ((double)dt_ergo * 1e-3);
        double ns_ergo   = (double)dt_ergo / (double)NUM_EVENTS;

        printf("N = %-6u | %-22s | %7.1f M/s | %6.2f ns | ",
               N, "Naive Window Rescan", mops_naive, ns_naive);
        print_bar(mops_naive, 800.0);
        printf("\n");

        printf("N = %-6u | %-22s | %7.1f M/s | %6.2f ns | ",
               N, "libergodev O(1) Ring", mops_ergo, ns_ergo);
        print_bar(mops_ergo, 800.0);
        printf(" (%.1fx faster)\n", mops_ergo / mops_naive);
        printf("---------------------------------------------------------------------------------------\n");
    }

    if (sink == 123456789.0) printf("%f\n", sink);
    free(stream);
    return 0;
}