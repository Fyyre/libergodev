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
 * to build the word histogram and Fiedler energy before scoring. */
static double naive_rescan_and_score(
    ergodev_ctx_t *scratch_ctx,
    const uint8_t *stream,
    size_t pos,
    uint32_t win_len
) {
    memset(scratch_ctx->counts, 0, sizeof(scratch_ctx->counts));
    uint32_t num_pairs = win_len - 1U;
    const uint8_t *win = stream + pos - win_len + 1U;

    double e_sum = 0.0;
    for (uint32_t i = 0; i < num_pairs; ++i) {
        uint32_t w = (uint32_t)win[i] * 3U + (uint32_t)win[i + 1U];
        scratch_ctx->counts[w]++;
        e_sum += scratch_ctx->jump_cost[w];
    }
    scratch_ctx->sym_filled = win_len;
    scratch_ctx->energy_sum = e_sum;

    double def = ergodev_defect(scratch_ctx);
    if (def > 0.25) {
        return ergodev_score_sinkhorn(scratch_ctx, 15, NULL, NULL);
    }
    return def;
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

    uint8_t *stream = (uint8_t *)malloc(NUM_EVENTS);
    if (!stream) return 1;

    /* Generate stationary baseline stream + periodic one-sided anomaly bursts */
    uint32_t rng = 0xC0FFEEU;
    uint8_t state = 1U;
    for (size_t i = 0; i < NUM_EVENTS; ++i) {
        rng ^= rng << 13;
        rng ^= rng >> 17;
        rng ^= rng << 5;
        if ((i % 1000000ULL) > 995000ULL) {
            /* Anomaly burst: confined to one side {0, 1} */
            stream[i] = (uint8_t)(rng & 1U);
        } else {
            /* Approximate baseline transitions */
            uint32_t r = rng % 10U;
            if (state == 1U) {
                state = (r < 4U) ? 0U : ((r < 8U) ? 2U : 1U);
            } else {
                state = (r < 8U) ? 1U : ((r == 8U) ? 0U : 2U);
            }
            stream[i] = state;
        }
    }

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

    const uint32_t windows[] = { 128U, 512U, 2048U, 8192U, 16384U };
    const size_t num_w = sizeof(windows) / sizeof(windows[0]);

    printf("=======================================================================================\n");
    printf("  LIBERGODEV STREAMING BENCHMARK: 10,000,000 Events (Query Stride = %u ticks)\n", EVAL_STRIDE);
    printf("=======================================================================================\n\n");
    printf("%-10s | %-22s | %-12s | %-10s | %s\n",
           "Window N", "Method", "Throughput", "Latency", "Relative Speed");
    printf("---------------------------------------------------------------------------------------\n");

    volatile double sink = 0.0;
    ergodev_ctx_t ctx;
    ergodev_ctx_t scratch_ctx;

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
                if (def > 0.25) {
                    acc_ergo += ergodev_score_sinkhorn(&ctx, 15, NULL, NULL);
                } else {
                    acc_ergo += def;
                }
            }
        }
        uint64_t dt_ergo = now_ns() - t1;
        sink += acc_ergo;

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