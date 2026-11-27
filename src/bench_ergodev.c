#define _POSIX_C_SOURCE 199309L
#include "ergodev.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

/*
 * Timed path: integer ring push, plus the Lemma 8.1 upper bound every
 * 16 symbols. That bound is one pass over the word counts. Certified
 * Sinkhorn stays off this loop; the alternation check below still
 * requires a converged coupling.
 */

static uint64_t now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

static double naive_upper(
    ergodev_ctx_t *scratch_ctx,
    const uint8_t *stream,
    size_t pos,
    uint32_t win_len
) {
    uint32_t num_pairs = win_len - 1U;
    const uint8_t *win = stream + pos - win_len + 1U;
    uint32_t i;
    memset(ergodev_counts(scratch_ctx), 0, (size_t)scratch_ctx->num_words * sizeof(uint32_t));
    for (i = 0; i < num_pairs; ++i) {
        uint32_t w = (uint32_t)win[i] * 3U + (uint32_t)win[i + 1U];
        ergodev_counts(scratch_ctx)[w]++;
    }
    scratch_ctx->sym_filled = win_len;
    scratch_ctx->active_words = num_pairs;
    return ergodev_es_upper(scratch_ctx);
}

static void print_bar(double mops, double max_mops) {
    int width = 28;
    int filled = (int)((mops / max_mops) * width);
    int i;
    if (filled < 1 && mops > 0.0) filled = 1;
    if (filled > width) filled = width;
    printf("[");
    for (i = 0; i < width; ++i) putchar(i < filled ? '#' : ' ');
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
    ergodev_ctx_t *ctx;
    ergodev_scratch_t scratch;
    double *store;
    size_t bytes, nd;
    uint8_t *stream;
    uint32_t rng = 0xC0FFEEU;
    uint8_t state = 1U;
    size_t i;
    int a, b;
    const uint32_t windows[] = { 128U, 512U, 2048U, 8192U, 16384U };
    const size_t num_w = sizeof(windows) / sizeof(windows[0]);
    volatile double sink = 0.0;

    for (a = 0; a < 9; ++a) {
        for (b = 0; b < 9; ++b) {
            dist_mat[a * 9 + b] = fabs((double)(a / 3 - b / 3)) + 0.5 * fabs((double)(a % 3 - b % 3));
        }
    }

    bytes = ergodev_ctx_bytes(9U, 3U, 129U - 2U + 1U, 0);
    ctx = (ergodev_ctx_t *)malloc(bytes);
    nd = ergodev_scratch_doubles(9U);
    store = (double *)malloc(nd * sizeof(double));
    if (!ctx || !store) return EXIT_FAILURE;
    if (ergodev_init(ctx, bytes, 3U, 2U, 129U, 0.5, 1.0, nu_2, jump_cost, dist_mat) != 0) return EXIT_FAILURE;
    if (ergodev_scratch_bind(&scratch, 9U, store, nd, dist_mat, 0.02) != 0) return EXIT_FAILURE;

    {
        ergodev_score_t sc;
        uint8_t alt_buf[129];
        int rc;
        for (i = 0; i < 129U; ++i) {
            alt_buf[i] = (uint8_t)((i & 1U) ? 2U : 0U);
            ergodev_push_fast(ctx, alt_buf[i]);
        }
        rc = ergodev_score(ctx, &scratch, 400U, 1e-5, &sc);
        if (rc != 0 || fabs(sc.es - 0.80) > 1e-3) {
            fprintf(stderr, "FAIL: alternation certified score rc=%d es=%.4f resid=%g\n",
                    rc, sc.es, sc.residual);
            return EXIT_FAILURE;
        }
        if (!(ergodev_es_upper(ctx) + 1e-9 >= sc.es)) {
            fprintf(stderr, "FAIL: Lemma 8.1 upper %.4f is below E_S %.4f\n",
                    ergodev_es_upper(ctx), sc.es);
            return EXIT_FAILURE;
        }
    }

    stream = (uint8_t *)malloc(NUM_EVENTS);
    if (!stream) return EXIT_FAILURE;
    for (i = 0; i < NUM_EVENTS; ++i) {
        uint32_t phase, r;
        rng ^= rng << 13;
        rng ^= rng >> 17;
        rng ^= rng << 5;
        phase = (uint32_t)(i % 1000000ULL);
        if (phase > 997500ULL) stream[i] = (uint8_t)((i & 1ULL) ? 2U : 0U);
        else if (phase > 995000ULL) stream[i] = (uint8_t)(rng & 1U);
        else {
            r = rng % 10U;
            if (state == 1U) state = (r < 4U) ? 0U : ((r < 8U) ? 2U : 1U);
            else state = (r < 8U) ? 1U : ((r == 8U) ? 0U : 2U);
            stream[i] = state;
        }
    }

    printf("=======================================================================================\n");
    printf("  LIBERGODEV STREAMING BENCHMARK: %zu Events (Query Stride = %u ticks)\n",
           NUM_EVENTS, EVAL_STRIDE);
    printf("  Timed path: O(1) ring + Lemma 8.1 upper bound. Certified Sinkhorn is off this loop.\n");
    printf("=======================================================================================\n\n");
    printf("%-10s | %-22s | %-12s | %-10s | %s\n",
           "Window N", "Method", "Throughput", "Latency", "Relative Speed");
    printf("---------------------------------------------------------------------------------------\n");

    for (i = 0; i < num_w; ++i) {
        uint32_t N = windows[i];
        uint32_t nwin = N - 2U + 1U;
        size_t nb = ergodev_ctx_bytes(9U, 3U, nwin, 0);
        ergodev_ctx_t *c_naive = (ergodev_ctx_t *)malloc(nb);
        ergodev_ctx_t *c_fast = (ergodev_ctx_t *)malloc(nb);
        double acc_naive = 0.0, acc_ergo = 0.0;
        uint64_t t0, dt_naive, dt_ergo;
        double mops_naive, ns_naive, mops_ergo, ns_ergo;
        size_t ev;
        if (!c_naive || !c_fast) return EXIT_FAILURE;
        if (ergodev_init(c_naive, nb, 3U, 2U, N, 0.5, 1.0, nu_2, jump_cost, dist_mat) != 0) return EXIT_FAILURE;
        if (ergodev_init(c_fast, nb, 3U, 2U, N, 0.5, 1.0, nu_2, jump_cost, dist_mat) != 0) return EXIT_FAILURE;

        t0 = now_ns();
        for (ev = N - 1U; ev < NUM_EVENTS; ev += EVAL_STRIDE) {
            acc_naive += naive_upper(c_naive, stream, ev, N);
        }
        dt_naive = now_ns() - t0;

        t0 = now_ns();
        for (ev = 0; ev < NUM_EVENTS; ++ev) {
            int ready = ergodev_push_fast(c_fast, stream[ev]);
            if (ready && (((ev - (N - 1U)) & (EVAL_STRIDE - 1U)) == 0U)) {
                acc_ergo += ergodev_es_upper(c_fast);
            }
        }
        dt_ergo = now_ns() - t0;
        sink += acc_naive + acc_ergo;
        if (fabs(acc_naive - acc_ergo) > 1e-6 * (1.0 + fabs(acc_naive))) {
            fprintf(stderr, "FAIL [N=%u]: naive %.6f vs streaming %.6f\n", N, acc_naive, acc_ergo);
            return EXIT_FAILURE;
        }

        mops_naive = (double)NUM_EVENTS / ((double)dt_naive * 1e-3);
        ns_naive = (double)dt_naive / (double)NUM_EVENTS;
        mops_ergo = (double)NUM_EVENTS / ((double)dt_ergo * 1e-3);
        ns_ergo = (double)dt_ergo / (double)NUM_EVENTS;
        printf("N = %-6u | %-22s | %7.1f M/s | %6.2f ns | ", N, "Naive Window Rescan", mops_naive, ns_naive);
        print_bar(mops_naive, 800.0);
        printf("\n");
        printf("N = %-6u | %-22s | %7.1f M/s | %6.2f ns | ", N, "libergodev O(1) Ring", mops_ergo, ns_ergo);
        print_bar(mops_ergo, 800.0);
        printf(" (%.1fx faster)\n", mops_ergo / mops_naive);
        printf("---------------------------------------------------------------------------------------\n");
        free(c_naive);
        free(c_fast);
    }

    if (sink == 123456789.0) printf("%f\n", sink);
    free(stream);
    free(store);
    free(ctx);
    return 0;
}
