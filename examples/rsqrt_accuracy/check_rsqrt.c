#define _ISOC11_SOURCE
#include <immintrin.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static inline __m512d newton_iter(__m512d a_half, __m512d y) {
    __m512d y2  = _mm512_mul_pd(y, y);
    __m512d c15 = _mm512_set1_pd(1.5);
    __m512d fac = _mm512_fnmadd_pd(a_half, y2, c15);
    return _mm512_mul_pd(y, fac);
}

int main(void) {
    const double log10_lo = -10.0;
    const double log10_hi = 10.0;
    const int per_dec = 200;
    const int N = (int)((log10_hi - log10_lo) * per_dec) + 1;

    double *rs = (double*)aligned_alloc(64, sizeof(double) * ((N + 7) & ~7));
    double *e0 = (double*)aligned_alloc(64, sizeof(double) * N);
    double *e1 = (double*)aligned_alloc(64, sizeof(double) * N);
    double *e2 = (double*)aligned_alloc(64, sizeof(double) * N);

    for (int i = 0; i < N; i++) {
        rs[i] = pow(10.0, log10_lo + (double)i * (log10_hi - log10_lo) / (double)(N - 1));
    }
    for (int i = N; i < ((N + 7) & ~7); i++) rs[i] = 1.0;

    const __m512d half = _mm512_set1_pd(0.5);

    for (int i = 0; i < N; i += 8) {
        __m512d a = _mm512_loadu_pd(&rs[i]);
        __m512d y0 = _mm512_rsqrt14_pd(a);
        __m512d a_half = _mm512_mul_pd(half, a);
        __m512d y1 = newton_iter(a_half, y0);
        __m512d y2 = newton_iter(a_half, y1);

        double y0a[8], y1a[8], y2a[8];
        _mm512_storeu_pd(y0a, y0);
        _mm512_storeu_pd(y1a, y1);
        _mm512_storeu_pd(y2a, y2);

        for (int k = 0; k < 8 && (i + k) < N; k++) {
            double a_s = rs[i + k];
            double true_v = 1.0 / sqrt(a_s);
            e0[i + k] = fabs(y0a[k] / true_v - 1.0);
            e1[i + k] = fabs(y1a[k] / true_v - 1.0);
            e2[i + k] = fabs(y2a[k] / true_v - 1.0);
        }
    }

    printf("# a   rel_err_y0   rel_err_y1   rel_err_y2\n");
    for (int i = 0; i < N; i++) {
        printf("%.17e %.17e %.17e %.17e\n", rs[i], e0[i], e1[i], e2[i]);
    }

    free(rs); free(e0); free(e1); free(e2);
    return 0;
}
