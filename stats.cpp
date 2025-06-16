/**
 * @file stats.cpp
 * @brief Реализация статистических функций и тестов NIST для анализа качества генераторов случайных чисел.
 */

#include "stats.h"
#include <cmath>
#include <cstdlib>

/**
 * @brief Вычисляет среднее значение массива.
 */
double mean(const uint32_t *data, int n) {
    uint64_t sum = 0;
    for (int i = 0; i < n; ++i) sum += data[i];
    return (double)sum / n;
}

/**
 * @brief Вычисляет стандартное отклонение массива.
 */
double stdev(const uint32_t *data, int n, double m) {
    double acc = 0.0;
    for (int i = 0; i < n; ++i)
        acc += (data[i] - m) * (data[i] - m);
    return sqrt(acc / n);
}

/**
 * @brief Возвращает коэффициент вариации.
 */
double coeff_var(double m, double sd) {
    return m == 0.0 ? 0.0 : sd / m;
}

/**
 * @brief Вычисляет критерий хи-квадрат по частотам попадания в корзины.
 */
double chi_squared(const uint32_t *data, int n, int bins, unsigned long long int max_val) {
    unsigned long int *freq = (unsigned long int *)calloc(bins, sizeof(unsigned long int));
    for (int i = 0; i < n; ++i) {
        int idx = ((uint64_t)data[i] * bins) / max_val;
        if (idx >= bins) idx = bins - 1;
        freq[idx]++;
    }

    double expected = (double)n / bins;
    double chi2 = 0.0;
    for (int i = 0; i < bins; ++i) {
        double diff = freq[i] - expected;
        chi2 += diff * diff / expected;
    }

    free(freq);
    return chi2;
}

/**
 * @brief Подсчитывает количество единичных битов в 32-битном числе.
 */
int popcount32(uint32_t x) {
    int count = 0;
    while (x) {
        count += x & 1u;
        x >>= 1;
    }
    return count;
}

/**
 * @brief NIST Monobit Test — проверяет, приблизительно ли равное количество 0 и 1.
 */
int nist_monobit(const uint32_t *w, size_t len) {
    int64_t ones = 0;
    for (size_t i = 0; i < len; ++i) ones += popcount32(w[i]);
    int n = (int)len * 32;
    double s = fabs(2.0 * ones - n) / sqrt((double)n);
    return erfc(s / sqrt(2.0)) >= 0.01;
}

/**
 * @brief NIST Block Frequency Test — проверяет равномерность битов в блоках размера M.
 */
int nist_block_frequency(const uint32_t *w, size_t len, size_t M) {
    size_t nBits = len * 32;
    size_t nBlocks = nBits / M;
    if (nBlocks < 20) return 0;

    double chi = 0.0;
    size_t bit = 0;
    for (size_t b = 0; b < nBlocks; ++b) {
        int ones = 0;
        for (size_t i = 0; i < M; ++i, ++bit) {
            ones += (w[bit / 32] >> (bit % 32)) & 1u;
        }
        double pi = (double)ones / (double)M;
        chi += (pi - 0.5) * (pi - 0.5);
    }
    chi *= 4.0 * M;
    double p = erfc(sqrt(chi / 2.0) / sqrt(nBlocks / 2.0));
    return p >= 0.01;
}

/**
 * @brief NIST Runs Test — проверяет, не слишком ли часто переключаются 0 и 1.
 */
int nist_runs(const uint32_t *w, size_t len) {
    size_t n = len * 32;
    int64_t ones = 0;
    int prev = w[0] & 1;
    ones += prev;
    int runsCnt = 1;

    for (size_t i = 1; i < n; ++i) {
        int bit = (w[i / 32] >> (i % 32)) & 1;
        ones += bit;
        if (bit != prev) {
            ++runsCnt;
            prev = bit;
        }
    }

    double pi = (double)ones / (double)n;
    if (fabs(pi - 0.5) > 2.0 / sqrt((double)n)) return 0;

    double expRuns = 2.0 * n * pi * (1.0 - pi);
    double z = fabs(runsCnt - expRuns) / (2.0 * sqrt(2.0 * n) * pi * (1.0 - pi));
    return erfc(z) >= 0.01;
}

/**
 * @brief NIST Cumulative Sums Test — проверяет смещения от нуля при суммировании битов.
 */
int nist_cumulative_sums(const uint32_t *w, size_t len) {
    int64_t s = 0;
    int64_t zmax = 0;
    size_t n = len * 32;

    for (size_t i = 0; i < n; ++i) {
        int bit = (w[i / 32] >> (i % 32)) & 1;
        s += bit ? 1 : -1;
        if (s > zmax) zmax = s;
        else if (-s > zmax) zmax = -s;
    }

    if (zmax == 0) return 0;

    double p = 1.0;
    int start = (int)((-((double)n / zmax) + 1.0) / 4.0);
    int end   = (int)(((double)n / zmax - 1.0) / 4.0);

    for (int k = start; k <= end; ++k) {
        double a = (4.0 * k + 1.0) * zmax / sqrt(2.0 * n);
        double b = (4.0 * k - 1.0) * zmax / sqrt(2.0 * n);
        p -= erfc(a) - erfc(b);
    }

    return p >= 0.01;
}

/**
 * @brief NIST Serial Test (2-й порядок) — проверяет частоты повторяющихся битовых шаблонов.
 */
int nist_serial2(const uint32_t *w, size_t len) {
    uint64_t c1[2] = {0, 0};   ///< Частоты одиночных битов (0, 1)
    uint64_t c2[4] = {0, 0, 0, 0}; ///< Частоты пар битов (00, 01, 10, 11)
    int first = w[0] & 1;
    int prev = -1;
    size_t n = len * 32;

    for (size_t i = 0; i < n; ++i) {
        int b = (w[i / 32] >> (i % 32)) & 1;
        ++c1[b];
        if (prev != -1) ++c2[(prev << 1) | b];
        prev = b;
    }
    ++c2[(prev << 1) | first];

    double dn = (double)(c1[0] + c1[1]);
    double psi1 = (c1[0]*c1[0] + c1[1]*c1[1]) * 2.0 / dn - dn;

    double psi2 = 0.0;
    for (int i = 0; i < 4; ++i) psi2 += c2[i] * c2[i];
    psi2 = psi2 * 4.0 / dn - dn;

    double diff = fabs(psi2 - psi1);
    return erfc(diff / (2.0 * sqrt(2.0 * dn))) >= 0.01;
}
