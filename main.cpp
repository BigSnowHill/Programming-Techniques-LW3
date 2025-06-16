#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <ctime>
#include "generators.h"
#include "stats.h"

/**
 * @brief Главная функция программы.
 *
 * Проводит статистический анализ последовательностей чисел, сгенерированных разными ГСЧ (LCG, XORShift, MWC),
 * с использованием тестов, таких как среднее значение, стандартное отклонение, коэффициент вариации,
 * критерий хи-квадрат и тесты из NIST.
 *
 * @return 0 при успешном завершении программы.
 */
int main() {
    /// Массив размеров выборок, которые будут протестированы.
    const int sample_sizes[] = {1000, 2000, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 55000, 60000,
    70000, 75000, 80000, 85000, 90000, 100000};

    /// Количество повторов тестирования для одной размерности.
    const int num_samples = 10;
    int sample_size;

    /// Диапазон значений для оценки распределения.
    const unsigned long long int range = (1ul << 32);

    /// Количество интервалов (bins) для гистограммы в тесте хи-квадрат.
    const int bins = 1000;

    // Переменные для хранения результатов тестов.
    double m = 0, s = 0, cv = 0, chi2 = 0, monobit = 0, block_frequency = 0, runs = 0, cumulative_sums = 0, serial2 = 0;
    double m_i, s_i, cv_i;

    /// Временные метки для измерения времени выполнения.
    clock_t t1, t2, t3, t4;

    /// Заголовок таблицы результатов.
    printf("Generator type  |       Mean     |      STDdev     |   CV   |     chi2      | monobit | block freq |  runs  | cumulative sums  | serial2 |   time\n");

    /**
     * @brief Тестирование генератора LCG (Linear Congruential Generator).
     */
    LCG lcg(1234);
    for (int ss = 0; ss < 20; ss++){
        t1 = clock();
        for (int i = 0; i < num_samples; ++i) {
            sample_size = sample_sizes[ss];
            uint32_t buffer[sample_size];

            // Генерация последовательности.
            for (int j = 0; j < sample_size; ++j)
                buffer[j] = lcg.next();

            t3 = clock();

            // Расчет статистик.
            m_i = mean(buffer, sample_size);
            s_i = stdev(buffer, sample_size, m_i);
            cv_i = coeff_var(m_i, s_i);

            m += m_i; s += s_i; cv += cv_i;
            chi2 += chi_squared(buffer, sample_size, bins, range);
            monobit += nist_monobit(buffer, sample_size);
            block_frequency += nist_block_frequency(buffer, sample_size, 128);
            runs += nist_runs(buffer, sample_size);
            cumulative_sums += nist_cumulative_sums(buffer, sample_size);
            serial2 += nist_serial2(buffer, sample_size);
            t4 = clock();
        }
        t2 = clock();

        /// Вывод результатов для LCG.
        printf("LCG      %-7d| %.2f  |  %.2f  | %.3f  | %-12.2f  |  %.2f   |    %.2f    |  %.2f  |      %.2f        |  %.2f   | %-6.2f ms\n",
                sample_sizes[ss], m / num_samples, s / num_samples, cv / num_samples, chi2 / num_samples,
                monobit / num_samples, block_frequency / num_samples, runs / num_samples, cumulative_sums / num_samples, serial2 / num_samples,
                1000.0*(t2 - t1 - (t4 - t3)) / CLOCKS_PER_SEC);

        // Сброс накопленных значений.
        m = 0; s = 0; cv = 0; chi2 = 0; monobit = 0; block_frequency = 0; runs = 0; cumulative_sums = 0; serial2 = 0;
    }

    /**
     * @brief Тестирование генератора XORShift32.
     */
    XORShift32 xor32(9876);
    for (int ss = 0; ss < 20; ss++){
        t1 = clock();
        for (int i = 0; i < num_samples; ++i) {
            sample_size = sample_sizes[ss];
            uint32_t buffer[sample_size];

            for (int j = 0; j < sample_size; ++j)
                buffer[j] = xor32.next();

            t3 = clock();

            m_i = mean(buffer, sample_size);
            s_i = stdev(buffer, sample_size, m_i);
            cv_i = coeff_var(m_i, s_i);

            m += m_i; s += s_i; cv += cv_i;
            chi2 += chi_squared(buffer, sample_size, bins, range);
            monobit += nist_monobit(buffer, sample_size);
            block_frequency += nist_block_frequency(buffer, sample_size, 128);
            runs += nist_runs(buffer, sample_size);
            cumulative_sums += nist_cumulative_sums(buffer, sample_size);
            serial2 += nist_serial2(buffer, sample_size);
            t4 = clock();
        }
        t2 = clock();

        /// Вывод результатов для XORShift32.
        printf("XORShift %-7d| %.2f  |  %.2f  | %.3f  | %-12.2f  |  %-.2f   |    %.2f    |  %.2f  |      %.2f        |  %.2f   | %-6.2f ms\n",
                sample_sizes[ss], m / num_samples, s / num_samples, cv / num_samples, chi2 / num_samples,
                monobit / num_samples, block_frequency / num_samples, runs / num_samples, cumulative_sums / num_samples, serial2 / num_samples,
                1000.0*(t2 - t1 - (t4 - t3)) / CLOCKS_PER_SEC);

        m = 0; s = 0; cv = 0; chi2 = 0; monobit = 0; block_frequency = 0; runs = 0; cumulative_sums = 0; serial2 = 0;
    }

    /**
     * @brief Тестирование генератора MWC (Multiply-With-Carry).
     */
    MWC mwc(13579);
    for (int ss = 0; ss < 20; ss++){
        t1 = clock();
        for (int i = 0; i < num_samples; ++i) {
            sample_size = sample_sizes[ss];
            uint32_t buffer[sample_size];

            for (int j = 0; j < sample_size; ++j)
                buffer[j] = mwc.next();

            t3 = clock();

            m_i = mean(buffer, sample_size);
            s_i = stdev(buffer, sample_size, m_i);
            cv_i = coeff_var(m_i, s_i);

            m += m_i; s += s_i; cv += cv_i;
            chi2 += chi_squared(buffer, sample_size, bins, range);
            monobit += nist_monobit(buffer, sample_size);
            block_frequency += nist_block_frequency(buffer, sample_size, 128);
            runs += nist_runs(buffer, sample_size);
            cumulative_sums += nist_cumulative_sums(buffer, sample_size);
            serial2 += nist_serial2(buffer, sample_size);
            t4 = clock();
        }
        t2 = clock();

        /// Вывод результатов для MWC.
        printf("MWC      %-7d| %.2f  |  %.2f  | %.3f  | %-12.2f  |  %.2f   |    %.2f    |  %.2f  |      %.2f        |  %.2f   | %-6.2f ms\n",
                sample_sizes[ss], m / num_samples, s / num_samples, cv / num_samples, chi2 / num_samples,
                monobit / num_samples, block_frequency / num_samples, runs / num_samples, cumulative_sums / num_samples, serial2 / num_samples,
                1000.0*(t2 - t1 - (t4 - t3)) / CLOCKS_PER_SEC);

        m = 0; s = 0; cv = 0; chi2 = 0; monobit = 0; block_frequency = 0; runs = 0; cumulative_sums = 0; serial2 = 0;
    }

    return 0;
}
