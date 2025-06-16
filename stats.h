/**
 * @file stats.h
 * @brief Заголовочный файл с функциями статистического анализа для оценки качества генераторов случайных чисел.
 */

#ifndef STATS_H
#define STATS_H

#include <cstdint>
#include <cstddef>

/**
 * @brief Вычисляет среднее значение выборки.
 * @param data Указатель на массив данных.
 * @param n Размер выборки.
 * @return Среднее арифметическое значений.
 */
double mean(const uint32_t *data, int n);

/**
 * @brief Вычисляет стандартное отклонение выборки.
 * @param data Указатель на массив данных.
 * @param n Размер выборки.
 * @param m Среднее значение выборки.
 * @return Стандартное отклонение.
 */
double stdev(const uint32_t *data, int n, double m);

/**
 * @brief Вычисляет коэффициент вариации.
 * @param m Среднее значение.
 * @param sd Стандартное отклонение.
 * @return Коэффициент вариации (sd / m).
 */
double coeff_var(double m, double sd);

/**
 * @brief Вычисляет значение критерия хи-квадрат для равномерности распределения.
 * @param data Указатель на массив данных.
 * @param n Размер выборки.
 * @param bins Количество интервалов (корзин).
 * @param max_val Максимально возможное значение случайной величины.
 * @return Значение критерия хи-квадрат.
 */
double chi_squared(const uint32_t *data, int n, int bins, unsigned long long int max_val);

/**
 * @brief Выполняет тест Моно-бита (NIST STS).
 * @param w Указатель на массив данных.
 * @param len Количество элементов в массиве.
 * @return Результат теста (1 — успешно, 0 — неуспешно).
 */
int nist_monobit(const uint32_t *w, size_t len);

/**
 * @brief Выполняет тест частот блоков (Block Frequency Test, NIST STS).
 * @param w Указатель на массив данных.
 * @param len Количество элементов.
 * @param M Размер блока в битах.
 * @return Результат теста (1 — успешно, 0 — неуспешно).
 */
int nist_block_frequency(const uint32_t *w, size_t len, size_t M);

/**
 * @brief Выполняет тест на количество последовательностей одинаковых битов (Runs Test, NIST STS).
 * @param w Указатель на массив данных.
 * @param len Количество элементов.
 * @return Результат теста (1 — успешно, 0 — неуспешно).
 */
int nist_runs(const uint32_t *w, size_t len);

/**
 * @brief Выполняет тест кумулятивных сумм (Cumulative Sums Test, NIST STS).
 * @param w Указатель на массив данных.
 * @param len Количество элементов.
 * @return Результат теста (1 — успешно, 0 — неуспешно).
 */
int nist_cumulative_sums(const uint32_t *w, size_t len);

/**
 * @brief Выполняет тест на серийность второго порядка (Serial Test, NIST STS).
 * @param w Указатель на массив данных.
 * @param len Количество элементов.
 * @return Результат теста (1 — успешно, 0 — неуспешно).
 */
int nist_serial2(const uint32_t *w, size_t len);

#endif // STATS_H
