/**
 * @file generators.h
 * @brief Заголовочный файл с реализациями различных генераторов псевдослучайных чисел.
 */

#ifndef GENERATORS_H
#define GENERATORS_H

#include <cstdint>

/**
 * @brief Линейный конгруэнтный генератор (LCG).
 *
 * Использует формулу: X_{n+1} = a * X_n + c (mod 2^32)
 */
struct LCG {
    uint32_t state; ///< Текущее состояние генератора.
    static constexpr uint32_t a = 1664525u; ///< Множитель.
    static constexpr uint32_t c = 1013904223u; ///< Прибавка.

    /**
     * @brief Конструктор генератора.
     * @param s Начальное значение состояния (по умолчанию 1).
     */
    LCG(uint32_t s = 1u) : state(s) {}

    /**
     * @brief Генерация следующего псевдослучайного числа.
     * @return Следующее значение в последовательности.
     */
    inline uint32_t next() { return state = a * state + c; }
};

/**
 * @brief Генератор XORShift32.
 *
 * Быстрый генератор на основе побитовых операций XOR и сдвигов.
 */
struct XORShift32 {
    uint32_t state; ///< Текущее состояние генератора.

    /**
     * @brief Конструктор генератора.
     * @param s Начальное значение состояния (по умолчанию 2463534242).
     */
    XORShift32(uint32_t s = 2463534242u) : state(s) {}

    /**
     * @brief Генерация следующего псевдослучайного числа.
     * @return Следующее значение в последовательности.
     */
    inline uint32_t next() {
        uint32_t x = state;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        return state = x;
    }
};

/**
 * @brief Генератор Multiply-With-Carry (MWC).
 *
 * Использует умножение и перенос для генерации случайных чисел.
 */
struct MWC {
    uint32_t state; ///< Основное состояние генератора.
    uint32_t carry; ///< Перенос при умножении.
    static constexpr uint32_t a = 4294957665u; ///< Множитель.

    /**
     * @brief Конструктор генератора.
     * @param s Начальное значение состояния в виде 64-битного числа,
     *          где старшие 32 бита — перенос, младшие — состояние.
     */
    MWC(uint64_t s = 88172645463325252ull)
        : state(uint32_t(s)), carry(uint32_t(s >> 32)) {}

    /**
     * @brief Генерация следующего псевдослучайного числа.
     * @return Следующее значение в последовательности.
     */
    inline uint32_t next() {
        uint64_t p = uint64_t(a) * state + carry;
        state = uint32_t(p);
        carry = uint32_t(p >> 32);
        return state;
    }
};

#endif // GENERATORS_H
