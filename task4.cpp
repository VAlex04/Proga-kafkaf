#include <iostream>
#include <cmath>
#include <TRandom.h>
#include <TMath.h>

// Функция-макрос для оценки числа π методом Бюффона
void task4() {
    // Параметры иглы и расстояния между линиями
    double distance_between_lines = 1.0; // Расстояние между параллельными линиями (D)
    double needle_length = 1.0; // Длина иглы (L)

    // Количество бросков иглы
    long long num_drops = 1000000; // N_drops
    long long num_crossings = 0; // N_crossings

    // Цикл бросков иглы
    for (long long i = 0; i < num_drops; ++i) {
        // Случайное положение центра иглы между 0 и D/2
        double center_position = (distance_between_lines / 2.0) * gRandom->Rndm();

        // Случайный угол между 0 и π/2
        double angle = (TMath::Pi() / 2.0) * gRandom->Rndm();

        // Проверка пересечения линии
        if (center_position <= (needle_length / 2.0) * TMath::Sin(angle)) {
            num_crossings++;
        }
    }

    // Проверка, чтобы избежать деления на ноль
    if (num_crossings == 0) {
        std::cerr << "Нет пересечений для оценки числа π. Попробуйте увеличить количество бросков." << std::endl;
        return;
    }

    // Оценка числа π по формуле Буффона
    double pi_estimate = (2.0 * needle_length * num_drops) / (distance_between_lines * num_crossings);

    // Вывод результатов
    std::cout << "Оценка π: " << pi_estimate << std::endl;
    std::cout << "Фактическое π: " << TMath::Pi() << std::endl;
    std::cout << "Относительная ошибка: " << std::abs(pi_estimate - TMath::Pi()) / TMath::Pi() * 100.0 << " %" << std::endl;
}
