// task3star.C
#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMarker.h>
#include <iostream>

// Константы
const Double_t hbar2_over_2m_star = 3.80998; // (ħ²)/(2м) в eV·Å²
const Double_t V0_star = 0.5;                // Глубина потенциальной ямы в эВ

// Волновая функция ψ_star(x, a)
Double_t psi_star(Double_t *x, Double_t *par)
{
    Double_t a = par[0]; // Параметр 'a'
    Double_t A = par[1]; // Нормировочный множитель 'A'
    Double_t abs_x = fabs(x[0]);

    // Вычисление ψ(x)
    Double_t numerator = a * x[0] * TMath::Exp(-abs_x / a);
    Double_t denominator = x[0] * x[0] + a * a;
    Double_t psi = A * numerator / denominator;
    return psi;
}

// Модуль квадрата волновой функции |ψ(x)|²
Double_t psi_star_squared(Double_t *x, Double_t *par)
{
    Double_t psi = psi_star(x, par);
    return psi * psi;
}

// Первая производная волновой функции ψ'_star(x, a)
Double_t psi_star_derivative(Double_t *x, Double_t *par)
{
    Double_t a = par[0];
    Double_t A = par[1];
    Double_t abs_x = fabs(x[0]);
    Double_t sgn = (x[0] >= 0) ? 1.0 : -1.0;

    // Вычисление компонентов
    Double_t exp_part = TMath::Exp(-abs_x / a);
    Double_t numerator = a * x[0];
    Double_t denominator = x[0] * x[0] + a * a;

    // Вычисление производной
    Double_t dpsi_dx = A * exp_part * (
        (a * (x[0] * x[0] - a * a) / (denominator * denominator))
        - (sgn * numerator / (a * denominator))
    );
    return dpsi_dx;
}

// Интеграл для нормировки
Double_t norm_star_integrand(Double_t *x, Double_t *par)
{
    return psi_star_squared(x, par);
}

// Интеграл для кинетической энергии
Double_t kinetic_star_integrand(Double_t *x, Double_t *par)
{
    Double_t dpsi_dx = psi_star_derivative(x, par);
    return dpsi_dx * dpsi_dx;
}

// Интеграл для потенциальной энергии
Double_t potential_star_integrand(Double_t *x, Double_t *par)
{
    return psi_star_squared(x, par);
}

// Функция энергии E_star(a)
Double_t E_star_func(Double_t *a_arr, Double_t *par)
{
    Double_t a = a_arr[0];

    // Устанавливаем параметры для интегралов
    Double_t par_int[2] = {a, 1.0}; // A = 1.0 временно

    // Интегрирование для нормировки
    TF1 *norm_int = new TF1("norm_int_star", norm_star_integrand, 0, 100, 2);
    norm_int->SetParameters(par_int);
    Double_t norm = 2.0 * norm_int->Integral(0, 100);

    // Вычисление нормировочного множителя A(a)
    Double_t A = 1.0 / TMath::Sqrt(norm);
    par_int[1] = A;

    // Интегрирование для кинетической энергии
    TF1 *kinetic_int = new TF1("kinetic_int_star", kinetic_star_integrand, 0, 100, 2);
    kinetic_int->SetParameters(par_int);
    Double_t T = hbar2_over_2m_star * 2.0 * kinetic_int->Integral(0, 100);

    // Интегрирование для потенциальной энергии
    TF1 *potential_int = new TF1("potential_int_star", potential_star_integrand, 0, 10, 2);
    potential_int->SetParameters(par_int);
    Double_t V = -V0_star * 2.0 * potential_int->Integral(0, 10);

    // Полная энергия
    Double_t E_total = T + V;

    // Освобождаем память
    delete norm_int;
    delete kinetic_int;
    delete potential_int;

    return E_total;
}

void task3star()
{
    // Задаём желаемый диапазон по оси X для графика волновой функции
    Double_t x_min = -30.0; // Минимальное значение X
    Double_t x_max = 30.0;  // Максимальное значение X

    // Создание функции энергии E_star(a)
    TF1 *E = new TF1("E_star", E_star_func, 0.1, 50.0, 0);
    E->SetNpx(1000); // Увеличиваем количество точек для более гладкого графика

    // Создание полотна для E(a)
    TCanvas *c2 = new TCanvas("c2_star", "Energy Expectation Value E(a)", 800, 600);
    E->SetTitle("Energy Expectation Value E(a);a (Å);E(a) (eV)");
    E->SetLineColor(kBlue);
    E->Draw();

    // Поиск минимума энергии и соответствующего 'a' в полном диапазоне
    Double_t a_min = E->GetMinimumX(0.1, 50.0);
    Double_t E_min = E->Eval(a_min);
    std::cout << "Минимальная энергия E(a) при a = " << a_min << " Å составляет E = " << E_min << " эВ" << std::endl;

    // Добавление точки минимума на график
    TMarker *min_marker = new TMarker(a_min, E_min, 20);
    min_marker->SetMarkerColor(kRed);
    min_marker->SetMarkerSize(1.5);
    min_marker->Draw("same");

    // Вычисление нормировочного множителя A(a_min)
    Double_t par_int[2] = {a_min, 1.0}; // A = 1.0 временно
    TF1 *norm_int = new TF1("norm_int_star", norm_star_integrand, 0, 100, 2);
    norm_int->SetParameters(par_int);
    Double_t norm = 2.0 * norm_int->Integral(0, 100);
    Double_t A_min = 1.0 / TMath::Sqrt(norm);
    std::cout << "Нормировочный множитель A = " << A_min << std::endl;

    // Создание волновой функции ψ(x) с оптимальным 'a' и 'A'
    TF1 *psi = new TF1("psi_star", psi_star, x_min, x_max, 2);
    psi->SetParameters(a_min, A_min); // Установка оптимального 'a' и 'A'
    psi->SetLineColor(kGreen + 2);
    psi->SetLineWidth(2);
    psi->SetTitle("Wavefunction ψ(x);x (Å);ψ(x)");

    // Создание полотна для ψ(x)
    TCanvas *c1 = new TCanvas("c1_star", "Wavefunction ψ(x)", 800, 600);
    psi->Draw();

    // Автоматическое масштабирование осей для волновой функции
    // Обновляем холст, чтобы получить доступ к гистограмме
    c1->Update();
    TH1 *hist = psi->GetHistogram();

    // Получаем минимальные и максимальные значения по Y
    Double_t y_min = hist->GetMinimum();
    Double_t y_max = hist->GetMaximum();

    // Устанавливаем новые пределы по оси Y с некоторым запасом
    hist->GetYaxis()->SetRangeUser(1.2 * y_min, 1.2 * y_max);

    // Устанавливаем пределы по оси X в соответствии с заданными x_min и x_max
    hist->GetXaxis()->SetRangeUser(x_min, x_max);

    // Обновление графика после изменения пределов
    c1->Modified();
    c1->Update();

    // Отображение графика энергии
    c2->Update();
}
