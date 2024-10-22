#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>

void variational_method() {
    // Константы
    const double hbar = 1.054571817e-34; // Дж·с
    const double me = 9.10938356e-31;    // кг
    const double eV = 1.602176634e-19;   // Дж
    const double angstrom = 1e-10;       // м

    // Параметры потенциальной ямы
    const double U0 = -0.5 * eV; // Потенциал в Дж
    const double L = 10 * angstrom; // Половина ширины ямы в метрах

    // Функция потенциала U(x)
    auto U = [&](double x) {
        if (fabs(x) < L) {
            return U0;
        } else {
            return 0.0;
        }
    };

    // Функция для расчета энергии в зависимости от параметра 'a'
    TF1 *EnergyFunc = new TF1("EnergyFunc", [&](double *a, double *) {
        double a_param = a[0];

        // Нормировочный множитель
        double norm = pow(2 / TMath::Pi(), 0.25) / sqrt(a_param);

        // Кинетическая энергия ⟨T⟩
        double T_integrand = [&](double *x, double *) {
            double psi = norm * exp(-x[0] * x[0] / (a_param * a_param));
            double d2psi_dx2 = norm * (4 * x[0] * x[0] - 2 * a_param * a_param) / (a_param * a_param * a_param * a_param) * exp(-x[0] * x[0] / (a_param * a_param));
            return - (hbar * hbar) / (2 * me) * psi * d2psi_dx2;
        };

        // Потенциальная энергия ⟨V⟩
        double V_integrand = [&](double *x, double *) {
            double psi = norm * exp(-x[0] * x[0] / (a_param * a_param));
            return psi * psi * U(x[0]);
        };

        // Пределы интегрирования (достаточно большие для экспоненты)
        double xmin = -5 * a_param;
        double xmax = 5 * a_param;

        // Интегрирование
        TF1 *TFunc = new TF1("TFunc", T_integrand, xmin, xmax, 0);
        double T = TFunc->Integral(xmin, xmax);

        TF1 *VFunc = new TF1("VFunc", V_integrand, xmin, xmax, 0);
        double V = VFunc->Integral(xmin, xmax);

        delete TFunc;
        delete VFunc;

        double Energy = T + V;
        return Energy;
    }, 0.1 * angstrom, 10 * angstrom, 0);

    // Поиск минимальной энергии и соответствующего 'a'
    double a_min = EnergyFunc->GetMinimumX(0.1 * angstrom, 10 * angstrom);

    std::cout << "Минимальная энергия: " << EnergyFunc->Eval(a_min) / eV << " эВ" << std::endl;
    std::cout << "Оптимальное значение параметра a: " << a_min / angstrom << " Ангстрем" << std::endl;

    // Построение волновой функции
    TCanvas *c1 = new TCanvas("c1", "Wave Function", 800, 600);

    TF1 *Psi = new TF1("Psi", [&](double *x, double *) {
        double psi = pow(2 / TMath::Pi(), 0.25) / sqrt(a_min) * exp(-x[0] * x[0] / (a_min * a_min));
        return psi;
    }, -2 * L, 2 * L, 0);

    Psi->SetTitle("Волновая функция основного состояния; x (м); Ψ(x)");
    Psi->SetLineColor(kBlue);
    Psi->Draw();

    c1->SaveAs("wave_function.png");
}
