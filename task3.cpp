// task3.C
#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>

// Константы
const Double_t hbar2_over_2m = 3.80998; // (ħ²) / (2m) в eV·Å²

// Функция для расчёта энергии <Ψ|H|Ψ>
Double_t E_func(Double_t *a, Double_t *par)
{
    // Кинетическая энергия: <T> = ħ² / (2m a²)
    Double_t kinetic = hbar2_over_2m / (a[0] * a[0]);

    // Потенциальная энергия: <U> = -0.5 * erf(10 * sqrt(2) / a)
    Double_t potential = -0.5 * TMath::Erf(10.0 * TMath::Sqrt(2.0) / a[0]);

    // Общая энергия
    Double_t total_energy = kinetic + potential;

    return total_energy;
}

// Функция волновой функции ψ(x) с параметром a
Double_t psi_func(Double_t *x, Double_t *par)
{
    Double_t a_min = par[0]; // Оптимальный параметр 'a'
    Double_t N = pow(2.0 / TMath::Pi(), 0.25) / sqrt(a_min);
    Double_t psi = N * exp(-x[0] * x[0] / (a_min * a_min));
    return psi;
}

void task3()
{
    // Создаём функцию энергии E(a)
    TF1 *E = new TF1("E", E_func, 0.1, 50.0, 0);
    E->SetNpx(1000); // Увеличиваем количество точек для точности

    // Создаём холст для графика E(a)
    TCanvas *c2 = new TCanvas("c2", "Energy Expectation Value E(a)", 800, 600);
    E->SetTitle("Energy Expectation Value E(a); a (Å); E(a) (eV)");
    E->SetLineColor(kBlue);
    E->Draw();

    // Находим минимум энергии и соответствующее значение 'a'
    Double_t a_min = E->GetMinimumX(0.1, 50.0);
    Double_t E_min = E->Eval(a_min);
    std::cout << "Минимальная энергия E(a) при a = " << a_min << " Å составляет E = " << E_min << " eV" << std::endl;

    // Отмечаем минимум на графике
    TMarker *marker = new TMarker(a_min, E_min, 20);
    marker->SetMarkerColor(kRed);
    marker->SetMarkerSize(1.5);
    marker->Draw("same");

    // Создаём функцию волновой функции ψ(x) с оптимальным 'a'
    TF1 *psi = new TF1("psi", psi_func, -20.0, 20.0, 1);
    psi->SetParameter(0, a_min); // Устанавливаем оптимальный 'a'
    psi->SetLineColor(kGreen+2);
    psi->SetTitle("Ground State Wavefunction ψ(x); x (Å); ψ(x)");
    psi->SetNpx(1000); // Увеличиваем количество точек для гладкости

    // Создаём холст для графика ψ(x)
    TCanvas *c1 = new TCanvas("c1", "Wavefunction ψ(x)", 800, 600);
    psi->Draw();

    // Сохраняем холсты в файлы (опционально)
    // c2->SaveAs("Energy_vs_a.png");
    // c1->SaveAs("Wavefunction.png");
}
