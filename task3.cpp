// task3.C
#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <iostream>

// Константы
const Double_t hbar2_over_2m = 3.80998; // (ħ²)/(2m) в eV·Å²
const Double_t V0 = 0.5;                // Глубина потенциальной ямы в эВ

// Энергия ожидания как функция a
Double_t E_func(Double_t *a, Double_t *par)
{
    // Кинетическая энергия: (ħ²)/(2m a²)
    Double_t kinetic = hbar2_over_2m / (a[0] * a[0]);

    // Потенциальная энергия: -V0 * Erf(√2 * 10 / a)
    Double_t potential = -V0 * TMath::Erf(TMath::Sqrt(2.0) * 10.0 / a[0]);

    Double_t total_energy = kinetic + potential;
    return total_energy;
}

// Волновая функция ψ(x) с параметром a
Double_t psi_func(Double_t *x, Double_t *par)
{
    Double_t a_min = par[0]; // Оптимальное значение 'a'
    Double_t N = TMath::Power(2.0 / TMath::Pi(), 0.25) / TMath::Sqrt(a_min);
    Double_t psi = N * TMath::Exp(-x[0] * x[0] / (a_min * a_min));
    return psi;
}

void task3()
{
    // Создание функции энергии E(a)
    TF1 *E = new TF1("E", E_func, 0.1, 50.0, 0);
    E->SetNpx(1000); // Увеличиваем количество точек для более гладкого графика

    // Создание полотна для E(a)
    TCanvas *c2 = new TCanvas("c2", "Energy Expectation Value E(a)", 800, 600);
    E->SetTitle("Energy Expectation Value E(a);a (Å);E(a) (eV)");
    E->SetLineColor(kBlue);
    E->Draw();

    // Поиск минимума энергии и соответствующего 'a'
    Double_t a_min = E->GetMinimumX(0.1, 50.0);
    Double_t E_min = E->Eval(a_min);
    std::cout << "Минимальная энергия E(a) при a = " << a_min << " Å составляет E = " << E_min << " эВ" << std::endl;

    // Добавление точки минимума на график
    TMarker *min_marker = new TMarker(a_min, E_min, 20);
    min_marker->SetMarkerColor(kRed);
    min_marker->SetMarkerSize(1.5);
    min_marker->Draw("same");

    // Создание волновой функции ψ(x) с оптимальным 'a'
    TF1 *psi = new TF1("psi", psi_func, -20.0, 20.0, 1);
    psi->SetParameter(0, a_min); // Установка оптимального 'a'
    psi->SetLineColor(kGreen + 2);
    psi->SetLineWidth(2);
    psi->SetTitle("Ground State Wavefunction ψ(x);x (Å);ψ(x)");
    
    // Создание полотна для ψ(x)
    TCanvas *c1 = new TCanvas("c1", "Wavefunction ψ(x)", 800, 600);
    psi->Draw();

    // Дополнительные улучшения графика ψ(x)
    psi->GetXaxis()->SetRangeUser(-20, 20);
    psi->GetYaxis()->SetRangeUser(0, psi->GetMaximum() * 1.2);
    
    // Отображение графиков
    c1->Update();
    c2->Update();
}
