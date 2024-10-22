#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMarker.h>
#include <iostream>

// Константы
const Double_t hbar2_over_2m_graph = 3.80998; // (ħ²)/(2m) в eV·Å²
const Double_t L_graph = 10.0;                // Полуширина потенциальной ямы в Å

// Энергия ожидания как функция a для заданного V0
Double_t E_graph_func(Double_t *a_arr, Double_t *par)
{
    Double_t a = a_arr[0];
    Double_t V0 = par[0];

    // Кинетическая энергия: (ħ²)/(2m a²)
    Double_t kinetic = hbar2_over_2m_graph / (a * a);

    // Потенциальная энергия: -V0 * Erf(√2 * L / a)
    Double_t potential = -V0 * TMath::Erf(TMath::Sqrt(2.0) * L_graph / a);

    Double_t total_energy = kinetic + potential;
    return total_energy;
}

// Вычисление ⟨x²⟩ для заданного a
Double_t x2_expectation(Double_t a)
{
    // ⟨x²⟩ = (a²) / 2
    Double_t x2 = (a * a) / 2.0;
    return x2;
}

void task3graph()
{
    // Диапазон значений V0
    const Int_t N = 50;        // Количество точек
    Double_t V0_min = 0.1;     // Минимальное значение V0 в эВ
    Double_t V0_max = 5.0;     // Максимальное значение V0 в эВ
    Double_t dV0 = (V0_max - V0_min) / (N - 1);

    // Массивы для хранения данных
    Double_t V0_values[N];
    Double_t x2_values[N];

    // Цикл по значениям V0
    for (Int_t i = 0; i < N; i++)
    {
        Double_t V0 = V0_min + i * dV0;
        V0_values[i] = V0;

        // Создание функции энергии E(a) для текущего V0
        TF1 *E = new TF1("E_graph", E_graph_func, 0.1, 50.0, 1);
        E->SetParameter(0, V0);
        E->SetNpx(1000);

        // Поиск оптимального a
        Double_t a_min = E->GetMinimumX(0.1, 50.0);
        Double_t E_min = E->Eval(a_min);

        // Вычисление ⟨x²⟩
        Double_t x2 = x2_expectation(a_min);
        x2_values[i] = x2;

        // Освобождение памяти
        delete E;
    }

    // Создание графика
    TGraph *gr = new TGraph(N, V0_values, x2_values);
    gr->SetTitle("<x^{2}> vs V_{0};V_{0} (eV);<x^{2}> (Å^{2})");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.8);
    gr->SetLineColor(kBlue);

    // Отображение графика
    TCanvas *c1 = new TCanvas("c1_graph", "<x^{2}> vs V_{0}", 800, 600);
    gr->Draw("APL");

    // Дополнительные настройки осей
    gr->GetXaxis()->SetLimits(V0_min, V0_max);
    gr->GetYaxis()->SetRangeUser(0, x2_values[N - 1] * 1.2);

    // Обновление холста
    c1->Update();
}
