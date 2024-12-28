#include <iostream>
#include <fstream>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TStyle.h>

// ----------------------------------------------------------------------------
// Функция для фитинга: Метод "P" с экспоненциальным спадом и гауссианом
// ----------------------------------------------------------------------------
double fitFunction(double *x, double *par) {
    // Экспоненциальный спад
    double expDecay = par[0] * TMath::Exp(-par[1] * x[0]);

    // Гауссиан с заданным средним и sigma
    double gaussian = par[2] * TMath::Gaus(x[0], 5.0, 1.0, true); // Среднее фиксировано на 5.0, sigma фиксирована на 1.0

    // Линейный хвост
    double linearTail = par[3] * (x[0] - 5.0); // Линейная коррекция с учетом смещения

    // Итоговая модель
    return expDecay + gaussian + linearTail;
}

// ----------------------------------------------------------------------------
// Функция для заполнения гистограммы
// ----------------------------------------------------------------------------
void FillHistogram(TH1D *hist, const std::string &filename) {
    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        std::cerr << "Ошибка открытия файла: " << filename << std::endl;
        exit(-1);
    }

    double value;
    while (inFile >> value) {
        hist->Fill(value);
    }
    inFile.close();
}

// ----------------------------------------------------------------------------
// Основная функция для выполнения задачи 12
// ----------------------------------------------------------------------------
void task12() {
    // Настройка стиля ROOT
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(1111);

    // Создание гистограммы
    TH1D *hist = new TH1D("hist", "Distribution data from task10Nov.dat;mm;yields / 0.1 mm", 100, 0, 10);
    hist->Sumw2();

    // Заполнение гистограммы данными из файла
    FillHistogram(hist, "task10Nov.dat");

    // Создание функции фитирования
    TF1 *fitFunc = new TF1("fitFunc", fitFunction, 0, 10, 4);

    // Установка начальных значений параметров
    fitFunc->SetParameter(0, 28);   // Амплитуда экспоненты
    fitFunc->SetParameter(1, 0.65);  // Коэффициент экспоненты
    fitFunc->SetParameter(2, 54.0); // Амплитуда гауссиана
    fitFunc->SetParameter(3, 1.0);  // Наклон линейного хвоста

    // Установка ограничений параметров
    fitFunc->SetParLimits(0, 25, 30);
    fitFunc->SetParLimits(1, 0.5, 0.7);
    fitFunc->SetParLimits(2, 52.0, 56.0);
    fitFunc->SetParLimits(3, 0.7, 1.5);

    // Выполнение фитирования методом "P"
    TFitResultPtr fitRes = hist->Fit(fitFunc, "PS");

    // Создание холста для отображения результата
    TCanvas *c1 = new TCanvas("c1", "Diagnostic Fit", 1200, 600);
    hist->SetLineColor(kBlue + 2);
    hist->SetLineWidth(2);
    hist->Draw("HIST");

    // Добавление функции фита на график
    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineWidth(2);
    fitFunc->Draw("SAME");

    // Сохранение графика
    c1->SaveAs("task12_fit_result_P.png");

    // Вывод параметров в консоль
    std::cout << "Chi2/NDF = " << fitRes->Chi2() / fitRes->Ndf() << std::endl;
    for (int i = 0; i < 4; ++i) {
        std::cout << "Parameter " << i << ": " << fitFunc->GetParameter(i)
                  << " +/- " << fitFunc->GetParError(i) << std::endl;
    }
}
