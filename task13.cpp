#include <TROOT.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>

// ----------------------------------------------------------------------------
// Функция резолюции: две гауссианы с амплитудами
// ----------------------------------------------------------------------------
double Resolution(double x, double A1, double sigma1, double A2, double mu2, double sigma2)
{
   // Первая гауссиана: основное разрешение
   double g1 = A1 * TMath::Exp(-0.5 * (x / sigma1) * (x / sigma1))
             / (TMath::Sqrt(2 * TMath::Pi()) * sigma1);

   // Вторая гауссиана: добавляет шум и сдвиг mu2
   double dx = x - mu2;
   double g2 = A2 * TMath::Exp(-0.5 * (dx / sigma2) * (dx / sigma2))
             / (TMath::Sqrt(2 * TMath::Pi()) * sigma2);

   return g1 + g2; // Сумма двух гауссиан
}

// ----------------------------------------------------------------------------
// Функция свёртки: Брейт-Вигнер * Резолюция + Фон
// ----------------------------------------------------------------------------
double BWConvolution(double *x, double *par)
{
   // 1) Линейный фон
   double bg = par[0] + par[1] * x[0];

   // 2) Амплитуда сигнала
   double A = par[2];

   // 3) Параметры резолюции
   double A1 = par[3];
   double sigma1 = par[4];
   double A2 = par[5];
   double mu2 = par[6];
   double sigma2 = par[7];

   // 4) Численное интегрирование
   double tMin = 3.0;
   double tMax = 3.2;
   int Nstep = 200;
   double step = (tMax - tMin) / Nstep;

   double sum = 0.0;
   for (int i = 0; i < Nstep; i++)
   {
      double t = tMin + (i + 0.5) * step;

      // Брейт-Вигнер
      double bwVal = TMath::BreitWigner(t, 3.0969, 0.000093);

      // Резолюция
      double rVal = Resolution(x[0] - t, A1, sigma1, A2, mu2, sigma2);

      sum += bwVal * rVal;
   }
   sum *= step;

   return bg + A * sum;
}

// ----------------------------------------------------------------------------
// Основная функция: анализ данных и фитинг
// ----------------------------------------------------------------------------
void task13()
{
    // Настройки стиля ROOT
    gStyle->SetOptStat(111); // Показывать только Entries, Mean, Std Dev
    gStyle->SetOptFit(1);    // Показывать Chi2/NDF и параметры фита

    // Удаляем старую канву, если она существует
    if (gROOT->FindObject("c1")) delete gROOT->FindObject("c1");

    // 1) Открытие файла с данными
    std::ifstream fin("m3piJPSI_cut.dat");
    if (!fin.is_open())
    {
        std::cerr << "Нет файла m3piJPSI_cut.dat" << std::endl;
        return;
    }

    // 2) Создание гистограммы
    TH1F *hMass = new TH1F("hMass", "J/psi -> 3#pi; m (GeV/c^{2}); Events", 400, 3.0, 3.2);
    double massVal;
    while (fin >> massVal)
    {
        if (!fin.good())
            break;
        hMass->Fill(massVal); // Заполнение гистограммы данными
    }
    fin.close();

    // 3) Создание функции для фита
    TF1 *fFit = new TF1("fFit", BWConvolution, 3.0, 3.2, 8);
    fFit->SetNpx(800); // Увеличиваем число точек для сглаживания графика

    // Установка имён параметров
    fFit->SetParName(0, "p0_bg");    // Константа фона
    fFit->SetParName(1, "p1_bg");    // Наклон фона
    fFit->SetParName(2, "A_signal"); // Амплитуда сигнала
    fFit->SetParName(3, "A1");       // Амплитуда первой гауссианы
    fFit->SetParName(4, "sigma1");   // Ширина первой гауссианы
    fFit->SetParName(5, "A2");       // Амплитуда второй гауссианы
    fFit->SetParName(6, "mu2");      // Сдвиг второй гауссианы
    fFit->SetParName(7, "sigma2");   // Ширина второй гауссианы

    // Установка начальных параметров
    fFit->SetParameter(0,  0.0);    // Константа фона
    fFit->SetParameter(1,  0.0);    // Наклон фона
    fFit->SetParameter(2,  500.0);  // Амплитуда сигнала
    fFit->SetParameter(3,  1.0);    // Амплитуда первой гауссианы
    fFit->SetParameter(4,  0.005);  // Начальное значение sigma1 = 5 МэВ
    fFit->SetParLimits(4, 1e-6, 0); // Только нижний предел для sigma1
    fFit->SetParameter(5,  1.0);    // Амплитуда второй гауссианы
    fFit->SetParameter(6,  0.0);    // Сдвиг второй гауссианы
    fFit->SetParameter(7,  0.003);  // Начальное значение sigma2 = 3 МэВ
    fFit->SetParLimits(7, 1e-6, 0); // Только нижний предел для sigma2

    // Фитинг гистограммы
    hMass->Fit(fFit, "R");

    // 4) Рисование результата
    TCanvas *c1 = new TCanvas("c1", "J/psi -> 3pi with BW * Gaus conv", 900, 700);
    c1->SetLogy();
    hMass->Draw("E");

    fFit->SetLineColor(kRed);
    fFit->SetLineWidth(2);
    fFit->Draw("same");

    // Линейный фон
    TF1 *fBg = new TF1("fBg", "[0] + [1]*x", 3.0, 3.2);
    fBg->SetParameters(fFit->GetParameter(0), fFit->GetParameter(1));
    fBg->SetLineColor(kBlack);
    fBg->SetLineStyle(2);
    fBg->Draw("same");

    // Убираем стандартную легенду ROOT
    delete c1->GetListOfPrimitives()->FindObject("TPave");

    // Сохранение графика
    c1->SaveAs("task13_result_corrected.png");
}
