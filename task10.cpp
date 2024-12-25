#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TMinuit.h"
#include "TF1.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLegend.h"

static TH1D *h1 = nullptr; // Гистограмма для data_1.dat
static TH1D *h2 = nullptr; // Гистограмма для data_2.dat

// Функция для минимизации суммы χ²
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
    double chi2 = 0.0;

    // Модель: 
    // h1: fval = const + (amplitude/(sigma*sqrt(2*pi)))*exp(-0.5*((x-mean)^2/sigma^2))
    // h2: fval = const
    // amplitude в таком определении равна общему числу событий под гауссианой.

    for (int i = 1; i <= h1->GetNbinsX(); i++) {
        double x  = h1->GetBinCenter(i);
        double y  = h1->GetBinContent(i);
        double ey = (y > 0) ? std::sqrt(y) : std::sqrt(par[0]);

        double gauss_val = (par[1]/(par[3]*std::sqrt(2*M_PI))) * std::exp(-0.5*(x - par[2])*(x - par[2])/(par[3]*par[3]));
        double fval = par[0] + gauss_val;

        double diff = (y - fval)/ey;
        chi2 += diff*diff;
    }

    for (int i = 1; i <= h2->GetNbinsX(); i++) {
        double x  = h2->GetBinCenter(i);
        double y  = h2->GetBinContent(i);
        double ey = (y > 0) ? std::sqrt(y) : std::sqrt(par[0]);

        double fval = par[0]; // Только фон для h2
        double diff = (y - fval)/ey;
        chi2 += diff*diff;
    }

    f = chi2;
}

void task10() {
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(0);

    h1 = new TH1D("h1","Data 1",100,500,600);
    h2 = new TH1D("h2","Data 2",100,500,600);

    {
        std::ifstream in("data_1.dat");
        if(!in) {
            std::cerr << "Не удалось открыть data_1.dat" << std::endl;
            return;
        }
        double val;
        while(in >> val) {
            if(val>=500 && val<=600) h1->Fill(val);
        }
        in.close();
    }

    {
        std::ifstream in("data_2.dat");
        if(!in) {
            std::cerr << "Не удалось открыть data_2.dat" << std::endl;
            return;
        }
        double val;
        while(in >> val) {
            if(val>=500 && val<=600) h2->Fill(val);
        }
        in.close();
    }

    TMinuit *gMinuit = new TMinuit(4);
    gMinuit->SetFCN(fcn);

    double p0=5.0;   // const
    double p1=40.0;  // amplitude (общее число событий под гауссианой)
    double p2=550.0; // mean
    double p3=10.0;  // sigma

    double step=0.1;
    gMinuit->DefineParameter(0,"const", p0, step, 0, 1e6);
    gMinuit->DefineParameter(1,"amplitude", p1, step, 0, 1e6);
    gMinuit->DefineParameter(2,"mean", p2, step, 500, 600);
    gMinuit->DefineParameter(3,"sigma", p3, step, 0.1, 100);

    gMinuit->Command("MIGRAD");
    gMinuit->Command("HESSE");

    double par[4], err[4];
    for (int i=0;i<4;i++){
        gMinuit->GetParameter(i,par[i],err[i]);
    }

    double final_chi2;
    {
        Int_t npar=4;
        Double_t *gin=0;
        Int_t iflag=0;
        fcn(npar, gin, final_chi2, par, iflag);
    }

    int ndof = (h1->GetNbinsX()+h2->GetNbinsX()) - 4;
    double chi2_per_ndf = final_chi2/ndof;

    // Число событий под гауссианой = par[1], ошибка = err[1]
    double N_signal = par[1];
    double N_signal_err = err[1];

    // Проверка численной интеграции: 
    // Интегрируем гауссиану от -∞ до +∞ и должны получить N_signal.
    TF1 *gaussCheck = new TF1("gaussCheck","[0]/([2]*sqrt(2*pi))*exp(-0.5*((x-[1])*(x-[1]))/([2]*[2]))",-1e6,1e6);
    gaussCheck->SetParameters(par[1],par[2],par[3]);
    double numeric_integral = gaussCheck->Integral(-1e6,1e6);

    std::cout << "Результаты подгонки:" << std::endl;
    std::cout << "const = " << par[0] << " ± " << err[0] << std::endl;
    std::cout << "amplitude (N_signal) = " << par[1] << " ± " << err[1] << std::endl;
    std::cout << "mean = " << par[2] << " ± " << err[2] << std::endl;
    std::cout << "sigma = " << par[3] << " ± " << err[3] << std::endl;
    std::cout << "Chi2 = " << final_chi2 << " for " << ndof << " d.o.f => Chi2/ndf = " << chi2_per_ndf << std::endl;
    std::cout << "Число событий под гауссом (аналитически) = " << N_signal << " ± " << N_signal_err << std::endl;
    std::cout << "Число событий под гауссом (числ. интеграл) = " << numeric_integral << std::endl;

    TCanvas *c = new TCanvas("c","Fits",1200,600);
    c->Divide(2,1);

    c->cd(1);
    h1->SetMarkerStyle(20);
    h1->SetMarkerColor(kBlue);
    h1->SetTitle("Data 1 with fit; x; counts");
    h1->Draw("E");

    TF1 *f1 = new TF1("f1","[0] + ([1]/([3]*sqrt(2*pi)))*exp(-0.5*((x-[2])*(x-[2]))/([3]*[3]))",500,600);
    f1->SetParameters(par[0],par[1],par[2],par[3]);
    f1->SetLineColor(kRed);
    f1->Draw("SAME");

    {
        TLegend *leg = new TLegend(0.65,0.75,0.9,0.9);
        leg->AddEntry(h1,"Data 1","p");
        leg->AddEntry(f1,"Const+Gauss fit","l");
        leg->Draw();
    }

    c->cd(2);
    h2->SetMarkerStyle(20);
    h2->SetMarkerColor(kBlue);
    h2->SetTitle("Data 2 with fit; x; counts");
    h2->Draw("E");

    TF1 *f2 = new TF1("f2","[0]",500,600);
    f2->SetParameter(0,par[0]);
    f2->SetLineColor(kRed);
    f2->Draw("SAME");

    {
        TLegend *leg2 = new TLegend(0.65,0.75,0.9,0.9);
        leg2->AddEntry(h2,"Data 2","p");
        leg2->AddEntry(f2,"Const fit","l");
        leg2->Draw();
    }

    c->Update();
    c->SaveAs("fit_results.png");
}
