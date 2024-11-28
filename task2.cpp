#include <TCanvas.h>
#include <TLegend.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <string>

void readDataFromFileToHist2D(std::string filename, TH2D *hist)
{
    std::ifstream ifile(filename, std::ios::in);
    if (ifile.is_open())
    {
        Double_t x, y;
        while (ifile >> x >> y)
            hist->Fill(x, y);
    } else 
    {
        std::cerr << "Возникли проблемы при попытке прочтения файла " 
                  << filename << "." << std::endl;
        exit(-1);
    }
}

Int_t task2()
{
    constexpr Double_t nbinX{50}, nbinY{50};
    constexpr Double_t x_min{50}, x_max{200};
    constexpr Double_t y_min{50}, y_max{200};


    TH2D *experiment_hist = new TH2D("experiment_hist",
                                     "Experiment (style = surf3)",
                                     nbinX, x_min, x_max,
                                     nbinY, y_min, y_max);
    
    TH2D *simulation_hist = new TH2D("simulation_hist",
                                     "Simulation (sytle = colz)",
                                     nbinX, x_min, x_max,
                                     nbinY, y_min, y_max);


    readDataFromFileToHist2D("fMC.dat",   simulation_hist);
    readDataFromFileToHist2D("fdata.dat", experiment_hist);

    TH1D *experiment_hist_px = experiment_hist->ProjectionX();
    TH1D *simulation_hist_px = simulation_hist->ProjectionX();

    TCanvas *canvas = new TCanvas("canvas", "Experiment vs Simulation", 1200, 900);
    canvas->Divide(2, 2);
    gStyle->SetOptStat(11111);

    canvas->cd(1);
    
    experiment_hist->SetXTitle("m_{#gamma_{1}#gamma_{2}}, MeV/c^{2}");
    experiment_hist->SetYTitle("m_{#gamma_{3}#gamma_{4}}, MeV/c^{2}");

    experiment_hist->SetTitleSize(0.04, "X");
    experiment_hist->SetTitleSize(0.04, "Y");

    experiment_hist->SetTitleOffset(1.8, "X");
    experiment_hist->SetTitleOffset(1.8, "Y");

    experiment_hist->SetNdivisions(5, "X");
    experiment_hist->SetNdivisions(5, "Y");
    experiment_hist->SetNdivisions(10, "Z");

    experiment_hist->SetLineWidth(2);

    experiment_hist->Draw("SURF3");


    canvas->cd(2);

    simulation_hist->SetXTitle("m_{#gamma_{1}#gamma_{2}}, MeV/c^{2}");
    simulation_hist->SetYTitle("m_{#gamma_{3}#gamma_{4}}, MeV/c^{2}");

    simulation_hist->SetTitleSize(0.04, "X");
    simulation_hist->SetTitleSize(0.04, "Y");

    simulation_hist->SetTitleOffset(1.0, "X");
    simulation_hist->SetTitleOffset(0.9, "Y");

    simulation_hist->SetNdivisions(5, "X");
    simulation_hist->SetNdivisions(5, "Y");

    simulation_hist->Draw("COLZ");


    canvas->cd(3);

    Int_t simulation_max_bin_content = 
        simulation_hist_px->GetBinContent(simulation_hist_px->GetMaximumBin());
    Int_t experiment_max_bin_content = 
        experiment_hist_px->GetBinContent(experiment_hist_px->GetMaximumBin());
    Double_t scale3 = (Double_t)experiment_max_bin_content / simulation_max_bin_content;

    simulation_hist_px->SetTitle("Experiment vs Simulation");
    simulation_hist_px->SetYTitle("yields");
    simulation_hist_px->SetXTitle("m_{#gamma_{1}#gamma_{2}}, MeV/c^{2}");
    simulation_hist_px->Scale(scale3);
    simulation_hist_px->SetLineColor(kRed);
    simulation_hist_px->SetTitleOffset(1.1, "X");
    simulation_hist_px->SetTitleOffset(0.9, "Y");
    simulation_hist_px->SetNdivisions(5, "X");
    simulation_hist_px->SetLineWidth(2);
    experiment_hist_px->SetLineWidth(2);

    simulation_hist_px->Draw("HIST");
    experiment_hist_px->Draw("HIST SAME");

    canvas->cd(4);
    TH1D *background_hist = (TH1D*)experiment_hist_px->Clone("background_hist");
    background_hist->Add(simulation_hist_px, -1);

    TString sum = ""; sum.Form("%.2f", (Double_t)background_hist->GetSum());
    TString title = "Estimated background, N = GetSum() = " + sum;

    background_hist->SetTitle(title);
    background_hist->SetYTitle("yields");
    background_hist->SetXTitle("m_{#gamma_{1}#gamma_{2}}, MeV/c^{2}");
    background_hist->SetTitleOffset(1.1, "X"); 
    background_hist->SetTitleOffset(0.9, "Y"); 
    background_hist->SetNdivisions(5, "X");
    background_hist->Draw();

    canvas->SaveAs("Experiment_vs_Simulation.png");

    return 0;
}
