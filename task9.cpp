#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TMath.h>
#include <iostream>
#include <vector>

void task9() {
    // Открытие файла results.root
    TFile *resultsFile = new TFile("results.root", "UPDATE");
    if (!resultsFile || resultsFile->IsZombie()) {
        std::cerr << "Ошибка открытия файла results.root!" << std::endl;
        return;
    }

    // Чтение дерева Pi0CandidatesTree
    TTree *candTree = (TTree *)resultsFile->Get("Pi0CandidatesTree");
    if (!candTree) {
        std::cerr << "Дерево Pi0CandidatesTree не найдено в results.root!" << std::endl;
        resultsFile->Close();
        return;
    }

    // Переменные для чтения данных
    Float_t theta_pi0, phi_pi0;
    candTree->SetBranchAddress("theta_pi0", &theta_pi0);
    candTree->SetBranchAddress("phi_pi0", &phi_pi0);

    // Настройка бинов для графиков
    const int nbins_theta = 18; // 0-180 градусов, шаг 10 градусов
    const int nbins_phi = 36;   // 0-360 градусов, шаг 10 градусов
    std::vector<double> thetaCenters(nbins_theta), phiCenters(nbins_phi);
    std::vector<double> thetaCounts(nbins_theta, 0), phiCounts(nbins_phi, 0);
    std::vector<double> thetaErrors(nbins_theta, 0), phiErrors(nbins_phi, 0);

    for (int i = 0; i < nbins_theta; ++i)
        thetaCenters[i] = 10.0 * i + 5.0; // Центр бина для theta
    for (int i = 0; i < nbins_phi; ++i)
        phiCenters[i] = 10.0 * i + 5.0;   // Центр бина для phi

    // Подсчет количества кандидатов по углам
    Long64_t nEntries = candTree->GetEntries();
    // Переменные для подсчета общего количества событий в графиках
    Long64_t totalThetaEvents = 0, totalPhiEvents = 0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        candTree->GetEntry(i);

        // Нормализация угла phi в диапазоне [0, 2*pi)
        if (phi_pi0 < 0) phi_pi0 += 2 * TMath::Pi();

        int thetaBin = static_cast<int>(theta_pi0 * 180.0 / TMath::Pi() / 10); // Перевод в градусы
        int phiBin = static_cast<int>(phi_pi0 * 180.0 / TMath::Pi() / 10);

        if (thetaBin >= 0 && thetaBin < nbins_theta) thetaCounts[thetaBin]++;
        if (phiBin >= 0 && phiBin < nbins_phi) phiCounts[phiBin]++;
    }

    // Подсчет общего количества событий после распределения
    for (int i = 0; i < nbins_theta; ++i) totalThetaEvents += thetaCounts[i];
    for (int i = 0; i < nbins_phi; ++i) totalPhiEvents += phiCounts[i];

    // Расчет ошибок (кв. корень из количества событий)
    for (int i = 0; i < nbins_theta; ++i)
        thetaErrors[i] = TMath::Sqrt(thetaCounts[i]);
    for (int i = 0; i < nbins_phi; ++i)
        phiErrors[i] = TMath::Sqrt(phiCounts[i]);

    // Создание графиков
    TCanvas *canvas = new TCanvas("canvas_task9", "Task 9 Results", 1200, 600);
    canvas->Divide(2, 1);

    TGraphErrors *thetaGraph = new TGraphErrors(nbins_theta, thetaCenters.data(), thetaCounts.data(), nullptr, thetaErrors.data());
    thetaGraph->SetTitle("#pi^{0} Candidates vs Polar Angle;Polar Angle (deg);# of Candidates");
    thetaGraph->SetMarkerStyle(20);
    thetaGraph->SetMarkerColor(kBlue);

    TGraphErrors *phiGraph = new TGraphErrors(nbins_phi, phiCenters.data(), phiCounts.data(), nullptr, phiErrors.data());
    phiGraph->SetTitle("#pi^{0} Candidates vs Azimuthal Angle;Azimuthal Angle (deg);# of Candidates");
    phiGraph->SetMarkerStyle(20);
    phiGraph->SetMarkerColor(kRed);

    // Рисование графиков
    canvas->cd(1);
    thetaGraph->Draw("AP");

    canvas->cd(2);
    phiGraph->Draw("AP");

    // Сохранение результатов
    TDirectory *task9Dir = (TDirectory *)resultsFile->Get("Task9");
    if (!task9Dir) {
        task9Dir = resultsFile->mkdir("Task9");
    }
    task9Dir->cd();
    thetaGraph->Write("Graph_theta");
    phiGraph->Write("Graph_phi");
    canvas->Write("Task9_Canvas");

    std::cout << "Total events in theta graph: " << totalThetaEvents << std::endl;
    std::cout << "Total events in phi graph: " << totalPhiEvents << std::endl;

    // Закрытие файла
    resultsFile->Close();
    std::cout << "Task 9 analysis complete. Results saved in results.root under Task9 directory." << std::endl;
}
