#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>

void task8() {
    // Открываем файл с данными
    TFile *file = TFile::Open("newroot.root");
    if (!file || file->IsZombie()) {
        std::cout << "Failed to open file newroot.root" << std::endl;
        return;
    }

    // Получаем дерево MyTree
    TTree *tree = (TTree*)file->Get("MyTree");
    if (!tree) {
        std::cout << "Failed to find tree MyTree in the file" << std::endl;
        return;
    }

    // Объявляем переменные для чтения данных
    const Int_t maxnph = 100; // Максимальное количество фотонов в событии
    Int_t nph;
    Float_t eph[maxnph];
    Float_t thetaph[maxnph];
    Float_t phiph[maxnph];

    // Устанавливаем адреса ветвей
    tree->SetBranchAddress("nph", &nph);
    tree->SetBranchAddress("eph", eph);
    tree->SetBranchAddress("thetaph", thetaph);
    tree->SetBranchAddress("phiph", phiph);

    // Создаем гистограммы
    // Важно: гистограммы остаются такими же, но теперь мы будем заполнять h_invmass только для событий с 2 кандидатами
    TH1F *h_invmass = new TH1F("h_invmass", "Invariant Mass of #pi^{0} Candidates; M_{#gamma#gamma} [GeV/c^{2}]; Number of Candidates", 100, 0, 0.3);
    TH1F *h_angle = new TH1F("h_angle", "Angle Between Photon Pairs; Angle [rad]; Number of Pairs", 100, 0, TMath::Pi());

    // Переменная для подсчёта общего количества кандидатов (по всем событиям)
    // Это общее число найденных кандидатов во всех событиях, не обязательно используется для условия,
    // но пусть останется для статистики.
    Long64_t total_candidates = 0;

    // Цикл по событиям
    Long64_t nentries = tree->GetEntries();
    std::cout << "Total number of events: " << nentries << std::endl;

    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        // Вектор для хранения инвариантных масс кандидатов текущего события
        std::vector<Float_t> inv_masses;

        // Цикл по всем парам фотонов (комбинации без повторений)
        for (Int_t j = 0; j < nph - 1; j++) {
            for (Int_t k = j + 1; k < nph; k++) {
                // Энергии фотонов
                Float_t E1 = eph[j];
                Float_t E2 = eph[k];

                // Углы фотонов
                Float_t theta1 = thetaph[j];
                Float_t phi1 = phiph[j];

                Float_t theta2 = thetaph[k];
                Float_t phi2 = phiph[k];

                // Расчет компонент импульса фотонов
                TVector3 p1;
                p1.SetMagThetaPhi(E1, theta1, phi1);
                TVector3 p2;
                p2.SetMagThetaPhi(E2, theta2, phi2);

                // Создание четырехмерных векторов фотонов
                TLorentzVector photon1(p1, E1);
                TLorentzVector photon2(p2, E2);

                // Расчет инвариантной массы пары фотонов
                TLorentzVector pi0_candidate = photon1 + photon2;
                Float_t inv_mass = pi0_candidate.M();

                // Расчет угла между фотонами и заполнение гистограммы углов для всех пар
                Float_t angle = photon1.Angle(photon2.Vect());
                h_angle->Fill(angle);

                // Проверка условия на инвариантную массу
                if (inv_mass >= 0.1 && inv_mass <= 0.2) {
                    // Сохраняем инвариантную массу кандидата в вектор
                    inv_masses.push_back(inv_mass);
                }
            }
        }

        // Теперь, после перебора всех пар в данном событии,
        // проверяем, сколько кандидатов мы нашли
        if (inv_masses.size() == 2) {
            // Если ровно два кандидата, заполняем гистограмму инвариантной массы
            for (auto mass : inv_masses) {
                h_invmass->Fill(mass);
                total_candidates++; // считаем общее количество кандидатов во всем анализе
            }
        }
        // Если кандидатов не 2, то мы игнорируем это событие в плане инвариантной массы
    }

    // Вывод общего количества кандидатов
    std::cout << "Total number of pi0 candidates (from all events with exactly 2 candidates): " << total_candidates << std::endl;

    // Сохраняем гистограммы в файл
    TFile *outfile = new TFile("results.root", "RECREATE");
    h_invmass->Write();
    h_angle->Write();
    outfile->Close();

    std::cout << "Analysis complete. Results saved in file results.root" << std::endl;

    // Отображение гистограмм
    TCanvas *c1 = new TCanvas("c1", "Invariant Mass", 800, 600);
    h_invmass->Draw();

    TCanvas *c2 = new TCanvas("c2", "Angle Between Photons", 800, 600);
    h_angle->Draw();

    // Тест угла, чтобы убедиться, что Angle() работает корректно
    {
        TVector3 p1(1.0, 0.0, 0.0);
        TVector3 p2(0.0, 1.0, 0.0);
        TLorentzVector photon1(p1, 1.0);
        TLorentzVector photon2(p2, 1.0);

        Float_t angle = photon1.Angle(photon2.Vect());
        std::cout << "Angle test: " << angle << " rad, expected ~" << TMath::Pi()/2 << " rad" << std::endl;
    }
}
