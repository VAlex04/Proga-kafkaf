#include <iostream>
#include <cmath>
#include <TRandom3.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TF1.h>

void task6() {
    // Проверяем и удаляем существующие объекты, чтобы избежать предупреждений
    if (gROOT->FindObject("h_theta_KS_lab")) delete gROOT->FindObject("h_theta_KS_lab");
    if (gROOT->FindObject("h_phi_KS_lab")) delete gROOT->FindObject("h_phi_KS_lab");
    if (gROOT->FindObject("h_theta_pi_lab")) delete gROOT->FindObject("h_theta_pi_lab");
    if (gROOT->FindObject("h_phi_pi_lab")) delete gROOT->FindObject("h_phi_pi_lab");
    if (gROOT->FindObject("h_theta_pi_KS")) delete gROOT->FindObject("h_theta_pi_KS");
    if (gROOT->FindObject("h_phi_pi_KS")) delete gROOT->FindObject("h_phi_pi_KS");

    // Устанавливаем стиль ROOT
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1); // Отображение результатов подгонки на графике

    // Инициализация генератора случайных чисел
    TRandom3 rand(0);

    // Константы
    const Double_t E_CM = 1020.0; // Полная энергия в системе ЦМ (МэВ)
    const Double_t m_KS = 497.611; // Масса K_S (МэВ/c^2)
    const Double_t m_pi = 139.570; // Масса пиона (МэВ/c^2)
    const Double_t c = 29.9792458; // Скорость света в см/нс
    const Double_t tau_KS = 0.0895; // Время жизни K_S в нс

    // Параметры детектора
    const Double_t detector_radius = 30.0; // Радиус детектора (см)
    const Double_t detector_length = 50.0; // Длина детектора (см)
    const Double_t min_pion_momentum = 40.0; // Минимальный импульс пиона (МэВ/c)

    // Гистограммы для результатов
    // 1. Распределения по углам для K_S в лабораторной системе
    TH1D* h_theta_KS_lab = new TH1D("h_theta_KS_lab", "", 100, 0, TMath::Pi());
    TH1D* h_phi_KS_lab = new TH1D("h_phi_KS_lab", "", 100, -TMath::Pi(), TMath::Pi());

    // 2. Распределения по углам для пионов в лабораторной системе
    TH1D* h_theta_pi_lab = new TH1D("h_theta_pi_lab", "", 100, 0, TMath::Pi());
    TH1D* h_phi_pi_lab = new TH1D("h_phi_pi_lab", "", 100, -TMath::Pi(), TMath::Pi());

    // 3. Распределения по углам для пионов в системе K_S
    TH1D* h_theta_pi_KS = new TH1D("h_theta_pi_KS", "", 100, 0, TMath::Pi());
    TH1D* h_phi_pi_KS = new TH1D("h_phi_pi_KS", "", 100, -TMath::Pi(), TMath::Pi());

    // Счетчики событий
    Int_t N_total = 100000; // Общее число сгенерированных событий
    Int_t N_registered = 0; // Число зарегистрированных событий

    // Создаем функцию распределения для θ_K
    TF1 *theta_dist = new TF1("theta_dist", "sin(x)*sin(x)*sin(x)", 0, TMath::Pi());

    // Основной цикл моделирования
    for (Int_t i = 0; i < N_total; ++i) {
        // 1. Генерация углов для K_S в лабораторной системе
        Double_t theta_K = theta_dist->GetRandom();
        Double_t phi_K = rand.Rndm() * 2 * TMath::Pi() - TMath::Pi(); // От -pi до pi

        Double_t cos_theta_K = cos(theta_K);
        Double_t sin_theta_K = sin(theta_K);

        // Заполняем гистограммы углов для K_S в лабораторной системе
        h_theta_KS_lab->Fill(theta_K);
        h_phi_KS_lab->Fill(phi_K);

        // 2. Импульс K_S в лабораторной системе
        Double_t E_K = E_CM / 2.0; // Энергия K_S
        Double_t p_K = sqrt(E_K * E_K - m_KS * m_KS); // Импульс K_S

        // Вектор импульса K_S
        Double_t p_Kx = p_K * sin_theta_K * cos(phi_K);
        Double_t p_Ky = p_K * sin_theta_K * sin(phi_K);
        Double_t p_Kz = p_K * cos_theta_K;
        TLorentzVector KS_lab(p_Kx, p_Ky, p_Kz, E_K);

        // 3. Генерация распада K_S -> pi+ pi- в системе покоя K_S
        // Генерируем изотропные направления для пионов
        Double_t cos_theta_pi = rand.Rndm() * 2 - 1;
        Double_t sin_theta_pi = sqrt(1 - cos_theta_pi * cos_theta_pi);
        Double_t theta_pi = acos(cos_theta_pi);
        Double_t phi_pi = rand.Rndm() * 2 * TMath::Pi() - TMath::Pi();

        // Заполняем гистограммы углов для пионов в системе K_S
        h_theta_pi_KS->Fill(theta_pi);
        h_phi_pi_KS->Fill(phi_pi);

        // Импульс пиона в системе покоя K_S
        Double_t E_pi = m_KS / 2.0; // Энергия пиона
        Double_t p_pi = sqrt(E_pi * E_pi - m_pi * m_pi); // Импульс пиона

        // Вектора импульсов пионов в системе K_S
        Double_t p_pion_x = p_pi * sin_theta_pi * cos(phi_pi);
        Double_t p_pion_y = p_pi * sin_theta_pi * sin(phi_pi);
        Double_t p_pion_z = p_pi * cos_theta_pi;
        TLorentzVector pi_plus(p_pion_x, p_pion_y, p_pion_z, E_pi);
        TLorentzVector pi_minus(-p_pion_x, -p_pion_y, -p_pion_z, E_pi);

        // 4. Преобразование Лоренца в лабораторную систему

        // Направление движения K_S в лабораторной системе
        TVector3 K_direction = KS_lab.Vect().Unit();

        double old_angle = pi_plus.Z() / pi_plus.P();
        // Поворачиваем импульсы пионов так, чтобы ось z совпала с направлением K_S
        pi_plus.RotateUz(K_direction);
        pi_minus.RotateUz(K_direction);
        double new_angle = pi_plus.Vect().Dot(K_direction) / pi_plus.P();
        std::cout << old_angle - new_angle << std::endl;

        // Скорость K_S в лабораторной системе
        TVector3 beta_KS = KS_lab.BoostVector();

        // Выполняем буст пионов в лабораторную систему
        pi_plus.Boost(beta_KS);
        pi_minus.Boost(beta_KS);

        // Заполняем гистограммы углов для пионов в лабораторной системе
        h_theta_pi_lab->Fill(pi_plus.Theta());
        h_phi_pi_lab->Fill(pi_plus.Phi());
        h_theta_pi_lab->Fill(pi_minus.Theta());
        h_phi_pi_lab->Fill(pi_minus.Phi());

        // 5. Симуляция длины пробега K_S
        Double_t beta_KS_mag = beta_KS.Mag();
        Double_t gamma_KS = KS_lab.Gamma();

        // Средняя длина пробега K_S в см
        Double_t lambda = beta_KS_mag * gamma_KS * c * tau_KS;

        // Генерируем длину пробега K_S
        Double_t r_decay = rand.Rndm();
        Double_t L_decay = -lambda * log(r_decay);

        // Положение распада K_S
        TVector3 decay_vertex = K_direction * L_decay;

        // 6. Проверка регистрации пионов детектором
        // Предполагаем, что пионы начинают движение из точки распада K_S

        // Траектории пионов
        TVector3 pi_plus_dir = pi_plus.Vect().Unit();
        TVector3 pi_minus_dir = pi_minus.Vect().Unit();

        // Расстояние до пересечения с цилиндрической поверхностью детектора
        Double_t t_plus = (detector_radius - decay_vertex.Perp()) / (pi_plus_dir.Perp());
        Double_t t_minus = (detector_radius - decay_vertex.Perp()) / (pi_minus_dir.Perp());

        // Положение пересечения с цилиндром
        TVector3 pi_plus_pos = decay_vertex + pi_plus_dir * t_plus;
        TVector3 pi_minus_pos = decay_vertex + pi_minus_dir * t_minus;

        // Проверяем, находятся ли точки пересечения внутри детектора по оси z
        Bool_t registered_plus = (fabs(pi_plus_pos.Z()) <= detector_length / 2.0) && (pi_plus.P() >= min_pion_momentum);
        Bool_t registered_minus = (fabs(pi_minus_pos.Z()) <= detector_length / 2.0) && (pi_minus.P() >= min_pion_momentum);

        // Если оба пиона зарегистрированы, увеличиваем счетчик
        if (registered_plus && registered_minus) {
            N_registered++;
        }
    }

    // 7. Расчет эффективности
    Double_t efficiency = (Double_t)N_registered / (Double_t)N_total;
    Double_t error = sqrt(efficiency * (1 - efficiency) / N_total);

    // Вывод результатов
    std::cout << "Total events: " << N_total << std::endl;
    std::cout << "Registered events: " << N_registered << std::endl;
    std::cout << "Reconstruction efficiency: " << efficiency * 100 << " ± " << error * 100 << " %" << std::endl;

    // 8. Отрисовка гистограмм и подгонка

    // Канвас для K_S углов в лабораторной системе
    TCanvas* c1 = new TCanvas("c1", "K_S angles in Lab System", 1200, 600);
    c1->Divide(2, 1);
    c1->cd(1);
    h_theta_KS_lab->GetXaxis()->SetTitle("#theta_{K_{S}} (rad)");
    h_theta_KS_lab->GetYaxis()->SetTitle("Number of events");
    h_theta_KS_lab->Draw();

    // Подгонка h_theta_KS_lab с функцией sin^3(θ)
    TF1 *sin3_theta_KS = new TF1("sin3_theta_KS", "[0]*pow(sin(x),3)", 0, TMath::Pi());
    sin3_theta_KS->SetParameter(0, h_theta_KS_lab->GetMaximum());
    h_theta_KS_lab->Fit(sin3_theta_KS, "R");
    sin3_theta_KS->SetLineColor(kRed);
    sin3_theta_KS->Draw("Same");

    c1->cd(2);
    h_phi_KS_lab->GetXaxis()->SetTitle("#phi_{K_{S}} (rad)");
    h_phi_KS_lab->GetYaxis()->SetTitle("Number of events");
    h_phi_KS_lab->Draw();

    // Подгонка h_phi_KS_lab с константой
    TF1 *const_phi_KS = new TF1("const_phi_KS", "[0]", -TMath::Pi(), TMath::Pi());
    const_phi_KS->SetParameter(0, h_phi_KS_lab->GetMaximum());
    h_phi_KS_lab->Fit(const_phi_KS, "R");
    const_phi_KS->SetLineColor(kRed);
    const_phi_KS->Draw("Same");

    // Канвас для углов пионов в лабораторной системе
    TCanvas* c2 = new TCanvas("c2", "Pion angles in Lab System", 1200, 600);
    c2->Divide(2, 1);
    c2->cd(1);
    h_theta_pi_lab->GetXaxis()->SetTitle("#theta_{#pi} (rad)");
    h_theta_pi_lab->GetYaxis()->SetTitle("Number of events");
    h_theta_pi_lab->Draw();

    // Подгонка h_theta_pi_lab с функцией sin(θ)
    TF1 *sin_theta_pi_lab = new TF1("sin_theta_pi_lab", "[0]*sin(x)", 0, TMath::Pi());
    sin_theta_pi_lab->SetParameter(0, h_theta_pi_lab->GetMaximum());
    h_theta_pi_lab->Fit(sin_theta_pi_lab, "R");
    sin_theta_pi_lab->SetLineColor(kRed);
    sin_theta_pi_lab->Draw("Same");

    c2->cd(2);
    h_phi_pi_lab->GetXaxis()->SetTitle("#phi_{#pi} (rad)");
    h_phi_pi_lab->GetYaxis()->SetTitle("Number of events");
    h_phi_pi_lab->Draw();

    // Подгонка h_phi_pi_lab с константой
    TF1 *const_phi_pi_lab = new TF1("const_phi_pi_lab", "[0]", -TMath::Pi(), TMath::Pi());
    const_phi_pi_lab->SetParameter(0, h_phi_pi_lab->GetMaximum());
    h_phi_pi_lab->Fit(const_phi_pi_lab, "R");
    const_phi_pi_lab->SetLineColor(kRed);
    const_phi_pi_lab->Draw("Same");

    // Канвас для углов пионов в системе покоя K_S
    TCanvas* c3 = new TCanvas("c3", "Pion angles in K_S rest frame", 1200, 600);
    c3->Divide(2, 1);
    c3->cd(1);
    h_theta_pi_KS->GetXaxis()->SetTitle("#theta_{#pi} (rad)");
    h_theta_pi_KS->GetYaxis()->SetTitle("Number of events");
    h_theta_pi_KS->Draw();

    // Подгонка h_theta_pi_KS с функцией sin(θ)
    TF1 *sin_theta_pi_KS = new TF1("sin_theta_pi_KS", "[0]*sin(x)", 0, TMath::Pi());
    sin_theta_pi_KS->SetParameter(0, h_theta_pi_KS->GetMaximum());
    h_theta_pi_KS->Fit(sin_theta_pi_KS, "R");
    sin_theta_pi_KS->SetLineColor(kRed);
    sin_theta_pi_KS->Draw("Same");

    c3->cd(2);
    h_phi_pi_KS->GetXaxis()->SetTitle("#phi_{#pi} (rad)");
    h_phi_pi_KS->GetYaxis()->SetTitle("Number of events");
    h_phi_pi_KS->Draw();

    // Подгонка h_phi_pi_KS с константой
    TF1 *const_phi_pi_KS = new TF1("const_phi_pi_KS", "[0]", -TMath::Pi(), TMath::Pi());
    const_phi_pi_KS->SetParameter(0, h_phi_pi_KS->GetMaximum());
    h_phi_pi_KS->Fit(const_phi_pi_KS, "R");
    const_phi_pi_KS->SetLineColor(kRed);
    const_phi_pi_KS->Draw("Same");

    // Сохраняем гистограммы в файлы
    c1->SaveAs("KS_angles_lab.png");
    c2->SaveAs("pion_angles_lab.png");
    c3->SaveAs("pion_angles_KS.png");
}
