#include <TFile.h>
 #include <TTree.h>
 #include <TRandom3.h>
 #include <TMath.h>
 #include <iostream>
 
 // --- Глобальные константы ---
 static const double mTau  = 1.77686;   // масса tau, ГэВ
 static const double Ebeam = 5.29;       // энергия e-/e+ пучка, ГэВ
 // Вычисление импульса tau: p_tau = sqrt(E_beam^2 - mTau^2)
 static const double pTau  = TMath::Sqrt(Ebeam*Ebeam - mTau*mTau);
 // Энергия tau: E_tau = sqrt(p_tau^2 + mTau^2)
 static const double E_tau = TMath::Sqrt(pTau*pTau + mTau*mTau);
 static const double pi    = TMath::Pi();
 
 // --- Прототипы функций ---
 double Weight(double cosTh, const double *xiMinus, const double *xiPlus);
 void generateTauTau(int N=100000);
 
 // --- Главная функция-макрос для ROOT ---
 void generator_tautau()
 {
    // Генерация 100k событий
    generateTauTau(100000);
 }
 
 //------------------------------------------------------------------------------
 // Функция для вычисления веса с учетом полного множителя
 // Реализует:
 //   w(θ) = [α²/(64 E_τ²)] β_τ [D₀(θ) + D_ij(θ) ξ⁻_i ξ⁺_j]
 // где D₀(θ) = 1 + cos²θ + (1/γ²) sin²θ,
 // Но я взял только D_zz ибо не знаю, все надо учитывать или только вдоль z куда они направлены: D_zz = cos²θ - sin²θ/γ².
 double Weight(double cosTh, const double *xiMinus, const double *xiPlus)
 {
    double sinTh2 = 1.0 - cosTh*cosTh;  // sin²θ = 1 - cos²θ
    double gamma  = E_tau / mTau;        // гамма-фактор τ
    // Спин-независимая часть:
    double D0 = (1.0 + cosTh*cosTh) + (sinTh2 / (gamma*gamma));
    // Модельная зависимость для D_ij: только компонента D_zz:
    double D_zz = cosTh*cosTh - (sinTh2/(gamma*gamma));
    // Извлечение z-компоненты векторов поляризации:
    double xiMz = xiMinus[2];
    double xiPz = xiPlus[2];
    double polTerm = D_zz * (xiMz * xiPz);
    double baseWeight = D0 + polTerm; // вес без нормировочного множителя
 
    // Теперь включение общего множителя:
    double beta_tau = pTau / E_tau; // β_τ = p_τ / E_τ
    double alpha_const = 1.0/137.0;   // α ≈ 1/137
    double kappa = (alpha_const * alpha_const) / (64 * E_tau * E_tau) * beta_tau;
    double w = kappa * baseWeight;
    return (w > 0.0) ? w : 0.0;  // отсечение отрицательных значений, ибо в реальности не может быть такого
 }
 
 //------------------------------------------------------------------------------
 // Генерация событий e+ e- -> tau+ tau-
 void generateTauTau(int N)
 {
    // 1) Открытие выходного ROOT-файла и создание TTree
    TFile *fout = new TFile("tauPolar.root", "RECREATE");
    TTree *tree = new TTree("tree","tau events with polarization");
 
    // 2) Объявление переменных для записи 4-векторов и поляризаций
    Float_t pxMinus, pyMinus, pzMinus, Eminus;
    Float_t pxPlus,  pyPlus,  pzPlus,  Eplus;
    Float_t polMinus[3], polPlus[3];
 
    // Создание ветви TTree для этих переменных
    tree->Branch("pxMinus",&pxMinus,"pxMinus/F");
    tree->Branch("pyMinus",&pyMinus,"pyMinus/F");
    tree->Branch("pzMinus",&pzMinus,"pzMinus/F");
    tree->Branch("Eminus", &Eminus, "Eminus/F");
 
    tree->Branch("pxPlus",&pxPlus,"pxPlus/F");
    tree->Branch("pyPlus",&pyPlus,"pyPlus/F");
    tree->Branch("pzPlus",&pzPlus,"pzPlus/F");
    tree->Branch("Eplus", &Eplus, "Eplus/F");
 
    tree->Branch("polMinus", polMinus, "polMinus[3]/F");
    tree->Branch("polPlus",  polPlus,  "polPlus[3]/F");
 
    // 3) Инициализация генератора случайных чисел TRandom3
    TRandom3 rand(0);  // seed=0: генерация на основе системного времени
 
    // 4) Задание (как пример) постоянных векторов поляризации для tau- и tau+
    //    Допустим: tau- поляризован вдоль +z, tau+ вдоль -z.
    double xiM[3] = {0.0, 0.0, +1.0};
    double xiP[3] = {0.0, 0.0, -1.0};
 
    // 5) Нахождение максимального значения веса (w_max) для метода acceptance-rejection - один из методов Монте-Карло
    double wmax = 0.0;
    for(int i=0; i<1000; i++){
       double cth  = -1.0 + 2.0*(i/999.0);
       double wval = Weight(cth, xiM, xiP);
       if(wval > wmax) wmax = wval;
    }
 
    // 6) Генерация событий до получения N принятых событий
    int nAccepted = 0;
    while(nAccepted < N)
    {
       // 6a) Случайная генерация cosTheta и phi
       double cth = rand.Uniform(-1.0,1.0);
       double phi = rand.Uniform(0, 2.0*pi);
 
       // 6b) Вычисление веса для сгенерированного угла
       double wval = Weight(cth, xiM, xiP);
 
       // 6c) Применение метода acceptance-rejection:
       double r = rand.Uniform(0, wmax);
       if(r > wval) continue;  // если случайное число больше, событие отклоняется
 
       // 6d) Если событие принято, идет вычисление 3-импульс tau- и tau+
       double sth = TMath::Sqrt(1.0 - cth*cth);
       double px  = pTau * sth * TMath::Cos(phi);
       double py  = pTau * sth * TMath::Sin(phi);
       double pz  = pTau * cth;
 
       // Присваивание компоненты 4-векторов:
       pxMinus = px;  pyMinus = py;  pzMinus = pz;   Eminus = E_tau;
       pxPlus  = -px; pyPlus  = -py; pzPlus  = -pz;  Eplus  = E_tau;
 
       // Записывание векторов поляризации (как заданы)
       polMinus[0] = xiM[0];  polMinus[1] = xiM[1];  polMinus[2] = xiM[2];
       polPlus[0]  = xiP[0];  polPlus[1]  = xiP[1];  polPlus[2]  = xiP[2];
 
       // 6e) Записывание события в TTree
       tree->Fill();
       nAccepted++;
    }
 
    // 7) Сохранение дерева и закрытие файла
    tree->Write();
    fout->Close();
 
    std::cout << "Done. Generated " << N << " events into tauPolar.root\n"
              << "Max weight was " << wmax << "\n";
 }
 
