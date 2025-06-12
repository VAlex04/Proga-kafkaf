#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TCanvas.h>

#include <cmath>
#include <iostream>
using namespace std;

// ---------- константы ----------
constexpr double PI   = M_PI;
constexpr double mTau = 1.77686;      // GeV
constexpr double mMu  = 0.105658;     // GeV
constexpr double mNu  = 0.0;          // GeV  (масса обеих ν)
constexpr int    Nevt = 100000;       // принятых событий

// ---------- util: λ(a,b,c) ----------
inline double lambda(double a,double b,double c)
{
    return a*a + b*b + c*c - 2.0*(a*b + a*c + b*c);
}

// ---------- dΓ события ----------
double DifferentialWidth(const TLorentzVector &tau,
                         const TLorentzVector &nuTau,
                         const TLorentzVector &mu,
                         const TLorentzVector &nuMu)
{
    constexpr double GF = 1.1663787e-5;               // GeV⁻²
    const double twoPI5 = std::pow(2.0*PI,5);         // (2π)^5

    // |M|² для чистого V-A  (см. Jadach–Waś 1991)
    const double Msq = 64.0 * (tau * nuMu) * (mu * nuTau);

    // Инвариантная масса промежуточного объекта X
    const double M2 = (mu + nuMu).M();

    // --- импульсы в вершинах через общую λ-формулу ---
    // τ → X + ντ
    const double pX  = std::sqrt(lambda(mTau*mTau, M2*M2, mNu*mNu)) / (2.0*mTau);
    // X → μ + νμ
    const double pMu = std::sqrt(lambda(M2*M2,  mMu*mMu, mNu*mNu))   / (2.0*M2);

    // Якобиан фазового объёма (безразмерный)
    const double jac = pX * pMu / (4.0 * mTau * M2);

    // Общий коэффициент
    const double pref = GF*GF / ( twoPI5 * 2.0 * mTau );

    return pref * Msq * jac;          // GeV
}

// ================================================================
//  главный код
// ================================================================
void generator_lepton()
{
    TRandom3 rng(0);

    // --- выходной файл, дерево, гистограмма ----------------------
    TFile *fout = TFile::Open("tau3body.root","RECREATE");
    TTree *T    = new TTree("T","tau->mu nu nu");

    TLorentzVector tauLV, muLV, nuMuLV, nuTauLV;
    T->Branch("tau",    &tauLV);
    T->Branch("mu",     &muLV);
    T->Branch("nu_mu",  &nuMuLV);
    T->Branch("nu_tau", &nuTauLV);

    TH1D *hX = new TH1D("hX",
        "Lepton spectrum; x = 2E_{#mu}/m_{#tau}; Events", 100, 0.0, 1.0);

    // ---------- 1) поиск Wmax ----------
    double Wmax = 0.0;
    for(int i=0;i<5000;++i){
        // τ -> X + ντ  (X масса M2)
        double x5 = rng.Uniform();
        double M2 = std::sqrt( mMu*mMu +
                      (mTau*mTau - mMu*mMu)*x5 );         // равномерно по фаз. объёму
        double pX = std::sqrt(lambda(mTau*mTau, M2*M2, mNu*mNu)) / (2.0*mTau);

        double c1  = rng.Uniform(-1,1);
        double s1  = std::sqrt(1-c1*c1);
        double p1  = rng.Uniform(0,2*PI);

        tauLV.SetPxPyPzE(0,0,0,mTau);
        nuTauLV.SetPxPyPzE(pX*s1*std::cos(p1), pX*s1*std::sin(p1), pX*c1, pX);

        TLorentzVector XLV = tauLV - nuTauLV;              // 4-импульс X

        double E_mu = (M2*M2 + mMu*mMu - mNu*mNu)/(2.0*M2);
        double pMu  = std::sqrt(lambda(M2*M2, mMu*mMu, mNu*mNu)) / (2.0*M2);

        double c2 = rng.Uniform(-1,1);
        double s2 = std::sqrt(1-c2*c2);
        double p2 = rng.Uniform(0,2*PI);

        TLorentzVector muX ( pMu*s2*std::cos(p2),  pMu*s2*std::sin(p2),  pMu*c2, E_mu);
        TLorentzVector nuMuX(-muX.Px(),           -muX.Py(),            -muX.Pz(),
                             M2-E_mu);             // mν=0, энергия = |p|

        muLV    = muX;
        nuMuLV  = nuMuX;
        muLV   .Boost(XLV.BoostVector());
        nuMuLV .Boost(XLV.BoostVector());

        double w = DifferentialWidth(tauLV, nuTauLV, muLV, nuMuLV);
        if(w > Wmax) Wmax = w;
    }
    Wmax *= 1.1;

    // ---------- 2) основная генерация ----------
    int        accepted = 0;
    long long  trials   = 0;
    while(accepted < Nevt){
        ++trials;
        // τ → X + ντ
        double x5 = rng.Uniform();
        double M2 = std::sqrt( mMu*mMu +
                      (mTau*mTau - mMu*mMu)*x5 );
        double pX = std::sqrt(lambda(mTau*mTau, M2*M2, mNu*mNu)) / (2.0*mTau);

        double c1 = rng.Uniform(-1,1);
        double s1 = std::sqrt(1-c1*c1);
        double p1 = rng.Uniform(0,2*PI);

        tauLV.SetPxPyPzE(0,0,0,mTau);
        nuTauLV.SetPxPyPzE(pX*s1*std::cos(p1), pX*s1*std::sin(p1), pX*c1, pX);

        TLorentzVector XLV = tauLV - nuTauLV;

        // X → μ + νμ
        double E_mu = (M2*M2 + mMu*mMu - mNu*mNu)/(2.0*M2);
        double pMu  = std::sqrt(lambda(M2*M2, mMu*mMu, mNu*mNu)) / (2.0*M2);

        double c2 = rng.Uniform(-1,1);
        double s2 = std::sqrt(1-c2*c2);
        double p2 = rng.Uniform(0,2*PI);

        TLorentzVector muX ( pMu*s2*std::cos(p2),  pMu*s2*std::sin(p2),  pMu*c2, E_mu);
        TLorentzVector nuMuX(-muX.Px(),           -muX.Py(),           -muX.Pz(),
                             M2-E_mu);

        muLV   = muX;
        nuMuLV = nuMuX;
        muLV   .Boost(XLV.BoostVector());
        nuMuLV .Boost(XLV.BoostVector());

        double w = DifferentialWidth(tauLV, nuTauLV, muLV, nuMuLV);
        if(rng.Uniform()*Wmax < w){
            T->Fill();
            hX->Fill( 2.0*muLV.E()/mTau );
            ++accepted;
        }
    }
    std::cout << "Generated " << accepted << " events in " << trials << " trials\n";

    // запись
    hX->SetDirectory(nullptr);
    T->Write();  hX->Write();  fout->Close();

    // быстрый просмотр
    TCanvas *c = new TCanvas("c","Lepton spectrum",800,600);
    hX->Draw("hist");
}
