#include <TFile.h>
#include <TTree.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH1D.h>
#include <TF1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TCanvas.h>
#include <TColor.h>
#include <iostream>
#include <iostream>
#include <string>
#include <string>
#include <cstdlib>


// ---------- физконстанты ----------
// ---------- физконстанты ----------
static const double mTau  = 1.77686;   // GeV
static const double mTau  = 1.77686;   // GeV
static const double Ebeam = 5.29;      // GeV  ( √s /2 )
static const double Ebeam = 5.29;      // GeV  ( √s /2 )
static const double pTau  = TMath::Sqrt(Ebeam*Ebeam - mTau*mTau);
static const double pTau  = TMath::Sqrt(Ebeam*Ebeam - mTau*mTau);
static const double E_tau = TMath::Sqrt(pTau*pTau + mTau*mTau);
static const double E_tau = TMath::Sqrt(pTau*pTau + mTau*mTau);
static const double gammaTau = E_tau / mTau;
static const double gammaTau = E_tau / mTau;
static const double betaTau  = pTau  / E_tau;
static const double betaTau  = pTau  / E_tau;
static const double alphaEM  = 1./137.035999;
static const double alphaEM  = 1./137.035999;
static const double kNorm   = (alphaEM*alphaEM)/(64.*E_tau*E_tau)*betaTau; // общий множитель
static const double kNorm   = (alphaEM*alphaEM)/(64.*E_tau*E_tau)*betaTau; // общий множитель


struct PolVec { double x{}, y{}, z{}; };
struct PolVec { double x{}, y{}, z{}; };


// ---------- вспом. функции для матрицы D ----------
// ---------- вспом. функции для матрицы D ----------
static inline double D0(double ct){
static inline double D0(double ct){
    double st2 = 1. - ct*ct;
    double st2 = 1. - ct*ct;
    return 1. + ct*ct + st2/(gammaTau*gammaTau);
    return 1. + ct*ct + st2/(gammaTau*gammaTau);
}
}
static inline double Dij(int i,int j,double ct){
static inline double Dij(int i,int j,double ct){
    double st2 = 1.-ct*ct;
    double st2 = 1.-ct*ct;
    double b2 = betaTau*betaTau;
    double b2 = betaTau*betaTau;
    double s2t = 2.*TMath::Sqrt(st2)*ct; // sin2θ = 2 sinθ cosθ
    double s2t = 2.*TMath::Sqrt(st2)*ct; // sin2θ = 2 sinθ cosθ
    if(i==0&&j==0) return (1.+1./(gammaTau*gammaTau))*st2;
    if(i==0&&j==0) return (1.+1./(gammaTau*gammaTau))*st2;
    if(i==1&&j==1) return -b2*st2;
    if(i==1&&j==1) return -b2*st2;
    if(i==2&&j==2) return 1.+ct*ct - st2/(gammaTau*gammaTau);
    if(i==2&&j==2) return 1.+ct*ct - st2/(gammaTau*gammaTau);
    if(i==0&&j==2) return s2t/gammaTau;
    if(i==0&&j==2) return s2t/gammaTau;
    if(i==2&&j==0) return s2t/gammaTau;
    if(i==2&&j==0) return s2t/gammaTau;
    return 0.;
    return 0.;
}
}


static void randomSpin(TRandom3 &rng, PolVec &v)
{
    double c = rng.Uniform(-1,1);
    double s = TMath::Sqrt(1 - c*c);
    double ph = rng.Uniform(0, 2*TMath::Pi());
    v.x = s*TMath::Cos(ph);
    v.y = s*TMath::Sin(ph);
    v.z = c;
}

// ---------- полный вес ----------
// ---------- полный вес ----------
double Weight(double ct, const PolVec &xiM, const PolVec &xiP)
double Weight(double ct, const PolVec &xiM, const PolVec &xiP)
{
{
    double w = D0(ct);
    double w = D0(ct);
    const double xm[3]={xiM.x,xiM.y,xiM.z};
    const double xm[3]={xiM.x,xiM.y,xiM.z};
    const double xp[3]={xiP.x,xiP.y,xiP.z};
    const double xp[3]={xiP.x,xiP.y,xiP.z};
    for(int i=0;i<3;++i)
    for(int i=0;i<3;++i)
        for(int j=0;j<3;++j)
        for(int j=0;j<3;++j)
            w += Dij(i,j,ct)*xm[i]*xp[j];
            w += Dij(i,j,ct)*xm[i]*xp[j];
    return kNorm * w;             // ≥0 для любых ξ
    return kNorm * w;             // ≥0 для любых ξ
}
}


// ---------- генерация τ‑пар ----------
// ---------- генерация τ‑пар ----------
void makeSample(Long64_t N, int spinMode, const char *outName)
void makeSample(Long64_t N, int spinMode, const char *outName)
{
{
    TFile  fout(outName,"RECREATE");
    TFile  fout(outName,"RECREATE");
    TTree  tree("tree","e+e- -> tau+tau-");
    TTree  tree("tree","e+e- -> tau+tau-");
    TH1D  hCos("hCos","cos#theta;cos#theta;Events",50,-1,1);
    TH1D  hCos("hCos","cos#theta;cos#theta;Events",50,-1,1);
    hCos.SetDirectory(&fout);            // keep histogram in the output file


    // ветки дерева
    // ветки дерева
    Float_t pxM,pyM,pzM,eM, pxP,pyP,pzP,eP;
    Float_t pxM,pyM,pzM,eM, pxP,pyP,pzP,eP;
    tree.Branch("pxMinus",&pxM,"pxMinus/F");
    tree.Branch("pxMinus",&pxM,"pxMinus/F");
    tree.Branch("pyMinus",&pyM,"pyMinus/F");
    tree.Branch("pyMinus",&pyM,"pyMinus/F");
    tree.Branch("pzMinus",&pzM,"pzMinus/F");
    tree.Branch("pzMinus",&pzM,"pzMinus/F");
    tree.Branch("Eminus", &eM ,"Eminus/F");
    tree.Branch("Eminus", &eM ,"Eminus/F");
    tree.Branch("pxPlus", &pxP,"pxPlus/F");
    tree.Branch("pxPlus", &pxP,"pxPlus/F");
    tree.Branch("pyPlus", &pyP,"pyPlus/F");
    tree.Branch("pyPlus", &pyP,"pyPlus/F");
    tree.Branch("pzPlus", &pzP,"pzPlus/F");
    tree.Branch("pzPlus", &pzP,"pzPlus/F");
    tree.Branch("Eplus",  &eP ,"Eplus/F");
    tree.Branch("Eplus",  &eP ,"Eplus/F");


    TRandom3 rng(0);
    TRandom3 rng(0);


    // предварительный поиск wMax
    // предварительный поиск wMax
    double wMax=0.;
    double wMax=0.;
    for(int i=0;i<1000;++i){
    for(int i=0;i<1000;++i){
        double ct = -1.+2.*i/999.;
        double ct = -1.+2.*i/999.;
        PolVec z0{0,0,0};
        PolVec z0{0,0,0};
        double w = Weight(ct,z0,z0);
        double w = Weight(ct,z0,z0);
        if(w>wMax) wMax=w;
        if(w>wMax) wMax=w;
    }
    }
    // небольшое увеличение запаса
    // небольшое увеличение запаса
    wMax*=1.2;
    wMax*=1.2;

    auto randomSpin=[&](PolVec &v){
        double c=rng.Uniform(-1,1); double s=TMath::Sqrt(1-c*c);
        double ph=rng.Uniform(0,2*TMath::Pi());
        v.x=s*TMath::Cos(ph); v.y=s*TMath::Sin(ph); v.z=c; };

    Long64_t accepted=0;
    Long64_t accepted=0;
    while(accepted<N){
    while(accepted<N){
        double ct  = rng.Uniform(-1,1);
        double ct  = rng.Uniform(-1,1);
        double phi = rng.Uniform(0,2*TMath::Pi());
        double phi = rng.Uniform(0,2*TMath::Pi());


        PolVec xiM,xiP;
        PolVec xiM,xiP;
        if(spinMode==0){ xiM={0,0,0}; xiP={0,0,0}; }
        if(spinMode==0){ xiM={0,0,0}; xiP={0,0,0}; }
        else if(spinMode==1){ xiM={0,0, +1}; xiP={0,0,-1}; }
        else if(spinMode==1){ xiM={0,0, +1}; xiP={0,0,-1}; }
        else { randomSpin(xiM); randomSpin(xiP);} // spinMode==2
        else { randomSpin(rng, xiM); randomSpin(rng, xiP);} // spinMode==2


        double w = Weight(ct,xiM,xiP);
        double w = Weight(ct,xiM,xiP);
        if(rng.Uniform(0,wMax)>w) continue;
        if(rng.Uniform(0,wMax)>w) continue;


        double st = TMath::Sqrt(1-ct*ct);
        double st = TMath::Sqrt(1-ct*ct);
        double px = pTau*st*TMath::Cos(phi);
        double px = pTau*st*TMath::Cos(phi);
        double py = pTau*st*TMath::Sin(phi);
        double py = pTau*st*TMath::Sin(phi);
        double pz = pTau*ct;
        double pz = pTau*ct;


        pxM=px; pyM=py; pzM=pz; eM=E_tau;
        pxM=px; pyM=py; pzM=pz; eM=E_tau;
        pxP=-px;pyP=-py;pzP=-pz;eP=E_tau;
        pxP=-px;pyP=-py;pzP=-pz;eP=E_tau;


        tree.Fill();
        tree.Fill();
        hCos.Fill(ct);
        hCos.Fill(ct);
        ++accepted;
        ++accepted;
    }
    }


    tree.Write();
    tree.Write();
    hCos.Write();
    hCos.Write();
    fout.Close();
    fout.Close();
    std::cout<<"Generated "<<accepted<<" events → "<<outName<<" (spinMode="<<spinMode<<")\n";
    std::cout<<"Generated "<<accepted<<" events → "<<outName<<" (spinMode="<<spinMode<<")\n";
}
}


// ---------- интерфейс для ROOT ----------
// ---------- функция подгонки гистограммы ----------
void plotCos(const char *file = "tauPolar.root")
void fitCos(const char *file = "tauPolar.root")
{
{
    TFile f(file,"READ");
    TFile f(file,"READ");
    TH1D *h = dynamic_cast<TH1D*>(f.Get("hCos"));
    TH1D *h = dynamic_cast<TH1D*>(f.Get("hCos"));
    if(!h){ std::cerr<<"hCos not found in "<<file<<"\n"; return; }
    if(!h){ std::cerr<<"hCos not found in "<<file<<"\n"; return; }


    // функция N·(1+x²)
    TF1 fitFunc("fitFunc","[0]*(1 + x*x)",-1.,1.);
    TF1 *f1 = new TF1("f1","[0]*(1 + x*x)",-1.,1.);
    fitFunc.SetParameter(0,h->Integral());
    f1->SetParameter(0,h->Integral());
    h->Fit(&fitFunc,"RQL");          // тихий (Q) лайклихуд-фит (L) в заданном (R) диапазоне
    h->Fit(f1,"RQL");          // тихий (Q) лайклихуд-фит (L) в заданном (R) диапазоне


    TCanvas *c = new TCanvas("cCos","cos theta",600,450);
    TCanvas c("cCos","cos theta",600,450);
    h->SetMarkerStyle(20);
    h->SetMinimum(0.);               // ось OY начинается с нуля
    h->Draw("E1");
    h->Draw("HIST");                // отображаем гистограмму ломаной линией
    f1->Draw("same");
    fitFunc.SetLineColor(kRed);
    c->SaveAs("fig_cosTheta.pdf");
    fitFunc.Draw("same");
    c.SaveAs("изображение.png");


    std::cout<<"\nFit result:"
    std::cout<<"\nFit result:"
             <<"\n  N          = "<<f1->GetParameter(0)
             <<"\n  N          = "<<fitFunc.GetParameter(0)
             <<"\n  chi2/ndf   = "<<f1->GetChisquare()/f1->GetNDF()
             <<"\n  chi2/ndf   = "<<fitFunc.GetChisquare()/fitFunc.GetNDF()
             <<std::endl;
             <<std::endl;
}
}


void generator_tautau(const char *mode="make", Long64_t N=300000, int spinMode=0, const char *fname="tauPolar.root")
void generator_tautau(const char *mode="make", Long64_t N=300000, int spinMode=0, const char *fname="tauPolar.root")
{
{
    std::string m(mode);
    std::string m(mode);
    if(m=="make")      makeSample(N,spinMode,fname);
    if(m=="make")      makeSample(N,spinMode,fname);
    else if(m=="plot") plotCos(fname);
    else if(m=="plot" || m=="fit") fitCos(fname);
    else std::cerr<<"Unknown mode: "<<mode<<" (use make/plot)\n";
    else std::cerr<<"Unknown mode: "<<mode<<" (use make/fit)\n";
}

#ifndef __CLING__
int main(int argc, char **argv)
{
    const char *mode = (argc>1) ? argv[1] : "make";
    Long64_t N = (argc>2) ? std::atoll(argv[2]) : 300000;
    int spinMode = (argc>3) ? std::atoi(argv[3]) : 0;
    const char *fname = (argc>4) ? argv[4] : "tauPolar.root";
    generator_tautau(mode,N,spinMode,fname);
    return 0;
}
}
#endif
