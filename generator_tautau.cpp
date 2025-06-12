#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include <string>
#include <cstdlib>

// ---------- физконстанты ----------
static const double mTau  = 1.77686;   // GeV
static const double Ebeam = 5.29;      // GeV  ( √s /2 )
static const double pTau  = TMath::Sqrt(Ebeam*Ebeam - mTau*mTau);
static const double E_tau = TMath::Sqrt(pTau*pTau + mTau*mTau);
static const double gammaTau = E_tau / mTau;
static const double betaTau  = pTau  / E_tau;
static const double alphaEM  = 1./137.035999;
static const double kNorm   = (alphaEM*alphaEM)/(64.*E_tau*E_tau)*betaTau; // общий множитель

struct PolVec { double x{}, y{}, z{}; };

// ---------- вспом. функции для матрицы D ----------
static inline double D0(double ct){
    double st2 = 1. - ct*ct;
    return 1. + ct*ct + st2/(gammaTau*gammaTau);
}
static inline double Dij(int i,int j,double ct){
    double st2 = 1.-ct*ct;
    double b2 = betaTau*betaTau;
    double s2t = 2.*TMath::Sqrt(st2)*ct; // sin2θ = 2 sinθ cosθ
    if(i==0&&j==0) return (1.+1./(gammaTau*gammaTau))*st2;
    if(i==1&&j==1) return -b2*st2;
    if(i==2&&j==2) return 1.+ct*ct - st2/(gammaTau*gammaTau);
    if(i==0&&j==2) return s2t/gammaTau;
    if(i==2&&j==0) return s2t/gammaTau;
    return 0.;
}

// ---------- полный вес ----------
double Weight(double ct, const PolVec &xiM, const PolVec &xiP)
{
    double w = D0(ct);
    const double xm[3]={xiM.x,xiM.y,xiM.z};
    const double xp[3]={xiP.x,xiP.y,xiP.z};
    for(int i=0;i<3;++i)
        for(int j=0;j<3;++j)
            w += Dij(i,j,ct)*xm[i]*xp[j];
    return kNorm * w;             // ≥0 для любых ξ
}

// ---------- генерация τ‑пар ----------
void makeSample(Long64_t N, int spinMode, const char *outName)
{
    TFile  fout(outName,"RECREATE");
    TTree  tree("tree","e+e- -> tau+tau-");
    TH1D  hCos("hCos","cos#theta;cos#theta;Events",50,-1,1);
    hCos.SetDirectory(&fout);            // keep histogram in the output file

    // ветки дерева
    Float_t pxM,pyM,pzM,eM, pxP,pyP,pzP,eP;
    tree.Branch("pxMinus",&pxM,"pxMinus/F");
    tree.Branch("pyMinus",&pyM,"pyMinus/F");
    tree.Branch("pzMinus",&pzM,"pzMinus/F");
    tree.Branch("Eminus", &eM ,"Eminus/F");
    tree.Branch("pxPlus", &pxP,"pxPlus/F");
    tree.Branch("pyPlus", &pyP,"pyPlus/F");
    tree.Branch("pzPlus", &pzP,"pzPlus/F");
    tree.Branch("Eplus",  &eP ,"Eplus/F");

    TRandom3 rng(0);

    // предварительный поиск wMax
    double wMax=0.;
    for(int i=0;i<1000;++i){
        double ct = -1.+2.*i/999.;
        PolVec z0{0,0,0};
        double w = Weight(ct,z0,z0);
        if(w>wMax) wMax=w;
    }
    // небольшое увеличение запаса
    wMax*=1.2;

@@ -111,40 +113,52 @@ void makeSample(Long64_t N, int spinMode, const char *outName)
        ++accepted;
    }

    tree.Write();
    hCos.Write();
    fout.Close();
    std::cout<<"Generated "<<accepted<<" events → "<<outName<<" (spinMode="<<spinMode<<")\n";
}

// ---------- интерфейс для ROOT ----------
void plotCos(const char *file = "tauPolar.root")
{
    TFile f(file,"READ");
    TH1D *h = dynamic_cast<TH1D*>(f.Get("hCos"));
    if(!h){ std::cerr<<"hCos not found in "<<file<<"\n"; return; }

    // функция N·(1+x²)
    TF1 *f1 = new TF1("f1","[0]*(1 + x*x)",-1.,1.);
    f1->SetParameter(0,h->Integral());
    h->Fit(f1,"RQL");          // тихий (Q) лайклихуд-фит (L) в заданном (R) диапазоне

    TCanvas *c = new TCanvas("cCos","cos theta",600,450);
    h->SetMarkerStyle(20);
    h->Draw("E1");
    f1->Draw("same");
    c->SaveAs("fig_cosTheta.pdf");
    c->SaveAs("изображение.png");

    std::cout<<"\nFit result:"
             <<"\n  N          = "<<f1->GetParameter(0)
             <<"\n  chi2/ndf   = "<<f1->GetChisquare()/f1->GetNDF()
             <<std::endl;
}

void generator_tautau(const char *mode="make", Long64_t N=300000, int spinMode=0, const char *fname="tauPolar.root")
{
    std::string m(mode);
    if(m=="make")      makeSample(N,spinMode,fname);
    else if(m=="plot") plotCos(fname);
    else std::cerr<<"Unknown mode: "<<mode<<" (use make/plot)\n";
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
#endif
