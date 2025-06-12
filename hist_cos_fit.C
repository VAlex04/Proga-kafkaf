//==============================================================
//  hist_cos_fit.C      (ROOT 6.x, C++11)
//  usage:  root -l -q hist_cos_fit.C
//==============================================================
void hist_cos_fit(const char* fname = "tauPolar.root")
{
    // ---------- открываем файл и гистограмму ----------
    TFile f(fname,"READ");
    TH1D* h = dynamic_cast<TH1D*>(f.Get<TH1D>("hCos"));
    if(!h){ std::cerr<<"hCos not found in "<<fname<<"\n"; return; }

    // ---------- глобальные опции ----------
    gStyle->SetOptStat(1111);   // Name, Entries, Mean, RMS
    gStyle->SetOptFit (1111);   // χ²/ndf, Prob, p0±σ

    // ---------- создаём канвас ----------
    TCanvas* c = new TCanvas("cCos","cos#theta",700,500);

    // ---------- рисуем гистограмму ----------
    h->SetMinimum(0);           // ось Y начинается с 0
    h->SetLineWidth(2);
    h->Draw("HIST");            // ступеньки без маркеров

    // ---------- функция фита ----------
    TF1* f1 = new TF1("f1","[0]*(1 + x*x)",-1.,1.);
    f1->SetParameter(0,h->Integral());   // начальное N
    f1->SetLineColor(kRed);
    f1->SetLineWidth(2);

    // ---------- фитируем (R – диапазон функции, Q – тихо) ----------
    h->Fit(f1,"RQL");

    // ---------- возвращаем ступеньки поверх возможной перерисовки ----------
    h->Draw("HIST SAME");
    gPad->Update();

    // ---------- при желании двигаем stats-блок ----------
    if(auto st = (TPaveStats*)h->FindObject("stats")){
        st->SetX1NDC(0.65); st->SetX2NDC(0.88);
        st->SetY1NDC(0.60); st->SetY2NDC(0.84);
    }
    gPad->Modified(); gPad->Update();

    // ---------- выводим параметры в консоль ----------
    std::cout << "\n=== Fit result ===============================\n"
              << "  N         = " << f1->GetParameter(0)
              << "  ± "          << f1->GetParError(0)     << '\n'
              << "  chi2/ndf  = " << f1->GetChisquare()/f1->GetNDF()
              << " / "            << f1->GetNDF()          << '\n'
              << "  Prob      = " << f1->GetProb()         << '\n'
              << "=============================================\n";

    // ---------- сохраняем картинку ----------
    c->SaveAs("fig_cosTheta.pdf");
}
//==============================================================
