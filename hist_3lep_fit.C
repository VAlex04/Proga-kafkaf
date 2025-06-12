//==============================================================
//  hist_3lep_fit.C   —  фит распределения hX из tau3body.root
//==============================================================
void hist_3lep_fit(const char* fname="tau3body.root")
{
    // 1) Открываем файл и достаём hX
    TFile f(fname,"READ");
    TH1D* h = dynamic_cast<TH1D*>(f.Get<TH1D>("hX"));
    if(!h){ std::cerr<<"hX not found in "<<fname<<"\n"; return; }

    // 2) Считаем x0 = m_l/Emax
    const double m_tau = 1.77686;        // GeV
    const double m_l   = 0.105658;       // GeV (muon)
    const double Emax  = m_tau/2*(1 + m_l*m_l/(m_tau*m_tau));
    const double x0    = m_l/Emax;

    // 3) Настройки stats и fit
    gStyle->SetOptStat(1111);  // Name, Entries, Mean, RMS
    gStyle->SetOptFit (1111);  // χ²/ndf, Prob, p0±σ, p1±σ

    // 4) Canvas + гистограмма
    TCanvas* c = new TCanvas("c3","3ℓ τ-spectrum",700,500);
    h->SetMinimum(0);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitle("x = 2E_{ℓ}/m_{τ}");
    h->GetYaxis()->SetTitle("Events");
    h->Draw("HIST");

    // 5) Теоретическая функция Мишеля с ρ и η
    //    F(x) = N * sqrt(x^2 - x0^2)
    //           * [ x(1-x) + (2/9)*ρ*(4x^2-3x-x0^2) + η*x0*(1-x) ]
    TString form = TString::Format(
        "[0]*sqrt(x*x - %f*%f)"
        "*( x*(1-x) + 2./9*[1]*(4*x*x - 3*x - %f*%f) + [2]*%f*(1-x) )",
        x0,x0, x0,x0, x0
    );
    TF1* fitf = new TF1("fitf", form, x0, 1.0);
    fitf->SetParName(0,"N");
    fitf->SetParName(1,"rho");
    fitf->SetParName(2,"eta");
    fitf->SetParameter(0, h->Integral());
    fitf->SetParameter(1, 0.75);  // старт ρ
    fitf->SetParameter(2, 0.01);  // старт η
    fitf->SetLineColor(kRed);
    fitf->SetLineWidth(2);

    // 6) Фит без перерисовки гистограммы
    h->Fit(fitf,"R0");
    fitf->Draw("same");
    gPad->Update();

    // 7) Подвинуть stats, чтобы не заслонял кривую
    if(auto st=(TPaveStats*)h->FindObject("stats")){
        st->SetX1NDC(0.60); st->SetX2NDC(0.88);
        st->SetY1NDC(0.65); st->SetY2NDC(0.88);
        gPad->Modified(); gPad->Update();
    }

    // 8) Вывод в консоль
    std::cout << "\n=== 3-lepton Michel fit =========================\n"
              << "  rho      = " << fitf->GetParameter(1)
              << " ± "        << fitf->GetParError(1) << "\n"
              << "  eta      = " << fitf->GetParameter(2)
              << " ± "        << fitf->GetParError(2) << "\n"
              << "  chi2/ndf = " << fitf->GetChisquare()/fitf->GetNDF()
              << " / "        << fitf->GetNDF() << "\n"
              << "  Prob     = " << fitf->GetProb() << "\n"
              << "===============================================\n";

    // 9) Сохраняем PDF
    c->SaveAs("fig_3lepSpectrum_fit.pdf");
}
//==============================================================
