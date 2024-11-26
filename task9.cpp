#include <iostream>
#include <vector>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <TDirectory.h>

void task9() {
    // Open the input file from Task 7
    TFile *inputFile = new TFile("newroot.root", "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cout << "Error opening file newroot.root" << std::endl;
        return;
    }

    // Retrieve the tree MyTree
    TTree *tree = (TTree*)inputFile->Get("MyTree");
    if (!tree) {
        std::cout << "Error retrieving tree MyTree" << std::endl;
        inputFile->Close();
        return;
    }

    // Set branch addresses
    const Int_t maxph = 100; // Maximum number of photons per event
    Int_t nph;
    Float_t eph[maxph];
    Float_t thetaph[maxph];
    Float_t phiph[maxph];

    tree->SetBranchAddress("nph", &nph);
    tree->SetBranchAddress("eph", eph);
    tree->SetBranchAddress("thetaph", thetaph);
    tree->SetBranchAddress("phiph", phiph);

    // Create histograms
    TH1F *hAzimuthal = new TH1F("hAzimuthal", "Number of #pi^{0} Candidates vs Azimuthal Angle;Azimuthal Angle (degrees);Number of #pi^{0} Candidates", 36, 0, 360);
    TH1F *hPolar = new TH1F("hPolar", "Number of #pi^{0} Candidates vs Polar Angle;Polar Angle (degrees);Number of #pi^{0} Candidates", 18, 0, 180);
    TH1F *hNumPions = new TH1F("hNumPions", "Distribution of Number of #pi^{0} Candidates per Event;Number of #pi^{0} Candidates;Events", 10, 0, 10);

    Int_t nEntries = tree->GetEntries();

    // Vector for variance calculation
    std::vector<Int_t> numPionsPerEvent;

    // Loop over events
    for (Int_t entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        Int_t numPionsInEvent = 0;

        // Loop over photon pairs
        for (Int_t i = 0; i < nph; ++i) {
            for (Int_t j = i + 1; j < nph; ++j) {
                // Photon energies and angles
                Float_t E1 = eph[i];
                Float_t theta1 = thetaph[i] * TMath::DegToRad();
                Float_t phi1 = phiph[i] * TMath::DegToRad();

                Float_t E2 = eph[j];
                Float_t theta2 = thetaph[j] * TMath::DegToRad();
                Float_t phi2 = phiph[j] * TMath::DegToRad();

                // Photon momentum components
                Float_t px1 = E1 * sin(theta1) * cos(phi1);
                Float_t py1 = E1 * sin(theta1) * sin(phi1);
                Float_t pz1 = E1 * cos(theta1);

                Float_t px2 = E2 * sin(theta2) * cos(phi2);
                Float_t py2 = E2 * sin(theta2) * sin(phi2);
                Float_t pz2 = E2 * cos(theta2);

                // Total energy and momentum
                Float_t E_sum = E1 + E2;
                Float_t px_sum = px1 + px2;
                Float_t py_sum = py1 + py2;
                Float_t pz_sum = pz1 + pz2;

                // Invariant mass of photon pair
                Float_t mass2 = E_sum * E_sum - (px_sum * px_sum + py_sum * py_sum + pz_sum * pz_sum);
                if (mass2 < 0) mass2 = 0;
                Float_t mass = sqrt(mass2);

                // Check if the pair is a candidate for π^0
                if (mass >= 0.1 && mass <= 0.2) {
                    numPionsInEvent++;

                    // Calculate angles of the π^0 candidate
                    Float_t p_sum = sqrt(px_sum * px_sum + py_sum * py_sum + pz_sum * pz_sum);
                    Float_t theta_pi0 = acos(pz_sum / p_sum) * TMath::RadToDeg();
                    Float_t phi_pi0 = atan2(py_sum, px_sum) * TMath::RadToDeg();
                    if (phi_pi0 < 0) phi_pi0 += 360;

                    hAzimuthal->Fill(phi_pi0);
                    hPolar->Fill(theta_pi0);
                }
            }
        }

        hNumPions->Fill(numPionsInEvent);
        numPionsPerEvent.push_back(numPionsInEvent);
    }

    // Calculate variance
    Double_t sum = 0;
    Double_t sumSq = 0;
    Int_t nEvents = numPionsPerEvent.size();

    for (Int_t i = 0; i < nEvents; ++i) {
        sum += numPionsPerEvent[i];
        sumSq += numPionsPerEvent[i] * numPionsPerEvent[i];
    }

    Double_t mean = sum / nEvents;
    Double_t variance = sumSq / nEvents - mean * mean;

    std::cout << "Distribution of number of π^0 candidates per event:" << std::endl;
    std::cout << "Mean = " << mean << std::endl;
    std::cout << "Variance = " << variance << std::endl;

    // Draw histograms with errors
    TCanvas *c1 = new TCanvas("c1", "Azimuthal Angle", 800, 600);
    hAzimuthal->SetMarkerStyle(21); // Square markers
    hAzimuthal->Draw("E1 P");

    // Add legend
    TLegend *leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg1->AddEntry(hAzimuthal, "#pi^{0} Candidates", "lep");
    leg1->Draw();

    c1->SaveAs("AzimuthalAngle.png");

    TCanvas *c2 = new TCanvas("c2", "Polar Angle", 800, 600);
    hPolar->SetMarkerStyle(21); // Square markers
    hPolar->Draw("E1 P");

    // Add legend
    TLegend *leg2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg2->AddEntry(hPolar, "#pi^{0} Candidates", "lep");
    leg2->Draw();

    c2->SaveAs("PolarAngle.png");

    // Before closing files, set histograms' directory to memory
    hAzimuthal->SetDirectory(0);
    hPolar->SetDirectory(0);
    hNumPions->SetDirectory(0);

    // Save histograms in a new directory in results.root
    TFile *outputFile = new TFile("results.root", "UPDATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cout << "Error opening file results.root" << std::endl;
        inputFile->Close();
        return;
    }

    // Check if the directory already exists
    TDirectory *dir = outputFile->GetDirectory("Task9Results");
    if (!dir) {
        dir = outputFile->mkdir("Task9Results");
    }
    outputFile->cd("Task9Results");

    hAzimuthal->Write("", TObject::kOverwrite);
    hPolar->Write("", TObject::kOverwrite);
    hNumPions->Write("", TObject::kOverwrite);

    // Close files properly
    //outputFile->Close();
    //inputFile->Close();

    // Keep canvases accessible after macro execution
    c1->Draw();
    c2->Draw();
}
