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
    // Open the data file
    TFile *file = TFile::Open("newroot.root");
    if (!file || file->IsZombie()) {
        std::cout << "Failed to open file newroot.root" << std::endl;
        return;
    }

    // Get the tree MyTree
    TTree *tree = (TTree*)file->Get("MyTree");
    if (!tree) {
        std::cout << "Failed to find tree MyTree in the file" << std::endl;
        return;
    }

    // Declare variables for reading data
    const Int_t maxnph = 100; // Maximum number of photons in an event
    Int_t nph;
    Float_t eph[maxnph];
    Float_t thetaph[maxnph];
    Float_t phiph[maxnph];

    // Set branch addresses
    tree->SetBranchAddress("nph", &nph);
    tree->SetBranchAddress("eph", eph);
    tree->SetBranchAddress("thetaph", thetaph);
    tree->SetBranchAddress("phiph", phiph);

    // Create histograms
    TH1F *h_invmass = new TH1F("h_invmass", "Invariant Mass of #pi^{0} Candidates; M_{#gamma#gamma} [GeV/c^{2}]; Number of Candidates", 100, 0, 0.3);
    TH1F *h_angle = new TH1F("h_angle", "Angle Between Photon Pairs; Angle [rad]; Number of Events", 100, 0, TMath::Pi());

    // Variable to count the total number of candidates
    Long64_t total_candidates = 0;

    // Loop over events
    Long64_t nentries = tree->GetEntries();
    std::cout << "Total number of events: " << nentries << std::endl;

    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        // Loop over all photon pairs
        for (Int_t j = 0; j < nph - 1; j++) {
            for (Int_t k = j + 1; k < nph; k++) {
                // Photon energies
                Float_t E1 = eph[j];
                Float_t E2 = eph[k];

                // Photon angles
                Float_t theta1 = thetaph[j];
                Float_t phi1 = phiph[j];

                Float_t theta2 = thetaph[k];
                Float_t phi2 = phiph[k];

                // Calculate momentum components of photons
                TVector3 p1;
                p1.SetMagThetaPhi(E1, theta1, phi1);
                TVector3 p2;
                p2.SetMagThetaPhi(E2, theta2, phi2);

                // Create four-momentum vectors of photons
                TLorentzVector photon1(p1, E1);
                TLorentzVector photon2(p2, E2);

                // Calculate invariant mass of the photon pair
                TLorentzVector pi0_candidate = photon1 + photon2;
                Float_t inv_mass = pi0_candidate.M();

                // Calculate angle between photons and fill the histogram
                Float_t angle = photon1.Angle(photon2.Vect());
                h_angle->Fill(angle);

                // Check if invariant mass is within the desired range
                if (inv_mass >= 0.1 && inv_mass <= 0.2) {
                    // Fill the invariant mass histogram
                    h_invmass->Fill(inv_mass);
                    // Increment the candidate counter
                    total_candidates++;
                }
            }
        }
    }

    // Output the total number of candidates
    std::cout << "Total number of pi0 candidates: " << total_candidates << std::endl;

    // Save histograms to a file
    TFile *outfile = new TFile("results.root", "RECREATE");
    h_invmass->Write();
    h_angle->Write();
    outfile->Close();

    std::cout << "Analysis complete. Results saved in file results.root" << std::endl;

    // Display histograms
    TCanvas *c1 = new TCanvas("c1", "Invariant Mass", 800, 600);
    h_invmass->Draw();

    TCanvas *c2 = new TCanvas("c2", "Angle Between Photons", 800, 600);
    h_angle->Draw();
}
