#include <iostream>
#include <cmath>
#include "TRandom3.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"

void task5()
{
    // Setting the number of sampling points
    const int N = 10000; // You can increase this number for higher accuracy

    // Initialize random number generator
    TRandom3 rand(0);

    // ------------------------------
    // 1. Inverse Transform Method for the distribution 2xe^{-x^2}
    // ------------------------------

    // Check and remove existing histogram to avoid warnings
    if (gDirectory->FindObject("hist_inverse"))
        gDirectory->Remove(gDirectory->FindObject("hist_inverse"));

    // Create histogram for the simulated distribution
    TH1F *hist_inverse = new TH1F("hist_inverse", "Distribution 2xe^{-x^{2}} (Inverse Transform Method)", 100, 0, 3);

    for (int i = 0; i < N; i++)
    {
        // Generate a uniformly distributed random number u in [0,1)
        double u = rand.Rndm();

        // Use the inverse function: x = sqrt(-ln(u))
        double x = sqrt(-log(u));

        hist_inverse->Fill(x);
    }

    // Check and remove existing Canvas to avoid warnings
    if (gROOT->FindObject("c1"))
        delete gROOT->FindObject("c1");

    // Draw the histogram
    TCanvas *c1 = new TCanvas("c1", "Inverse Transform Method", 800, 600);
    hist_inverse->GetXaxis()->SetTitle("x");
    hist_inverse->GetYaxis()->SetTitle("Frequency");
    hist_inverse->Draw();

    // ------------------------------
    // 2. Computing the integral of e^{-x} x^{3} from 0 to 1 using various methods
    // ------------------------------

    // Define the function to be integrated
    TF1 *func_exact = new TF1("func_exact", "exp(-x)*pow(x,3)", 0, 1);

    // Compute the exact value of the integral numerically
    double exact_value = func_exact->Integral(0, 1, 1e-10); // High-precision numerical integration

    std::cout << "Exact value of the integral (numerically): " << exact_value << std::endl;

    // Arrays to store results
    const int num_methods = 4;
    double results[num_methods];
    double errors[num_methods];

    // Integration limits
    double a = 0.0;
    double b = 1.0;

    // ------------------------------
    // Method 1: Monte Carlo Integration
    // ------------------------------

    double sum_mc = 0.0;
    for (int i = 0; i < N; i++)
    {
        double x = a + (b - a) * rand.Rndm();
        sum_mc += exp(-x) * pow(x, 3);
    }
    double integral_mc = (b - a) * sum_mc / N;
    results[0] = integral_mc;
    errors[0] = fabs(integral_mc - exact_value);

    // ------------------------------
    // Method 2: Neumann's Method (Rejection Sampling)
    // ------------------------------

    // Determine the maximum value of the function on the interval [0,1]
    double fmax = func_exact->GetMaximum();

    int count_under_curve = 0;
    int total_points = N;

    for (int i = 0; i < total_points; i++)
    {
        double x = a + (b - a) * rand.Rndm();
        double y = fmax * rand.Rndm();

        if (y <= exp(-x) * pow(x, 3))
        {
            count_under_curve++;
        }
    }

    double integral_neyman = fmax * (b - a) * count_under_curve / total_points;
    results[1] = integral_neyman;
    errors[1] = fabs(integral_neyman - exact_value);

    // ------------------------------
    // Method 3: Main Part Integration
    // ------------------------------

    // Split the integral into two sections
    double c = 0.5;
    double integral_main = 0.0;
    double sum_main = 0.0;

    for (int i = 0; i < N / 2; i++)
    {
        double x = a + (c - a) * rand.Rndm();
        sum_main += exp(-x) * pow(x, 3);
    }
    integral_main += (c - a) * sum_main / (N / 2);

    sum_main = 0.0;
    for (int i = 0; i < N / 2; i++)
    {
        double x = c + (b - c) * rand.Rndm();
        sum_main += exp(-x) * pow(x, 3);
    }
    integral_main += (b - c) * sum_main / (N / 2);

    results[2] = integral_main;
    errors[2] = fabs(integral_main - exact_value);

    // ------------------------------
    // Method 4: Rectangular Integration
    // ------------------------------

    double sum_rect = 0.0;
    double dx = (b - a) / N;

    for (int i = 0; i < N; i++)
    {
        double x = a + i * dx + dx / 2.0; // Center of the rectangle
        sum_rect += exp(-x) * pow(x, 3) * dx;
    }

    double integral_rect = sum_rect;
    results[3] = integral_rect;
    errors[3] = fabs(integral_rect - exact_value);

    // ------------------------------
    // 3. Comparing the results with the exact value
    // ------------------------------

    std::cout << "Method\t\tResult\t\t\tError" << std::endl;
    std::cout << "Monte Carlo\t" << results[0] << "\t" << errors[0] << std::endl;
    std::cout << "Neumann\t\t" << results[1] << "\t" << errors[1] << std::endl;
    std::cout << "Main Part\t" << results[2] << "\t" << errors[2] << std::endl;
    std::cout << "Rectangles\t" << results[3] << "\t" << errors[3] << std::endl;

    // ------------------------------
    // 4. Determining the error of each method and plotting histograms
    // ------------------------------

    const int num_experiments = 1000;
    TH1F *hist_methods[num_methods];
    const char *method_names[num_methods] = {"Monte Carlo", "Neumann", "Main Part", "Rectangles"};

    // Remove existing histograms to avoid warnings
    for (int m = 0; m < num_methods; m++)
    {
        if (gDirectory->FindObject(Form("hist_%d", m)))
            gDirectory->Remove(gDirectory->FindObject(Form("hist_%d", m)));
    }

    for (int m = 0; m < num_methods; m++)
    {
        hist_methods[m] = new TH1F(Form("hist_%d", m), method_names[m], 100, exact_value - 0.005, exact_value + 0.005);
    }

    for (int n = 0; n < num_experiments; n++)
    {
        // Monte Carlo
        sum_mc = 0.0;
        for (int i = 0; i < N; i++)
        {
            double x = a + (b - a) * rand.Rndm();
            sum_mc += exp(-x) * pow(x, 3);
        }
        integral_mc = (b - a) * sum_mc / N;
        hist_methods[0]->Fill(integral_mc);

        // Neumann
        count_under_curve = 0;
        for (int i = 0; i < total_points; i++)
        {
            double x = a + (b - a) * rand.Rndm();
            double y = fmax * rand.Rndm();

            if (y <= exp(-x) * pow(x, 3))
            {
                count_under_curve++;
            }
        }
        integral_neyman = fmax * (b - a) * count_under_curve / total_points;
        hist_methods[1]->Fill(integral_neyman);

        // Main Part
        integral_main = 0.0;
        sum_main = 0.0;
        for (int i = 0; i < N / 2; i++)
        {
            double x = a + (c - a) * rand.Rndm();
            sum_main += exp(-x) * pow(x, 3);
        }
        integral_main += (c - a) * sum_main / (N / 2);

        sum_main = 0.0;
        for (int i = 0; i < N / 2; i++)
        {
            double x = c + (b - c) * rand.Rndm();
            sum_main += exp(-x) * pow(x, 3);
        }
        integral_main += (b - c) * sum_main / (N / 2);

        hist_methods[2]->Fill(integral_main);

        // Rectangles
        sum_rect = 0.0;
        for (int i = 0; i < N; i++)
        {
            double x = a + i * dx + dx / 2.0;
            sum_rect += exp(-x) * pow(x, 3) * dx;
        }
        integral_rect = sum_rect;
        hist_methods[3]->Fill(integral_rect);
    }

    // Check and remove existing Canvas to avoid warnings
    if (gROOT->FindObject("c2"))
        delete gROOT->FindObject("c2");

    // Draw the histograms
    TCanvas *c2 = new TCanvas("c2", "Method Comparison", 1200, 800);
    c2->Divide(2, 2);

    for (int m = 0; m < num_methods; m++)
    {
        c2->cd(m + 1);
        hist_methods[m]->Draw();
        hist_methods[m]->Fit("gaus");

        hist_methods[m]->SetLineColor(kBlue);
        hist_methods[m]->SetFillColor(kCyan);
        hist_methods[m]->SetFillStyle(3001);
        hist_methods[m]->GetXaxis()->SetTitle("Integral Value");
        hist_methods[m]->GetYaxis()->SetTitle("Frequency");
    }

    // ------------------------------
    // 5. Determining the best method
    // ------------------------------

    // Output statistics for each method
    std::cout << std::endl << "Statistics for methods:" << std::endl;
    for (int m = 0; m < num_methods; m++)
    {
        double mean = hist_methods[m]->GetMean();
        double stddev = hist_methods[m]->GetRMS();
        std::cout << method_names[m] << ":\tMean = " << mean << "\tStd Dev = " << stddev << std::endl;
    }

    // Note: We are not manually deleting canvases and histograms to avoid segmentation faults.
    // ROOT will handle the cleanup when the session ends.
}
