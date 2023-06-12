//This code is going to plot and then extrapolate the mPMT behaviour at different R and rayff
#include "../include/truth_alpha.hh"
#include "../include/findBin.hh"
#include "../include/readTxtFile.hh"

#include <iostream>
#include "../chisq/chisq.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TFile.h"
// #include "TRandom3.h"
#include "TSpline.h"
#include "TLine.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"



void plot_Scattering_Absoption_ExtrapolateR(){
    float R = 120.;

    std::unique_ptr<TFile> myFile( TFile::Open("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/reference_root/reference_bin24/SimpleSpline/results_SimpleSpline_Abs_Scat_bin24_theta0.31_phi1.02_R120.00.root") );
    std::unique_ptr<TF1> best_fit_scat = std::unique_ptr<TF1>(myFile->Get("best_fit_scat"));
    std::unique_ptr<TF1> best_fit_abs(myFile->Get<TF1>("best_fit_abs"));
    std::unique_ptr<TGraphErrors> data_distribution(myFile->Get<TGraphErrors>("data_distribution"));
    std::unique_ptr<TGraphErrors> data_scat_distribution(myFile->Get<TGraphErrors>("data_scat_distribution"));

// // //     TFile f("","recreate");
//     TGraphErrors *data_distribution = (TGraphErrors*)f.Get("data_distribution");
//     TGraphErrors *data_scat_distribution = (TGraphErrors*)f.Get("data_scat_distribution");
//     TF1 *best_fit_scat = (TF1*)f.Get("best_fit_scat");
//     TF1 *best_fit_abs = (TF1*)f.Get("best_fit_abs");
    std::cout << "oh oui" << std::endl;

    TCanvas *c = new TCanvas("c", "c", 800,800);
    c->Draw();
    TPad *p1 = new TPad("p1","p1",0.,0.35,1.0,1.);
    p1->SetBorderSize(2);

    c->cd(0);
    TPad *p2 = new TPad("p2","p2",0.,0.01,1.0,0.377);
    p2->SetBorderSize(0);
    p2->Draw();
    p1->Draw();

    std::cout << "oh oui 2" << std::endl;



//     data_distribution->SetMarkerColor(kBlue);
//     data_distribution->SetLineColor(kBlue);
//     data_distribution->SetFillColor(kBlue);
//     data_scat_distribution->SetMarkerColor(kRed);
//     data_scat_distribution->SetLineColor(kRed);
//     data_scat_distribution->SetFillColor(kRed);
//     best_fit_abs->SetLineColor(kBlue);

    std::cout << "oh oui 3 " << std::endl;
    //The charge is in bin 24 and collected with test data
    std::vector<double> test_Q = {0,0,0,0,0};   //empty array for initialisation
    if (R == 20)  test_Q = {321.48, 325.56, 324.30,325.31,327.45 };    //(1)
    if (R == 10)  test_Q = {304.91, 314.10, 305.54, 325.36, 305.67};   //(0)
    if (R == 20)  test_Q = {321.48, 325.56, 324.30,325.31,327.45 };    //(1)
    if (R == 40)  test_Q = {314.23, 322.06, 325.50, 325.24, 325.16};   //(2)
    if (R == 80)  test_Q = {292.42, 298.51,297.05, 306.54, 299.17};    //(3)
    if (R == 120)  test_Q = {297.54, 293.80, 305.80, 295.85, 300.33};  //(4)
    if (R == 140)  test_Q = {285.01, 303.03, 298.70, 302.91, 301.76 }; //(5)
    if (R == 160)  test_Q = {282.40, 294.01, 298.86, 298.14, 303.84};  //(6)
    if (R == 180)  test_Q = {279.21, 284.30, 289.52, 298.80, 305.94};  //(7)
    if (R == 210)  test_Q = {265.02, 291.86, 290.76,299.71,290.61 };   //(8)
    if (R == 250)  test_Q = {260.07, 281.21, 282.66, 290.70, 291.62 }; //(9)

    //at these scatteing attenuation lengths
    std::vector<double> test_rayff = {0.08322,0.13888, 0.19444, 0.24999, 0.2222};
    std::vector<double> test_scattering_length;
    for (double r : test_rayff) test_scattering_length.push_back(truth_alpha(401.9, 10e10, r));
    // 	std::cout << test_scattering_length[1] << std::endl;
    std::cout << "oh oui 4 " << std::endl;
    double x, y, y_err;
    double data, pred, compa;
    double total_y_test = 0;
    TGraphErrors * test = new TGraphErrors;
    TGraphErrors * compa_test = new TGraphErrors;

    std::cout << "oh oui 5" << std::endl;
    for (int l=0; l < test_scattering_length.size(); l++)
    {
        x = (double) test_scattering_length[l];
        data = (double) test_Q[l];
        y_err = TMath::Sqrt(data * 70 * (1 - data/1000)) / (data * 70) * data;
        std::cout << "oh oui 6" << std::endl;

        pred = best_fit_scat->Eval(x);
        std::cout << "oh oui 8" << std::endl;
        compa = (data - pred)/data;
        std::cout << "oh oui 7" << std::endl;
        total_y_test += ((data- pred) * (data - pred))/pred;
        std::cout << x << " " << y << std::endl;
        test->SetPoint(l,x,data);
        test->SetPointError(l, 0, y_err);
        compa_test->SetPoint(l, x, compa);
        compa_test->SetPointError(l, 0, y_err/pred);
    }

    TGraphErrors * compa_graph = new TGraphErrors;
    std::cout << "oh oui" << std::endl;
    double total_y = 0;
    for (int i=0; i<data_scat_distribution->GetN(); i++){
        std::cout << i << std::endl;
        x =(double) data_scat_distribution->GetX()[i];
        y = data_scat_distribution->GetY()[i];
        double yerr = data_scat_distribution->GetErrorY(i);
        double y_pred = best_fit_scat->Eval(x);
        std::cout << y-y_pred << " x " << x << " "<< y_pred<< std::endl;
        if (y_pred == 0) y_pred = 1.;
        yerr = yerr/y;
        y = (y - y_pred)/y_pred;
        total_y += y * y * y_pred; //this is the chi2
        compa_graph->SetPoint(i,x,y);
        compa_graph->SetPointError(i,0,yerr);
    }
    TLine * zeros = new TLine(0.0, 0.0, x + 100, 0.0 );
    p2->cd(0);
    compa_graph->Draw("AE*");
    zeros->Draw("same");
    std::cout << "oh oui" << std::endl;
    zeros->SetLineColor(kRed);
    std::cout << "oh oui" << std::endl;
    TLegend *legend2 = new TLegend(0.35,0.65,0.86,0.82);
    legend2->AddEntry(compa_graph,Form("Reference data\n(y-y_pred)^2/y_pred per point: %.2e", total_y/data_scat_distribution->GetN()),"f");
    legend2->AddEntry(compa_test,Form("Test data\n(y-y_pred)^2/y_pred per point: %.2e", total_y_test/test_rayff.size()),"f");
    compa_graph->GetXaxis()->SetTitle("Attenuation length (cm)");
    compa_graph->GetXaxis()->SetLabelSize(0.055);
    compa_graph->GetXaxis()->SetTitleSize(0.05);
    compa_graph->GetXaxis()->SetTitleOffset(0.95);
    compa_graph->GetYaxis()->SetTitle("(data - fit)/fit");
    compa_graph->GetYaxis()->SetTitleOffset(0.7);
    compa_graph->GetYaxis()->SetLabelSize(0.055);
    compa_graph->GetYaxis()->SetTitleSize(0.06);
    std::cout << "oh oui" << std::endl;
    compa_graph->SetFillColor(kRed);
    compa_test->SetFillColor(kBlack);
    compa_test->SetMarkerColor(kBlack);
    compa_test->SetLineColor(kBlack);
    compa_test->Draw("sameE*");
    compa_graph->SetMarkerColor(kRed);
    compa_graph->SetLineColor(kRed);
    legend2->Draw("same");
    std::cout << "oh oui" << std::endl;

    p1->cd(0);
    std::cout << "oh oui" << std::endl;
    TLegend *legend = new TLegend(0.5,0.55,0.66,0.72);
//     legend->AddEntry(data_scat_distribution,"Scattering","f");
//     legend->AddEntry(data_distribution,"Absorption","f");
//     // 	legend->AddEntry(test,"Test data points","f");
    data_scat_distribution->Draw("AE*");
    test->SetFillColor(kBlack);
    test->SetMarkerColor(kBlack);
    test->SetLineColor(kBlack);
    std::cout << "oh oui" << std::endl;
    test->Draw("sameE*");
    best_fit_scat->Draw("same");
    data_distribution->Draw("sameE*");
    best_fit_abs->Draw("same");
    legend->Draw("same");
    c->SaveAs(Form("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/Pictures/Bin24_Scattering_FitQuality_15nodes/Bin24_15nodes_R%.2f.pdf", R));




}

int main() {
    plot_Scattering_Absoption_ExtrapolateR();
    return 0;
}

