//This code is going to plot and then extrapolate the mPMT behaviour at different R and rayff

#include "./include/truth_alpha.hh"




//This is the fitter part - copied and adapted from the cern home page: https://root.cern/doc/master/classTGraph2D.html
void graph2dfit(TGraph2D *dt)
{
	const char* fimpName = "./Maps/test_maps_oneBin/all_test_files_bin356.txt";
	std::ifstream in(fimpName);
	double temp;
	int count = 0;
//The reference datasets to use for comparision
	std::vector<double> x_vector, y_vector, z_vector;
	while ((in >> temp)) {
		if (count %3 == 0) {
// 			double r = std::strtod(in, NULL);

			x_vector.push_back(truth_alpha(401.9,10e10,temp));
			std::cout << truth_alpha(401.9,10e10,temp) << " ";
// 			std::cout << truth_alpha(401.9,10e10,temp) << " " << temp << std::endl;
		}//scat len
		if (count %3 == 1) {
			y_vector.push_back(temp); //R
		}
		if (count %3 == 2){
			z_vector.push_back(temp); //Q per 1000 photons
		}
		count+=1;
// 		first rayff
// 		then R
// 		then charge per 1000
// 		hist->Fill(x);
	}


	//dt is the thing we are going to fit for
	gStyle->SetOptStat(0);
	gStyle->SetOptFit();

	auto cq = new TCanvas("c","Graph2D example",0,0,600,800);
	cq->Divide(2,2);

	double rnd, x, y, z;
	double e = 100;
	int nd = 400;
	int np = 10000;

	TRandom r;
	double fl = 400;
	//1000*(([0]*sin(x)/x)*([1]*sin(y)/y))+200
	auto f2 = new TF2("f2","[0] + (1/(-x+[2])) * [1]",
					  0,6000, 0, 250);
	f2->SetParameters(1,1,1,1, 1);
// 	auto dt = new TGraph2D()

	double hr = 0.5;
	auto h1 = new TH1D("h1",
					   "Difference between Original splitline{function and Function}{with noise}",
					100, -hr, hr);
	auto h2 = new TH1D("h2",
					   "Test - Delaunay interpolation",
					50, -hr/5, hr/5);
	auto h3 = new TH1D("h3",
					   "Test - Minuit fit : p0 + (1/(-x + p1)) * p2",
					50, -hr, hr);

// 	f2->SetParameters(0.5,1.5);
	//use this function to fit
	dt->Fit(f2);
	auto fit2 = (TF2*)dt->FindObject("f2");

	f2->SetParameters(300,4000, -129, 45, 67);
	auto f3 = new TGraph2D();

	//here we need to compare with our true data points!
	for (int N=0; N<x_vector.size(); N++) {
		// 		f2->GetRandom2(x,y);
		double x = x_vector[N];
		double y = y_vector[N];
		double z_true = z_vector[N];
		f3->SetPoint(N, x, y, z_true);
		// Generate a random number in [-e,e]
		//rnd = 2*r.Rndm()*e-e;
		//this is the function (without plus noise)
		z = f2->Eval(x,y); //*(1+rnd);

		h1->Fill((z_true-z)/z_true);
		//this is a simple interpolation
		z = dt->Interpolate(x,y);
		h2->Fill((z_true-z)/z_true);

		//this is the function fited
		z = fit2->Eval(x,y);
		h3->Fill((z_true-z)/z_true);
	}




// 	c->cd(1);
// 	f2->SetTitle("Original function with Graph2D points on top");
// // 	f2->SetMaximum(zmax);
// 	gStyle->SetHistTopMargin(0);
// 	f2->Draw("surf1");
// 	f3->Draw("same p");
// 	dt->Draw("same p0");

	cq->cd(1);
	dt->SetMargin(0.1);
	dt->SetFillColor(36);
	dt->SetTitle("Delaunay interp.");
// 	dt->SetTitleSize(0.03);
	dt->Draw("surf1");
	f3->Draw("same p0");
// 	dt->Draw("same p0");

	cq->cd(3);
	fit2->SetTitle("Minuit fit result\n p0 + (1/(-x + p1)) * p2");

	fit2->Draw("surf1");
	f3->Draw("same p0");
// 	dt->Draw("same p0");

	h1->SetFillColor(47);
	h2->SetFillColor(38);
	h3->SetFillColor(29);

// 	c->cd(2); h1->Fit("gaus","Q") ; h1->Draw();
// 	h2->SetTitleSize(1);
	h2->GetXaxis()->SetTitle("True-Pred/True charge");
	h3->GetXaxis()->SetTitle("True-Pred/True charge");
	cq->cd(2); h2->Fit("gaus","Q") ; h2->Draw();
	cq->cd(4); h3->Fit("gaus","Q") ; h3->Draw();
	cq->cd();
}


void plot_Scattering_Absoption_ExtrapolateR(){
	std::vector<float> R_list = {10., 20., 40., 80., 120., 140., 160., 180., 210., 250.};

	TCanvas *c = new TCanvas("c", "c",0,0,600,400);
	c->cd(0);
	auto f2 = new TGraph2D();
	gStyle->SetPalette(1);
	int n = 0;
	double x, y;
	std::vector<double> test_Q;
	std::vector<double> x_true, y_true, z_true;
	for (int ri=0; ri < R_list.size(); ri++)
	{
		TFile f(Form("./reference_root/reference_bin_all/5nodesRef/results_5nodesRef_Abs_Scat_bin356_theta0.90_phi3.04_R%.2f.root", R_list[ri]));

		TGraphErrors *data_distribution = (TGraphErrors*)f.Get("data_distribution");
		TGraphErrors *data_scat_distribution = (TGraphErrors*)f.Get("data_scat_distribution");
		TF1 * best_fit_scat = (TF1*)f.Get("best_fit_scat");
		TF1 * best_fit_abs = (TF1*)f.Get("best_fit_abs");

		for (int i=0; i<data_scat_distribution->GetN(); i++){
			x =(double) data_scat_distribution->GetX()[i];
			y = (double) data_scat_distribution->GetY()[i];
			f2->SetPoint(n, x, R_list[ri], y);
// 			x_true.push_back(x);
// 			y_true.push_back(R_list[ri]);
// 			z_true.push_back(y);
			n +=1;
		}

		//comment all of the following to only look at "proper reference datasets" - not the test ones
		double R = R_list[ri];
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
		std::vector<double> test_rayff = {0.08322,0.13888, 0.19444, 0.24999, 0.2222};
		std::vector<double> test_scattering_length;
		for (double r : test_rayff) test_scattering_length.push_back(truth_alpha(401.9, 10e10, r));


		for (int i = 0; i < test_Q.size(); i++){
// 			std::cout << i+ri*5 << std::endl;
// 			f2->SetPoint(n, test_scattering_length[i], R,  test_Q[i] );
			x_true.push_back(test_scattering_length[i]);
			y_true.push_back(R);
			z_true.push_back(test_Q[i]);
// 			std::cout << i << "!" << test_Q[i]<< std::endl;
			n += 1;
		}
// 		if (ri) f2->Draw("same");
	}
	f2->GetXaxis()->SetTitle("Scattering length(cm)");
	f2->GetYaxis()->SetTitle("mPMT-source distance R (cm)");
	f2->GetZaxis()->SetTitle("Charge Collected per 1000 photons");
	f2->GetXaxis()->SetTitleOffset(3.6);
// 	f2->GetYaxis()->SetLabelSize(0.055);
	f2->GetXaxis()->SetTitleSize(0.03);
	f2->GetYaxis()->SetTitleOffset(3.6);
	f2->GetYaxis()->SetTitleSize(0.03);
	f2->GetZaxis()->SetTitleOffset(1.6);
	f2->GetZaxis()->SetTitleSize(0.03);
	f2->SetTitle("Bin 24");
	//remove "Surf1" to see the individual points
	f2->Draw("Surf1");
// 	f2->Close();
	graph2dfit(f2);

}


void plot_Scattering_Absoption_ExtrapolateR_safe(){




	float R = 140.;

	TFile f(Form("reference_root/reference_bin24/SimpleSpline/results_SimpleSpline_Abs_Scat_bin24_theta0.31_phi1.02_R%.2f.root", R));


	TGraphErrors *data_distribution = (TGraphErrors*)f.Get("data_distribution");
	TGraphErrors *data_scat_distribution = (TGraphErrors*)f.Get("data_scat_distribution");
	TF1 * best_fit_scat = (TF1*)f.Get("best_fit_scat");
	TF1 * best_fit_abs = (TF1*)f.Get("best_fit_abs");



/*
	TFile myFile = (TFile) TFile::Open("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/reference_root/reference_bin24/SimpleSpline/results_SimpleSpline_Abs_Scat_bin24_theta0.31_phi1.02_R120.00.root");
	TF1 best_fit_scat(myFile.Get<TF1>("best_fit_scat"));
	TF1 best_fit_abs(myFile.Get<TF1>("best_fit_abs"));
	TGraphErrors data_distribution(myFile.Get<TGraphErrors>("data_distribution"));
	TGraphErrors data_scat_distribution(myFile.Get<TGraphErrors>("data_scat_distribution"));*/


	TCanvas *c = new TCanvas("c", "c", 800,800);
	c->Draw();
	TPad *p1 = new TPad("p1","p1",0.,0.35,1.0,1.);
	p1->SetBorderSize(2);

	c->cd(0);
	TPad *p2 = new TPad("p2","p2",0.,0.01,1.0,0.377);
	p2->SetBorderSize(0);
	p2->Draw();
	p1->Draw();


	data_distribution->SetMarkerColor(kBlue);
	data_distribution->SetLineColor(kBlue);
	data_distribution->SetFillColor(kBlue);
	data_scat_distribution->SetMarkerColor(kRed);
	data_scat_distribution->SetLineColor(kRed);
	data_scat_distribution->SetFillColor(kRed);
	best_fit_abs->SetLineColor(kBlue);
// 	The charge is in bin 24 and collected with test data
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

	double x, y, y_err;
	double data, pred, compa;
	double total_y_test = 0;
	TGraphErrors * test = new TGraphErrors;
	TGraphErrors * compa_test = new TGraphErrors;


	for (int l=0; l < test_scattering_length.size(); l++)
	{
		x = (double) test_scattering_length[l];
		data = (double) test_Q[l];
		y_err = TMath::Sqrt(data * 70 * (1 - data/1000)) / (data * 70) * data;
		std::cout << x << " " << y_err << std::endl;
		pred = best_fit_scat->Eval(x);
		compa = (data - pred)/data;
		total_y_test += ((data- pred) * (data - pred))/pred;
		std::cout << x << " " << y << std::endl;
		test->SetPoint(l,x,data);
		test->SetPointError(l, 0, y_err);
		compa_test->SetPoint(l, x, compa);
		compa_test->SetPointError(l, 0, y_err/pred);
	}

	TGraphErrors * compa_graph = new TGraphErrors;

	double total_y = 0;
	for (int i=0; i<data_scat_distribution->GetN(); i++){
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
	zeros->SetLineColor(kRed);
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
	compa_graph->SetFillColor(kRed);
	compa_test->SetFillColor(kBlack);
	compa_test->SetMarkerColor(kBlack);
	compa_test->SetLineColor(kBlack);
	compa_test->Draw("sameE*");
	compa_graph->SetMarkerColor(kRed);
	compa_graph->SetLineColor(kRed);
	legend2->Draw("same");


	p1->cd(0);
	TLegend *legend = new TLegend(0.5,0.55,0.66,0.72);
	legend->AddEntry(data_scat_distribution,"Scattering","f");
// 	legend->AddEntry(data_distribution,"Absorption","f");
	legend->AddEntry(test,"Test data points","f");
	data_scat_distribution->Draw("AE*");
	test->SetFillColor(kBlack);
	test->SetMarkerColor(kBlack);
	test->SetLineColor(kBlack);
	test->Draw("sameE*");
	best_fit_scat->Draw("same");
	data_distribution->Draw("sameE*");
	best_fit_abs->Draw("same");
	legend->Draw("same");
	c->SaveAs(Form("./Pictures/Bin24_Scattering_FitQuality_15nodes/Bin24_15nodes_R%.2f.pdf", R));




}

