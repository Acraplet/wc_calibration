//After having imported a given reference root file we can run this macro to plot the behaviour of the charge as 
//the attenuation or scattering length is increased

void plot_Scattering_Absoption(){
	data_distribution->SetMarkerColor(kBlue);
	data_distribution->SetLineColor(kBlue);
	data_distribution->SetFillColor(kBlue);
	data_scat_distribution->SetMarkerColor(kRed);
	data_scat_distribution->SetLineColor(kRed);
	data_scat_distribution->SetFillColor(kRed);
	best_fit_abs->SetLineColor(kBlue);
	TLegend *legend = new TLegend(0.55,0.65,0.76,0.82);
	legend->AddEntry(data_scat_distribution,"Scattering","f");
	legend->AddEntry(data_distribution,"Absorption","f");
	data_scat_distribution->Draw("AE*");
	best_fit_scat->Draw("same");
	data_distribution->Draw("sameE*");
	best_fit_abs->Draw("same");
	legend->Draw("same");
}
