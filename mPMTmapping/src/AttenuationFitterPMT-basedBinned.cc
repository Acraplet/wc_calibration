//This code is extracting the absorption length corresponding to one or more maps that have the same abwff
//and ideally rayff = inf
//For this is is looking up in Maps the reference number of pe collected per photon for this given bin
//
//Then the chi is minimised between the observed data and y_pred = Fake_Data_Spline(scattering length pred)
//expected scattering length. So far this is position dependant but it will have to be upgraded to a binned approach
#include "../include/findBin.hh"
#include "../include/readReferenceFiles.hh"
#include "../include/readTxtFileWithTime.hh"
#include "../include/truth_alpha.hh"
#include <iostream>
#include "../chisq/chisq.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1.h"
#include "TGraphErrors.h"
//#include "../include/toml/toml_helper.h"

void HelpMessage()
{
    std::cout   << "USAGE: "
                << "AttenuationFitterBinned" << "\nOPTIONS:\n"
                << "-b : base of the set of files you want to use \n"
                << "-c : specific files in that set that you want to read together (added to the base)\n"
                << "-n : Number of bins max \n"
                << "-o : output file (one _withText and one _withoutText) to which we will append the information and outcome of this specific fit - useful for scanning \n"
                << "-g : initial guess for the scattering length \n"
                << "-Q : charge threshold for a bin to be included in the fit \n";
}



int main(int argc, char **argv){

    std::vector<double> list_Q, list_i,list_A, list_R;
    //std::vector<std::vector<double>> list_ynodes, list_xnodes;
    //char* R_test;
    //char* Q_test;
    //char* theta_test;
    //double  theta_test_num, Q_test_num;
    //char* phi_test;
    double noAttenuation_pred; //not actually useful - just ot keep track and compare the
    //pred vs truth for absorption
    //The config number is a rough way to output the distances that we have looked out 
    //useful for run comparisions
    long int config_number = 0.0;
    double w = 0.0;
    
    //TODO: change for a config file where the max number of bins is an entry
    int Q_thresh = 0; //minimum number of photons in a given test bin to be included in the fit
    int base = 740; //this is the base of the file we have

    std::string configuration = "0123"; //these are the R configurations we are looking at 
    
    double initGuessAbs = 1500.;
    double trueScat;
    double R;

    int nBins = 20; //in the new scheme we have 20 bins and for each bins we have 1 ref value per PMT
    int nPMTs = 19; //in the new scheme we have 20 bins and for each bins we have 1 ref value per PMT
    std::string output_file = "out";
    std::string config_file;
    
    char option;
    while((option = getopt(argc, argv, "x:m:c:b:Q:i:hn:g:o:")) != -1)
    {
	switch(option)
        {
		case 'c':
			configuration = optarg;
			break;

     	        case 'n':
			nBins = std::stoi(optarg);
			break;

		case 'b':
			base = std::stoi(optarg);
			break;

		case 'Q':
			Q_thresh = std::stoi(optarg);
			break;

		case 'g':
			initGuessAbs = std::stof(optarg);
			break;

		case 'o':
			output_file = optarg;
			break;

		case 'h':
			HelpMessage();
			return 0;
			break;


		default:
			return 0;

		}
    }
    //Setup for a .toml config - not working yet
   // auto const &card_toml = toml_h::parse_card(config_file);
   // int const &nBins = toml_h::find(card_toml, "nBins");
   // int const &bin_min = toml_h::find(card_toml, "min_bin");
   // int const &Q_thresh = toml_h::find(card_toml, "Q_tresh");
   // double const &spline_min = toml_h::find(card_toml, "spline_min");
   // double const &spline_max = toml_h::find(card_toml, "spline_max");
   // double const &spline_increment = toml_h::find(card_toml, "spline_increment");
   // double const &initGuessAbs = toml_h::find(card_toml, "initGuessAbs");

    std::vector<int> list_files;
    for (int c=0; c<=configuration.size()-1; c++){
	    list_files.push_back(base + int(configuration[c]) - 48); //char to int conversion needs -48
    }

    std::cout << "The total number of files to read together is: " << configuration.size() << std::endl;
    //for (int f=0; f<= list_files.size()-1; f++){ //argc-1
	//This is the test file we are going to extract the scattering length out of
	char * filename = "test_hadded/ReferenceAttenuation_ALLhits_PMT-basedBin_R80.00.txt"; //Form("./Maps/maps_txtFiles/mPMT_map_ID%d.txt") //,list_files[f]);
	//in these double we are storing the total charge, one entry per bin
	std::vector<std::vector<double>> total_bin_charge;
	
	std::cout << "Hello" << std::endl;
	double total_bin_photons[nBins];
	//initialisation of the bin-dependant counters
	for (int k=0; k<=nBins; k++){
		//std::cout << k << std::endl;
		total_bin_photons[k] = 0.;
		std::vector<double> a;
		for (int p=0; p<=nBins-1; p++){
			//std::cout << p << std::endl;
			a.push_back(0.);
		        //total_bin_charge[k][p] = 0.;
		}	
		total_bin_charge.push_back(a);
	}
	std::cout << "Hello" << std::endl;
	//This is reading in the test file 
	std::vector<DataWithTime> test_positions = readTxtFileWithTime(filename);
	std::cout << test_positions.size() << std::endl;
	
	for (int i = 0; i < test_positions.size(); i++){
		//read each source position in the file one by one
		DataWithTime pos =  test_positions[i];
		//Then we check which bin it belongs to - now we already have it 
		//Bin closestBin = findPMTBin(pos.theta, pos.phi, pos.time);
		//std::cout << "The closest bin is " << pos.bin << " and the hit PMT " << pos.mPMT_pmt<< std::endl;
		//add to the charge the fractionnal charge collected at this position
		//we do not use direct charge because the comparision has to be made
		//with respect to 1000 photons for now
		//std::cout << "Q " << pos.Q << " existing " << total_bin_charge[int(pos.bin)][int(pos.mPMT_pmt)] << std::endl;
		total_bin_charge[int(pos.bin)][int(pos.mPMT_pmt)] = total_bin_charge[int(pos.bin)][int(pos.mPMT_pmt)] + float(pos.Q);
		total_bin_photons[int(pos.bin)] += 1;
		R = pos.R;
		std::cout << "Bin " << int(pos.bin) << " PMT " << int(pos.mPMT_pmt) << " charge " << pos.Q <<std::endl;
	}

	for (int binTarget=0; binTarget < nBins; binTarget++){
		for (int PMTTarget=0; PMTTarget < nPMTs; PMTTarget++){
		//I am not sure we need any of the default, just normalise properly 
	 	//total_bin_charge[binTarget][PMTTarget] = total_bin_charge[binTarget][PMTTarget] / total_bin_photons[binTarget] * 1000;
		if (total_bin_charge[binTarget][PMTTarget]>= Q_thresh and total_bin_photons[binTarget]!=0 ){
			std::cout << "Total charge " << total_bin_charge[binTarget][PMTTarget] <<  "Bin " << binTarget << " pmt " << PMTTarget << std::endl;
//			int n = 0;
//			double x, y;
			double ref_info; //empty variable to read the file with
			int count = 0;
			//Now fitting multiple bins together

			//need to get the reference amplitude
			//const char* fimpName = Form("Absorption_PMT-basedBin%i_R%.2f.txt", binTarget, test_positions[0].R);
			//std::ifstream in(fimpName);
			//The reference datasets to use for comparision
			//std::vector<double> x_vector, y_vector, z_vector;
			
			//readAbsorption returns the fractionnal number of photons in that PMT for the given bin
			double ref = readAbsorptionRef(binTarget, PMTTarget, R);
			list_A.push_back(ref * total_bin_photons[binTarget]);
			list_R.push_back(R);

//Old way of doing things			
/*			
			while ((in >> ref_info)) {
				if (count - PMTTarget == 0) {
					//std::cout << ref_info << " " << binTarget << " " << test_positions[0].R << std::endl;
					noAttenuation_pred = ref_info *  total_bin_photons[binTarget];
					list_A.push_back(noAttenuation_pred); //need to weight by the
					//number of photons we have simulated in that given bin! - for now all the files have
					//the same number of
                    // events and distance but eventually we'll have to chage that !!
					list_R.push_back(test_positions[0].R);
					break; //the max amplitude at the bin of interest
				}
				count += 1;
			}
*/

			//w is just for bookkepping of the configuration (i.e. which files we fit together)
			list_i.push_back(w);
			double Q = total_bin_charge[binTarget][PMTTarget];
			list_Q.push_back(Q);
 			
			std::cout << " bin " << binTarget << " Hit PMT " << PMTTarget << " excat attenuation length length =" << truth_alpha(401.9, test_positions[0].abwff,test_positions[0].rayff);
			std::cout << " abwff " << test_positions[0].abwff << " rayff " << test_positions[0].rayff << " R = " <<   test_positions[0].R << " charge collected: " ;
			std::cout << total_bin_charge[binTarget][PMTTarget] << " total number of photons sent in this bin " << total_bin_photons[binTarget] << " prediction without attenuation " << ref << std::endl;
		}//if we have more than Q_thresh hits in a given bin
	     }//run through all the PMTs
	}//run through all the bins
	//The true scattering length that we know from the file
	trueScat = truth_alpha(401.9, test_positions[0].abwff,test_positions[0].rayff);
	//Keeping track of the files we fitted together
	config_number = std::stoi(configuration);
	w+=1;
        //file_ref->Close();
    //} //finished reading all of the test maps
    
    //Here the fitting begins
    std::cout << "Final config number = " << config_number << std::endl;
    const int nPars = 1; //the only parameter we fit is scattering length
    
// FOR NOW NOT FITTING BECAUSE CHISQ IS NOT YET WORKING    
    Chisq *chi = new Chisq(nPars);
    //In this case list_i is a index storing, later we will have as many entries as our number of bins
    chi->setData(list_i, list_Q);
    chi->setRef(list_A, list_R);

    ROOT::Math::Functor functor(chi, &Chisq::fcn_abwff, nPars);
    ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    min->SetStrategy(3);
    min->SetFunction(functor);
    min->SetMaxFunctionCalls(10000);
    min->SetVariable(0, "attenuation_length", initGuessAbs, 0.01);

    min->Minimize();
    min->PrintResults();
    const double * res_scat = min->X();
    const double * err_scat = min->Errors();
    //to evaluate accuracy of the minimizer
    double pres_scat = min->MinValue();

    //Here output some numbers for easier analysis
    std::ofstream outfile;
    outfile.open(Form("%s_withText.txt", output_file.c_str()), std::ofstream::app);
    outfile << "True: " << trueScat << " config: " << config_number << " reco: " << res_scat[0] << " +/- " << err_scat[0];
    outfile << " Q_thresh " << Q_thresh << " initGuessAbs " << initGuessAbs << " FVAL " << pres_scat<< std::endl;
   
    //and without text for analysis and doing scans
    std::ofstream outfile_noText;
    outfile_noText.open(Form("%s_withoutText.txt", output_file.c_str()), std::ofstream::app);
    outfile_noText << " " << trueScat << " " << config_number << " " << res_scat[0] << " " << err_scat[0];
    outfile_noText << " " << Q_thresh << " " << initGuessAbs << " " << pres_scat << std::endl;

}


