//This code is extracting the absorption length corresponding to one or more maps that have the same abwff
//and ideally rayff = inf
//For this is is looking up in Maps the reference number of pe collected per photon for this given bin
//multiply that by the expected number of photons reaching that specific PMT
//Then the chi is minimised between the observed data and y_pred = nPhotonsReaching * perpendicularQE * exp(-R/absff_length)
//to get the expected absorption length. So far this is position dependant but it will have to be upgraded to a binned approach
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
    double noAttenuation_pred; //not actually useful - just to keep track and compare the
    //pred vs truth for absorption
    //The config number is a rough way to output the distances that we have looked out 
    //useful for run comparisions
    long int config_number = 0.0;
    double w = 0.0;
    
    //TODO: change for a config file where the max number of bins is an entry
    int Q_thresh = 0; //minimum number of photons in a given test bin to be included in the fit
    int base = 740; //this is the base of the file we have
    double min_R_distance = 0.1; //cm the minium distance between the source and a PMT to make it in the fit  
    double direct_hit_time = 0.15; //ns is the TOF corrected time after which the hits aren't direct anymore. 
    std::string configuration = "0123"; //these are the R configurations we are looking at 
    
    double initGuessAbs = 1500.;
    double trueScat;
    double R;

    int nBins = 20; //in the new scheme we have 20 bins and for each bins we have 1 ref value per PMT
    int nPMTs = 20; //in the new scheme we have 20 bins and for each bins we have 1 ref value per PMT
    int n_mPMTs = 106;
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

    std::vector<int> list_files;
    for (int c=0; c<=configuration.size()-1; c++){
	    list_files.push_back(base + int(configuration[c]) - 48); //char to int conversion needs -48
    }

    std::cout << "The total number of files to read together is: " << configuration.size() << std::endl;
    for (int f=0; f<= list_files.size()-1; f++){ //argc-1
	
	
	//We need to read in the data: the txt map with the number of PE detected in each mPMT-PMT couple
	char * filename_data = Form("./Maps/maps_txtFiles/mPMT_map_ID%i.txt",list_files[f]);
	
	//Need to read in the predicted number of photon reaching each PMT, that is from the solid angle
	char * filename_solidAngle = Form("./Maps/Expected_Number_Photons_FileID%i.txt",list_files[f]);

	//and also will need to read in the perpendicular effiency maximal, 
	//that is already done with readAbsorptionRef?
	

	//in these double we are storing the total charge, one entry per mPMT-PMT pair
	std::vector<std::vector<double>> total_bin_charge;
	std::vector<std::vector<double>> total_bin_photons;
	std::vector<std::vector<double>> R_array;
	
	//initialisation of the bin-dependant counters
	for (int k=0; k<=n_mPMTs; k++){
		//total_bin_charge[k] = 0.;
		std::vector<double> a;
		for (int p=0; p<=nPMTs; p++){
			a.push_back(0.);
		        //total_bin_charge[k][p] = 0.;
		}	
		total_bin_charge.push_back(a);
		total_bin_photons.push_back(a);
		R_array.push_back(a);
	}
	//This is reading in the test file 
	std::vector<DataWithTime> test_positions = readTxtFileWithTime(filename_data);
	
	//this is reading in the expected number of photons (in Q)
	std::vector<DataWithTime> expected_photons = readLBpredNbPhotons(filename_solidAngle);
	
	for (int i = 0; i < test_positions.size(); i++){
		//read each source position in the file one by one
		DataWithTime pos =  test_positions[i];
		//Then we check which bin it belongs to - now we already have it 
		// we have one R per PMT, we need to take that into account
		if (pos.R>=min_R_distance and pos.time < direct_hit_time){
			//we cannot trust the reconstructed bin, we need to assume that the 
			total_bin_charge[int(pos.mPMT)][int(pos.mPMT_pmt)] += double(pos.Q);
			R_array[int(pos.mPMT)][int(pos.mPMT_pmt)] = pos.R;
		}
	}

	//save the expected number of photons from the solid angle calculation
	for (int i = 0; i < expected_photons.size(); i++){
		DataWithTime pos =  expected_photons[i];
		if (R_array[int(pos.mPMT)][int(pos.mPMT_pmt)]>=min_R_distance){
			total_bin_photons[int(pos.mPMT)][int(pos.mPMT_pmt)] = pos.Q;
		}
	}


	for (int mPMTTarget=0; mPMTTarget < n_mPMTs; mPMTTarget++){
	    for (int PMTTarget=0; PMTTarget < nPMTs; PMTTarget++){
		if (total_bin_charge[mPMTTarget][PMTTarget]>= Q_thresh and total_bin_photons[mPMTTarget][PMTTarget]>Q_thresh ){
			double ref_info; //empty variable to read the file with
			int count = 0;
			//Now fitting multiple bins together
			//photons in that PMT for the given bin
			double QE_perp = 1; //getPerpendicularQE(PMTTarget);
			//std::cout << "The QE to perpendicular light is " << QE_perp  << " for PMT " << PMTTarget << std::endl;
			//arbitrary scaling
			//if (PMTTarget < 13) QE_perp = 1/1.43;
			//else if (PMTTarget == 19) QE_perp = 1/0.77;
			//else QE_perp = 1/0.91;
			
			if (QE_perp>=1e-3) {
				list_A.push_back(QE_perp * total_bin_photons[mPMTTarget][PMTTarget]);
				list_R.push_back(R_array[mPMTTarget][PMTTarget]);
				//w is just for bookkepping of the 
				//configuration (i.e. which files we fit together)
				list_i.push_back(w);
				//Q is the data, A is what we will modify to have the correct prediction
				list_Q.push_back(total_bin_charge[mPMTTarget][PMTTarget]);
				std::cout << "Predicted charge without arb. scaling and with t<0.15ns " << QE_perp * total_bin_photons[mPMTTarget][PMTTarget] << " and real charge " << total_bin_charge[mPMTTarget][PMTTarget] << " in PMT " << PMTTarget << " of mPMT " << mPMTTarget << std::endl;
 			
			//	std::cout << " bin " << binTarget << " Hit PMT " << PMTTarget << " excat attenuation length length =" << truth_alpha(401.9, test_positions[0].abwff,test_positions[0].rayff);
			//	std::cout << " abwff " << test_positions[0].abwff << " rayff " << test_positions[0].rayff << " R = " <<   test_positions[0].R << " charge collected: " ;
			//	std::cout << total_bin_charge[binTarget][PMTTarget] << " total number of photons sent in this bin " << total_bin_photons[binTarget] << " prediction without attenuation " << QE_perp << std::endl;
		         }
		}//if we have more than Q_thresh hits in a given bin
	     }//run through all the PMTs
	}//run through all the mPMTs
	//The true scattering length that we know from the file
	trueScat = truth_alpha(401.9, test_positions[0].abwff,test_positions[0].rayff);
	//Keeping track of the files we fitted together
	config_number = std::stoi(configuration);
	w+=1;
        //file_ref->Close();
    } //finished reading all of the test maps
    
    //Here the fitting begins
    //std::cout << "Final config number = " << config_number << std::endl;
    const int nPars = 1; //the only parameter we fit is scattering length
    
// FOR NOW NOT FITTING BECAUSE CHISQ IS NOT YET WORKING    
    Chisq *chi = new Chisq(nPars);
    //In this case list_i is a index storing, later we will have as
    //many entries as our number of bins
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


