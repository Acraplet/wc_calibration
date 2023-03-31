//This code is extracting the scattering length corresponding to one or more maps that have the same rayff
//and ideally abwff = inf
//For this is is looking up in reference_maps what the y values of the 5 nodes making up the charge vs scattering length fake-data spline. This number of 5 nodes can be changed but for this we need to re-run the Fitter. 
//Then the chi is minimised between the observed data and y_pred = Fake_Data_Spline(scattering length pred)
//expected scattering length. So far this is position dependant but it will have to be upgraded to a binned approach
#include "../include/findBin.hh"
#include "../include/readTxtFile.hh"
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


int main(int argc, char **argv){

    std::vector<double> list_Q, list_i;
    std::vector<std::vector<double>> list_ynodes, list_xnodes;
    char* R_test;
    char* Q_test;
    char* theta_test;
    double  theta_test_num, Q_test_num;
    char* phi_test;
    //The config number is a rough way to output the distances that we have looked out 
    //useful for run comparisions
    long int config_number = 0.0;
    double w = 0.0;
    //TODO: change for a config file where the max number of bins is an entry
    int nBins = 800;  //total max number of possible bins
    int Q_thresh = 10; //minimum number of photons in a given test bin to be included in the fit
    double spline_min = 925.;  //the precision and range of the 
    double spline_max = 6000.; //spline we are drawing from the reference 2D distribution
    double spline_increment = 100.;
    int bin_min = 10; //for now the first bins have been messed up so we ignore them but 
		      //eventually will change that
    double initGuessScat = 1500.;
    double trueScat;

    std::cout << "The total number of files to read together is: " << argc << std::endl;
    for (int f=1; f<= argc-1; f++){
	//This is the test file we are going to extract the scattering length out of
        char * filename = Form("./Maps/maps_txtFiles/mPMT_map_ID%s.txt",argv[f]);
        //in these double we are storing the total charge, one entry per bin
	double total_bin_charge[nBins];
        double total_bin_photons[nBins];
	//initialisation of the bin-dependant counters
	for (int k=0; k<=nBins; k++){
		total_bin_charge[k] = 0.;
	        total_bin_photons[k] = 0.;	

	}
	//This is reading in the test file 
        std::vector<Data> test_positions = readTxtFile(filename);
	std::cout << test_positions.size() << std::endl;
        for (long unsigned i = 0; i < test_positions.size(); i++){
            //read each source position in the file one by one
            Data pos =  test_positions[i];
            //Then we check which bin it belongs to
            Bin closestBin = findBin(pos.theta, pos.phi);
	    //std::cout << "The closest bin is " << closestBin.ID<< std::endl;
	    //add to the charge the fractionnal charge collected at this position
            //we do not use direct charge because the comparision has to be made
            //with respect to 1000 photons for now
            total_bin_charge[closestBin.ID] += float(pos.Q);
            total_bin_photons[closestBin.ID] += float(pos.nEvents);
        }

	for (int binTarget=bin_min; binTarget <= nBins; binTarget++){
        	total_bin_charge[binTarget] = total_bin_charge[binTarget] / total_bin_photons[binTarget] * 1000;
		if (total_bin_charge[binTarget]>= Q_thresh){
			std::cout << "Fitting bin " << binTarget << std::endl;
		        int n = 0;
		        double x, y;
		        double ref_info; //empty variable to read the file with
		        int count = 0;
		        std::vector<double> bufY;
			std::vector<double> bufX;
			
			//Now fitting multiple bins together
			//now actually read the 2D entries into a histogram and then extract the position of the nodes from the Delaunay interpolation 
			TGraph2D * f2 = new TGraph2D();
		        gStyle->SetPalette(1);
			//This is where the 2D map is saved
			const char* fimpName = Form("./Maps/2D_ref_maps/all_ref_files_bin%i.txt",binTarget);
		        std::ifstream in(fimpName);
		       //The reference datasets to use for comparision
		        std::vector<double> x_vector, y_vector, z_vector;

		        while ((in >> ref_info)) {
		                if (count %3 == 0) {
		                        x_vector.push_back(truth_alpha(401.9,10e10,ref_info));
		                }//scatering len
		                if (count %3 == 1) {
		                        y_vector.push_back(ref_info);
		                }//R i.e. mPMT-source distance
		                if (count %3 == 2){
		                        z_vector.push_back(ref_info);
		                }//Q per 1000 photons
		                count+=1;
		        }

			//Add points to the bin reference scattering histogram
		        for (long unsigned int N=0; N<x_vector.size(); N++) {
		                f2->SetPoint(N, x_vector[N], y_vector[N], z_vector[N]);
		        }
			
			//the x values of the splines to which we will minimise, the y is drawn from the interpolation
			for (int x_i = spline_min; x_i < spline_max; x_i += spline_increment){
				bufX.push_back(x_i);
			}
			//the y values of the spline are taken from the interpolation of the 2D surface	
			for (long unsigned int count=0; count < bufX.size(); count++){
		            bufY.push_back(f2->Interpolate(bufX[count],test_positions[0].R ));
		        }
			//Add the relevant spline coordinates to the vector of splines that we will fit together
		        list_xnodes.push_back(bufX);
		        list_ynodes.push_back(bufY);
			//w is just for bookkepping of the configuration (i.e. which files we fit together)
		        list_i.push_back(w);
		        list_Q.push_back(total_bin_charge[binTarget]);

        		std::cout << "File " << argv[f] << " bin " << binTarget << " excat scattering length =" << truth_alpha(401.9, 10e10,test_positions[0].rayff);
		        std::cout << " " << test_positions[0].rayff  << " R = " <<   test_positions[0].R << " charge for 1000 photons: " ;
			std::cout << total_bin_charge[binTarget]  << " config " << pow(10., w) << std::endl;
		}//if we have more than Q_thresh hits in a given bin
	}//run through all the bins
	//The true scattering length that we know from the file
        trueScat = truth_alpha(401.9, 10e10,test_positions[0].rayff);
	//Keeping track of the files we fitted together
        double temp = pow(10., w);
        temp *= std::stoi(argv[f]) % 10;
        config_number += (double) temp;
        w+=1;
        //file_ref->Close();
    } //finished reading all of the test maps
    //Here the fitting begins
    std::cout << "Final config number = " << config_number << std::endl;
    const int nPars = 1; //the only parameter we fit is scattering length
    Chisq *chi = new Chisq(nPars);
    //In this case list_i is a index storing, later we will have as many entries as our number of bins
    chi->setData(list_i, list_Q);
    chi->setRef_scat(list_xnodes, list_ynodes);
    ROOT::Math::Functor functor(chi, &Chisq::fitter_rayleigh, nPars);
    ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    min->SetStrategy(3);
    min->SetFunction(functor);
    min->SetMaxFunctionCalls(10000);
    min->SetVariable(0, "scattering_length", initGuessScat, 0.01);

    min->Minimize();
    min->PrintResults();
    const double * res_scat = min->X();
    const double * err_scat = min->Errors();


    //Here output some numbers for easier analysis
    std::ofstream outfile;
    outfile.open("ScatteringLengthEstimation_2Dinterpolation_withText.txt", std::ofstream::app);
    outfile << "True: " << trueScat << " config: " << config_number << " reco: " << res_scat[0] << " +/- " << err_scat[0];
    outfile << " Q_thresh " << Q_thresh << " initGuessScat " << initGuessScat << " spline_min " << spline_min << " spline_max " << spline_max;
    outfile << " spline_increment " << spline_increment << std::endl;
    outfile.close();
   
    //and without text for analysis
    //std::ofstream outfile;
    outfile.open("ScatteringLengthEstimation_2Dinterpolation_withoutText.txt", std::ofstream::app);
    outfile << " " << trueScat << " " << config_number << " " << res_scat[0] << " " << err_scat[0];
    outfile << " " << Q_thresh << " " << initGuessScat << " " << spline_min << " " << spline_max;
    outfile << " " << spline_increment << std::endl;
    outfile.close();

}




//         if (newfile.is_open()){   //checking whether the file is open
//             std::string tp;
//             while(getline(newfile, tp)){  //read data from file object and put it into string.
//                 //this is for each test source position (theta, phi)
//                 double w = 0.;
//                 char *ptr;
//                 //convert to s char the string of the line we are extracting
//                 char* character = std::strcpy(new char[tp.length() + 1], tp.c_str());
//                 ptr = std::strtok(character, " "); //split the string after the blanks
//                 // use while loop to check ptr is not null
//                 int i = 0;
//
//                 while (ptr != NULL)
//                     //loop over the characteristics of the given position
//                 {
//                     if (i==3){
//                         theta_test = ptr;
//                         std::string fs(ptr);
//                         theta_test_num=std::stof(fs);
//                     }
//                     if (i==4){
//                         phi_test = ptr;
//                     }
//                     if (i==5){
//                         R_test = ptr;
//                     }
//                     if (i==6){
//                         Q_test = ptr;
//                         std::string fs(ptr);
//                         Q_test_num=std::stof(fs);
//                     }
//                     ptr = std::strtok (NULL, " ");
//                     i +=1;
//                 }
//Now we have the coordinate of the position we are looking for
//Next: open the reference txt file for this position and from there extract the fitting
//information we are looking for here it is the y value corresponding to the nodes of the fake

//Need to fill the correct bin

//Get the closest bin to that given point


