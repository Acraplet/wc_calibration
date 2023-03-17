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
    long int config_number = 0.0;
    double w = 0.0;
    double true_attLen;

    std::cout << "The total number of files to read together is: " << argc << std::endl;
    for (int f=1; f<= argc-1; f++){
//         std::fstream newfile;
	//This is the test file we are going to extract the scattering length out of
//         newfile.open(,std::ios::in); //open a file to perform read operation using file object
//         std::cout << argv[f] << std::endl;
        char * filename = Form("./Maps/maps_txtFiles/mPMT_map_ID%s.txt",argv[f]);
        
	double total_bin_charge[800];
        double total_bin_photons[800];
	//initialisation of the bin-dependant counters
	for (int k=0; k<=800; k++){
		total_bin_charge[k] = 0.;
	        total_bin_photons[k] = 0.;	

	}

        //For the next step: not useful just yet as we only have ref for bin 24
//         std::vector<int> binID;
//         std::vector<double> total_charge;
        std::vector<Data> test_positions = readTxtFile(filename);
	std::cout << test_positions.size() << std::endl;
        for (int i = 0; i < test_positions.size(); i++){
            //read each position in the file one by one
            Data pos =  test_positions[i];
            //Then we check which bin it belongs to
            Bin closestBin = findBin(pos.theta, pos.phi);
	    //std::cout << "The closest bin is " << closestBin.ID<< std::endl;
	    //std::cout << "Made it here " << i << std::endl;
            //For now I am only looking at Bin 24 for my fit but later we will have a
            //vector of bins with eveything in there
//            if (closestBin.ID != binTarget) continue;
            //add to the charge the fractionnal charge collected at this position
            //we do not use direct charge because the comparision has to be made
            //with respect to 1000 photons for now
            total_bin_charge[closestBin.ID] += float(pos.Q);
            total_bin_photons[closestBin.ID] += float(pos.nEvents);
        }

	for (int binTarget=10; binTarget <= 800; binTarget++){
        	total_bin_charge[binTarget] = total_bin_charge[binTarget] / total_bin_photons[binTarget] * 1000;
		if (total_bin_charge[binTarget]>= 10){
			std::cout << "Fitting bin " << binTarget << std::endl;
		//Now fitting multiple bins together
	//now actually read the 2D entries into a histogram and then extract the position of the nodes from the Delaunay interpolation 
			TGraph2D * f2 = new TGraph2D();
		        gStyle->SetPalette(1);
		        int n = 0;
		        double x, y;
			
			const char* fimpName = Form("./Maps/2D_ref_maps/all_ref_files_bin%i.txt",binTarget);
		        std::ifstream in(fimpName);
		        double ref_info;
		        int count = 0;
		       //The reference datasets to use for comparision
		        std::vector<double> x_vector, y_vector, z_vector;
		
		        while ((in >> ref_info)) {
		                if (count %3 == 0) {
		                        x_vector.push_back(truth_alpha(401.9,10e10,ref_info));
		                }//scat len
		                if (count %3 == 1) {
		                        y_vector.push_back(ref_info); //R
		                }
		                if (count %3 == 2){
		                        z_vector.push_back(ref_info); //Q per 1000 photons
		                }
		                count+=1;
		        }
		
		        for (int N=0; N<x_vector.size(); N++) {
		                double x = x_vector[N];
		                double y = y_vector[N];
		                double z_true = z_vector[N];
		                f2->SetPoint(N, x, y, z_true);
		        }
			
		        std::vector<double> bufY;
			std::vector<double> bufX;
			//the x values of the splines to which we will minimise, the y is drawn from the interpolation
			for (int x_i = 925; x_i < 6000; x_i += 100){
				bufX.push_back(x_i);
			}
		
			for (int count=0; count < bufX.size(); count++){
		//             std::cout << count << std::endl;
		            //std::cout << bufX[count] << " " << test_positions[0].R << " " << f2->Interpolate(bufX[count],test_positions[0].R) << std::endl;
		            bufY.push_back(f2->Interpolate(bufX[count],test_positions[0].R ));
		        }
		        list_xnodes.push_back(bufX);
		        list_ynodes.push_back(bufY);
		        list_i.push_back(w);
		        list_Q.push_back(total_bin_charge[binTarget]);
        		std::cout << "File " << argv[f] << " bin " << binTarget << " excat scattering length =" << truth_alpha(401.9, 10e10,test_positions[0].rayff) << " " << test_positions[0].rayff  << " R = " <<   test_positions[0].R << " charge for 1000 photons: " << total_bin_charge[binTarget]  << " config " << pow(10., w) << std::endl;
		}//if we have more than 100 hits in a given bin
	}//run through all the bins

	//for (int check=0; check < bufY.size(); check++ ){
	//	std::cout << "Test charge "<< total_bin_charge << " " << " Nodes y, x "<< bufY[check] << " " << bufX[check] <<std::endl;
	//}

        true_attLen = truth_alpha(401.9, 10e10,test_positions[0].rayff);

        double temp = pow(10., w);
        std::cout << temp << std::endl;
        temp *= std::stoi(argv[f]) % 10;
        std::cout << temp << std::endl;
        std::cout << config_number << std::endl;
        config_number += (double) temp;


        w+=1;
        //file_ref->Close();
    } //finished reading all of the test maps
    //Here the fitting begins

    std::cout << "Config number = " << config_number << std::endl;
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
    min->SetVariable(0, "scattering_length", 1500, 0.01);

    min->Minimize();
    min->PrintResults();
    const double * res_scat = min->X();
    const double * err_scat = min->Errors();


    	//Here output some numbers for easier analysis
    std::ofstream outfile;
    outfile.open("Bin24_ScatteringLengthEstimation_SimpleSpline_testRuns_2Dbin83.txt", std::ofstream::app);
    outfile << true_attLen << " " << config_number << " " << res_scat[0] << " " << err_scat[0] << std::endl;

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


