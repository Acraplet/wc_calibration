#include "./findBin.hh"
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

//this is a code that reads the reference root file and spits out how many PE are detected per photon sent perpendi
//cular to the PMT direction.  

int main(int argc, char **argv){
    float R;
    int bin;
    char c;
    char * filename = NULL;

    while( (c = getopt(argc,argv,"R:b:f:")) != -1 ){//input in c the argument (-f etc...) and in optarg the next argument. When the above test becomes -1, it means it fails to find a new argument.
        switch(c){
            case 'b':
                bin = std::stoi(optarg);
                break;
        }
        switch(c){
            case 'R':
                R = std::stof(optarg);
                break;
        }
        switch(c){
            case 'f':
                filename = optarg;
                break;
	}
	//might need to add an input of the  number of events 
    }

	if (filename == NULL) {
                std::cout << "Call the filename with the -f option" << std::endl;
        }


	       //	= Form("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/test_hadded/ReferenceAttenuation_ALLhitPMT_PMT-basedBin_R%.2f.txt", R);
	std::vector<DataWithTime> test_positions = readTxtFileWithTime(filename);

	int nBins = 20;
	int nPMTs = 19;
	std::vector<std::vector<float>> total_bin_charge;
	std::vector<int> total_bin_photons; //should be input instead of calculated
    	//double total_bin_photons[nBins];
    	//initialisation of the bin-dependant counters
	for (int l =0; l<nBins; l++){
		std::vector<float> bufDouble;
		std::vector<int> buf;
	       	for (int k = 1; k<=nPMTs+1; k++){
		       	buf.push_back(0);
		       	bufDouble.push_back(0.);
		}
		total_bin_charge.push_back(bufDouble);
		total_bin_photons = buf;
	}
	for (int i = 0; i < test_positions.size(); i++){
		DataWithTime refPoint = test_positions[i];
		//double time = refPoint.time;
		R = refPoint.R;
		double mPMT_pmt = refPoint.mPMT_pmt;
		//we need to read the source position even if there is nothing in there, that also means recon
		//structing the bin, for the expected number of photons
		//std::cout << "Bin " << refPoint.bin << " PMT " << refPoint.mPMT_pmt << total_bin_charge[int(refPoint.bin)][int(refPoint.mPMT_pmt)] << " " << total_bin_photons[int(refPoint.bin)][int(refPoint.mPMT_pmt)] << std::endl;
		//I am looking at the number of photons that are falling inside a bin and how many of them are going in each PMT otherwise it makes no sense, we could have 1 photon that happens to fall in PMT i but is in bin j and that would be a 100% efficiency - non-sense.
		total_bin_charge[int(refPoint.bin)][int(refPoint.mPMT_pmt)] = total_bin_charge[int(refPoint.bin)][int(refPoint.mPMT_pmt)] + float(refPoint.Q);
        	total_bin_photons[int(refPoint.bin)] = total_bin_photons[int(refPoint.bin)] + int(refPoint.nPhotons);
	}



	char *name = Form("./Maps/AbsorptionTable_PMT-basedBins_wholemPMT.txt");
        std::ofstream outfile;
        outfile.open(name, std::ofstream::app);

	//for now we are only simulating one quarter and then fill in the symetrical PMT
	//we ahve checked, sub-percent PMT to PMT efficiency difference 
	//int quarter[20] = {0, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 13, 14, 13, 14, 13, 14, 19}; 
	for (int bin=0; bin<nBins; bin++){
		outfile << R << " " << bin << " ";
		for (int pmt=0; pmt<=nPMTs; pmt++){
			float ratio = 0.;
			
			if (total_bin_photons[bin] == 0 or total_bin_charge[bin][pmt] == 0.){
				outfile << ratio << " ";
			}
		//in the case were we didn't simulate it -> want to use the matching one
			//if (bin == pmt and total_bin_photons[bin] == 0.) {
			//	ratio = total_bin_charge[quarter[bin]][quarter[pmt]]/total_bin_photons[quarter[bin]];
			//}
			else {
				ratio = total_bin_charge[bin][pmt]/total_bin_photons[bin];
				outfile << ratio << " ";
			}
			std::cout << "bin " << bin << " pmt " << pmt  << " ratio : " << ratio << std::endl;
            		//outfile << "PMT " << i << " recieved " <<  total_bin_charge[i] << " of the " << total_bin_photons << " photons that were sent to bin " << bin << std::endl;
		}
		outfile << std::endl;
	}
        outfile.close();

}

