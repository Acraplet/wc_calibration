#include "./findBin.hh"
#include "./readTxtFileWithTime.hh"
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
    float R;
    int bin;
    char c;
    char * filename;

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
    }




	       //	= Form("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/test_hadded/ReferenceAttenuation_ALLhitPMT_PMT-basedBin_R%.2f.txt", R);
	std::vector<DataWithTime> test_positions = readTxtFileWithTime(filename);

	int nBins = 20;
	int nPMTs = 19;
	std::vector<std::vector<float>> total_bin_charge;
	std::vector<std::vector<int>> total_bin_photons;
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
		total_bin_photons.push_back(buf);
	}
	for (int i = 0; i < test_positions.size(); i++){
		DataWithTime refPoint = test_positions[i];
		//double time = refPoint.time;
		R = refPoint.R;
		double mPMT_pmt = refPoint.mPMT_pmt;
		//std::cout << "Bin " << refPoint.bin << " PMT " << refPoint.mPMT_pmt << total_bin_charge[int(refPoint.bin)][int(refPoint.mPMT_pmt)] << " " << total_bin_photons[int(refPoint.bin)][int(refPoint.mPMT_pmt)] << std::endl;
		total_bin_charge[int(refPoint.bin)][int(refPoint.mPMT_pmt)] = total_bin_charge[int(refPoint.bin)][int(refPoint.mPMT_pmt)] + float(refPoint.Q);
        	total_bin_photons[int(refPoint.bin)][int(refPoint.mPMT_pmt)] = total_bin_photons[int(refPoint.bin)][int(refPoint.mPMT_pmt)] + 1;
	}



	char *name = Form("./Maps/AbsorptionTable_PMT-basedBins.txt");
        std::ofstream outfile;
        outfile.open(name, std::ofstream::app);


	for (int bin=0; bin<nBins; bin++){
		outfile << R << " " << bin << " ";
		for (int pmt=0; pmt<=nPMTs; pmt++){
			if (total_bin_photons[bin][pmt] == 0 or total_bin_charge[bin][pmt] == 0.) outfile << 0 << " ";
			else {
				float ratio = total_bin_charge[bin][pmt]/total_bin_photons[bin][pmt];
				outfile << ratio << " ";
			}
			std::cout << "bin " << bin << " pmt " << pmt << " ratio " << total_bin_charge[bin][pmt] << " " << total_bin_photons[bin][pmt] << std::endl;
            		//outfile << "PMT " << i << " recieved " <<  total_bin_charge[i] << " of the " << total_bin_photons << " photons that were sent to bin " << bin << std::endl;
		}
		outfile << std::endl;
	}
        outfile.close();

}
