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

    while( (c = getopt(argc,argv,"R:b:")) != -1 ){//input in c the argument (-f etc...) and in optarg the next argument. When the above test becomes -1, it means it fails to find a new argument.
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
    }



	char * filename = Form("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/test_hadded/ReferenceAttenuation_ALLhitPMT_PMT-basedBin_R%.2f.txt", R);
	std::vector<DataWithTime> test_positions = readTxtFileWithTime(filename);

	int nBins = 19;
	double total_bin_charge[nBins];
    	double total_bin_photons = 0;
    	//initialisation of the bin-dependant counters
    	for (int k=0; k<=nBins; k++){
        	total_bin_charge[k] = 0.;
    	}
	for (int i = 0; i < test_positions.size(); i++){
		DataWithTime refPoint = test_positions[i];
		double time = refPoint.time;
		double mPMT_pmt = refPoint.mPMT_pmt;
		total_bin_charge[int(refPoint.mPMT_pmt)] += float(refPoint.Q);
        	total_bin_photons += 1;
	}

	char *name = Form("AbsorptionTable_PMT-basedBins.txt");
        std::ofstream outfile;
        outfile.open(name, std::ofstream::app);
	outfile << R << " ";
	for (int i=0; i<=nBins; i++){
            	//outfile << "PMT " << i << " recieved " <<  total_bin_charge[i] << " of the " << total_bin_photons << " photons that were sent to bin " << bin << std::endl;
		outfile << total_bin_charge[i]/total_bin_photons << " ";
	}
	outfile << std::endl;
        outfile.close();

}

