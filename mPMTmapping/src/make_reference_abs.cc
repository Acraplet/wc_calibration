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
    int nBins = 800;
    char * filename = Form("./Maps/maps_txtFiles/mPMT_map_ID%s.txt",argv[1]);
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
    for (int i = 0; i < test_positions.size(); i++){
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


        //need to extrapolate the map to the other four quarters
        //note that if more than one exact quarter has been simulated we will have an average of the
        //different overlapping bins
        closestBin = findBin(pos.theta, TMath::Pi() - pos.phi);
        total_bin_charge[closestBin.ID] += float(pos.Q);
        total_bin_photons[closestBin.ID] += float(pos.nEvents);

        closestBin = findBin(pos.theta, TMath::Pi() + pos.phi);
        total_bin_charge[closestBin.ID] += float(pos.Q);
        total_bin_photons[closestBin.ID] += float(pos.nEvents);

        closestBin = findBin(pos.theta, -pos.phi);
        total_bin_charge[closestBin.ID] += float(pos.Q);
        total_bin_photons[closestBin.ID] += float(pos.nEvents);

    }

    for (int i=0; i<=nBins; i++){
        if (total_bin_photons[i]!=0){
            std::cout << total_bin_charge[i]/total_bin_photons[i] << std::endl;
        }
        else {std::cout << 0 << std::endl;};
    }
    return 0;
}
