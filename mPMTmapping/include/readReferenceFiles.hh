//This is a code that returns the reference values of number of hits in that PMT with respect to 
//number of photons in that bin 
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <TROOT.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TGraph.h>
#include <stdexcept>
#include <TVector3.h>
#include <TRandom3.h>

double readAbsorptionRef(int bin, int PMT, double R){
	//this is where lies the table of reference mPMT with R bin PMT1 PMT2 ....
	std::fstream position_file;
	char* filename = "Maps/AbsorptionTable_PMT-basedBins.txt";
	
	position_file.open(filename, std::ios::in);
	if (position_file.is_open()){   //checking whether the file is open
        	std::string tp_ref;
		while(getline(position_file, tp_ref)){
		    char *ptr_ref;
		    char* character_ref = std::strcpy(new char[tp_ref.length() + 1], tp_ref.c_str());
	            ptr_ref = std::strtok(character_ref, " "); //split the string after the blanks
		    //the first entry is the R - it needs to be correct
		    std::string fs(ptr_ref);
	            double R_ref=std::stod(fs);
		    //no point readinf the rest of the line if it is not the correct R
		    if (R_ref == R){
			    //the next entry is the bin number which also needs to be correct before we go in
			    ptr_ref = std::strtok (NULL, " ");
			    fs = ptr_ref; //need to convert of a std::string
			    int bin_ref = std::stoi(fs);
			    if (bin_ref == bin){
				std::cout << "This is bin " << bin_ref << " at distance " << R << std::endl;
				//now we read the rest of the lines
			    	ptr_ref = std::strtok (NULL, " ");
		    		int j = 0; //the first entry is PMT 1
			    	while (ptr_ref != NULL){
					fs = ptr_ref;
					if (j == PMT){
						std::cout << "fractional value " << std::stod(fs) << std::endl;
					       	return std::stod(fs);
					}
			    		ptr_ref = std::strtok (NULL, " ");
					j += 1;
			    	}
			    }
		    }
		}
	}
	//if we haven't found the correct entry in the reference file - trow an exception
	throw std::invalid_argument( Form("Couldn't find a reference entry in file %s corresponding to bin %i and PMT %i with R %.2f", filename, bin, PMT, R) );
	

}
