//This is a cpp code reading the WCSim files counting the number of photons hitting the mPMT from all the angles

//This is to store for each source position the max number of charge when there is no refraction and no scattering
//it is important to only use this code with the correct reference file as these max charge are then used as the 
//amplitude for the absorption fits

//here the input file is the WCSim output root file directly 
// it is input with the -f argument

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
#include <TVector3.h>
#include <TRandom3.h>

#include "ColorOutput.hh"

const std::string ERR = color::RED_STR + "[ERROR]: " + color::RESET_STR;
const std::string WAR = color::GREEN_STR + "[WARNING]: " + color::RESET_STR;

void HelpMessage()
{
    std::cout << "USAGE: "
    << "make_mPMTmap" << "\nOPTIONS:\n"
    << "-f : Input file\n" << std::endl;
}

int main(int argc, char **argv){
    char * filename=NULL;
    char c;

    while( (c = getopt(argc,argv,"f:o:b:s:e:l:r:p:c:w:hdtv")) != -1 ){
	    //input in c the argument (-f etc...) and in optarg the next argument. 
	    //When the above test becomes -1, it means it fails to find a new argument.
        switch(c){
            case 'f':
                filename = optarg;
                break;
        }
    }

    if (filename==NULL){
        std::cout << ERR << "Error, no input file: " << std::endl;
        HelpMessage();
        return -1;
    }

    std::string x, y, z, t, p, R; 
    //the cartesian and polar coordinate, careful R is from the surface of the mPMT dome and not from its centre of the mPMT
    int str; //this will help break down the string
    //To read the characteristics of the thing
    std::string T;
    std::string abwff;
    std::string rayff;
    std::stringstream X(filename);
    while(std::getline(X, T, '_')){
        if (!T.find("Absff")){
            str = T.find("Absff");
            abwff = T.erase(0,str+5);
        }
        if (!T.find("Rayff")){
            str = T.find("Rayff");
            rayff = T.erase(0,str+5);
        }
        if (!T.find("x")){
            str = T.find("x");
            x = T.erase(0,str+1);
        }
        if (!T.find("y")){
            str = T.find("y");
            y = T.erase(0,str+1);
        }
        if (!T.find("z")){
            str = T.find("z");
            z = T.erase(0,str+1);
        }
        if (!T.find("t")){
            str = T.find("t");
            t = T.erase(0,str+1);
        }
        if (!T.find("p")){
            str = T.find("p");
            p = T.erase(0,str+1);
        }
        if (!T.find("R")){
            str = T.find("R");

            R = T.erase(0,str+1);
        }
    }
    std::cout << "Source position : " << x << " " << y << " " << z << " " << t << " " << p << " " << R << std::endl;
    std::cout << "Absorption : " << abwff << " Scattering : " << rayff << std::endl;
    
    //In case we made a mistake in selecting the reference code
    if (abwff <= 10 or rayff <= 10){
	    std::cout << "Please use files that have a infinite attenuation and scattering length as reference" << std::endl;
	    throw;
    }
    //so far we are looking at the raw Cherenkov hits instead of the digitised ones
    TFile *infile = new TFile(filename, "READ");
    TTree *events = (TTree*)infile->Get("CherenkovHits");

    int event;
    int NHits;
    int mPMT;
    int mPMT_pmt;
    int nHitsTot = 0;

    events->SetBranchAddress("Event", &event);
    events->SetBranchAddress("NHits", &NHits);
    events->SetBranchAddress("mPMT", &mPMT);
    events->SetBranchAddress("mPMT_pmt", &mPMT_pmt);

    int n_entries = events->GetEntries();
    //Read the root file
    for(int i=0; i<n_entries; i++){
        events->GetEntry(i);
	//looking at the 58th mPMT but WCSim bug means that it is split over two mPMT number...
	if (mPMT == 58){
		if (mPMT_pmt != 19){
        		nHitsTot += NHits;
		}
	}
	if (mPMT == 59){
                if (mPMT_pmt == 19){
                        nHitsTot += NHits;
                }
        }
    }
    //print out some values
    std::cout << "Total number of hits : "<< nHitsTot << std::endl;

    const char* theta = t.c_str();
    const char* phi = p.c_str();
    const char* dist = R.c_str();

    //This is the file onto which we output
    std::string onePosition_outfile;
    if(std::getenv("WCCALIB")) onePosition_outfile = Form("%s/mPMTmapping/Maps/maps_Reference/OnePosition_theta%s_phi%s_R%s.txt", std::getenv("WCCALIB"), theta, phi, dist);
    else onePosition_outfile = Form("./Maps/maps_Reference/OnePosition_theta%s_phi%s_R%s.txt", theta, phi, dist);
    std::cout << onePosition_outfile << std::endl;
    std::ofstream onePosition;
    onePosition.open(onePosition_outfile, std::ofstream::app);
    onePosition << nHitsTot << " " << n_entries << std::endl;
}
