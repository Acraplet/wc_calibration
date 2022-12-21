//This is a cpp code reading the WCSim files counting the number of photons hitting the mPMT from all the angles
//It stores the information in a .txt file that has  x << " " << y << " " << z << " " << t << " " << p << " " << R << " " << nHitsTot << " " << n_entries where nEntries is the total number of photons that were simulated
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
// #include "../include/WCSimRootEvent.hh"

//Need to:
/*
 - open all the data with a given configurations
 - read the output giving out: R, theta, phi, (source x,y,z + dir x,y,z)
 - then compile everything into a 2D map
*/
const std::string ERR = color::RED_STR + "[ERROR]: " + color::RESET_STR;
const std::string WAR = color::GREEN_STR + "[WARNING]: " + color::RESET_STR;

void HelpMessage()
{
    std::cout << "USAGE: "
    << "make_mPMTmap" << "\nOPTIONS:\n"
    << "-f : Input file\n" << std::endl
    << "-o : Output file\n" << std::endl;
}

int main(int argc, char **argv){
    char * filename=NULL;
    char * outfilename=NULL;
    char c;

    while( (c = getopt(argc,argv,"f:o:b:s:e:l:r:p:c:w:hdtv")) != -1 ){//input in c the argument (-f etc...) and in optarg the next argument. When the above test becomes -1, it means it fails to find a new argument.
        switch(c){
            case 'f':
                filename = optarg;
                break;
        }
        switch(c){
            case 'o':
                outfilename = optarg;
                break;
        }
    }

    if (filename==NULL){
        std::cout << ERR << "Error, no input file: " << std::endl;
        HelpMessage();
        return -1;
    }

    if (outfilename==NULL){
        std::cout << WAR << "Warning, no output file given, using test.txt " << std::endl;
//         HelpMessage();
    }

    std::string x, y, z, t, p, R; //the cartesian and polar coordinate, careful R is from the surface of the mPMT dome and not from its centre of the mPMT
    int str; //this will help break down the string
    //To read the characteristics of the thing
    std::string T;
    std::stringstream X(filename);
    while(std::getline(X, T, '_')){
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

    for(int i=0; i<n_entries; i++){
        events->GetEntry(i);
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

    std::cout << "Total number of hits : "<< nHitsTot << std::endl;
    std::ofstream outfile;
    outfile.open(outfilename, std::ofstream::app);
    outfile <<  x << " " << y << " " << z << " " << t << " " << p << " " << R << " " << nHitsTot << " " << n_entries << std::endl;
}//End of main


