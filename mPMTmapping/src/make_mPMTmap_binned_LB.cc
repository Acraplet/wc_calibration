//This is a cpp code reading the WCSim files counting the number of photons hitting the mPMT from all the angles
//It stores the information in a .txt file that has  x << " " << y << " " << z << " " << t << " " << p << " " << R << " " << nHitsTot << " " << n_entries and bins where nEntries is the total number of photons that were simulated
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
#include "findBin.hh"

//#include "ColorOutput.hh"
// #include "../include/WCSimRootEvent.hh"

//Need to:
/*
 - open all the data with a given configurations
 - read the output giving out: R, theta, phi, (source x,y,z + dir x,y,z)
 - then compile everything into a 2D map
*/
//const std::string ERR = color::RED_STR + "[ERROR]: " + color::RESET_STR;
//const std::string WAR = color::GREEN_STR + "[WARNING]: " + color::RESET_STR;

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
        std::cout << "Error, no input file: " << std::endl;
        HelpMessage();
        return -1;
    }

    if (outfilename==NULL){
    	outfilename = Form("./Maps/output.txt");
        std::cout << "Warning, no output file given, using " << outfilename << std::endl;
//         HelpMessage();
    }

    std::string x, y, z, t, p, R_buf, abwff, rayff; //the cartesian and polar coordinate, careful R is from the surface of the mPMT dome and not from its centre of the mPMT
    int str; //this will help break down the string
    //To read the characteristics of the thing
    float R;
    std::string T;
    std::stringstream X(filename);
    while(std::getline(X, T, '_')){
        if (!T.find("Absff")){
            str = T.find("Absff");
            abwff = T.erase(0,str+5);
        }
        if (!T.find("Rayff")){
            str = T.find("Rayff");
            rayff = T.erase(0,str+5);
            // Deletes str+5 characters from index number 0
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
  //      if (!T.find("t")){
  //          str = T.find("t");
  //          t = T.erase(0,str+1);
  //      }
  //      if (!T.find("p")){
  //          str = T.find("p");
  //          p = T.erase(0,str+1);
  //      }
        if (!T.find("R")){
            str = T.find("R");
            R_buf = T.erase(0,str+1);
	    R = std::stof(R_buf);
        }
    }

     //we need to read the position and direction of the mPMTs to calculate the incoming theta and phi angles

    std::cout << "Source position : " << x << " " << y << " " << z << " " << t << " " << p << " " << R << " " << abwff << " " << rayff <<  std::endl;

    TFile *infile = new TFile(filename, "READ");
    TTree *events = (TTree*)infile->Get("CherenkovHits");
    TTree *source_position = (TTree*)infile->Get("EventInfo");

    int event;
    int NHits;
    int mPMT;
    int mPMT_pmt;
    
    float PMT_x;
    float PMT_y;
    float PMT_z;    

    int nHitsTot = 0;
    float Time;

    int run;
    double source_x;
    double source_y;
    double source_z;

    events->SetBranchAddress("Event", &event);
    events->SetBranchAddress("NHits", &NHits);
    events->SetBranchAddress("mPMT", &mPMT);
    events->SetBranchAddress("mPMT_pmt", &mPMT_pmt);
    events->SetBranchAddress("Time", &Time);

    events->SetBranchAddress("PMT_x", &PMT_x);
    events->SetBranchAddress("PMT_y", &PMT_y);
    events->SetBranchAddress("PMT_z", &PMT_z);

    source_position->SetBranchAddress("Run", &run);
    source_position->SetBranchAddress("Vertex_x", &source_x);
    source_position->SetBranchAddress("Vertex_y", &source_y);
    source_position->SetBranchAddress("Vertex_z", &source_z);

    //MPMT test = getMPMT(65);
    //std::cout << test.ID << " is ID" << std::endl;

    int n_entries = events->GetEntries();

    std::ofstream outfile; 
    outfile.open(outfilename, std::ofstream::app);
    for(int i=0; i<n_entries; i++){
	std::cout << i << " nentries " << n_entries << std::endl;
        events->GetEntry(i);
	source_position->GetEntry(i);
	//double[3] source_xyz = {source_x, source_y, source_z};
	//double[3] PMT_xyz = {PMT_x, PMT_y, PMT_z};
	//std::cout << NHits << " " << source_z << " " << PMT_x  << std::endl;
    //this is for the _flat issue
    	//if (NHits!=0) {
//	if (mPMT == 0 and mPMT_pmt != 19){
//			std::cout << "mPMT: " << mPMT << std::endl; 
//        		nHitsTot += NHits;
//			//need to find the theta and phi of each of the source positions
//			std::vector<double> theta_phi = findThetaPhi(source_x, source_y, source_z);
//		        Bin bin = findPMTBin(theta_phi[0], theta_phi[1]);		
//			//std::cout << mPMT << " PMT " << mPMT_pmt << " Q " << NHits << " T " << Time  << " x " << source_x << " y " << source_y << " z " << source_z << " theta " << theta_phi[0] << " phi " << theta_phi[1] << " prefered bin " << bin.ID << " bin difference " << bin.ID-mPMT_pmt << std::endl;
//    			//std::cout << R << std::endl;	
//			double TOF_Time = TOFCorrect(Time, R);
////			char *name = Form("ReferenceAttenuation_ALLhits_PMT-basedBin%i_R%.2f.txt", bin.ID, R);
////			std::ofstream outfile; 
////			outfile.open(name, std::ofstream::app);
//   		 	outfile << mPMT << " " << mPMT_pmt << " " <<  source_x << " " << source_y << " " << source_z << " " << theta_phi[0] << " " << theta_phi[1] << " " << R << " " << NHits << " " << TOF_Time << " " << abwff << " " << rayff << " " << bin.ID << std::endl;
////			outfile.close();
//			//delete [] name;
//	}
//	
//	
//	else if (mPMT == 0 and mPMT_pmt == 19){
//        		nHitsTot += NHits;
//			//need to find the theta and phi of each of the source positions
//			std::vector<double> theta_phi = findThetaPhi(source_x, source_y, source_z);
//		        Bin bin = findPMTBin(theta_phi[0], theta_phi[1]);		
//			double TOF_Time = TOFCorrect(Time, R);
//			//std::cout << mPMT << " PMT " << mPMT_pmt << " Q " << NHits << " T " << Time  << " x " << source_x << " y " << source_y << " z " << source_z << " prefered bin " << bin.ID << " bin difference " << bin.ID-mPMT_pmt << std::endl;
//    			
//   		 	outfile << mPMT << " " << mPMT_pmt << " " <<  source_x << " " << source_y << " " << source_z << " " << theta_phi[0] << " " << theta_phi[1] << " " << R << " " << NHits << " " << TOF_Time << " " << abwff << " " << rayff << " " << bin.ID << std::endl;
//    //            }
//	}
	//if (mPMT == 58 or mPMT == 59){
	//this is in the case where we have different mPMT than the one at the bottom
		R = TMath::Sqrt((source_x - PMT_x) * (source_x - PMT_x) + (source_y - PMT_y) * (source_y - PMT_y) +  (source_z - PMT_z) * (source_z - PMT_z));
		if (mPMT_pmt == 19 and NHits>0){ //because of the issue in numbering, this one is always central one
		//we need to calculate the PMT-source distance
			double TOF_Time = TOFCorrect_usingXYZ(Time, source_x, source_y, source_z, PMT_x, PMT_y, PMT_z);
			std::vector<double> theta_phi = findThetaPhi_usingXYZ(source_x, source_y, source_z, mPMT+1);
			Bin bin = findPMTBin(theta_phi[0], theta_phi[1],mPMT);

			outfile << mPMT << " " << mPMT_pmt << " " <<  source_x << " " << source_y << " " << source_z << " " << theta_phi[0] << " " << theta_phi[1] << " " << R << " " << NHits << " " << TOF_Time << " " << abwff << " " << rayff << " " << bin.ID << std::endl;
	    	}
		if (mPMT_pmt <= 18 and NHits>0){ //because of the issue in numbering, this one is always central one
			//we need to calculate the PMT-source distance
			//double TOF_Time = 0; 
			double TOF_Time = TOFCorrect_usingXYZ(Time, source_x, source_y, source_z, PMT_x, PMT_y, PMT_z);
			std::vector<double> theta_phi = findThetaPhi_usingXYZ(source_x, source_y, source_z, mPMT);
			Bin bin = findPMTBin(theta_phi[0], theta_phi[1], mPMT);
			std::cout << "mPMT " << mPMT << " mPMT_pmt " << mPMT_pmt<<" Theta " << theta_phi[0] << " phi " << theta_phi[1] << " reco_bin " << bin.ID << std::endl;

			outfile << mPMT << " " << mPMT_pmt << " " <<  source_x << " " << source_y << " " << source_z << " " << theta_phi[0] << " " << theta_phi[1] << " " << R << " " << NHits << " " << TOF_Time << " " << abwff << " " << rayff << " " << bin.ID << std::endl;
		}

    //}
    }
    outfile.close();

    //std::cout << "Total number of hits : "<< nHitsTot << std::endl;
    //std::ofstream outfile;
    //outfile.open(outfilename, std::ofstream::app);
    //outfile <<  x << " " << y << " " << z << " " << t << " " << p << " " << R << " " << nHitsTot << " " << n_entries << " " << abwff << " " << rayff << std::endl;
}//End of main


