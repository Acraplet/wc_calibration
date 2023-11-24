//This is a cpp code reading the WCSim files and calculating the solid angle of each PMT, 
//saving that into a text file 
//and mutlipy this value by the total number of photons 
//simulated to predict the total number of photons that will reach each PMT 
//This assumes a uniformly emitting sphere, TODO: update code to be able to factor in the actualy source profile
//Also add a condition to predict that no photons leave from the shaded region of the LB (will have to calculate 
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <string>
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
#include <stdexcept>

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

    while( (c = getopt(argc,argv,"f:o:b:s:e:l:r:p:c:w:hdtv")) != -1 ){//input in c the argument (-f etc...) 
	//and in optarg the next argument. When the above test becomes -1, it means it fails to find a new argument.
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
    	outfilename = Form("./Maps/predictionNumberOfPhotons.txt");
        std::cout << "Warning, no output file given, using " << outfilename << std::endl;
//         HelpMessage();
    }

    //technically we do not need the details of the run but it doesn't hurt to keep it 
    std::string x, y, z, t, p, R_buf, abwff, rayff; 
    //the cartesian and polar coordinate, careful R is from the surface of the mPMT dome 
    //and not from its centre of the mPMT
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
    TTree *geometry = (TTree*)infile->Get("Geometry");

    int flag_int = -9999;
    double flag = -9999.;	

    int nPMT = 2014; 

    int event=flag_int;
    //std::vector<int> mPMT;
    Int_t mPMT[nPMT];
    Int_t mPMT_pmt[nPMT]; 
    Double_t PMT_x[nPMT];
    Double_t PMT_y[nPMT];
    Double_t PMT_z[nPMT];    
    Double_t dirPMT_x[nPMT];
    Double_t dirPMT_y[nPMT];
    Double_t dirPMT_z[nPMT];    

    int run;
    double source_x=flag;
    double source_y=flag;
    double source_z=flag;

    //the theta span of the shadowed region at the top of the LB
    double thetaDarkCone = 160 * TMath::Pi() / 180; 
    std::vector<double> LBaxis = {0, 1, 0};


    double solidAngle;

    events->SetBranchAddress("Event", &event);
    //events->SetBranchAddress("NHits", &NHits);
    //issue the geometry branch have very different sizes
    geometry->SetBranchAddress("mPMT", &mPMT);
    geometry->SetBranchAddress("mPMT_pmt", &mPMT_pmt);
//    //events->SetBranchAddress("Time", &Time);
//
//    //now we are looking at all the PMTs, not only the ones that were hit 
    geometry->SetBranchAddress("x", &PMT_x);
    geometry->SetBranchAddress("y", &PMT_y);
    geometry->SetBranchAddress("z", &PMT_z);
    geometry->SetBranchAddress("direction_x", &dirPMT_x);
    geometry->SetBranchAddress("direction_y", &dirPMT_y);
    geometry->SetBranchAddress("direction_z", &dirPMT_z);

    source_position->SetBranchAddress("Run", &run);
    source_position->SetBranchAddress("Vertex_x", &source_x);
    source_position->SetBranchAddress("Vertex_y", &source_y);
    source_position->SetBranchAddress("Vertex_z", &source_z);

    //we are looking at each PMT
    int n_entries = geometry->GetEntries();
    
    //the total number of photons that was sent is
    int n_photons_tot = events->GetEntries();
    //the number of photns per unit solid angle that is not in the dark cone
    double solidAngleEmisson = 2 * TMath::Pi() * (1 - TMath::Cos(TMath::Pi() - thetaDarkCone));
    double densityPhotons = n_photons_tot / (solidAngleEmisson);

    source_position->GetEntry(0); //technically the source doesn't move, we can keep it as 1, eventually, 
    //if we have more than one source position in the same root file, we might need to have an extra for loop
    //over the source positions 
    geometry->GetEntry(0);
    std::ofstream outfile; 
    outfile.open(outfilename, std::ofstream::app);
    for(int i=0; i<nPMT; i++){
	//std::cout << i << " number of PMTs " << n_entries << " total number of photons " << n_photons_tot<< std::endl;
    	//std::cout << mPMT[35] << std::endl;
        //geometry->GetEntry(0);
	//calculate the absolute distance between the source and the PMT
	//std::cout << "PMTxyz" << PMT_x[i] << " " << PMT_y[i] << " "<<PMT_z[i] << " " << std::endl;
	//std::cout << "sourcexyz" << source_x << " " << source_y << " "<<source_z << " " << std::endl;
	R = TMath::Sqrt((source_x - PMT_x[i]) * (source_x - PMT_x[i]) + (source_y - PMT_y[i]) * (source_y - PMT_y[i]) +  (source_z - PMT_z[i]) * (source_z - PMT_z[i]));
	//std::cout << "aLL IS OK SO FAR" << std::endl;
	//prepare the command for calculating (using a python code) the solid angle correcponding to this source 
	//position and PMT info 
	std::string command = "python ./src/elliptical_solid_angle.py " + std::to_string(source_x) + " " + std::to_string(source_y) + " " + std::to_string(source_z) + " " + std::to_string(PMT_x[i]) + " " + std::to_string(PMT_y[i]) + " " + std::to_string(PMT_z[i]) + " " + std::to_string(dirPMT_x[i]) + " " + std::to_string(dirPMT_y[i]) + " " + std::to_string(dirPMT_z[i]);
	//now calculate it 
	FILE* pythonProcess = popen(command.c_str(), "r");

	//read out the solid angle
	if (!pythonProcess) { //in case we failed to open it
        	perror("popen");
        	return -1;
    	}
	char buffer[128];
    	while (!feof(pythonProcess)) {
       		if (fgets(buffer, 128, pythonProcess) != NULL) {
            	// std::cout << buffer;  //reading out the output of the python code, i.e. the solid angle
                //information
                	std::string str_buffer(buffer);
            		try {
                		solidAngle = std::stod(str_buffer);
                		//std::cout << "Received double: " << solidAngle << std::endl;
            		} 
			catch (const std::invalid_argument& e) {
                		std::cerr << "Invalid double value: " << str_buffer;
            		}
            	}
    	}
	//calculate the theta to the PMT from the centre of the source to check if we are inside the 'dark cone'
	std::vector<double> sourcePMTdir = {PMT_x[i]-source_x, PMT_y[i]-source_y, PMT_z[i]-source_z};
	double norm = 0;
	double dotProduct = 0;
	for (int j=0; j<3; j++) norm += sourcePMTdir[j]*sourcePMTdir[j];
	for (int j=0; j<3; j++) sourcePMTdir[j] = sourcePMTdir[j]/TMath::Sqrt(norm);
	for (int j=0; j<3; j++) dotProduct += sourcePMTdir[j] * LBaxis[j];
	double thetaToPMT = TMath::ACos(dotProduct); 
	int multiplicativeFactor = 1.; 
	if (thetaToPMT < thetaDarkCone) multiplicativeFactor = 0;
	double expectedNbPhotonsPerPMT = densityPhotons * multiplicativeFactor * solidAngle;	
	//TODO: add the phi calculation so we can apply a corrective factor for teh soruce anisotropy
	double phiToPMT = flag; //placeholder for now

	std::cout << "Source position: [" << source_x << ", " << source_y << ", " << source_z << "], mPMT: " << mPMT[i] << " mPMT_pmt " << mPMT_pmt[i] << " PMT x y z = " <<  PMT_x[i] << " " << PMT_y[i] << " "<<PMT_z[i] << " solid angle " << solidAngle << " angle to the PMT " << thetaToPMT << std::endl;
	outfile << mPMT[i] << " " << mPMT_pmt[i] << " " <<  solidAngle << " " << expectedNbPhotonsPerPMT << " " << source_x << " " << source_y << " " << source_z << " " <<  PMT_x[i] << " " << PMT_y[i] << " "<<PMT_z[i] << " " << thetaToPMT << " " << phiToPMT << " " << n_photons_tot << std::endl;

	pclose(pythonProcess);
//			outfile << mPMT << " " << mPMT_pmt << " " <<  source_x << " " << source_y << " " << source_z << " " << theta_phi[0] << " " << theta_phi[1] << " " << R << " " << NHits << " " << TOF_Time << " " << abwff << " " << rayff << " " << bin.ID << std::endl;
    }

    //}
    //}
    outfile.close();

    //std::cout << "Total number of hits : "<< nHitsTot << std::endl;
    //std::ofstream outfile;
    //outfile.open(outfilename, std::ofstream::app);
    //outfile <<  x << " " << y << " " << z << " " << t << " " << p << " " << R << " " << nHitsTot << " " << n_entries << " " << abwff << " " << rayff << std::endl;
}//End of main


