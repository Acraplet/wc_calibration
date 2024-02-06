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
#include <stdexcept>
#include <TRandom3.h>
#include "findBin.hh"
#include <stdexcept>

double FresnelTransmission(double theta, double n1 = 1.333, double n2 = 1.50){
	//Function to calculate the fraction of light that is transmitted through the 
	//theta should be given in radiants
	//acrylic surface depending on the incoming angle, n1 is ref index of water, n2 of acrylic
	//we approximate that the position of the PMT centre is the angle for all of the PMT
	//
		
	double sqrtFractor = TMath::Sqrt(1 - ((n1/n2 * TMath::Sin(theta)) * (n1/n2 * TMath::Sin(theta))) )  ;
	
	double Rs_num = n1 * TMath::Cos(theta) - n2 * sqrtFractor; 
	double Rs_denom = n1 * TMath::Cos(theta) + n2 * sqrtFractor; 
	
	double FresnelRs = TMath::Abs( Rs_num / Rs_denom) * TMath::Abs( Rs_num / Rs_denom);

	double Rp_num = n1 * sqrtFractor - n2 * TMath::Cos(theta);
        double Rp_denom = n1 * sqrtFractor + n2 * TMath::Cos(theta);
	std::cout << "sqrtFractor " << sqrtFractor << std::endl;
	std::cout << "Rp_num " << Rp_num << std::endl;
	std::cout << "Rp_denom " << Rp_denom << std::endl;

        double FresnelRp = TMath::Abs( Rp_num / Rp_denom) * TMath::Abs( Rp_num / Rp_denom);
	std::cout << "FresnelRp " << FresnelRp << std::endl;
	std::cout << "FresnelRs " << FresnelRs << std::endl;

	double meanReflection = 0.5 * (FresnelRs + FresnelRp);
	std::cout << " meanReflection " << meanReflection << std::endl;
	return 1 - meanReflection; //this is the total power that is transmitted	
}	


std::vector<double> GetSolidAngle(double source_x, double source_y, double source_z, double surfaceX, double surfaceY, double surfaceZ, double dirPMT_x, double dirPMT_y, double dirPMT_z, double effectivePMTradius, double defaultMultiplicativeFactor = 1., std::vector<double> LBaxis = {0, 1, 0}, double thetaDarkCone = 40.){
	std::string command = "python ./src/elliptical_solid_angle.py " + std::to_string(source_x) + " " + std::to_string(source_y) + " " + std::to_string(source_z) + " " + std::to_string(surfaceX) + " " + std::to_string(surfaceY) + " " + std::to_string(surfaceZ) + " " + std::to_string(dirPMT_x) + " " + std::to_string(dirPMT_y) + " " + std::to_string(dirPMT_z) + " " + std::to_string(effectivePMTradius);
	double solidAngle;
        std::vector<double> results;
       	FILE* pythonProcess = popen(command.c_str(), "r");
        //read out the solid angle
        if (!pythonProcess) { //in case we failed to open it
                perror("popen");
                return results;
        }
        char buffer[128];
        while (!feof(pythonProcess)) {
        if (fgets(buffer, 128, pythonProcess) != NULL) {
            std::string str_buffer(buffer);
             try {
                 solidAngle = std::stod(str_buffer);
             }
             catch (const std::invalid_argument& e) {
             	std::cerr << "Invalid double value: " << str_buffer;
             }
        }
        }
        //calculate the theta to the PMT from the centre of the source to check if we are inside the 'dark cone'
        std::vector<double> sourcePMTdir = {surfaceX-source_x, surfaceY-source_y, surfaceZ-source_z};
        std::vector<double> PMTdir = {dirPMT_x, dirPMT_y, dirPMT_z};
        double norm = 0;
        double dotProduct = 0;
        double fromPMTdotProduct = 0;
        for (int j=0; j<3; j++) norm += sourcePMTdir[j]*sourcePMTdir[j];
        for (int j=0; j<3; j++) sourcePMTdir[j] = sourcePMTdir[j]/TMath::Sqrt(norm);
        for (int j=0; j<3; j++) dotProduct += sourcePMTdir[j] * LBaxis[j];
        double thetaToPMT = TMath::ACos(dotProduct);
	for (int j=0; j<3; j++) fromPMTdotProduct += sourcePMTdir[j] * PMTdir[j];
	double thetaFromPMT = TMath::Pi() - TMath::ACos(fromPMTdotProduct);
	//this is the anlge w,r,t the PMT normal 
	//(to make the Fresnel reflection calculation)
        double multiplicatorFactor = defaultMultiplicativeFactor;

        if (thetaToPMT < thetaDarkCone){
		 multiplicatorFactor = 0.;
	}
	else{
		if (defaultMultiplicativeFactor >= 1e-4) {
			//Need to calculate the transmittance (Fresnel), which depends
			// on the incoming angle, only if we did expect some light (reflectors, to be cleaned up)
			multiplicatorFactor = FresnelTransmission(thetaFromPMT);
			std::cout << "Yay" << std::endl;
		}
	}

	if (defaultMultiplicativeFactor >= 1e-4) {	
		std::cout << "The angle between the source and the PMT is " << thetaToPMT << " radians (" << thetaToPMT * 180 / TMath::Pi() << " degrees) and the source has an angular half opening of  "<< thetaDarkCone  << " radians. The photon incoming angle w.r.t to the PMT surface is " << thetaFromPMT * 180/TMath::Pi()<< " degrees. The multiplicative factor applied is  " <<  multiplicatorFactor << ", the default multiplicative factor was:  " << defaultMultiplicativeFactor << "\n" << std::endl;
	}
	
	results.push_back(solidAngle);
	results.push_back(multiplicatorFactor);
	results.push_back(thetaToPMT);
	results.push_back(thetaFromPMT);
	pclose(pythonProcess);
	return results;
}

void HelpMessage()
{
    std::cout << "USAGE: "
    << "make_mPMTmap" << "\nOPTIONS:\n"
    << "-f : Input file\n" << std::endl
    << "-o : Output file\n" << std::endl
    << "-r : effective PMT radius\n" << std::endl
    << "-c : LB cone angle\n" << std::endl;
}

int main(int argc, char **argv){
    char * filename=NULL;
    char * outfilename=NULL;
    char c;
    double effectivePMTradius = 4.75; //cm 
    double coneSize = NULL; //degrees 
    double trueRadius = 4.1; //cm
    double multiplicatorReflector = 0.5;
    double multiplicatorPMTsurface = 1;

    while( (c = getopt(argc,argv,"f:o:b:s:e:l:n:m:r:p:c:w:hdtv")) != -1 ){//input in c the argument (-f etc...) 
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

	switch(c){
            case 'e':
                effectivePMTradius = std::stod(optarg);
                break;
	}
	switch(c){
            case 'r':
                trueRadius = std::stod(optarg);
                break;
        }
	switch(c){
            case 'm':
                multiplicatorPMTsurface = std::stod(optarg);
                break;
        }
	switch(c){
            case 'n':
                multiplicatorReflector = std::stod(optarg);
                break;
        }
	
	switch(c){
            case 'c':
                coneSize = std::stod(optarg);
                break;
        }

    }

    if (filename==NULL){
        std::cout << "Error, no input file: " << std::endl;
        HelpMessage();
	throw std::invalid_argument( "Missing the input file");
    }

    if (outfilename==NULL){
    	outfilename = Form("./Maps/predictionNumberOfPhotons.txt");
        std::cout << "Warning, no output file given, using " << outfilename << std::endl;
//         HelpMessage();
    }

    if (coneSize==NULL){
        std::cout << "Warning, no coneSize file given, please use the -c argument" << std::endl;
        HelpMessage();
	throw std::invalid_argument( "Missing the cone size");
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

    std::cout << "Starting calculation of solid angle, for a source located at " << x << " " << y << " " << z <<  " with an opening angle of " << coneSize << " degrees and a PMT effective radius of " << effectivePMTradius << " cm." << " The multiplicator factors are " << multiplicatorPMTsurface << " for the surface and " << multiplicatorReflector <<  " for the reflector."<< std::endl; 
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


    double thetaDarkCone = (180-coneSize) * TMath::Pi() / 180; 
    std::vector<double> LBaxis = {0, 1, 0};

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
    double PMTcentreToSurfaceOffset = 4.; 
    for(int i=0; i<nPMT; i++){ //change back for 0 and nPMT
	double surfaceX = PMT_x[i] + PMTcentreToSurfaceOffset * dirPMT_x[i];
	double surfaceY = PMT_y[i] + PMTcentreToSurfaceOffset * dirPMT_y[i];
	double surfaceZ = PMT_z[i] + PMTcentreToSurfaceOffset * dirPMT_z[i];

     	std::vector<double> surface_results  = GetSolidAngle(source_x, source_y, source_z, surfaceX, surfaceY, surfaceZ, dirPMT_x[i], dirPMT_y[i],  dirPMT_z[i], trueRadius, multiplicatorPMTsurface, LBaxis , thetaDarkCone);

	double solidAngle = surface_results[0];
	double thetaToPMT = surface_results[2];
	double thetaFromPMT = surface_results[3];
	double outputMultiplicatorPMTsurface = surface_results[1]; 
     	
	std::vector<double> reflector_results  = GetSolidAngle(source_x, source_y, source_z, surfaceX, surfaceY, surfaceZ, dirPMT_x[i], dirPMT_y[i],  dirPMT_z[i], effectivePMTradius, multiplicatorReflector, LBaxis , thetaDarkCone);
	
	double solidAngleReflector = reflector_results[0] - solidAngle; 
	double phiToPMT = flag;
	double outputMultiplicatorReflector = reflector_results[1];

	double expectedNbPhotonsPerPMTsolidAngle = (solidAngle * outputMultiplicatorPMTsurface + solidAngleReflector * outputMultiplicatorReflector) * densityPhotons;

	//Here, calculate the angular response correction that is necessary
	//WCTE factors, obtained from the simlation, sorry for the hardcoded factors
	//TODO: improve this 
	double angularResponseFactor = 1.0;
	double theta = thetaFromPMT* 180./TMath::Pi();
	if (theta < 40 and theta >= 1){
		std::cout << "yay" << std::endl;
		std::cout << "yay, angle is "  << theta << " degrees." << std::endl;
		angularResponseFactor = theta * 7.67e-3 - 3.93e-2;
		angularResponseFactor = 1/(1+angularResponseFactor);
		std::cout << "The AR factor is " << angularResponseFactor << std::endl;
	}
	else if (theta >= 40){
		std::cout << "yay, angle is "  << theta << " degrees." << std::endl;
		angularResponseFactor = theta * 15.67e-3 - 14.93e-2;
		angularResponseFactor = 1/(1+angularResponseFactor);
		std::cout << "The AR factor is " << angularResponseFactor << std::endl;
	} 
	double expectedNbPhotonsPerPMT = expectedNbPhotonsPerPMTsolidAngle * angularResponseFactor;

	outfile << mPMT[i] << " " << mPMT_pmt[i] << " " <<  solidAngle << " " << expectedNbPhotonsPerPMT << " " << source_x << " " << source_y << " " << source_z << " " <<  surfaceX << " " << surfaceY << " "<<surfaceZ << " " << thetaToPMT << " " << phiToPMT << " " << n_photons_tot << std::endl;
    	std::cout << "Source position: [" << source_x << ", " << source_y << ", " << source_z << "], mPMT: " << mPMT[i] << " mPMT_pmt " << mPMT_pmt[i] << " PMT x y z = " <<  PMT_x[i] << " " << PMT_y[i] << " "<<PMT_z[i] << " solid angle " << solidAngle << " angle to the PMT " << thetaToPMT << " multiplicatorSurface = " << outputMultiplicatorPMTsurface << " total expected nB photons = " << expectedNbPhotonsPerPMT << ". \n \n " << std::endl;

	}

    //}
    //}
    outfile.close();

    //std::cout << "Total number of hits : "<< nHitsTot << std::endl;
    //std::ofstream outfile;
    //outfile.open(outfilename, std::ofstream::app);
    //outfile <<  x << " " << y << " " << z << " " << t << " " << p << " " << R << " " << nHitsTot << " " << n_entries << " " << abwff << " " << rayff << std::endl;
}//End of main


