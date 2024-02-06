//This is a code that given any text file returns a set of vectors will all of the informations
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
struct fileData {
    float x, y, z, theta, phi, R, abwff, rayff, Q, time, bin, nPhotons, mPMT, mPMT_pmt;
};

typedef struct fileData DataWithTime;

std::vector<DataWithTime> readTxtFileWithTime(char* filename){
    std::fstream position_file;
    std::cout << "Trying " << std::endl;
    //this is the the vector of 'line' Data where one line is one simulated source position (or content of a bin)
    std::vector<DataWithTime> list_all;
    position_file.open(filename, std::ios::in);
    if (position_file.is_open()){   //checking whether the file is open
        std::string tp_ref;
        while(getline(position_file, tp_ref)){ //Each reference point
	    //std::cout << "ok" << std::endl;
            char *ptr_ref;
            //convert to s char the string of the line we are extracting
            char* character_ref = std::strcpy(new char[tp_ref.length() + 1], tp_ref.c_str());
            ptr_ref = std::strtok(character_ref, " "); //split the string after the blanks
            int j=0;
            float x, y, z, t, p, R, abwff, rayff, time, bin, mPMT, mPMT_pmt;
            float Q_ref, nPhotons;

            while (ptr_ref != NULL)
            {//loop over the characteristics of the given position
                if (j==0){
                    std::string fs(ptr_ref);
                    mPMT=std::stof(fs);
                }
                if (j==1){
                    std::string fs(ptr_ref);
                    mPMT_pmt=std::stof(fs);
                }
                if (j==2){
                    std::string fs(ptr_ref);
                    x=std::stof(fs);
                }
                if (j==3){
                    std::string fs(ptr_ref);
                    y=std::stof(fs);
                }
                if (j==4){
                    std::string fs(ptr_ref);
                    z=std::stof(fs);
                }
                if (j==5){
                    std::string fs(ptr_ref);
                    t=std::stof(fs);
                }
                if (j==6){
                    std::string fs(ptr_ref);
                    p=std::stof(fs);
                }
                if (j==7){
                    std::string fs(ptr_ref);
                    R=std::stof(fs);
                }
		//if (j==8){
                //    std::string fs(ptr_ref);
                //    nPhotons=std::stof(fs);
                //}

                if (j==8){
                    std::string fs(ptr_ref);
                    Q_ref=std::stof(fs);
                }
                if (j==9){
                    std::string fs(ptr_ref);
                    time=std::stof(fs);
                }
                if (j==10){
                    std::string fs(ptr_ref);
                    abwff=std::stof(fs);
                }
                if (j==11){
                    std::string fs(ptr_ref);
                    rayff=std::stof(fs);
                }
                if (j==12){
                    std::string fs(ptr_ref);
                    bin=std::stof(fs);
                }
                ptr_ref = std::strtok (NULL, " ");
                j +=1;
            }//characteristic of the given line
            DataWithTime pos;
            pos.x = x;
            pos.y = y;
            pos.z = z;
            pos.time = time;
            pos.theta = t;
            pos.phi = p;
            pos.R = R;
	    pos.nPhotons = nPhotons;
            pos.Q = Q_ref;
            pos.bin = bin;
            pos.mPMT = mPMT;
            pos.mPMT_pmt = mPMT_pmt;
            pos.abwff = abwff;
            pos.rayff = rayff;
            list_all.push_back(pos);
	    std::cout << "The data that we are reading in has " << Q_ref << " PE collected in the PMT " << mPMT_pmt << " of mPMT " << mPMT << " the source - PMT distance is " << R << " and the TOF corrected time for this hit is " << time << " ns." << std::endl;
        }//now have read all the lines
        position_file.close();
    }//now we finished reading the file
    return list_all;
}//now we finished the function


std::vector<DataWithTime> readLBpredNbPhotons(char* filename){
    std::fstream position_file;
    //this is the the vector of 'line' Data where one line is one simulated source position (or content of a bin)
    std::vector<DataWithTime> list_all;
    position_file.open(filename, std::ios::in);
    if (position_file.is_open()){   //checking whether the file is open
        std::string tp_ref;
        while(getline(position_file, tp_ref)){ //Each mPMT-PMT pair
            char *ptr_ref;
            //convert to s char the string of the line we are extracting
            char* character_ref = std::strcpy(new char[tp_ref.length() + 1], tp_ref.c_str());
            ptr_ref = std::strtok(character_ref, " "); //split the string after the blanks
            int j=0;
            float x_PMT, y_PMT, z_PMT, t, p, x_source, y_source, z_source, mPMT, mPMT_pmt;
            float nPhotons, solid_angle, pred_nb_photons;

            while (ptr_ref != NULL)
            {//loop over the characteristics of the given position
                if (j==0){
                    std::string fs(ptr_ref);
                    mPMT=std::stof(fs);
                }
                if (j==1){
                    std::string fs(ptr_ref);
                    mPMT_pmt=std::stof(fs);
                }
		if (j==2){
                    std::string fs(ptr_ref);
                    solid_angle=std::stof(fs);
                }
		if (j==3){
                    std::string fs(ptr_ref);
                    pred_nb_photons=std::stof(fs);
                }
                if (j==4){
                    std::string fs(ptr_ref);
                    x_source=std::stof(fs);
                }
                if (j==5){
                    std::string fs(ptr_ref);
                    y_source=std::stof(fs);
                }
                if (j==6){
                    std::string fs(ptr_ref);
                    z_source=std::stof(fs);
                }
                if (j==7){
                    std::string fs(ptr_ref);
                    x_PMT=std::stof(fs);
                }
                if (j==8){
                    std::string fs(ptr_ref);
                    y_PMT=std::stof(fs);
                }
                if (j==9){
                    std::string fs(ptr_ref);
                    z_PMT=std::stof(fs);
                }
                if (j==10){
                    std::string fs(ptr_ref);
                    t=std::stof(fs);
                }
                if (j==11){
                    std::string fs(ptr_ref);
                    p=std::stof(fs);
                }
		if (j==12){
                    std::string fs(ptr_ref);
                    nPhotons=std::stof(fs);
                }
                ptr_ref = std::strtok (NULL, " ");
                j +=1;
            }//characteristic of the given line
            DataWithTime pos;
            pos.x = x_PMT;
            pos.y = y_PMT;
            pos.z = z_PMT;
            pos.time = -9999.;
            pos.theta = t;
            pos.phi = p;
            pos.R = TMath::Sqrt((x_PMT-x_source) * (x_PMT-x_source) + (y_PMT-y_source) * (y_PMT-y_source) + (z_PMT-z_source) * (z_PMT-z_source));
	    pos.nPhotons = nPhotons;
            pos.Q = pred_nb_photons;
            pos.bin = -9999.;
            pos.mPMT = mPMT;
            pos.mPMT_pmt = mPMT_pmt;
            pos.abwff = -9999.;
            pos.rayff = -9999.;
            list_all.push_back(pos);
	    if (mPMT >= 50 and mPMT <=60){ 
	    	std::cout << "The number of photons: " << pos.Q << " photon reaches the PMT " << mPMT_pmt << " of mPMT " << mPMT << " the source - PMT distance is " << pos.R << " and the total number of photons sent was " << nPhotons << std::endl;
	    }
        }//now have read all the lines
        position_file.close();
    }//now we finished reading the file
    return list_all;
}//now we finished the function

