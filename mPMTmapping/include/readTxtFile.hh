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
    float x, y, z, theta, phi, R, abwff, rayff, Q, nEvents;
};

typedef struct fileData Data;

std::vector<Data> readTxtFile(char* filename){
    std::fstream position_file;
    //this is the the vector of 'line' Data where one line is one simulated source position (or content of a bin)
    std::vector<Data> list_all;
    position_file.open(filename, std::ios::in);
    if (position_file.is_open()){   //checking whether the file is open
        std::string tp_ref;
        while(getline(position_file, tp_ref)){ //Each reference point
            char *ptr_ref;
            //convert to s char the string of the line we are extracting
            char* character_ref = std::strcpy(new char[tp_ref.length() + 1], tp_ref.c_str());
            ptr_ref = std::strtok(character_ref, " "); //split the string after the blanks
            int j=0;
            float x, y, z, t, p, R, abwff, rayff;
            float Q_ref, nevents;

            while (ptr_ref != NULL)
            {//loop over the characteristics of the given position
                if (j==0){
                    std::string fs(ptr_ref);
                    x=std::stof(fs);
                }
                if (j==1){
                    std::string fs(ptr_ref);
                    y=std::stof(fs);
                }
                if (j==2){
                    std::string fs(ptr_ref);
                    z=std::stof(fs);
                }
                if (j==3){
                    std::string fs(ptr_ref);
                    t=std::stof(fs);
                }
                if (j==4){
                    std::string fs(ptr_ref);
                    p=std::stof(fs);
                }
                if (j==5){
                    std::string fs(ptr_ref);
                    R=std::stof(fs);
                }
                if (j==6){
                    std::string fs(ptr_ref);
                    Q_ref=std::stof(fs);
                }
                if (j==7){
                    std::string fs(ptr_ref);
                    nevents=std::stof(fs);
                }
                if (j==8){
                    std::string fs(ptr_ref);
                    abwff=std::stof(fs);
                }
                if (j==9){
                    std::string fs(ptr_ref);
                    rayff=std::stof(fs);
                }
                ptr_ref = std::strtok (NULL, " ");
                j +=1;
            }//characteristic of the given line
            Data pos;
            pos.x = x;
            pos.y = y;
            pos.z = z;
            pos.theta = t;
            pos.phi = p;
            pos.R = R;
            pos.Q = Q_ref;
            pos.nEvents = nevents;
            pos.abwff = abwff;
            pos.rayff = rayff;
            list_all.push_back(pos);
        }//now have read all the lines
        position_file.close();
    }//now we finished reading the file
    return list_all;
}//now we finished the function
