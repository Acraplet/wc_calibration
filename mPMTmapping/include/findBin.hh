//This is a small code to return the bin associated with a given theta and phi of a source position
//we return a struct wich has all of the coordinates of the cloest bin
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
#include <cstring>

struct binInfo {
    int ID;
    double x, y, z, theta, phi;
};

typedef struct binInfo Bin;

Bin findBin(double theta_pos, double phi_pos)
{
    std::vector<double> xbins, ybins, zbins, thetabins, phibins;
    //Read the file holding the bin information so we can then measure the distance to the given source pos
    std::fstream bins_file;
    bins_file.open("./uniform_top_bins_withBinNumber.txt");
    double x_bin, y_bin, z_bin, theta_bin, phi_bin;//, dist;
    int minID = -1; //id of the closest bin
    Bin bin; //this is the struct where we store the bin information
    if (bins_file.is_open()){   //checking whether the file is open
        std::string tp_ref;
        //the variable you're fitting in
        while(getline(bins_file, tp_ref)){ //each reference point
            char *ptr_ref;
            //convert to char the string of the line we are extracting
            char* character_ref = std::strcpy(new char[tp_ref.length() + 1], tp_ref.c_str());
            ptr_ref = std::strtok(character_ref, " "); //split the string after the blanks
            int j=0;
            while (ptr_ref != NULL)
            {//loop over the characteristics of the given position
                if (j==0){
                    std::string fs(ptr_ref);
                    x_bin=double(std::stof(fs));
                }
                if (j==1){
                    std::string fs(ptr_ref);
                    y_bin=double(std::stof(fs));
                }
                if (j==2){
                    std::string fs(ptr_ref);
                    z_bin=double(std::stof(fs));
                }
                if (j==3){
                    std::string fs(ptr_ref);
                    theta_bin=double(std::stof(fs));
                }
                if (j==4){
                    std::string fs(ptr_ref);
                    phi_bin=double(std::stof(fs));
                }
                ptr_ref = std::strtok (NULL, " ");
                j +=1;

            }//characteristic of the given bin
            xbins.push_back(x_bin);
            ybins.push_back(y_bin);
            zbins.push_back(z_bin);
            thetabins.push_back(theta_bin);
            phibins.push_back(phi_bin);
        }//end of reading the bin positions
        //Now we need to calculate the distance - run through the vector to find the minimum
        double minBinDist = 1e10;
        double Dx, Dy, dist;
        for(long unsigned int i=0; i<=xbins.size(); i++) {
            Dx = thetabins[i] - theta_pos;
            Dy = phibins[i] - phi_pos;
            //limit computing time: look at the square of the angular distance as our reference measure
            dist = Dx * Dx + Dy * Dy;
            if (dist<minBinDist){
                minBinDist = dist;
                minID = i;
            }
        }
        bins_file.close();
        bin.ID = minID;
        bin.theta = thetabins[minID];
        bin.phi = phibins[minID];
        bin.x = xbins[minID];
        bin.y = ybins[minID];
        bin.z = zbins[minID];
    }
    return bin;
}
