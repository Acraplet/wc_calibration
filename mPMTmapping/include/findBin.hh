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
    	    delete [] character_ref ;
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

Bin findPMTBin(double theta_pos, double phi_pos, double hitTime=0., double R_dome=34.2, double l_PMT = 5.30, double indirectLightCut = 1.)
{
//R_dome is the radius of the mPMT dome and l_PMT is the arc length spanned on the mPMT donme by the PMT radius (projection) to defined whether the hit is incoming directly or not
    std::vector<float> xbins, ybins, zbins, thetabins, phibins;
    std::vector<int> numberbins;
    //Read the file holding the bin information so we can then measure the distance to the given source pos
    std::fstream bins_file;
    bins_file.open("./PMT-basedBins.txt");
    double x_bin, y_bin, z_bin, theta_bin, phi_bin;//, dist
    int totalNbBins = 19;
    int minID = 0; //The default value is the 20th bin, 'left over' bin
    int bin_nb;
    Bin bin; //this is the struct where we store the bin information
    if (bins_file.is_open()){   //checking whether the file is open
        std::string tp_ref;
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
                    bin_nb=double(std::stof(fs));
		    //std::cout << bin_nb << std::endl;
                }
                if (j==1){
                    std::string fs(ptr_ref);
                    x_bin=double(std::stof(fs));
                }
                if (j==2){
                    std::string fs(ptr_ref);
                    y_bin=double(std::stof(fs));
                }
                if (j==3){
                    std::string fs(ptr_ref);
                    z_bin=double(std::stof(fs));
                }
                if (j==4){
                    std::string fs(ptr_ref);
                    theta_bin=double(std::stof(fs));
                }
                if (j==5){
                    std::string fs(ptr_ref);
                    phi_bin=double(std::stof(fs));
		    //std::cout << phi_bin << std::endl;
                }
                ptr_ref = std::strtok (NULL, " ");
                j +=1;

            }//characteristic of the given bin
        xbins.push_back(x_bin);
        ybins.push_back(y_bin);
        zbins.push_back(z_bin);
        thetabins.push_back(theta_bin);
        phibins.push_back(phi_bin);
        numberbins.push_back(bin_nb);
    	delete [] character_ref;
        }//end of reading the bin positions
        //Now we need to calculate the distance - run through the vector to find the minimum
        double minBinDist = 1e10;
        double Dx, Dy, dist;
        for(long unsigned int i=0; i<=xbins.size(); i++) {
            Dx = TMath::Abs(thetabins[i] - theta_pos);
            Dy = TMath::Abs(phibins[i] - phi_pos);

            double mean_theta = (thetabins[i] + theta_pos)/2;
	    
	    //the max distance in phi is 180deg -> reverse the direction of rotation if we are larger
	    if (Dy>TMath::Pi()) Dy = 2 * TMath::Pi() - Dy;
	    
	    //the theta=0 OR pi is apex and so is on the same line as all of the thetas
            if (thetabins[i]==0. or theta_pos==0 or theta_pos == TMath::Pi() or thetabins[i]== TMath::Pi() ) Dy = 0;

            if (Dy == 0 ) dist = Dx * R_dome;
            else dist = R_dome * TMath::Sqrt(Dx*Dx + (TMath::Sin(mean_theta) * Dy) * (TMath::Sin(mean_theta) * Dy));

            if (dist<minBinDist){
                minBinDist = dist;
		//only put in the bin if we are in the radius of view
                if (minBinDist<l_PMT) {
		 	//std::cout << numberbins[i] << " " << minBinDist << std::endl;      
			minID = numberbins[i];
		}
            }
        }
        bins_file.close();
	//now we give another ID to the bins with later hits //just offset by 20 bins
	//if (hitTime>indirectLightCut) bin.ID = minID+totalNbBins+1;
	bin.ID = minID;
        bin.theta = thetabins[minID];
        bin.phi = phibins[minID];
        bin.x = xbins[minID];
        bin.y = ybins[minID];
        bin.z = zbins[minID];
    }
    return bin;

}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

std::vector<double> findThetaPhi(double source_x, double source_y, double source_z, double R=10, double mPMT = 58, double R_dome=34.2, double l_PMT = 5.30){
	//Sometimes we only have the carthesian coordinates of the source and we need to project the position on the mPMT dome. For now we are looking at the mPMT 58 (at the bottom) only but later on we will actually first check which mPMT we are in and then calculate the angle w.r.t to that specific PMT direction. 
	double mPMTDir[3] = {0, -1, 0}; //this will be changed later
	
	//we need to substract the central position of the mPMT so we can do trigonometry accurately
	double mPMTCentre[3] = {0, -155.45, 0};
	double sourcePos[3] = {source_x-mPMTCentre[0], source_y-mPMTCentre[1], source_z-mPMTCentre[2]};
	
	//we also need to normalise the position of the source back onto the dome of the mPMT
	double norm = 0;
	for (int i = 0; i < 3; i++){
		norm += sourcePos[i] * sourcePos[i];
	}
	norm = TMath::Sqrt(norm);

	//now calculate the angle
	double dotProd = 0;
        for (int i = 0; i < 3; i++){
		dotProd += mPMTDir[i] * sourcePos[i]/norm;
	};
	double theta_source = TMath::ACos(-dotProd);
	double phi_source = TMath::ACos(sourcePos[0]/(norm * TMath::Sin(theta_source)));
	if (sourcePos[2] <= 0) phi_source = 2 * TMath::Pi() - phi_source ; 

	//std::cout << "theta " << theta_source << " phi " << phi_source << " (" << phi_source * 180/TMath::Pi() << " deg )" << std::endl;
	//std::cout << std::endl;
	//Bin bin = findPMTBin(theta_source, phi_source, R_dome, l_PMT);

	std::vector<double> theta_phi = {theta_source, phi_source};
	//std::cout << "Closest bin " << bin.ID << " Theta " << theta_source << " Phi " << phi_source << std::endl;
	return theta_phi;	


}

double TOFCorrect(double Time, double R, double water_velocity = 22.4047797){
	//units used are cm/ns, using https://refractiveindex.info/?shelf=3d&book=liquids&page=water
	return Time - R/water_velocity; 
	
	
	
};

