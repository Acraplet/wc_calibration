//This code is extracting the scattering length corresponding to one or more maps that have the same rayff
//and ideally abwff = inf
//For this is is looking up in reference_maps what the y values of the 5 nodes making up the charge vs scattering length fake-data spline. This number of 5 nodes can be changed but for this we need to re-run the Fitter. 
//Then the chi is minimised between the observed data and y_pred = Fake_Data_Spline(scattering length pred)
//expected scattering length. So far this is position dependant but it will have to be upgraded to a binned approach

#include "../include/truth_alpha.hh"
#include <iostream>
#include "../chisq/chisq.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include <cstring>

int main(int argc, char **argv){
    std::vector<double> list_Q, list_i;
    std::vector<std::vector<double>> list_ynodes, list_xnodes;
    char* R_test;
    char* Q_test;
    char* theta_test;
    double  theta_test_num, Q_test_num;
    char* phi_test;

    std::cout << "The total number of files to read together is: " << argc << std::endl;
    for (int f=1; f<= argc-1; f++){
        std::fstream newfile;
	//This is the test file we are going to extract the scattering length out of
        newfile.open(Form("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/Maps/maps_txtFiles/mPMT_map_ID%s.txt",argv[f]),std::ios::in); //open a file to perform read operation using file object
        std::cout << argv[f] << std::endl;
        if (newfile.is_open()){   //checking whether the file is open
            std::string tp;
            while(getline(newfile, tp)){  //read data from file object and put it into string.
                //this is for each test source position (theta, phi)
                double w = 0.;
                char *ptr;
                //convert to s char the string of the line we are extracting
                char* character = std::strcpy(new char[tp.length() + 1], tp.c_str());
                ptr = std::strtok(character, " "); //split the string after the blanks
                // use while loop to check ptr is not null
                int i = 0;

                while (ptr != NULL)
                    //loop over the characteristics of the given position
                {
                    if (i==3){
                        theta_test = ptr;
                        std::string fs(ptr);
                        theta_test_num=std::stof(fs);
                    }
                    if (i==4){
                        phi_test = ptr;
                    }
                    if (i==5){
                        R_test = ptr;
                    }
                    if (i==6){
                        Q_test = ptr;
                        std::string fs(ptr);
                        Q_test_num=std::stof(fs);
                    }
                    ptr = std::strtok (NULL, " ");
                    i +=1;
                }
                //Now we have the coordinate of the position we are looking for
                //Next: open the reference txt file for this position and from there extract the fitting 
		//information we are looking for here it is the y value corresponding to the nodes of the fake
		//data spline. Note that usually we have 5 nodes, we need to modify this value if we have another
		//number of nodes. The profile is stored in reference_root as Graph 
                std::string position_file = Form("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/reference_root/results_Abs_Scat_theta%s_phi%s_R%s.root", theta_test, phi_test, R_test);
                TFile *f = new TFile(position_file.c_str(), "READ");
                TGraphErrors *h = (TGraphErrors*)f->Get("Graph");
                std::vector<double> bufX, bufY;
                for (int count=0; count < 5; count++){
                    bufX.push_back(h->GetX()[count]);
                    bufY.push_back(h->GetY()[count]);
                }
                list_xnodes.push_back(bufX);
                list_ynodes.push_back(bufY);
                list_i.push_back(w);
                list_Q.push_back(Q_test_num);
                w+=1;
                f->Close();
            }
        } //finished reading all of the positions for a given tets map and stored the relevant info
    } //finished reading all of the test maps
    //Here the fitting begins
    const int nPars = 1; //the only parameter we fit is scattering length
    Chisq *chi = new Chisq(nPars);
    chi->setData(list_i, list_Q);
    chi->setRef_scat(list_xnodes, list_ynodes);
    ROOT::Math::Functor functor(chi, &Chisq::fitter_rayleigh, nPars);
    ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    min->SetStrategy(3);
    min->SetFunction(functor);
    min->SetMaxFunctionCalls(10000);
    min->SetVariable(0, "scattering_length", 5.0, 0.01);
    min->Minimize();
    min->PrintResults();
}
