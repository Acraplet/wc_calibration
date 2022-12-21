
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



int main(int argc, char **argv){
    //we need the true value of alpha
    //std::cout << "hello world" << std::endl;
    //std::cout << truth_alpha(401.9, 0.000486, 10e10) << std::endl;

    std::vector<double> list_Q, list_i;
    std::vector<std::vector<double>> list_ynodes, list_xnodes;
    char* R_test;
    char* Q_test;
    char* theta_test;
    double  theta_test_num, Q_test_num;
    char* phi_test;


    std::cout << "Thew total number of files to read together is: " << argc << std::endl;
    for (int f=1; f<= argc-1; f++){
        std::fstream newfile;
        newfile.open(Form("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/maps_txtFiles/mPMT_map_ID%s.txt",argv[f]),std::ios::in); //open a file to perform read operation using file object
        std::cout << argv[f] << std::endl;
        if (newfile.is_open()){   //checking whether the file is open
            std::string tp;
            while(getline(newfile, tp)){  //read data from file object and put it into string.
                //this is for each test source position
                //             std::cout << tp << "  " << std::endl;   //print the data of the string
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
//                         std::cout << "theta : " << ptr  << " " << i << std::endl; // print the string token
                        theta_test = ptr;
                        std::string fs(ptr);
                        theta_test_num=std::stof(fs);
                    }
                    if (i==4){
//                         std::cout << "phi : " << ptr  << " " << i << std::endl; // print the string token
                        phi_test = ptr;
                    }
                    if (i==5){
//                         std::cout << "R : " << ptr  << " " << i << std::endl; // print the string token
                        R_test = ptr;
                    }
                    if (i==6){
//                         std::cout << "Q : " << ptr  << " " << i << std::endl; // print the string token
                        Q_test = ptr;
                        std::string fs(ptr);
                        Q_test_num=std::stof(fs);
                    }
                    ptr = std::strtok (NULL, " ");
                    i +=1;
                }
                //Now we have the coordinate of the position we are looking for
                //Next: open the reference txt file for this position and from there extract the fitting information we are looking for
                //             std::fstream position_file;
                std::string position_file = Form("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/reference_root/results_Abs_Scat_theta%s_phi%s_R%s.root", theta_test, phi_test, R_test);
//                 std::cout<<  Form("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/reference_root/results_Abs_Scat_theta%s_phi%s_R%s.root", theta_test, phi_test, R_test) << std::endl;
                TFile *f = new TFile(position_file.c_str(), "READ");
                TGraphErrors *h = (TGraphErrors*)f->Get("Graph");
//                 std::cout << *h->GetY() << std::endl;
                std::vector<double> bufX, bufY;
                for (int count=0; count < 5; count++){
                    bufX.push_back(h->GetX()[count]);
                    bufY.push_back(h->GetY()[count]);
//                     std::cout << h->GetX()[count] << std::endl;
                }
//                 std::cout<< std::endl;
                list_xnodes.push_back(bufX);
                list_ynodes.push_back(bufY);
                list_i.push_back(w);
                list_Q.push_back(Q_test_num);
                //std::cout << std::endl;
                w+=1;
                //Here - extract the function (maybe the parameters directly?)
                //Maybe can save the coefficients for each source position for a given value?
                //Then can circle through the histogram - that's the best way I think - but then we need to save everything?
                //also the positions aren't in the same order? oh but they will be, that's okay
                f->Close();
            }
        } //finished reading all of the positions, now need to fit
    } //finished readin all of the test maps
    //Here the fitting begins
    const int nPars = 1; //the only parameter we fit is abwff
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
    //std::cout << min->Errors()[0] << std::endl;




}
