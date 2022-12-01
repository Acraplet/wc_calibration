#include "truth_alpha.hh"
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

    char* R_test;
    char* Q_test;
    char* theta_test;
    char* phi_test;
    double  theta_test_num, Q_test_num, R_test_num;
    std::vector<double> list_R;
    std::vector<double> list_A;
    std::vector<double> list_Q;
    std::vector<double> list_i;

    std::cout << "Thew total number of files to read together is: " << argc << std::endl;
    for (int f=1; f<= argc-1; f++){
        std::fstream newfile;
        newfile.open(Form("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/maps_txtFiles/mPMT_map_ID%s.txt",argv[f]),std::ios::in); //open a file to perform read operation using file object
        std::cout << argv[f] << std::endl;
        if (newfile.is_open()){   //checking whether the file is open
            std::string tp;
            while(getline(newfile, tp)){  //read data from file object and put it into string.
                //this is for each test source position
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
                        //std::cout << "theta : " << ptr  << " " << i << std::endl; // print the string token
                        theta_test = ptr;
                        std::string fs(ptr);
                        theta_test_num=std::stof(fs);
                    }
                    if (i==4){
                        //std::cout << "phi : " << ptr  << " " << i << std::endl; // print the string token
                        phi_test = ptr;
                    }
                    if (i==5){
                        //std::cout << "R : " << ptr  << " " << i << std::endl; // print the string token
                        R_test = ptr;
                        std::string fs(ptr);
                        R_test_num = std::stof(fs);
                        list_R.push_back(std::stof(fs));
                    }
                    if (i==6){
                        //std::cout << "Q : " << ptr  << " " << i << std::endl; // print the string token
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
                std::fstream position;
                position.open(Form("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/maps_Reference/OnePosition_theta%s_phi%s.txt", theta_test, phi_test),std::ios::in);
                if (position.is_open()){   //checking whether the file is open
                    std::string tp2;
//                     std::cout << "HUUU" << std::endl;
                    while(getline(position, tp2)){  //read data from file object and put it into string.
                        //this is for each test source position
                        char *ptr2;
                        //convert to s char the string of the line we are extracting
                        char* character2 = std::strcpy(new char[tp2.length() + 1], tp2.c_str());
                        ptr2 = std::strtok(character2, " "); //split the string after the blanks
                        // use while loop to check ptr is not null
                        int i2 = 0;
                        while (ptr2 != NULL)
                            //loop over the characteristics of the given position
                        {
                            if (i2==0){

                                std::string fs(ptr2);
                                list_A.push_back(std::stof(fs));
                                std::cout << "A: "<< std::stof(fs) << std::endl;
                            }
                            ptr2 = std::strtok (NULL, " ");
                            i2 +=1;
                        }
                    }
                }

                //std::cout << *h->GetY() << std::endl;
//                 list_A.push_back(*h->GetX());
//                 list_R.push_back(*h->GetY()); // we are not using the list of reference R anymore - we are using the true R
                list_i.push_back(w);
                list_Q.push_back(Q_test_num);
                std::cout << "Q: " << Q_test_num * 1/TMath::Exp(-R_test_num/60) << std::endl;
                std::cout << std::endl;
                //std::cout << std::endl;
                w+=1;
                //Here - extract the function (maybe the parameters directly?)
                //Maybe can save the coefficients for each source position for a given value?
                //also the positions aren't in the same order? oh but they will be, that's okay

            }
        } //finished reading all of the positions, now need to fit_output_xA_y
    } //finished reading all of the count positions
    //Here the fitting begins

//     std::cout << list_i.size() << " " << list_A.size() << std::endl;
    const int nPars = 1; //the only parameter we fit is abwff
    Chisq *chi = new Chisq(nPars);
    chi->setData(list_i, list_Q);
    chi->setRef(list_A, list_R);
    ROOT::Math::Functor functor(chi, &Chisq::fcn_abwff, nPars);
    ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    min->SetStrategy(3);
    min->SetFunction(functor);
    min->SetMaxFunctionCalls(10000);
    min->SetVariable(0, "attenuation_length", 5.0, 0.01);
    min->Minimize();
    min->PrintResults();
}



/* This is the old way of doing it where both the R and A were saved in a specific external .root file
 *
 * std::string position_file = Form("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/reference_root/results_theta%s_phi%s_R%s.root", theta_test, phi_test, R_test);
 *
 T File *f = new TFile(position*_file.c_str(), "READ");
 TGraphErrors *h = (TGraphErrors*)f->Get("fit_output_xA_yR");
 //std::cout << *h->GetY() << std::endl;
 list_A.push_back(*h->GetX());
 f->Close();
 */
