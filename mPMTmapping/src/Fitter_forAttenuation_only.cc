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
    std::cout << "hello world" << std::endl; 
    std::cout << truth_alpha(401.9, 0.000486, 10e10) << std::endl;

    std::fstream newfile;
    newfile.open(Form("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/maps_txtFiles/mPMT_map_ID%s.txt",argv[1]),std::ios::in); //open a file to perform read operation using file object
    std::vector<double> list_R, list_A;
    char* R_test;
    if (newfile.is_open()){   //checking whether the file is open
        std::string tp;
        while(getline(newfile, tp)){  //read data from file object and put it into string.
            //this is for each test source position
//             std::cout << tp << "  " << std::endl;   //print the data of the string
            char *ptr;
            //convert to s char the string of the line we are extracting
            char* character = std::strcpy(new char[tp.length() + 1], tp.c_str());
            ptr = std::strtok(character, " "); //split the string after the blanks
            // use while loop to check ptr is not null
            int i = 0;
            char* theta_test;
            double  theta_test_num;
            char* phi_test;


            while (ptr != NULL)
            //loop over the characteristics of the given position
            {
                if (i==3){
                    std::cout << "theta : " << ptr  << " " << i << std::endl; // print the string token
                    theta_test = ptr;
                    std::string fs(ptr);
                    theta_test_num=std::stof(fs);
                }
                if (i==4){
                    std::cout << "phi : " << ptr  << " " << i << std::endl; // print the string token
                    phi_test = ptr;
                }
                if (i==5){
                    std::cout << "R : " << ptr  << " " << i << std::endl; // print the string token
                    R_test = ptr;
                }
                ptr = std::strtok (NULL, " ");
                i +=1;
            }
            //Now we have the coordinate of the position we are looking for
            //Next: open the reference txt file for this position and from there extract the fitting information we are looking for
            std::fstream position_file;
            position_file.open(Form("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/maps_onePosition/OnePosition_theta%s_phi%s_R%s.txt", theta_test, phi_test, R_test),std::ios::in);

            if (position_file.is_open()){   //checking whether the file is open
                std::string tp_ref;
                std::vector<double> data_xval, err_xval; //the variable you're fitting in
                std::vector<double> data_yval, err_yval;
                while(getline(position_file, tp_ref)){ //each reference point
//                     std::cout << tp_ref << "  " << std::endl;
                    char *ptr_ref;
                    //convert to s char the string of the line we are extracting
                    char* character_ref = std::strcpy(new char[tp_ref.length() + 1], tp_ref.c_str());
                    ptr_ref = std::strtok(character_ref, " "); //split the string after the blanks
                    int j=0;
                    float abwff, rayff;
                    double  Q_ref;
                    while (ptr_ref != NULL)
                    {//loop over the characteristics of the given position
                        if (j==6){
                            //std::cout << "Q : " << ptr_ref  << std::endl; // print the string token
                            std::string fs(ptr_ref);
                            Q_ref=std::stof(fs);
                            //std::cout << Q_ref << std::endl;
                        }
                        if (j==8){
                            //std::cout << "abwff : " << ptr_ref  << std::endl; // print the string token
                            std::string fs(ptr_ref);
                            abwff=std::stof(fs);
//                             float  = (float) (ptr_ref)
                            //std::cout << truth_alpha(401.9, abwff, 10e10);
                        }
                        if (j==9){
                            //std::cout << "rayff : " << ptr_ref  << std::endl; // print the string token
                            std::string fs(ptr_ref);
                            rayff=std::stof(fs);
                            //std::cout << truth_alpha(401.9, 10e10, rayff) << std::endl;
                        }
                        ptr_ref = std::strtok (NULL, " ");
                        j +=1;
                    }//characteristic of the given position
//                     std::cout << abwff << " and "<< rayff<< std::endl;
                    //now append to the correct data points:
                    if (abwff <= 1e10 && rayff >= 1e10){
                        double alpha_abs = truth_alpha(401.9, abwff, rayff);
                        data_xval.push_back(alpha_abs);
                        //std::cout << "Value :" << alpha_abs << " " << abwff << " "<< Q_ref << std::endl;
                        data_yval.push_back(Q_ref);
                        err_xval.push_back(0);
                        err_yval.push_back(TMath::Sqrt(Q_ref * (1 - Q_ref/1000)));
//                         std::cout << std::endl;
                    }
                }//each reference point - now set up the fitting
                const int nPars = 2;

//                 std::cout << data_yval[4] << std::endl;

                Chisq *chi = new Chisq(nPars);
                chi->setData(data_xval, data_yval);
                ROOT::Math::Functor f(chi, &Chisq::fcn, nPars);

                ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
                min->SetStrategy(3);
                min->SetFunction(f);
                min->SetMaxFunctionCalls(10000);

                min->SetVariable(0, "amplitude", 0.0, 0.01);
                min->SetVariable(1, "radius", 1.0, 0.01);

                min->Minimize();
                min->PrintResults();
                std::cout << std::endl;
                const double * res = min->X();
                const double * err = min->Errors();
                std::cout << min->Status() << std::endl;
                if (min->Status() == 0) {
                    std::cout << "haye"<< std::endl;
                    list_A.push_back(res[0]);
                    list_R.push_back(res[1]);
                }
                std::cout << res[0] << std::endl;

                TF1 func = chi->getFunction(0, 220);
                TGraphErrors *data = new TGraphErrors(data_xval.size(), &data_xval[0], &data_yval[0], &err_xval[0], &err_yval[0]);
                TGraphErrors *fit_output;
                if (min->Status() == 0){
                    fit_output = new TGraphErrors(1, &res[0], &res[1], &err[0], &err[1]);
                }
                else {
                    auto a = res[0] * 0;
                    fit_output = new TGraphErrors(1, &a, &a, &a, &a);
                }

                TFile *outf = new TFile(Form("reference_root/results_theta%s_phi%s_R%s.root", theta_test, phi_test, R_test), "RECREATE");

                func.Write();
                data->SetTitle("Charge vs attenuation lenght (absorption only)");
                data->GetXaxis()->SetTitle("Attenuation lenght (absorption only) (cm)");
                data->GetYaxis()->SetTitle("Charge collected for 1000 photons");
                data->Write("data_distribution");
                fit_output->SetTitle("Exponential fit output");
                fit_output->GetXaxis()->SetTitle("Amplitude (nb of charge collected per 1000 photon)");
                fit_output->GetYaxis()->SetTitle("Distance R to the mPMT dome (cm)");
                fit_output->Write("fit_output_xA_yR");
                outf->Close();
            }//fiunished reading the reference file for the given position
        }
    newfile.close();   //close the file object.
    //now plot the histogram -
    TH1* h_R = new TH1D("h_R", Form("Histogram of Fitted R - R True: %s", R_test), 100,0.,100.);
    TH1* h_A = new TH1D("h_A", Form("Histogram of Fitted A - R True: %s", R_test), 100,200.,400.);
    for (int u=0; u<list_R.size(); u++){
        std::cout <<  list_R[u] << std::endl;
        h_R->Fill(list_R[u]);
        h_A->Fill(list_A[u]);

    }
    h_R->GetXaxis()->SetTitle("Fitted R(cm)");
    h_A->GetXaxis()->SetTitle("Fitted max charge per 1000 photon");
    TFile f(Form("histos_R%s.root", R_test), "RECREATE");
    h_A->Write();
    h_R->Write();
    f.Close();

    }
}
