#include "../include/truth_alpha.hh"
#include "../include/findBin.hh"
#include "../include/readTxtFile.hh"

#include <iostream>
#include "../chisq/chisq.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TFile.h"
// #include "TRandom3.h"
#include "TSpline.h"
#include "TGraph.h"
#include "TH1.h"
#include "TGraphErrors.h"

int main(int argc, char **argv){
    //we need the true value of alpha
    std::fstream newfile;
    //This is the configuration of which we are going to fit all of the source positions to get the reference
    //from the OnePositon files
    newfile.open(Form("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/Maps/maps_txtFiles/mPMT_map_ID%s.txt",argv[1]),std::ios::in); //open a file to perform read operation using file object
    std::vector<double> list_R, list_A;
    char* R_test;
    if (newfile.is_open()){   //checking whether the file is open
        std::string tp;
        while(getline(newfile, tp)){  //read data from file object and put it into string.
            //this is for each test source position
            char *ptr;
            //convert to s char the string of the line we are extracting
            char* character = std::strcpy(new char[tp.length() + 1], tp.c_str());
            ptr = std::strtok(character, " "); //split the string after the blanks
            // use while loop to check ptr is not null
            int i = 0;
            char* theta_test;
            double  theta_test_num, phi_test_num;
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
                    std::string fs(ptr);
                    phi_test_num=std::stof(fs);

                }
                if (i==5){
                    std::cout << "R : " << ptr  << " " << i << std::endl; // print the string token
                    R_test = ptr;
                }
                ptr = std::strtok (NULL, " ");
                i +=1;
            }
            //Now we have the coordinate of the position we are looking for
            //Next: open the reference txt file for this position and from 
            //there extract the fitting information we are looking for

            //Here the reference is a given bin - the closest one:
            Bin closestBin = findBin(theta_test_num, phi_test_num);
//             std::cout << "Bin " << closestBin.ID << std::endl;
            int minID = closestBin.ID;
            double theta_bin = closestBin.theta;
            double phi_bin = closestBin.phi;

            char * referenceFile = (char *)Form("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/Maps/maps_onePosition/OnePosition_theta%s_phi%s_R%s.txt", theta_test, phi_test, R_test);

//             std::cout << "Bin " << closestBin.ID << std::endl;

            //now we have all the reference positions from the correct bin
            std::vector<Data> referenceData = readTxtFile(referenceFile);

//             std::cout << "referenceData " << referenceData[0].abwff << std::endl;

            //We need to sum up the positions of the given bin at the given theta, phi combo for the final ref file
            std::vector<float> data_Q, data_nEvents, data_rayff, data_abwff;
            for (int pos = 0; pos<= referenceData.size(); pos++){
//                 std::cout << "Pos " << pos << std::endl;

                Data line = referenceData[pos];
                //need to circle over the reference that we already have to merge the positions instead of making new data points
                int check = 0;
                bool added = false;

                while (check < data_rayff.size() and added == false){
//                     std::cout << referenceData[pos].Q << " check " << check << std::endl;
                    if (data_rayff[check] == line.rayff and data_abwff[check] == line.abwff){
                        data_Q[check] += line.Q;
                        data_nEvents[check] += line.nEvents;
                        added == true;
                    }
                    check += 1;
                }
                //If it is the first time we encounter this rayff, abwff combo:
                if (added == false){
//                     std::cout << "New entry : rayff = " << line.rayff << " abwff = " << line.abwff << std::endl;
                    data_Q.push_back(line.Q);
                    data_nEvents.push_back(line.nEvents);
                    data_rayff.push_back(line.rayff);
                    data_abwff.push_back(line.abwff);
                }
            }//read a new entry of the reference dataset

//             std::cout << "Bin uum" << closestBin.ID << std::endl;

            std::vector<double> data_xval, err_xval, data_scat_xval, err_scat_xval;
            std::vector<double> data_yval, err_yval, data_scat_yval, err_scat_yval;
            //TO DO: update this to be the actual value of the test number of photons sent
            float  nEvents_test = 1000.0;

            //And now that the dataset is merged we can use it as a reference
            for (int i=0; i<=data_abwff.size(); i++){
                float abwff = data_abwff[i];
                float rayff = data_rayff[i];
                float Q = data_Q[i];
                float nEvents = data_nEvents[i];
                float Q_ref = Q/nEvents * nEvents_test;

                if (abwff <= 1e10 && rayff >= 1e10){
                    double alpha_abs = truth_alpha(401.9, abwff, rayff);
                    data_xval.push_back(alpha_abs);
                    data_yval.push_back(Q_ref);
                    err_xval.push_back(0);
                    //note: maybe have to change this 1000 by the actual true nEvents -> not difficult to do
                    err_yval.push_back(TMath::Sqrt(Q_ref * (1 - Q_ref/1000)));
                }
                //Scattering data points
                if (abwff >= 1e10 && rayff <= 1e10){
                    double alpha_scat = truth_alpha(401.9, abwff, rayff);
                    data_scat_xval.push_back(alpha_scat);
                    data_scat_yval.push_back(Q_ref);
                    err_scat_xval.push_back(0);
                    err_scat_yval.push_back(TMath::Sqrt(Q_ref * (1 - Q_ref/1000)));
                }
            }
            int nPars = 2; //Amplitude and R the distance from the mPMT
            //This is the absorption part
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
            const double * res = min->X();
            const double * err = min->Errors();
            if (min->Status() == 0) {
                list_A.push_back(res[0]);
                list_R.push_back(res[1]);
            }
            //put the information in root objects
            TF1 func = chi->getFunction(0, 220, "best_fit_abs");
            TGraphErrors *data = new TGraphErrors(data_xval.size(), &data_xval[0], &data_yval[0], &err_xval[0], &err_yval[0]);
            TGraphErrors *fit_output;
            if (min->Status() == 0){
                fit_output = new TGraphErrors(1, &res[0], &res[1], &err[0], &err[1]);
            }
            else {
                //save 0 if the fit didn't work
                auto a = res[0] * 0;
                fit_output = new TGraphErrors(1, &a, &a, &a, &a);
            }


            //This is the scattering part
            nPars = 5; // the y of the 5  nodes for the fake-data spline
            Chisq *chi_scat = new Chisq(nPars);
            chi_scat->setData(data_scat_xval, data_scat_yval);
            ROOT::Math::Functor f_scat(chi_scat, &Chisq::fcn_rayff, nPars);
            ROOT::Math::Minimizer *min_scat = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
            min_scat->SetStrategy(3);
            min_scat->SetFunction(f_scat);
            min_scat->SetMaxFunctionCalls(10000);
            //At the moment we are working with 5 nodes the y of these 5 nodes is what we fit for the scattering
            min_scat->SetVariable(0, "a", 0.0, 0.01);
            min_scat->SetVariable(1, "b", 1.0, 0.01);
            min_scat->SetVariable(2, "c", 0.0, 0.01);
            min_scat->SetVariable(3, "d", 1.0, 0.01);
            min_scat->SetVariable(4, "e", 1.0, 0.01);
            std::cout << "Results scat: " << std::endl;
            min_scat->Minimize();
//             min_scat->PrintResults();
            //Save as a TGraphError the charge collected at this given position as a function
            //of the attenuation and then scattering length
            const double * res_scat = min_scat->X();
            const double * err_scat = min_scat->Errors();
            double x_pos[5] = {5., 25., 75., 150., 220.};
            double x_err[5] = {0., 0., 0., 0., 0.};
            TGraphErrors *gr = new TGraphErrors(5,x_pos,res_scat,x_err,err_scat);
            gr->SetTitle("Scattering Spline nodes");
            TGraphErrors *data_scat = new TGraphErrors(data_scat_xval.size(), &data_scat_xval[0], &data_scat_yval[0],&err_scat_xval[0], &err_scat_yval[0]);
            //This is saving everything the fit output in the reference_root folder
            TFile *outf = new TFile(Form("reference_root/results_Abs_Scat_theta%s_phi%s_R%s.root", theta_test, phi_test, R_test), "RECREATE");

            TF1 func_scat = chi_scat->getFunction_rayff(0, 240, "best_fit_scat");
            func.SetTitle("absorption_best_fit");
            func.Write();
            func_scat.SetTitle("scattering_best_fit");
            func_scat.Write();
            gr->Write();
            data->SetTitle("Charge vs attenuation lenght");
            data->GetXaxis()->SetTitle("Attenuation lenght (cm)");
            data->GetYaxis()->SetTitle("Charge collected for 1000 photons");

            data_scat->SetTitle("Charge vs attenuation lenght");
            data_scat->GetXaxis()->SetTitle("Attenuation lenght (cm)");
            data_scat->GetYaxis()->SetTitle("Charge collected for 1000 photons");

            data->Write("data_distribution");
            data_scat->Write("data_scat_distribution");
            fit_output->SetTitle("Exponential fit output");
            fit_output->GetXaxis()->SetTitle("Amplitude (nb of charge collected per 1000 photon)");
            fit_output->GetYaxis()->SetTitle("Distance R to the mPMT dome (cm)");
            fit_output->Write("fit_output_xA_yR");
            outf->Close();
        }
        newfile.close();   //close the file object.
        //now plot the histogram - to check whether the fitted attenuation is consistent with
        //an exponential behaviour
        TH1* h_R = new TH1D("h_R", Form("Histogram of Fitted R - R True: %s", R_test), 100,0.,100.);
        TH1* h_A = new TH1D("h_A", Form("Histogram of Fitted A - R True: %s", R_test), 100,200.,400.);
        for (long unsigned int u=0; u<list_R.size(); u++){
            h_R->Fill(list_R[u]);
            h_A->Fill(list_A[u]);

        }
        h_R->GetXaxis()->SetTitle("Fitted R(cm)");
        h_A->GetXaxis()->SetTitle("Fitted max charge per 1000 photon");
        TFile f(Form("histogram_fit/histos_R%s.root", R_test), "RECREATE");
        h_A->Write();
        h_R->Write();
        f.Close();
    }
}
