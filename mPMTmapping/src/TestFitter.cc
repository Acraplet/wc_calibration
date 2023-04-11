//This code is extracting the absorption length corresponding to  one or more maps that have the same abwff
//For this is is looking up in reference_maps what the max value of the charge can be at each position 
//(which is also dependant of R) when there is no absorption or scattering
//The value of R is also stored. Then the chi is minimised between 
//the observed data and y_pred = A_max * exp(-R/absorption length pred) to give the expected absorption length
//So far this is position dependant but it will have to be upgraded to a binned approach
#include "truth_alpha.hh"
#include <iostream>
#include "chisq.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TSystemDirectory.h"
#include <cstring>

#include "ColorOutput.hh"

const std::string TAG = color::GREEN_STR + "[TestFitter]: " + color::RESET_STR;
const std::string ERR = color::RED_STR + "[ERROR]: " + color::RESET_STR;
const std::string WAR = color::RED_STR + "[WARNING]: " + color::RESET_STR;

using namespace std;

int main(int argc, char **argv){
    char* R_test;
    char* Q_test;
    char* theta_test;
    char* phi_test;
    double  theta_test_num, Q_test_num, R_test_num;
    std::vector<double> list_R;
    std::vector<double> list_A;
    std::vector<double> list_Q;
    std::vector<double> list_i;

    char * filename=NULL;
    char * reffilename=NULL;

    char c;
    while( (c = getopt(argc,argv,"f:r:h")) != -1 ){//input in c the argument (-f etc...) and in optarg the next argument. When the above test becomes -1, it means it fails to find a new argument.
        switch(c){
        case 'f':
            filename = optarg;
            break;
        case 'r':
            reffilename = optarg;
            break;
        case 'h':
            std::cout << TAG << "Print help message"<<std::endl;
            break;
        default:
            return 0;
        }
    }

    // Open the data file
    if (filename==NULL){
        std::cout << ERR << "Error, no input file: " << std::endl;
        return -1;
    }

    TFile* f = new TFile(filename);
    if (!f->IsOpen()){
        std::cout << ERR << "Error, could not open input file: " << filename << std::endl;
        return -1;
    }

    TTree* pmt_type0 = (TTree*)f->Get("pmt_type0");
    TTree* hitRate_pmtType0 = (TTree*)f->Get("hitRate_pmtType0");
    double R;
    pmt_type0->SetBranchAddress("R",&R);
    int nPMTs = pmt_type0->GetEntries();
    hitRate_pmtType0->Draw(Form("PMT_id>>(%i,0,%i)",nPMTs,nPMTs));
    TH1D* hData = (TH1D*)hitRate_pmtType0->GetHistogram();
    for (int i = 0 ; i<nPMTs; i++)
    {
        pmt_type0->GetEntry(i);
        list_i.push_back(i);
        list_R.push_back(R);
        list_Q.push_back(hData->GetBinContent(i+1));
    }

    // Open the ref file
    if (reffilename==NULL){
        std::cout << ERR << "Error, no ref file: " << std::endl;
        return -1;
    }

    TFile* fr = new TFile(reffilename);
    if (!fr->IsOpen()){
        std::cout << ERR << "Error, could not open ref file: " << reffilename << std::endl;
        return -1;
    }
    TTree* hitRate_pmtType0r = (TTree*)fr->Get("hitRate_pmtType0");
    hitRate_pmtType0r->Draw(Form("PMT_id>>(%i,0,%i)",nPMTs,nPMTs));
    TH1D* hRef = (TH1D*)hitRate_pmtType0r->GetHistogram();
    for (int i = 0 ; i<nPMTs; i++)
    {
        list_A.push_back(hRef->GetBinContent(i+1));
    }

    for (int i = 0 ; i<nPMTs; i++)
    {
        std::cout<<"i = "<<i<<" "<<list_R[i]<<" "<<list_Q[i]<<" "<<list_A[i]<<"\n";
    }

    //Here the fitting begins
    const int nPars = 1; //the only parameter we fit is absorption length
    Chisq *chi = new Chisq(nPars);
    chi->setData(list_i, list_Q);
    chi->setRef(list_A, list_R);
    chi->AddParameters(kAttenuation);
    ROOT::Math::Functor functor(chi, &Chisq::CalcChiSq, nPars);
    ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    min->SetStrategy(3);
    min->SetFunction(functor);
    min->SetMaxFunctionCalls(10000);
    min->SetVariable(0, "attenuation_length", 5.0, 0.01);

    std::cout << TAG << "Number of defined parameters: " << min->NDim() << std::endl
              << TAG << "Number of free parameters   : " << min->NFree() << std::endl
              << TAG << "Number of fixed parameters  : " << min->NDim() - min->NFree()
              << std::endl;

    min->Minimize();
    min->PrintResults();
    const int ndim        = min->NDim();
    const int nfree       = min->NFree();
    const double* par_val = min->X();
    const double* par_err = min->Errors();
    for (int i=0;i<nPars;i++)
    {
        std::cout<<"Variable "<<i<<" : "<<par_val[i]<<" +/- "<<par_err[i]<<"\n";
    }
}

