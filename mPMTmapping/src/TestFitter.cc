// This code reads the "data" generated by WCSIM_TreeConvert, which stores the photon hits per PMT,
// then reads the "reference map" (again generated by WCSIM_TreeConvert), and makes predictions based on the fit parameters.
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
#include "TH1.h"
#include <TMatrixDSym.h>
#include "TSystemDirectory.h"
#include <cstring>

#include "toml/toml_helper.h"
#include "ColorOutput.hh"

const std::string TAG = color::GREEN_STR + "[TestFitter]: " + color::RESET_STR;
const std::string ERR = color::RED_STR + "[ERROR]: " + color::RESET_STR;
const std::string WAR = color::RED_STR + "[WARNING]: " + color::RESET_STR;

using namespace std;

int main(int argc, char **argv){

    std::vector<double> list_R;
    std::vector<double> list_A;
    std::vector<double> list_Q;
    std::vector<double> list_i;

    char * filename=NULL;
    char * reffilename=NULL;
    char * splinefilename=NULL;

    std::string config_file;
    std::string outfilename;

    char c;
    while( (c = getopt(argc,argv,"f:r:s:c:o:h")) != -1 ){//input in c the argument (-f etc...) and in optarg the next argument. When the above test becomes -1, it means it fails to find a new argument.
        switch(c){
            case 'c':
                config_file = optarg;
                break;
            case 'f': // data
                filename = optarg;
                break;
            case 'o': // output file
                outfilename = optarg;
                break;
            case 'r': // reference map
                reffilename = optarg;
                break;
            case 's': // spline file for cathode parameters
                splinefilename = optarg;
                break;
            case 'h':
                std::cout << TAG    << "USAGE: "
                                    << "TestFitter" << "\nOPTIONS:\n"
                                    << "-c : config file\n"
                                    << "-f : Data file\n"
                                    << "-r : Reference map file\n"
                                    << "-s : Spline file for cathode parameters\n"
                                    << "-o : Output file\n"
                                    ;
                return 0;
            default:
                return 0;
        }
    }

    if (config_file.size()==0)
    {
        std::cout<< ERR << "No config file!" << std::endl;
        return -1;
    }
        if (outfilename.size()==0)
    {
        std::cout<< TAG << "Output file name default to: fitoutput.root" << std::endl;
        outfilename = "fitoutput.root";
    }
    auto const &card_toml = toml_h::parse_card(config_file);
    auto const &fitparameters_config = toml_h::find(card_toml, "fitparameters");

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
    // distance to source
    double R;
    pmt_type0->SetBranchAddress("R",&R);
    int nPMTs = pmt_type0->GetEntries();
    // hit rate histogram
    hitRate_pmtType0->Draw(Form("PMT_id>>(%i,0,%i)",nPMTs,nPMTs),"nPE");
    TH1D* hData = (TH1D*)hitRate_pmtType0->GetHistogram();
    // Fill data entries for fitter
    for (int i = 0 ; i<nPMTs; i++)
    {
        pmt_type0->GetEntry(i);
        list_i.push_back(i);
        list_R.push_back(R);
        list_Q.push_back(hData->GetBinContent(i+1));
    }
    f->Close();

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
    hitRate_pmtType0r->Draw(Form("PMT_id>>(%i,0,%i)",nPMTs,nPMTs),"nPE");
    TH1D* hRef = (TH1D*)hitRate_pmtType0r->GetHistogram();
    // Fill reference map
    for (int i = 0 ; i<nPMTs; i++)
    {
        list_A.push_back(hRef->GetBinContent(i+1));
    }
    fr->Close();

    // for (int i = 0 ; i<nPMTs; i++)
    // {
    //     std::cout<<"i = "<<i<<" "<<list_R[i]<<" "<<list_Q[i]<<" "<<list_A[i]<<"\n";
    // }

    //Here the fitting begins
    Chisq *chi = new Chisq();
    chi->setData(list_i, list_Q);
    chi->setRef(list_A, list_R);
    int nPars = 0; // number of fit parameters 
    // Add fit parameters
    std::vector<std::string> paramNames;
    std::vector<double> paramPriors;
    std::vector<double> paramSteps;
    std::vector<double> paramLows;
    std::vector<double> paramHighs;
    std::vector<bool> paramFixeds;
    for (auto const &name : toml_h::find<std::vector<std::string>>(fitparameters_config, "names"))
    {
        std::cout << TAG << "Parameter name: " << name << std::endl;
        auto const &ele = toml_h::find<toml::array>(fitparameters_config, name);
        auto paramtype = chi->GetParameterType(toml_h::find<std::string>(ele,0));
        //std::cout<<"paramtype = "<<paramtype<<std::endl;
        if (paramtype==kCathode)
        {
            std::cout << TAG << "Load spline file " << splinefilename << std::endl;
            chi->LoadCathodeSpline(splinefilename);
        }
        else if (paramtype==kInValid)
        {
            std::cout << WAR << "Invalid parameter type, skipping this" << std::endl;
            continue;
        }
        
        auto npar = toml_h::find<int>(ele,1);
        nPars += npar;

        chi->AddParameters(paramtype, npar);

        auto par_setup = toml_h::find<toml::array>(ele,2);
        for (auto const &par : par_setup)
        {
            paramNames.push_back(toml_h::find<std::string>(par,0));
            paramPriors.push_back(toml_h::find<double>(par,1));
            paramSteps.push_back(toml_h::find<double>(par,2));
            paramLows.push_back(toml_h::find<double>(par,3));
            paramHighs.push_back(toml_h::find<double>(par,4));
            paramFixeds.push_back(toml_h::find<bool>(par,5));

            //std::cout<<paramNames.back()<<" "<<paramPriors.back()<<" "<<paramSteps.back()<<" "<<paramLows.back()<<" "<<paramHighs.back()<<" "<<paramFixeds.back()<<"\n";
        }
    }
    chi->SetNPars(nPars);

    ROOT::Math::Functor functor(chi, &Chisq::CalcChiSq, nPars);
    ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    min->SetStrategy(3);
    min->SetFunction(functor);
    min->SetMaxFunctionCalls(1000000);
    min->SetPrintLevel(2);
    // Set parameter priors
    for (int i=0;i<nPars;i++)
    {
        min->SetVariable(i, paramNames[i].c_str(), paramPriors[i], paramSteps[i]);
        min->SetVariableLimits(i, paramLows[i], paramHighs[i]);
        if (paramFixeds[i]) min->FixVariable(i);
    }

    std::cout << TAG << "Number of defined parameters: " << min->NDim() << std::endl
              << TAG << "Number of free parameters   : " << min->NFree() << std::endl
              << TAG << "Number of fixed parameters  : " << min->NDim() - min->NFree()
              << std::endl;

    min->Minimize();
    min->PrintResults();
    const double* par_val = min->X();
    const double* par_err = min->Errors();
    const int ndim = nPars;
    double cov_array[ndim * ndim];
    min->GetCovMatrix(cov_array);
    TMatrixDSym cov_matrix;
    TMatrixDSym cor_matrix;
    cov_matrix.ResizeTo(ndim,ndim);
    cov_matrix = TMatrixDSym(ndim, cov_array);
    cor_matrix.ResizeTo(ndim,ndim);
    for(int r = 0; r < ndim; ++r)
    {
        for(int c = 0; c < ndim; ++c)
        {
            cor_matrix[r][c] = cov_matrix[r][c] / std::sqrt(cov_matrix[r][r] * cov_matrix[c][c]);
            if(std::isnan(cor_matrix[r][c]))
                cor_matrix[r][c] = 0;
        }
    }
    // for (int i=0;i<nPars;i++)
    // {
    //     std::cout<<"Variable "<<i<<" : "<<par_val[i]<<" +/- "<<par_err[i]<<"\n";
    // }

    TFile* fo = TFile::Open(outfilename.c_str(),"RECREATE");
    TH1D* hist_fit = new TH1D("fit_result","fit_result",nPars,0,nPars);
    for (int i=0;i<nPars;i++)
    {
        hist_fit->SetBinContent(i+1,par_val[i]);
        hist_fit->SetBinError(i+1,par_err[i]);
    }
    hist_fit->Write();
    cov_matrix.Write("res_cov_matrix");
    cor_matrix.Write("res_cor_matrix");
    fo->Close();
}

