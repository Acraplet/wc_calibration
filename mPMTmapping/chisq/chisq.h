#include <iostream>
#include <vector>
#include <map>
#include <memory>
#include "TMath.h"
#include "TF1.h"
#include "TFile.h"
#include "TH2.h"
#include "TH3.h"
#include "TVector3.h"
#include "TGraph.h"
//#include "BinManager.hh"

#ifdef _OPENMP
#include <omp.h>
#endif

enum ParameterType
{
    kNorm = 0,
    kAttenuation = 1,
    kRayleigh = 2,
    kCathode = 3,
    kSource = 4,
    kReflectivity = 5,
    kSourceCathodeReflectivity = 6,

    kInValid = 999,
};

class Chisq{
    public:
    Chisq(int npars=0);
    ~Chisq();
    void setData(std::vector<double> xin, std::vector<double> yin);
    void setRef(std::vector<double> Ain, std::vector<double> Rin);
    void setRef(std::vector<double> Ain, std::vector<double> Rin, std::vector<double> Cosths_in, TVector3 spin, TVector3 sdin, TVector3 sdxin, TVector3 sdyin/*, std::vector<TVector3> spdin, std::vector<double> sprofilein*/);
    void setRef_spline(std::vector<TF1*> Ref_splines_rayff);
    void setRef_scat(std::vector<std::vector<double>> spline_X_ref, std::vector<std::vector<double>> spline_Y_ref);
    void setPars(const double *pars);
    void setNodesX(std::vector<double> xx);
    double makePredictionX(double xval);
    double makePredictionX_abwff(int xval);
    double makePredictionX_rayleigh(int xval);
    double spline_4nodes(double xval);
    double getPoint(double* xval, double *parameters);
    double getPoint_rayff(double* xval, double *parameters);
    void makePrediction();
    void makePrediction_abwff();
    void makePrediction_rayff();
    void makePrediction_rayleigh();
    double fcn(const double *pars);
    double fcn_abwff(const double *pars);
    double fcn_rayff(const double *pars);
    double fitter_rayleigh(const double *pars);
    void print();
    TF1 getFunction(double xlow, double xhigh, const char* title);
    TF1 getFunction_rayff(double xlow, double xhigh, const char* title);
    void AddParameters(ParameterType kType, int nPars);
    ParameterType GetParameterType(std::string pname);
    double CalcChiSq(const double *pars);
    void LoadCathodeSpline(std::string fname);
    void LoadCathodeSpline(TH1F* h) { hit_template_f = h; }
    void LoadReflectivitySpline(TH1F* h) { hit_template_tof = h; }
    void LoadSourceBin(std::string fname, TGraph* source_profile_in);
    void SetNPars(int val) {pars.resize(val);}
    void PrintBins(TH1D* hPostfit, TH1D* hChi2);
    bool RemoveExtremeBins(double chi2_cut);
    void SetHitTemplate(TH3F* h);
    void SetHitTemplate(TH3F* h, TH3F* hf, TH3F* hft, TH3F* hnft, TH2F* hnbs);
    void SetNTreads(int n) {nthreads = n;}

    private:
    std::vector<double> x;
    std::vector<double> R;
    std::vector<double> A;
    std::vector<double> Cosths;
    std::vector<double> y;
    std::vector<TF1*> ray_spline;
    std::vector<std::vector<double>> spline_X;
    std::vector<std::vector<double>> spline_Y;
    std::vector<double> y_pred;
    std::vector<double> pars;
    int nNodes;
    std::vector<double> XNodes;
    std::vector<std::pair<ParameterType,int>> ParameterList;
    std::map<int, std::unique_ptr<TH3>> cathodeSpline;
    std::map<int, std::unique_ptr<TH1>> cathodeAngleSpline;
    std::map<int, std::unique_ptr<TVector3>> cathodeVSpline;
    std::unique_ptr<TH1> testSpline;
    //std::map<ParameterType, std::vector<int>> BinMap;
    //std::map<ParameterType, BinManager> BinManagerMap;
    TVector3 source_position;
    TVector3 source_direction;
    TVector3 source_xaxis;
    TVector3 source_yaxis;
    std::vector<std::vector<TVector3>> hit_template_direction;
    std::vector<std::vector<double>> hit_template_profile_weight;
    //std::vector<TVector3> SourcePMTDirection;
    //std::vector<double> SourceProfile;
    TGraph* source_profile;
    TH3F* hit_template;
    TH3F* hit_template__f;
    TH3F* hit_template__ft;
    TH3F* hit_template__nft;
    TH1F* hit_template_f;
    TH1F* hit_template_tof;
    TH2F* hit_template_NotBlacksheet;
    int nthreads;

};

