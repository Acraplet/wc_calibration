#include <iostream>
#include <vector>
#include "TMath.h"
#include "TF1.h"

enum ParameterType
{
    kNorm = 0,
    kAttenuation = 1,
    kRayleigh = 2,
    kCathode = 3,
};

class Chisq{
    public:
    Chisq(int npars);
    ~Chisq();
    void setData(std::vector<double> xin, std::vector<double> yin);
    void setRef(std::vector<double> Ain, std::vector<double> Rin);
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
    void AddParameters(ParameterType kType);
    double CalcChiSq(const double *pars);

    private:
    std::vector<double> x;
    std::vector<double> R;
    std::vector<double> A;
    std::vector<double> y;
    std::vector<TF1*> ray_spline;
    std::vector<std::vector<double>> spline_X;
    std::vector<std::vector<double>> spline_Y;
    std::vector<double> y_pred;
    std::vector<double> pars;
    int nNodes;
    std::vector<double> XNodes;
    std::vector<ParameterType> ParameterList;
};

