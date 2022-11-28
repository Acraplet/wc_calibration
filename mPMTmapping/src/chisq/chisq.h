#include <iostream>
#include <vector>
#include "TMath.h"
#include "TF1.h"

class Chisq{
    public:
    Chisq(int npars);
    ~Chisq();
    void setData(std::vector<double> xin, std::vector<double> yin);
    void setPars(const double *pars);
    double makePredictionX(double xval);
    double getPoint(double* xval, double *parameters);
    void makePrediction();
    double fcn(const double *pars);
    void print();
    TF1 getFunction(double xlow, double xhigh);

    private:
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> y_pred;
    std::vector<double> pars;

};

