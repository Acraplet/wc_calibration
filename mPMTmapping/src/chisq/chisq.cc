#include <iostream>
#include "chisq.h"

Chisq::Chisq(int npars){
    pars.resize(npars);
}

Chisq::~Chisq(){
}

void Chisq::setData(std::vector<double> xin, std::vector<double> yin){
    if(xin.size() != yin.size()){
        std::cerr<<"Chisq::setData(x, y) mismatch between size of x and y vectors"<<std::endl;
        throw;
    }
    x.resize(xin.size());
    y.resize(yin.size());
    y_pred.resize(xin.size());
    for(int i=0; i<xin.size(); i++){
        x[i] = xin[i];
        y[i] = yin[i];
    }
}

void Chisq::setPars(const double *parameters){
    for(int i=0; i<pars.size(); i++) pars[i] = parameters[i];

}
double Chisq::makePredictionX(double xval){
    //put your user defined function here
    return xval*pars[0] + pars[1] + pars[2]*TMath::Gaus(xval, pars[3], pars[4], true);


}
double Chisq::getPoint(double *xval, double *parameters){
    return makePredictionX(xval[0]);
}
void Chisq::makePrediction(){
    //put your user defined function in here
    for(int i=0; i<x.size(); i++){
        y_pred[i] = makePredictionX(x[i]);
    }
}

double Chisq::fcn(const double *parameters){
    setPars(parameters);
    makePrediction();


    //use a log-likelihood measure
    double ret_val = 0;
    for(int i=0; i<x.size(); i++){
        if(y_pred[i] < 1e-5) continue;
        ret_val += (y_pred[i] - y[i]) + y[i]*TMath::Log(y[i]/y_pred[i]);
    }
    return 2*ret_val; //-2Ln(L)

    //use a chisq measure:
//    double ret_val = 0;
//    for(int i=0; i<x.size(); i++){
//        ret_val += pow(y[i]-y_pred[i], 2) / sqrt(y[i]);
//    }
//    return ret_val;
}


void Chisq::print(){
    std::cout<<"test"<<std::endl;
}

TF1 Chisq::getFunction(double xlow, double xhigh){
    int npoints = 1000;
    TF1 func("best_fit", this, &Chisq::getPoint, xlow, xhigh, pars.size());
    func.SetNpx(npoints);
    return func;
}

