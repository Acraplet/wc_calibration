#include <iostream>
#include "chisq.h"
#include "TSpline.h"

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

void Chisq::setNodesX(std::vector<double> xx){
	nNodes = xx.size();
	for (double node : xx) XNodes.push_back(node);
}

void Chisq::setRef(std::vector<double> Ain, std::vector<double> Rin){
    if(Ain.size() != Rin.size()){
        std::cerr<<"Chisq::setReference(A, R) mismatch between size of A and R vectors"<<std::endl;
        throw;
    }
    A.resize(Ain.size());
    R.resize(Rin.size());
    //y_pred.resize(xin.size());
    for(int i=0; i<Ain.size(); i++){
        A[i] = Ain[i];
        R[i] = Rin[i];
    }
}

void Chisq::setRef_spline(std::vector<TF1*> Ref_splines_rayff){
//     /*std::cout*/ << "uwu"<<std::endl;
    ray_spline = Ref_splines_rayff;
//     std::cout << ray_spline[0]<<std::endl;
}

void Chisq::setRef_scat(std::vector<std::vector<double>> nodes_X, std::vector<std::vector<double>> nodes_Y){
     //This is storing in the reference nodes position and their positions
    //carefull! here nodes_X is a vecotr of vector of nodes {for each reference position}
    //TODO: change this so we can use reference files that have different number of nodes each time (a simple array with nNodes also being looked up as a function of i) - for now not nmecessary
     nNodes = nodes_X[0].size();
//      std::cout << nNodes << std::endl;
     //here we append the reference for each of our different distances R
     for(int i=0; i<nodes_X.size(); i++){
         spline_X.push_back(nodes_X[i]);
         spline_Y.push_back(nodes_Y[i]);
     }
}

void Chisq::setPars(const double *parameters){
    for(int i=0; i<pars.size(); i++) pars[i] = parameters[i];

}
double Chisq::makePredictionX(double xval){
    //put your user defined function here
    return pars[0] * TMath::Exp(- 1/xval * pars[1]); // + pars[2]; // + pars[2] + pars[3]- pars[4];

}

double Chisq::makePredictionX_abwff(int xval){
    return A[xval] * TMath::Exp(- 1/pars[0] * R[xval]); // + pars[2]; // + pars[2] + pars[3]- pars[4];

}


double Chisq::getPoint(double *xval, double *parameters){
    return makePredictionX(xval[0]);
}

double Chisq::getPoint_rayff(double *xval, double *parameters){
    return TMath::Max(spline_4nodes(xval[0]), 0.);
}

void Chisq::makePrediction(){
    //put your user defined function in here
    for(int i=0; i<x.size(); i++){
        y_pred[i] = makePredictionX(x[i]);
    }
}

void Chisq::makePrediction_abwff(){
    //put your user defined function in here
    for(int i=0; i<x.size(); i++){
        y_pred[i] = makePredictionX_abwff(i);
    }
}

void Chisq::makePrediction_rayleigh(){
    //put your user defined function in here
    for(int i=0; i<x.size(); i++){
        y_pred[i] = makePredictionX_rayleigh(i);
    }
}


void Chisq::makePrediction_rayff(){
    //put your user defined function in here
    for(int i=0; i<x.size(); i++){
        y_pred[i] = spline_4nodes(x[i]);
    }
}

double Chisq::fcn_abwff(const double *parameters){
    setPars(parameters);
    makePrediction_abwff();
    double ret_val = 0;
    for(int i=0; i<y.size(); i++){
	if (y[i]== 0. and y_pred[i] <= 10.e-5) continue;
	if (y_pred[i] <= 10e-5 ){
		ret_val += (y[i]-y_pred[i])*(y[i]-y_pred[i])/0.01;
	}
	else {
        	ret_val += (y[i]-y_pred[i])*(y[i]-y_pred[i])/y_pred[i];
	}
    }
    return 2*ret_val; //-2Ln(L)
}

double Chisq::fitter_rayleigh(const double *parameters){
    setPars(parameters);
    makePrediction_rayleigh();
    double ret_val = 0;
    for(int i=0; i<y.size(); i++){
        if (y[i]== 0. and y_pred[i] <= 10.e-5) continue;
        if (y_pred[i] <= 10e-5 ){
            ret_val += (y[i]-y_pred[i])*(y[i]-y_pred[i])/0.0001;
        }
        else {
            ret_val += (y[i]-y_pred[i])*(y[i]-y_pred[i])/y_pred[i];
        }
    }
    return 2*ret_val; //-2Ln(L)
}

double Chisq::fcn_rayff(const double *parameters){
    setPars(parameters);
    makePrediction_rayff();
    double ret_val = 0;
    for(int i=0; i<y.size(); i++){
        if (y[i]== 0. and y_pred[i] <= 10.e-5) continue;
        if (y_pred[i] <= 10e-5 ){
                ret_val += (y[i]-y_pred[i])*(y[i]-y_pred[i])/0.001;
        }
        else {
                ret_val += (y[i]-y_pred[i])*(y[i]-y_pred[i])/y_pred[i];
        }
    }

    return 2*ret_val; //-2Ln(L)
}


double Chisq::fcn(const double *parameters){
    setPars(parameters);
    makePrediction();
    //use a log-likelihood measure - not working so well...
    double ret_val = 0;
    for(int i=0; i<x.size(); i++){
	if (y[i]== 0. and y_pred[i] <= 10.e-5) continue;
        if (y_pred[i] <= 10e-5 ){
                ret_val += (y[i]-y_pred[i])*(y[i]-y_pred[i])/0.01;
        }
        else {
                ret_val += (y[i]-y_pred[i])*(y[i]-y_pred[i])/y_pred[i];
        }
    }
//std::cout << std::endl;
    return 2*ret_val; //-2Ln(L)
}

double Chisq::spline_4nodes(double xval)
{
    //int nNodes = pars.size();
    if (nNodes!= pars.size()){
	    //a bit of proofchechking that we have given the right number of nodes cooordinates
	    std::cout << "Number of nodes : " << nNodes << " is different to the number of fitted parameters : " << pars.size() << " which is an ISSUE " << std::endl;
    }
    double xx[nNodes], yy[nNodes];
    for (int i=0; i<=nNodes; i++){
            yy[i] = pars[i];
	    xx[i] = XNodes[i];
		
    }
    //These are the hard coded source positions: DONE: write a proper way to do it
//    xx[0] = 5;
//    xx[1] = 25;
//    xx[2] = 75;
//    xx[3] = 150;
//    xx[4] = 220;
    TSpline3 *spline3 = new TSpline3("Test",xx,yy,nNodes,"b1e1", (y[1]-y[0])/float(x[1]-x[0]), 0.);

    return spline3->Eval(xval);
}

double Chisq::makePredictionX_rayleigh(int i){
//     std::cout << "a" << std::endl;
    double yy[nNodes], xx[nNodes];
//     std::cout << spline_Y[i].size() << " " << nNodes<< std::endl;
    for (int k = 0; k < spline_X[i].size(); k++){
        yy[k] = spline_Y[i][k];
        xx[k] = spline_X[i][k];
//         std::cout << k << " "<<yy[k] << " "<< xx[k] <<std::endl;
    }
    //set the end point derivative to be 0 and the start point derivative to be the gradient betweeen the 
    //first two data points
    TSpline3 *spline3 = new TSpline3("Test",xx,yy, nNodes,"b1e1", (yy[1]-yy[0])/float(xx[1]-xx[0]), 0.);
    return spline3->Eval(pars[0]);

}

TF1 Chisq::getFunction(double xlow, double xhigh, const char* title = "best_fit"){
    int npoints = 1000;
    TF1 func(Form("%s",title), this, &Chisq::getPoint, xlow, xhigh, pars.size());
    func.SetNpx(npoints);
    return func;
}
TF1 Chisq::getFunction_rayff(double xlow, double xhigh, const char* title = "best_fit_scat"){
    int npoints = 1000;
    TF1 func(Form("%s",title), this, &Chisq::getPoint_rayff, xlow, xhigh, pars.size());
    func.SetNpx(npoints);
    return func;
}

