#include <iostream>
#include "chisq.h"
#include "TSpline3.h"

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
    //y_pred.resize(xin.size());
    std::cout << "uwu"<<std::endl;
    ray_spline = Ref_splines_rayff;
//     for(int i=0; i<Ref_splines_rayff.size(); i++){
// //         ray_spline[i] = &Ref_splines_rayff[i];
//          memcpy(ray_spline[i], Ref_splines_rayff[i], 1);
//     }
    std::cout << ray_spline[0]<<std::endl;
}

void Chisq::setRef_scat(std::vector<std::vector<double>> nodes_X, std::vector<std::vector<double>> nodes_Y){
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
    //return xval*pars[0] + pars[1] + pars[2]*TMath::Gaus(xval, pars[3], pars[4], true);
    return pars[0] * TMath::Exp(- 1/xval * pars[1]); // + pars[2]; // + pars[2] + pars[3]- pars[4];

}

double Chisq::makePredictionX_abwff(int xval){
    //std::cout << A[xval] << std::endl;
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

        //at position i we are looking up the value of the reference spline
        //correspionding to this scattering length where i in theta, phi, R pos
//         std::cout << ray_spline.size() << std::endl;
//         y_pred[i] = ray_spline[i]->Eval(pars[0]);
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
		//std::cout << "y_pred is 0 " << std::endl;
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
            //std::cout << "y_pred is 0 " << std::endl;
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
                //std::cout << "y_pred is 0 " << std::endl;
        }
        else {
                ret_val += (y[i]-y_pred[i])*(y[i]-y_pred[i])/y_pred[i];
        }
//         std::cout << "xval: " << x[i] << " y_pred: " << y_pred[i] << " y: " << y[i] << std::endl;
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
        //        std::cout << "y_pred is 0 " << std::endl;
        }
        else {
                ret_val += (y[i]-y_pred[i])*(y[i]-y_pred[i])/y_pred[i];
        }
    }
    std::cout << std::endl;
    return 2*ret_val; //-2Ln(L)
//    return ret_val;
}

double Chisq::spline_4nodes(double xval)
{
    int nNodes = pars.size();
    double xx[nNodes], yy[nNodes];
    for (int i=0; i<=nNodes; i++){
        //use this for equidistand nodes: doesn't work as well
           // xx[i] = 5. + i*(215./nNodes);
        //the fit parameters are the y of the nodes
            yy[i] = pars[i];
    }
    xx[0] = 5;
    xx[1] = 25;
    xx[2] = 75;
    xx[3] = 150;
    xx[4] = 220;

    TSpline3 *spline3 = new TSpline3("Test",xx,yy,nNodes,"b1e1", (y[1]-y[0])/5, 0.);

    return spline3->Eval(xval);

//     return pars[0] + pars[1]* xval + pars[2] * xval * xval + pars[3] * xval * xval * xval +  pars[4] * 1/xval;
}

double Chisq::makePredictionX_rayleigh(int i){
    double yy[5], xx[5];
    for (int k = 0; k < spline_X[i].size(); k++){
        yy[k] = spline_Y[i][k];
        xx[k] = spline_X[i][k];
    }
    //QUESTION: what of the dreivatives at the end?
    TSpline3 *spline3 = new TSpline3("Test",xx,yy, 5,"b1e1", (yy[1]-yy[0])/20.,0.);
    //Here pars 0 is the sctatering length guess
    //std::cout << "Best_guess :" << spline3->Eval(pars[0]) << " True : " << y[i] << std::endl;
    return spline3->Eval(pars[0]);

}

    //pars[1] * TMath::Exp(pars[2] / xval) + pars[4] * TMath::Exp(-pars[3] / xval);
    //TSpline3("my scattering spline", xval, const TF1* func, Int_t n, const char* opt = 0, Double_t valbeg = 0, Double_t valend = 0)

// pars[2] * xval + pars[3] * (xval*xval)
  //+  pars[4] * xval * xval * xval *xval;

   /*Fit parameters:
   par[0-3]=X of nodes (to be fixed in the fit!)
   par[4-7]=Y of nodes
   par[8-9]=first derivative at begin and end (to be fixed in the fit!)
   */
   //double xx = xval;

   //double xn[4] = { 10.0 , 40.0, 100.0, 100.0 };// { pars[0], pars[1], pars[2], pars[3] };
   //double yn[4] = { pars[0], pars[1], pars[2], pars[3] };

   //double b1 = 10.; //pars[8];
   //double e1 = 10.; //pars[9];

   //TSpline3 sp3("sp3", xn, yn, 4, "b1e1", b1, e1);
   //std::cout << "yn [0] : " << yn[0] << "sp3.Eval: " << sp3.Eval(xx) << std::endl;
   //return pars[0] * TMath::Exp( - 1/xval * pars[1]); // pars[2] * TMath::Exp(- xval * pars[3]);//sp3.Eval(xx);


void Chisq::print(){
    std::cout<<"test"<<std::endl;
}

TF1 Chisq::getFunction(double xlow, double xhigh, const char* title = "best_fit"){
    int npoints = 1000;
    TF1 func(Form("%s",title), this, &Chisq::getPoint, xlow, xhigh, pars.size());
    func.SetNpx(npoints);
    return func;
}
TF1 Chisq::getFunction_rayff(double xlow, double xhigh, const char* title = "best_fit"){
    int npoints = 1000;
    TF1 func(Form("%s",title), this, &Chisq::getPoint_rayff, xlow, xhigh, pars.size());
    func.SetNpx(npoints);
    return func;
}

