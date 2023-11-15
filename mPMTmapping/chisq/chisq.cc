#include <iostream>
#include "chisq.h"
#include <TSpline.h>

Chisq::Chisq(int npars){
    pars.resize(npars);
    ParameterList.clear();
    cathodeSpline.clear();
    cathodeAngleSpline.clear();
    cathodeVSpline.clear();

    nthreads = 1;
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

void Chisq::setRef(std::vector<double> Ain, std::vector<double> Rin, 
                   std::vector<double> Cosths_in, TVector3 spin, TVector3 sdin, TVector3 sdxin, TVector3 sdyin
                   /*,std::vector<TVector3> spdin, std::vector<double> sprofilein*/){
    setRef(Ain,Rin);
    Cosths.resize(Cosths_in.size());
    //SourcePMTDirection.resize(Cosths_in.size());
    //SourceProfile.resize(Cosths_in.size());
    for(int i=0; i<Cosths_in.size(); i++){
        Cosths[i] = Cosths_in[i];
        //SourcePMTDirection[i] = spdin[i].Unit();
        //SourceProfile[i] = sprofilein[i];
    }
    source_position = spin;
    source_direction = sdin.Unit();
    source_xaxis = sdxin.Unit();
    source_yaxis = sdyin.Unit();

    std::cout<<"source_direction = "<<source_direction.x()<<" "<<source_direction.y()<<" "<<source_direction.z()<<std::endl;
}

void Chisq::SetHitTemplate(TH3F* h) 
{
    hit_template=h;

    hit_template_direction.clear();
    for (int j=1;j<=hit_template->GetNbinsY();j++)
    {
        std::vector<TVector3> vDirs;
        for (int k=1;k<=hit_template->GetNbinsZ();k++)
        {
            double c = hit_template->GetYaxis()->GetBinCenter(j);
            double p = hit_template->GetZaxis()->GetBinCenter(k);
            TVector3 dir = c*source_direction+sqrt(1-c*c)*(cos(p*TMath::Pi()/180)*source_xaxis+sin(p*TMath::Pi()/180)*source_yaxis);
            vDirs.emplace_back(dir.Unit());
        }
        hit_template_direction.emplace_back(vDirs);
    }
}

void Chisq::SetHitTemplate(TH3F* h, TH3F* hf, TH3F* hft, TH3F* hnft, TH2F* hnbs) 
{
    SetHitTemplate(h);
    hit_template__f = hf;
    hit_template__ft = hft;
    hit_template__nft = hnft;
    hit_template_NotBlacksheet = hnbs;
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

void Chisq::AddParameters(ParameterType kType, int nPars)
{
    ParameterList.emplace_back(std::make_pair(kType, nPars));
}

ParameterType Chisq::GetParameterType(std::string pname)
{
    if (pname=="Norm")
        return kNorm;
    else if (pname=="Attenuation")
        return kAttenuation;
    else if (pname=="Cathode")
        return kCathode;
    else if (pname=="Source")
        return kSource;
    else if (pname=="Reflectivity")
        return kReflectivity;
    else if (pname=="SourceCathodeReflectivity")
        return kSourceCathodeReflectivity;
    else return kInValid;
}

void Chisq::LoadCathodeSpline(std::string fname)
{
    cathodeSpline.clear();
    cathodeAngleSpline.clear();
    cathodeVSpline.clear();

    TFile f(fname.c_str());
    if (!f.IsOpen()){
        std::cout << "Error, could not open cathode spline file: " << fname << std::endl;
        return;
    }

    // for (int i=0;i<x.size(); i++)
    // {
    //     TH3* h = (TH3*)f.Get(Form("CathodeSpline_PMT%i",i));
    //     if (h)
    //     {
    //         h->SetDirectory(nullptr);
    //         cathodeSpline[i] = std::unique_ptr<TH3>(h);
    //     }
    //     else
    //         std::cout << "Warning, could not find CathodeSpline_PMT" << i << std::endl;
    // }

    // for (int i=0;i<x.size(); i++)
    // {
    //     TH1* h = (TH1*)f.Get(Form("AngleSpline_PMT_%i",i));
    //     if (h)
    //     {
    //         h->SetDirectory(nullptr);
    //         cathodeAngleSpline[i] = std::unique_ptr<TH1>(h);
    //     }
    //     else
    //         std::cout << "Warning, could not find AngleSpline_PMT_" << i << std::endl;
    // }

    // TH1* h = (TH1*)f.Get("testSpline");
    // if (h)
    // {
    //     h->SetDirectory(nullptr);
    //     testSpline = std::unique_ptr<TH1>(h);
    // }
    // else
    //     std::cout << "Warning, could not find testSpline" << std::endl;

    for (int i=0;i<x.size(); i++)
    {
        TVector3* v = (TVector3*)f.Get(Form("v_spline_PMT_%i",i));
        if (v)
        {
            cathodeVSpline[i] = std::unique_ptr<TVector3>(v);
        }
        else
            std::cout << "Warning, could not find v_spline_PMT_" << i << std::endl;
    }

    f.Close();
}

void Chisq::LoadSourceBin(std::string fname, TGraph* source_profile_in)
{
    // BinManager bm(fname);
    // std::vector<int> binmap(x.size());
    // for (int i=0;i<x.size(); i++)
    // {
    //     binmap[i] = bm.GetBinIndex(std::vector<double>{Cosths[i]});
    // }
    // BinMap[kSource] = binmap;
    // BinManagerMap[kSource] = bm;

    source_profile = source_profile_in;
    hit_template_profile_weight.clear();
    for (int j=0;j<hit_template->GetNbinsY();j++)
    {   
        std::vector<double> wgt;
        for (int k=0;k<hit_template->GetNbinsZ();k++)
        {
            double c = source_direction.Dot(hit_template_direction[j][k]);
            wgt.emplace_back(source_profile->Eval(c));
        }
        hit_template_profile_weight.emplace_back(wgt);
    }

}

double Chisq::CalcChiSq(const double *pars)
{
    // Reset prediction
#pragma omp parallel for num_threads(nthreads)
    for(int i=0; i<x.size(); i++){
        y_pred[i] = A[i];
    }

    int p = 0;
    // Make prediction
    for (auto& k : ParameterList)
    {
        switch (k.first)
        {
            case kNorm: 
#pragma omp parallel for num_threads(nthreads)
                for(int i=0; i<x.size(); i++){
                    y_pred[i] *= pars[p];
                }
                p+=k.second;
                break;
            case kAttenuation:
                for(int i=0; i<x.size(); i++){
                    y_pred[i] *= TMath::Exp(- 1/pars[p] * R[i]); 
                }
                p+=k.second;
                break;
            case kCathode:
#pragma omp parallel for num_threads(nthreads)
                for(int i=0; i<x.size(); i++){
                    // if (cathodeSpline[i])
                    //     y_pred[i] *= cathodeSpline[i]->Interpolate(pars[p],pars[p+1],pars[p+2]); 
                    // if (cathodeAngleSpline[i])
                    // {
                    //     double wgt = 0;
                    //     for (int j=0;j<k.second;j++)
                    //     {
                    //         wgt += cathodeAngleSpline[i]->GetBinContent(j+1)*pars[p+j];
                    //     }
                    //     y_pred[i] *= wgt;
                    // }
                    // if (testSpline)
                    //     y_pred[i] *= (1.-testSpline->GetBinContent(i+1))*pars[p] + testSpline->GetBinContent(i+1)*pars[p+1] ;
                    // if (cathodeVSpline[i])
                    // {
                    //     y_pred[i] *= (1-cathodeVSpline[i]->x())*pars[p] + cathodeVSpline[i]->x()*pars[p+1]*(cathodeVSpline[i]->y()+cathodeVSpline[i]->z()*pars[p+2]);
                    // }
                    y_pred[i] *= pars[p]+(hit_template_f->GetBinContent(i+1))*(pars[p+1]-pars[p]);
                }
                p+=k.second;
                break;
            case kSource:
                {      
                    TVector3 source_direction_new = cos(pars[p]*TMath::Pi()/180)*source_direction+sin(pars[p]*TMath::Pi()/180)*(cos(pars[p+1])*source_xaxis+sin(pars[p+1])*source_yaxis);
                    std::vector<std::vector<double>> weight(hit_template->GetNbinsY(),std::vector<double>(hit_template->GetNbinsZ()));
                    for (int j=0;j<hit_template->GetNbinsY();j++)
                    {   
#pragma omp parallel for num_threads(nthreads)
                        for (int k=0;k<hit_template->GetNbinsZ();k++)
                        {
                            double c = source_direction_new.Dot(hit_template_direction[j][k]);
                            weight[j][k] = (TMath::Gaus(c,1,(1-cos(TMath::Pi()*pars[p+2]/180))/sqrt(2.*log(2)))/hit_template_profile_weight[j][k]);
                        }
                    }
#pragma omp parallel for num_threads(nthreads)
                    for(int i=0; i<x.size(); i++){
                        // double cosths_new = SourcePMTDirection[i].Dot(source_direction_new);
                        // //y_pred[i] *= pars[p+2+BinManagerMap[kSource].GetBinIndex(std::vector<double>{cosths_new})]*source_profile->Eval(cosths_new)/source_profile->Eval(Cosths[i]);
                        // y_pred[i] *= TMath::Gaus(cosths_new,1,(1-cos(TMath::Pi()*pars[p+2]/180))/sqrt(2.*log(2)))/source_profile->Eval(Cosths[i]);
                        double val = 0;
                        for (int j=0;j<hit_template->GetNbinsY();j++)
                            for (int k=0;k<hit_template->GetNbinsZ();k++)
                            {
                                val += hit_template->GetBinContent(i+1,j+1,k+1)*weight[j][k];//*(pars[p+3]+(hit_template__f->GetBinContent(i+1,j+1,k+1))*(pars[p+4]-pars[p+3]));
                            }
                        //std::cout<<"val = "<<val<<std::endl;
                        y_pred[i] *= val;
                    }
                    p+=k.second;
                    break;
                }
            case kReflectivity: 
#pragma omp parallel for num_threads(nthreads)
                for(int i=0; i<x.size(); i++){
                    y_pred[i] *= pars[p] + (1-pars[p])*hit_template_tof->GetBinContent(i+1);
                }
                p+=k.second;
                break;
            case kSourceCathodeReflectivity:
                {      
                    TVector3 source_direction_new = cos(pars[p]*TMath::Pi()/180)*source_direction+sin(pars[p]*TMath::Pi()/180)*(cos(pars[p+1])*source_xaxis+sin(pars[p+1])*source_yaxis);
                    std::vector<std::vector<double>> weight(hit_template->GetNbinsY(),std::vector<double>(hit_template->GetNbinsZ()));
                    for (int j=0;j<hit_template->GetNbinsY();j++)
                    {   
#pragma omp parallel for num_threads(nthreads)
                        for (int k=0;k<hit_template->GetNbinsZ();k++)
                        {
                            double c = source_direction_new.Dot(hit_template_direction[j][k]);
                            weight[j][k] = (TMath::Gaus(c,1,(1-cos(TMath::Pi()*pars[p+2]/180))/sqrt(2.*log(2)))/hit_template_profile_weight[j][k]);
                        }
                    }
#pragma omp parallel for num_threads(nthreads)
                    for(int i=0; i<x.size(); i++){
                        double val = 0;
                        for (int j=0;j<hit_template->GetNbinsY();j++)
                            for (int k=0;k<hit_template->GetNbinsZ();k++)
                            {
                                double pw = hit_template->GetBinContent(i+1,j+1,k+1)*weight[j][k];
                                double f = hit_template__f->GetBinContent(i+1,j+1,k+1);
                                double s1 = pars[p+3], s2 = pars[p+4];
                                if (hit_template_NotBlacksheet->GetBinContent(j+1,k+1)>0.5)
                                {
                                    val += pw*( (1-f)*s1 + f*s2 );
                                }
                                else
                                {
                                    double ft = hit_template__ft->GetBinContent(i+1,j+1,k+1);
                                    double nft = hit_template__nft->GetBinContent(i+1,j+1,k+1);
                                    double reflec = pars[p+5];
                                    val += pw*( (1-f)*s1*( nft + (1-nft)*reflec ) + f*s2*( ft + (1-ft)*reflec ) );
                                }
                                
                            }
                        //std::cout<<"val = "<<val<<std::endl;
                        y_pred[i] *= val;
                    }
                    p+=k.second;
                    break;
                }
            default:
                // do nothing
                p+=k.second;
        }
    }

    // Calculate chi2
    double ret_val = 0;
    for(int i=0; i<y.size(); i++){
        double chi2 = 0;
        if (y_pred[i]>0)
        {
            chi2 = 2*(y_pred[i]-y[i]);
            if (y[i]>0)
                chi2 += 2*y[i]*std::log(y[i]/y_pred[i]);
        }
        if (chi2>0) ret_val += chi2;
    }
    return ret_val; //-2Ln(L)
}

void Chisq::PrintBins(TH1D* hPostfit, TH1D* hChi2)
{
    //hPostfit = new TH1D("","",y.size(),0,y.size());
    //hChi2 = new TH1D("","",y.size(),0,y.size());
    for(int i=0; i<y.size(); i++){
        double chi2 = 0;
        if (y_pred[i]>0)
        {
            chi2 = 2*(y_pred[i]-y[i]);
            if (y[i]>0)
                chi2 += 2*y[i]*std::log(y[i]/y_pred[i]);
        }
        hPostfit->SetBinContent(i+1,y_pred[i]);
        hChi2->SetBinContent(i+1,chi2);
        //std::cout<<"Bin "<<i<<" PMT "<<x[i]<<" cosths = "<< Cosths[i] <<" nom = "<<A[i]<<" pred = "<<y_pred[i]<<" data = "<<y[i]<<" chi2 = "<<chi2<<std::endl;
    }
}

bool Chisq::RemoveExtremeBins(double chi2_cut)
{
    bool result = false;
    std::vector<int> remove_list;
    for(int i=0; i<y.size(); i++){
        double chi2 = 0;
        if (y_pred[i]>0)
        {
            chi2 = 2*(y_pred[i]-y[i]);
            if (y[i]>0)
                chi2 += 2*y[i]*std::log(y[i]/y_pred[i]);
        }
        if (chi2>chi2_cut)
        {
            std::cout<<"Remove Bin "<<i<<" PMT "<<x[i]<<" cosths = "<< Cosths[i] <<" nom = "<<A[i]<<" pred = "<<y_pred[i]<<" data = "<<y[i]<<" chi2 = "<<chi2<<std::endl;
            result = true;
            remove_list.push_back(i);
            A[i] = 0;
        }
    }
    // for (int i=remove_list.size()-1;i>=0;i--)
    // {
    //     x.erase(x.begin()+remove_list[i]);
    //     R.erase(R.begin()+remove_list[i]);
    //     A.erase(A.begin()+remove_list[i]);
    //     Cosths.erase(Cosths.begin()+remove_list[i]);
    //     y.erase(y.begin()+remove_list[i]);
    //     y_pred.erase(y_pred.begin()+remove_list[i]);
    //     SourcePMTDirection.erase(SourcePMTDirection.begin()+remove_list[i]);
    // }
    return result;
}