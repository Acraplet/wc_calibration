//This is a cpp code reading the WCSim files counting the number of photons hitting the mPMT from all the angles
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <TROOT.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TRandom3.h>

#include "ColorOutput.hh"
// #include "../include/WCSimRootEvent.hh"


//This is a work in progress - the python file needs to be translated into cpp
//Need to:
/*
 - open all the data with a given configurations
 - read the output giving out: R, theta, phi, (source x,y,z + dir x,y,z)
 - then compile everything into a 2D map
*/
const std::string ERR = color::RED_STR + "[ERROR]: " + color::RESET_STR;
const std::string WAR = color::GREEN_STR + "[WARNING]: " + color::RESET_STR;

void HelpMessage()
{
    std::cout << "USAGE: "
    << "make_mPMTmap" << "\nOPTIONS:\n"
    << "-f : Input file\n" << std::endl
    << "-o : Output file\n" << std::endl;
}

int main(int argc, char **argv){
    bool verbose=false;
    char * filename=NULL;
    char * outfilename=NULL;
    long int endEvent=0;
    long int startEvent=0;
    char c;

    while( (c = getopt(argc,argv,"f:o:b:s:e:l:r:p:c:w:hdtv")) != -1 ){//input in c the argument (-f etc...) and in optarg the next argument. When the above test becomes -1, it means it fails to find a new argument.
        switch(c){
            case 'f':
                filename = optarg;
                break;
        }
        switch(c){
            case 'o':
                outfilename = optarg;
                break;
        }
    }

    if (filename==NULL){
        std::cout << ERR << "Error, no input file: " << std::endl;
        HelpMessage();
        return -1;
    }

    if (outfilename==NULL){
        std::cout << WAR << "Warning, no output file given, using test.txt " << std::endl;
        outfilename = "test.txt";
//         HelpMessage();
    }

    std::string x, y, z, t, p, R; //the cartesian and polar coordinate, careful R is from the surface of the mPMT dome and not from its centre of the mPMT
    int str; //this will help break down the string
    //To read the characteristics of the thing
    std::string T;
    std::stringstream X(filename);
    while(std::getline(X, T, '_')){
//         std::cout << T << std::endl; // print split string
        if (!T.find("x")){
            str = T.find("x");
            x = T.erase(0,str+1);
        }
        if (!T.find("y")){
            str = T.find("y");
            y = T.erase(0,str+1);
        }
        if (!T.find("z")){
            str = T.find("z");
            z = T.erase(0,str+1);
        }
        if (!T.find("t")){
            str = T.find("t");
            t = T.erase(0,str+1);
        }
        if (!T.find("p")){
            str = T.find("p");
            p = T.erase(0,str+1);
        }
        if (!T.find("R")){
            str = T.find("R");
            R = T.erase(0,str+1);
        }
    }
    std::cout << "Source position : " << x << " " << y << " " << z << " " << t << " " << p << " " << R << std::endl;

//     int pos = inputFilename.find("fitoutput");
//     std::string str3 = inputFilename.substr(pos);
//     int pos2 = str3.find(".root");
//     std::string str4 = str3.erase(pos2);
//     const char* folder = str4.c_str();
/*
    TChain *tree = new TChain("wcsimT");
    tree->Add(filename);*/

    TFile *infile = new TFile(filename, "READ");
    TTree *events = (TTree*)infile->Get("CherenkovHits");

    int event;
//     int PMT_QTot;
    int NHits;
    int nHitsTot = 0;

    events->SetBranchAddress("Event", &event);
    events->SetBranchAddress("NHits", &NHits);
//     events->SetBranchAddress("PMT_QTot", &PMT_QTot);


    int n_entries = events->GetEntries();

    for(int i=0; i<n_entries; i++){
        events->GetEntry(i);
//         std::cout<<"event: "<<event<< " NHits " << NHits <<std::endl;
        nHitsTot += NHits;
    }

    std::cout << "Total number of hits : "<< nHitsTot << std::endl;
    std::ofstream outfile;
    outfile.open(outfilename, std::ofstream::app);
    outfile <<  x << " " << y << " " << z << " " << t << " " << p << " " << R << " " << nHitsTot << " " << n_entries << std::endl;

    //To have the output file! - obviously will need to adapt
//     std::ofstream outfile;
//     outfile.open(Form("Calibration_output/CleanFiles_v6/Angular_results_%s.txt", folder), std::ios_base::trunc);//
//
//     for (int i=1; i<h_res->GetEntries()+1 ; i++){
//         outfile << "mPMT_angular_" << i << ": "<< h_res->GetBinContent(i) << " +/- " << h_err->GetBinContent(i) << std::endl;
}//End of main



/*
TTree *geotree = (TTree*)file->Get("wcsimGeoT");

TFile * outfile = new TFile(outfilename,"RECREATE");
std::cout<<TAG<<"File "<<outfilename<<" is open for writing"<<std::endl;

TTree* hitNb = new TTree("hitNb","hitNb");
hitNb->Branch("nPE",&nPE)

double vtxpos[3];
wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
for (int i=0;i<3;i++) vtxpos[i]=wcsimrootevent->GetVtx(i); //this is the source position
*/


//////////////// The way inspiured by KM's code: way too convoluted



/*
//     std::string single_file_name = tree->GetFile()->GetName();
//     TFile *file = TFile::Open(single_file_name.c_str());

//Number of events
long int nevent = ((int)tree->GetEntries());
if(verbose) printf("nevent %ld\n",nevent);

//A place holder for the root event
WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
tree->SetBranchAddress("wcsimrootevent",&wcsimrootsuperevent);

WCSimRootTrigger* wcsimrootevent;
wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);

double nPE; //the number of photon-electrons collected

// Now loop over events
for (long int ev=startEvent; ev<nevent; ev++)
{
    delete wcsimrootsuperevent;
    wcsimrootsuperevent = 0;  // EXTREMELY IMPORTANT
    // Read the event from the tree into the WCSimRootEvent instance
    tree->GetEntry(ev);
    int ncherenkovhits     = wcsimrootevent->GetNcherenkovhits(); //the total number of cherenkov hits which can induce sub events later
    int nhits;
    nhits = ncherenkovhits;

    for (int i=0; i< nhits ; i++)
    {
        //if(verbose) cout << "Hit #" << i << endl;
        WCSimRootCherenkovHit *wcsimrootcherenkovhit;
        wcsimrootcherenkovhit = (WCSimRootCherenkovHit*) (wcsimrootevent->GetCherenkovHits())->At(i);
        int peForTube      = wcsimrootcherenkovhit->GetTotalPe(1);
        nPE = peForTube;

    } // End of loop over Cherenkov hits
} // End of loop over events
}//end of main
*/
