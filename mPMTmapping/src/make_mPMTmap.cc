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
#include "WCSimRootEvent.hh"


//This is a work in progress - the python file needs to be translated into cpp
//Need to:
/*
 - open all the data with a given configurations
 - read the output giving out: R, theta, phi, (source x,y,z + dir x,y,z)
 - then compile everything into a 2D map
*/
const std::string ERR = color::RED_STR + "[ERROR]: " + color::RESET_STR;

void HelpMessage()
{
    std::cout << "USAGE: "
    << "WCSIM_TreeConvert" << "\nOPTIONS:\n"
    << "-f : Input file\n" << std::endl;
}

int main(int argc, char **argv){
    bool verbose=false;
    char * filename=NULL;
    long int endEvent=0;
    char c;

    while( (c = getopt(argc,argv,"f:o:b:s:e:l:r:p:c:w:hdtv")) != -1 ){//input in c the argument (-f etc...) and in optarg the next argument. When the above test becomes -1, it means it fails to find a new argument.
        switch(c){
            case 'f':
                filename = optarg;
                break;
        }
    }

    if (filename==NULL){
        std::cout << ERR << "Error, no input file: " << std::endl;
        HelpMessage();
        return -1;
    }

    TChain *tree = new TChain("wcsimT");
    tree->Add(filename);

    std::string single_file_name = tree->GetFile()->GetName();
    TFile *file = TFile::Open(single_file_name.c_str());

    //Number of events
    long int nevent = ((int)tree->GetEntries());
    if(verbose) printf("nevent %ld\n",nevent);

    //A place holder for the root event
    WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
    tree->SetBranchAddress("wcsimrootevent",&wcsimrootsuperevent)

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
            nPE = peForTube

        } // End of loop over Cherenkov hits
    } // End of loop over events
}//end of main



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

