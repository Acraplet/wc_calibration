//This stores the bin position xyz, Rtp and the charge and total number of photons
//Any map that is called by this code should be a reference map as the information will be appended to the end - reads in the Maps/txtMap and places the events in the correct bin and in the correct folder
//for a given bin we are storing all the source positions that fit  it
//of the test "oneBin" files. These are intended as a profile for how he charge collected evolves with
//absorption and scattering length and from there be a look-up table for extracting the attenuation length 
//This is a file to sum directly the positions in the same bin

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
#include <../include/findBin.hh>
#include "ColorOutput.hh"
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

    if (outfilename!=NULL){
        std::cout << WAR << "Warning, output file given : we are working without giving names to the filename, it is automatically done w.r.t the bin that the source position lies in" << std::endl;
    }

    double x, y, z, t, p, R,nHitsTot, n_entries, abwff, rayff; //the cartesian and polar coordinate of the source position, careful R is from the surface of the mPMT dome and not from its centre of the mPMT
//     int str; //this will help break down the string
    //To read the characteristics of the file
//     std::string T;
//     float theta_pos, phi_pos;
//     std::string abwff;
//     int fileID;
//     std::string rayff;
//     std::stringstream X(filename);
    std::ifstream in(filename);
    double temp;
    int count = 0;

    while ((in >> temp)) {
        //every entry in the txt map test
        // 		std::cout << "aa" << count << std::endl;
        if (count %10 == 0) x = temp;
        if (count %10 == 1) y = temp;
        if (count %10 == 2) z = temp;
        if (count %10 == 3) t = temp;
        if (count %10 == 4) p = temp;
        if (count %10 == 5) R = temp;
        if (count %10 == 6) nHitsTot = temp; //number of p.e. collected in mPMT58
        if (count %10 == 7) n_entries = temp;
        if (count %10 == 8) abwff = temp;
        if (count %10 == 9) {
            //we have finished a line
            rayff = temp;
//             std::cout << "Total number of hits : "<< nHitsTot << std::endl;

            //Here we need to choose which bin we are into for then appending the given source position to the
            //correct bin reference text file
            Bin closestBin = findBin(t, p);
            int minID = closestBin.ID;
            double theta_bin = closestBin.theta;
            double phi_bin = closestBin.phi;

            //         std::cout << x << " " << y << " " << z << " closest to bin " << minID << " at " << xbins[minID]<< " " << ybins[minID] << " "<< zbins[minID] <<std::endl;

            /*onst char* theta = t.c_str();
            const char* phi = p.c_str();
            const char* dist = R.c_str();*/

            std::string onePosition_outfile = Form("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/Maps/test_maps_oneBin/OneBin_uwu_bin%d_theta%.2f_phi%.2f_R%.2f.txt", minID, theta_bin, phi_bin, R);
//            std::cout << onePosition_outfile << std::endl;
            //carefull, we are appending to the end of the existing file,
            //if this file is not a test file it will cause problems
            std::ofstream onePosition;
            onePosition.open(onePosition_outfile, std::ofstream::app);
            onePosition <<  x << " " << y << " " << z << " " << t << " " << p << " " << R << " " << nHitsTot << " " << n_entries << " " << abwff << " " << rayff << std::endl;
        }//only at the end of the line do we add it to the correct bin test
        count+=1;

    }

    std::cout << "File " << filename << " has been added to the test files" << std::endl;

   //moving it to looking at the text file directly



}//End of main


/*This is the older version of doing the reading - dig directly in the root file -> very slow
 *
 *
 w *hile(std::getline(X, T, '_')){
 if (!T.find("FileID")){
     str = T.find("FileID");
     p = T.erase(0,str+6);
     fileID = std::stoi(p);

     std::cout<< fileID << std::endl;
     }

     if (!T.find("Absff")){
         str = T.find("Absff");
         abwff = T.erase(0,str+5);
         }
         if (!T.find("Rayff")){
             str = T.find("Rayff");
             rayff = T.erase(0,str+5);
             // Deletes str+5 characters from index number 0
             }
             if (!T.find("x")){
                 str = T.find("x");
                 x = T.erase(0,str+1);
                 //need a float counter point for making distance calculations
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
                             theta_pos = std::stof(t);
                             }
                             if (!T.find("p")){
                                 str = T.find("p");
                                 p = T.erase(0,str+1);
                                 phi_pos =  std::stof(p);
                                 }
                                 if (!T.find("R")){
                                     str = T.find("R");
                                     R = T.erase(0,str+1);
                                     }
                                     }
                                     std::cout << "Source position : " << x << " " << y << " " << z << " " << t << " " << p << " " << R << std::endl;
                                     std::cout << "Absorption : " << abwff << " Scattering : " << rayff << std::endl;

                                     TFile *infile = new TFile(filename, "READ");
                                     TTree *events = (TTree*)infile->Get("CherenkovHits");

                                     int event;
                                     int NHits;
                                     int mPMT;
                                     int mPMT_pmt;
                                     int nHitsTot = 0;

                                     events->SetBranchAddress("Event", &event);
                                     events->SetBranchAddress("NHits", &NHits);
                                     events->SetBranchAddress("mPMT", &mPMT);
                                     events->SetBranchAddress("mPMT_pmt", &mPMT_pmt);

                                     int n_entries = events->GetEntries();
                                     for(int i=0; i<n_entries; i++){
                                         events->GetEntry(i);
                                         if (mPMT == 58){
                                             if (mPMT_pmt != 19){
                                                 nHitsTot += NHits;
                                                 }
                                                 }
                                                 if (mPMT == 59){
                                                     if (mPMT_pmt == 19){
                                                         nHitsTot += NHits;
                                                         }
                                                         }
                                                         }
 *
 *
 *
 *
 */
