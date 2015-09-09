#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSParameters.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"

#include "HGCSSCaloProperties.hh"

int main(int argc, char** argv){//main  

  if (argc < 7) {
    std::cout << " Usage: "
              << argv[0] << " <nEvts to process (0=all)>"
              << " <path to input files>"
              << " <name of input sim file>"
              << " <name of input reco file>"
              << " <full path to output file>"
              << " <number of si layers to consider: 1,2 or 3>"  
              << " <optional: debug (default=0)>"
              << std::endl;
    return 1;
  } 

 //////////////////////////////////////////////////////////
 //  //// Hardcoded config ////////////////////////////////////
 //    //////////////////////////////////////////////////////////
 //for HGCAL, true means only 12 FHCAL layers considered (24 are simulated)
 
 
  bool concept = true;

  bool selectEarlyDecays = true;

  double minX=-1700,maxX=1700;
  double minY=-1700,maxY=1700;
  double minZ=3170,maxZ=5070;

  unsigned nX=(maxX-minX)/10,nY=(maxY-minY)/10;
  unsigned nZ=maxZ-minZ;

  double FHcalEMCalib = 118;//40.4;//39.81;//38;
  double FHcalEMOffset = -209;//-3.9;//1.9;//-15;
  double BHcalEMCalib = 9.92;//40.4;//39.81;//38;
  double BHcalEMOffset = -5.1;//1.9;//-15;
  double HcalPionCalib = 0.92;//1/1.19;//0.901;//1./0.9;//1/0.846;
  double HcalPionOffset = 0;//-0.81;
  double BHcalSlope = 2.7;
  double G4BHcalSlope = 0.24;
  // choose a jet definition
  //   //double R = 0.5;
  //     //JetDefinition jet_def(antikt_algorithm, R);
  //////////////////////////////////////////////////////////
  //  //// End Hardcoded config ////////////////////////////////////
  //    //////////////////////////////////////////////////////////
  
  const unsigned pNevts = atoi(argv[1]);
  std::string filePath = argv[2];
  std::string simFileName = argv[3];
  std::string recoFileName = argv[4];

  std::string inFilePath = filePath+simFileName;

  std::string outPath = argv[5];
  unsigned nSiLayers = 2;
  nSiLayers = atoi(argv[6]);

  unsigned debug = 0;
  if (argc >7) debug = atoi(argv[7]);

  unsigned genEn;
  size_t end=outPath.find_last_of(".root");
  size_t start=outPath.find_last_of("e");
  std::istringstream(outPath.substr(start+1,end))>>genEn;

  bool isEM = false;

  if (inFilePath.find("e-")!=inFilePath.npos ||
      inFilePath.find("e+")!=inFilePath.npos) isEM = true;

  if (selectEarlyDecays && isEM) {
    selectEarlyDecays = false;
    HcalPionCalib = 1;
    HcalPionOffset = 0;
  }

  std::cout << " -- Input parameters: " << std::endl
            << " -- Input file path: " << filePath << std::endl
            << " -- Output file path: " << outPath << std::endl
            << " -- Generated energy: " << genEn << std::endl
            << " -- Requiring " << nSiLayers << " si layers." << std::endl
            << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;

  TRandom3 lRndm(1);
  std::cout << " -- Random number seed: " << lRndm.GetSeed() << std::endl;

  /////////////////////////////////////////////////////////////
  //  //input
  //    /////////////////////////////////////////////////////////////

  std::ostringstream input;
  input << filePath << "/" << simFileName;

  TFile *simFile = TFile::Open(input.str().c_str());

  if (!simFile) {
    std::cout << " -- Error, input file " << input.str() << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  else std::cout << " -- input file " << simFile->GetName() << " successfully opened." << std::endl;

  TTree *lSimTree = (TTree*)simFile->Get("HGCSSTree");
  if (!lSimTree){
    std::cout << " -- Error, tree HGCSSTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

   input.str("");
  input << filePath << "/" << recoFileName;

  TFile *recFile = TFile::Open(input.str().c_str());

  if (!recFile) {
    std::cout << " -- Error, input file " << input.str() << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  else std::cout << " -- input file " << recFile->GetName() << " successfully opened." << std::endl;

  TTree *lRecTree = (TTree*)recFile->Get("RecoTree");
  if (!lRecTree){
    std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  /////////////////////////////////////////////////////////////
  //  //Info
  //    /////////////////////////////////////////////////////////////
  
  HGCSSInfo * info=(HGCSSInfo*)simFile->Get("Info");
  const double cellSize = info->cellSize();
  const unsigned versionNumber = info->version();
  const unsigned model = info->model();

  //models 0,1 or 3.
  bool isTBsetup = (model != 2);
  bool isCaliceHcal = versionNumber==23;//inFilePath.find("version23")!=inFilePath.npos || inFilePath.find("version_23")!=inFilePath.npos;
  
  //extract input energy
  
  std::cout << " -- Version number is : " << versionNumber
            << ", model = " << model
            << ", cellSize = " << cellSize
            << std::endl;

  //initialise detector
  HGCSSDetector & myDetector = theDetector();
  
  myDetector.buildDetector(versionNumber,concept,isCaliceHcal);

  //initialise calibration class
  HGCSSCalibration mycalib(inFilePath);
  HGCSSDigitisation myDigitiser;
  myDigitiser.setRandomSeed(lRndm.GetSeed());

  const unsigned nLayers = myDetector.nLayers();
  const unsigned nSections = myDetector.nSections();

  std::cout << " -- N layers = " << nLayers << std::endl
            << " -- N sections = " << nSections << std::endl;

  HGCSSGeometryConversion geomConv5(inFilePath,model,cellSize);
     //assemble in 5*5 and 10*10 to fill maxE
   std::vector<unsigned> granularity5;
   granularity5.resize(nLayers,2);
   geomConv5.setGranularity(granularity5);
   geomConv5.initialiseHistos(false,"_5");
   HGCSSGeometryConversion geomConv10(inFilePath,model,cellSize);
   std::vector<unsigned> granularity10;
   granularity10.resize(nLayers,4);
   geomConv10.setGranularity(granularity10);
   geomConv10.initialiseHistos(false,"_10");
   HGCSSGeometryConversion geomConv15(inFilePath,model,cellSize);
   std::vector<unsigned> granularity15;
   granularity15.resize(nLayers,6);
   geomConv15.setGranularity(granularity15);
   geomConv15.initialiseHistos(false,"_15");
   HGCSSGeometryConversion geomConv2d5(inFilePath,model,cellSize);
   std::vector<unsigned> granularity2d5;
   granularity2d5.resize(nLayers,1);
   geomConv2d5.setGranularity(granularity2d5);
   geomConv2d5.initialiseHistos(false,"_2d5");

   //TFile *outputFile = TFile::Open(outPath.c_str(),"RECREATE");
   std::string histfilename=(outPath).c_str();
   TFile *outputFile = new TFile(histfilename.c_str(), "RECREATE");

  if (!outputFile) {
    std::cout << " -- Error, output file " << outPath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }

  std::cout << " -- 2-D histograms: " << std::endl
            << " -- X: " << nX << " " << minX << " " << maxX << std::endl
            << " -- Y: " << nY << " " << minY << " " << maxY << std::endl
            << " -- Z: " << nZ << " " << minZ << " " << maxZ << std::endl
    ;
  outputFile->cd();


  TH2F *p_ShowervsLayer = new TH2F("p_ShowervsLayer",";layer ; Shower",nLayers,0,nLayers,1000,0,5000); 
   
  HGCSSEvent * event = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  std::vector<HGCSSRecoHit> * rechitvec = 0;
  std::vector<HGCSSGenParticle> * genvec = 0;
 

  CaloProperties *caloProperties = new CaloProperties("V04");
 
  const unsigned nEvts = ((pNevts > lSimTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lSimTree->GetEntries()) : pNevts) ;

  std::vector<Float_t> eList={16};
  caloProperties->setEnergiesToScan(eList);
  caloProperties->characterizeCalo(); 
  //TGraph *sapta = (TGraph*)caloProperties->gr_showerMax;
  //sapta->Write("sapta");
  //for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

     
     //std::cout  << caloProperties->characterizeCalo() << std::endl;


  //}


  outputFile->Write();
  return 0;


}
