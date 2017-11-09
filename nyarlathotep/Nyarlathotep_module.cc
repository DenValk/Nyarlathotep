/*
 *Metrics Analysis program by Jason Stock of South Dakota School of Mines and Technology
 *September 2015.
 *
 *
 * This Version of the Nyarlathotep analysis uses the backtracker and is built for events with multiple particles. As such, its structure is different than previous analysis used before Dec 2015.
 *
 */

#ifndef Nyarlathotep_Module
#define Nyarlathotep_Module

// LArSoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h"

#include "larsim/MCCheater/PhotonBackTracker.h"
#include "larsim/MCCheater/BackTracker.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
//#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/Exception.h"
#include "nusimdata/SimulationBase/MCParticle.h"
//#include "art/Persistency/Common/PtrVector.h"
//#include "art/Utilities/Exception.h"
//#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <http://root.cern.ch/root/html532/ClassIndex.html>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"

// C++ Includes
#include <map>
#include <iomanip>
#include <sstream>
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

extern int global_event;
int global_event=0;
/////////////////////////////ANONYMOUS NAMESPACE
namespace {
  /*
  class hitContibutors{
    public:
      simb::MCParticle EVE();
      std::vector<double> EVEpos();
      std::vector<double> Pos();
      double Energy();
      double NPhotons();
      hitContributors::hitContibutors(simb::MCParticle part, std::vector<double> position, double energy, double nphotons);

    private:
      simb::MCParticle pEve;
      std::vector<double> pEVEpos;
      std::vector<double> pPos;
      double pEnergy;
      double pNPhotons;
  }
  hitContributors::hitContibutors(simb::MCParticle part, std::vector<double> position, double energy, double nphotons){
    
  }
  */


} // end local namespace





namespace Nyarlathotep {

  class Nyarlathotep : public art::EDAnalyzer
  {
  public:
 
    explicit Nyarlathotep(fhicl::ParameterSet const& parameterSet);

    virtual void beginJob();
    virtual void reconfigure(fhicl::ParameterSet const& parameterSet) override;
    virtual void analyze (const art::Event& event) override;
    virtual void endJob();
  

  private:

    int private_nEvt =0;
    int private_nPEsBins = 101;
    int private_stdBinWidth = 2; //Standard bin width in centimeters
    double private_stdBinWidthPE = 0.5; //Standard bin width for PE plots in centimeters
    double private_nPEsMin = -0.5;
    double private_nPEsMax = 100.5;

/*    double modulo (double a, double b);
    double modulo (double a, int b);
    int modulo (int a, int b);
    int signF (int a);
    int signF (double a);*/

    std::string private_trackLabel;
    std::string private_opHitLabel;
    std::string private_opFlashLabel;
    std::string private_hitLabel;
    std::string private_BTLabel;
    std::string private_PBTLabel;

    art::ServiceHandle<art::TFileService>        private_service_tfs;
    art::ServiceHandle<cheat::PhotonBackTracker> private_service_pbt;
    art::ServiceHandle<cheat::BackTracker>       private_service_bt;
    std::unique_ptr<art::TFileDirectory> generatorPlotsPtr;
    //art::TFileDirectory* generatorPlotsPtr;
    geo::GeometryCore const* private_provider_geom = lar::providerFrom<geo::Geometry>();
//    art::ServiceHandle<geo::Geometry> geom;

    Double_t private_xCryoMin, private_xCryoMax, private_yCryoMin, private_yCryoMax, private_zCryoMin, private_zCryoMax, private_xTpcMin, private_xTpcMax, private_yTpcMin, private_yTpcMax, private_zTpcMin, private_zTpcMax, private_edgeDelta;

    int private_nBinsX = 1;

    std::map<std::string, TH1D*> private_generatorHitIntegralsVX;
    std::map<std::string, TH1D*> private_generatorOpHitPEs;
    std::map<std::string, TH1D*> private_generatorHitIntegrals;
    std::map<std::string, TH1D*> private_generatorHitsPerEvt;

    TH1D* private_deltaTime;
    TH1D* private_trackLength;
    TH1D* fOpHitPEs;
    TH1D* fHitCharge;

    TH2D* fPeakAmpVXZ;
    TH2D* fIntegralVXZ;
    TH2D* fHitPEVXZ;

    //Unscalled Histograms
    TH1D* fPeakAmpVXhist;
    TH1D* fIntegralVXhist;
    TH1D* fHitPEVXhist;
    TH1D* fFlashPEsHist;

    //Scaled histograms
    TH1D* fPeakAmpScaledVXhist;
    TH1D* fIntegralScaledVXhist;
    TH1D* fHitPEScaledVXhist;

    //Calibration histograms Nov 2017
    TH1D* fEventCharge;
    TH1D* fEventLight;
    TH1D* fNHits_Charge;
    TH1D* fNHits_Light;

    TH2D* fZVT_Charge;
    TH2D* fEventChargeVLight;
    TH2D* fBinHit_Charge;
    TH2D* fBinHit_Charge_Buffer;
    TH2D* fBinHit_Light;
    TH2D* fBinHit_Light_Buffer;
    TH2D* fBinNHits_Charge;
    TH2D* fBinNHits_Light;
     
  }; // class Nyarlathotep


  
  Nyarlathotep::Nyarlathotep(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet) 
  {
    reconfigure(parameterSet);
  }

  //-----------------------------------------------------------------------
  void Nyarlathotep::beginJob()
  {
    std::cout<<"Begin Job\n";

    for(geo::CryostatID const& cID: private_provider_geom->IterateCryostatIDs()){
      Double_t* bounds = new double[6];
      private_provider_geom->CryostatBoundaries(bounds, cID);
      std::cout<< "Bounds xMax " << bounds[1] << "min " << bounds[0] <<"\n";
      private_xCryoMin = std::min(private_xCryoMin, bounds[0]);
      private_xCryoMax = std::max(private_xCryoMax, bounds[1]);
      private_yCryoMin = std::min(private_yCryoMin, bounds[2]);
      private_yCryoMax = std::max(private_yCryoMax, bounds[3]);
      private_zCryoMin = std::min(private_zCryoMin, bounds[4]);
      private_zCryoMax = std::max(private_zCryoMax, bounds[5]);
    }
    //Setting TPC boundaries
    for(geo::TPCGeo const& TPC: private_provider_geom->IterateTPCs()){
      private_xTpcMin = std::min(private_xTpcMin, TPC.MinX());
      private_xTpcMax = std::max(private_xTpcMax, TPC.MaxX());
      private_yTpcMin = std::min(private_yTpcMin, TPC.MinY());
      private_yTpcMax = std::max(private_yTpcMax, TPC.MaxY());
      private_zTpcMin = std::min(private_zTpcMin, TPC.MinZ());
      private_zTpcMax = std::max(private_zTpcMax, TPC.MaxZ());
    }
    private_nBinsX = 1;
    int private_nBinsXPE = 1;
    //int private_nBinsX = int((private_xCryoMax - private_xCryoMin)/10);
    {
      int nStdBins = int((private_xCryoMax - private_xCryoMin)/private_stdBinWidth);
      if(nStdBins<100){ private_nBinsX=100; }
      else{ private_nBinsX=nStdBins;}
      int nStdBinsPE = int((private_xCryoMax-private_xCryoMin)/private_stdBinWidthPE);
      if(nStdBinsPE<100){ private_nBinsXPE=100; }
      else{ private_nBinsXPE=nStdBinsPE;}

    }
    int private_nBinsZ = int((private_zCryoMax - private_zCryoMin)/private_stdBinWidth);
    
    //TDirectories
    
    art::TFileDirectory allParticleDir          = private_service_tfs->mkdir("allParticles", "Contains histograms with all hit information ");

    generatorPlotsPtr = 
      std::make_unique<art::TFileDirectory>
      (allParticleDir.mkdir("generators", 
                            "Contains histograms for each " 
                            "generator process"));
    //generatorPlotsPtr = new art::TFileDirectory(allParticleDir.mkdir("generators", "Contains histograms for each generator process"));

    art::TFileDirectory scaledAllParticlesDir   = allParticleDir.mkdir("scale_AllParticles", "Contains scaled histograms with all hit information ");

    art::TFileDirectory unscaledAllParticlesDir = allParticleDir.mkdir("unscaled_AllParticles", "Contains histograms with all hit information ");


    /*
    art::TFileDirectory scaled   =       private_service_tfs->mkdir("scaledHistograms", "Contains the histograms rescaled by energy deposited.");
    art::TFileDirectory unscaled =       private_service_tfs->mkdir("unscaledHistograms", "Contains the histograms without rescaling by energy deposited.");
    art::TFileDirectory unscaled2D =     private_service_tfs->mkdir("unscaled2D", "Contains the 2D XZ histograms without rescaling by energy deposited.");
    */


    private_deltaTime     =           private_service_tfs->make<TH1D>    ("histdeltaTime", 
        "Plot of Times between hits. Restricted to times > 500ns", 4500, 500, 5000);
    private_trackLength   =           allParticleDir.make<TH1D>("hist_trackLengths",
        "Histogram of Track Lengths",1000, 0.0, 100.0  );
    fOpHitPEs             =           allParticleDir.make<TH1D>("hist_OpHitPEs",
        "Histogram of PEs from each OpHit.",21, -0.5, 20.5  );
    fHitCharge            =           allParticleDir.make<TH1D>("hist_HitCharge",
        "Histogram of #int ADC Units from each Hit",5000, 0.0, 1000.0  );
    fPeakAmpVXZ           =           allParticleDir.make<TH2D>("hist_PeakAmpVXZ",
        "Histogram of HitPeakAmplitude vs Z Position (cm) and X Position (cm)",
        private_nBinsX, private_xCryoMin, private_xCryoMax, private_nBinsZ,private_zCryoMin,private_zCryoMax );
    fPeakAmpVXZ -> SetYTitle("Z Position (cm)");
    fPeakAmpVXZ -> SetXTitle("X Position (cm)");

    fIntegralVXZ          =           allParticleDir.make<TH2D>("hist_IntegralVXZ",
        "Histogram of HitIntegral vs Z Position (cm) and X Position (cm)",
        private_nBinsXPE, private_xCryoMin, private_xCryoMax,        private_nBinsZ, private_zCryoMin, private_zCryoMax        );
    fIntegralVXZ-> SetYTitle("Z Position (cm)");
    fIntegralVXZ-> SetZTitle("X Position (cm)");

    fHitPEVXZ             =           allParticleDir.make<TH2D>("hist_HitPEVXZ",
        "Histogram of PhotoElectrons Detected vs Z Position (cm) and X Position (cm)",
        private_nBinsXPE, private_xCryoMin, private_xCryoMax,        private_nBinsZ, private_zCryoMin, private_zCryoMax        );
    fHitPEVXZ-> SetYTitle("Z Position (cm)");
    fHitPEVXZ-> SetXTitle("X Position (cm)");

    //make scaled histograms
    fPeakAmpScaledVXhist         =    scaledAllParticlesDir.make<TH1D>("hist_PeakAmpScaledVXhist", 
        "Plot of Summed Hit Peak Amplitudes vs X position (Distance from APA in cm)", private_nBinsX, private_xCryoMin, private_xCryoMax);
    fIntegralScaledVXhist        =    scaledAllParticlesDir.make<TH1D>("hist_IntegralScaledVXhist", 
        "Plot of Summed Hit Integral Values vs X position (Distance from APA in cm)", private_nBinsX, private_xCryoMin, private_xCryoMax);
    fHitPEScaledVXhist           =    scaledAllParticlesDir.make<TH1D>("hist_HitPEScaledVXhist", 
        "Plot of the number of Detected Photo Electrons vs X position (Distance from APA in cm)", private_nBinsXPE, private_xCryoMin, private_xCryoMax);

    

    //Make unscaled Histograms
    fPeakAmpVXhist               =    unscaledAllParticlesDir.make<TH1D>("hist_PeakAmpVXhist", 
        "Plot of Summed Hit Peak Amplitudes vs X position (Distance from APA in cm)", private_nBinsX, private_xCryoMin, private_xCryoMax);
    fIntegralVXhist              =    unscaledAllParticlesDir.make<TH1D>("hist_IntegralVXhist", 
        "Plot of Summed Hit Integral Values vs X position (Distance from APA in cm)", private_nBinsX, private_xCryoMin, private_xCryoMax);
    fHitPEVXhist                 =    unscaledAllParticlesDir.make<TH1D>("hist_HitPEVXhist", 
        "Plot of the number of Detected Photo Electrons vs X position (Distance from APA in cm)", private_nBinsX, private_xCryoMin, private_xCryoMax);
    fFlashPEsHist                =    unscaledAllParticlesDir.make<TH1D>("hist_FlashPEsHist", 
        "Plot of the number of Detected Photo Electrons per Flash", private_nPEsBins, private_nPEsMin, private_nPEsMax);


    //Make Calibration histograms (2017).
    art::TFileDirectory calibDir  = allParticleDir.mkdir("CalibrationTestAnalysis","A directory for as of yet unclassified analysis code added for the Nov 2017 Calibration meetings.");
    fEventCharge                  = calibDir.make<TH1D>("hist_EventCharge", "Total Charge per event.", 5001, -0.5, 5000.5);
    fEventLight                   = calibDir.make<TH1D>("hist_EventLight", "Total PEs per event.", 501,-0.05,50.05);
    fNHits_Charge                 = calibDir.make<TH1D>("hist_NHits_Charge", "The number of charge hits in each event.", 101,-0.5,100.5);
    fNHits_Light                  = calibDir.make<TH1D>("hist_NHits_Light", "The number of optical hits in each event.", 41,-0.5,40.5);
    fZVT_Charge                   = calibDir.make<TH2D>("hist_ChargeVsZT",  "The charge collected vs Z and T", private_nBinsZ, private_zCryoMin, private_zCryoMax, 1000,0,5000);
    //Note, no light here. This is because the light should all be at 0 with respect to the TPC.
    fEventChargeVLight            = calibDir.make<TH2D>("hist_ChargeVLight", "The charge collected vs the light collected", 5001, -0.05, 5000.5, 501, -0.05, 50.05);
    fBinHit_Charge                = calibDir.make<TH2D>("hist_BinHit_Charge", "Each bin is filled only once per event when a hit is registered in that bin.", 
        private_nBinsZ, private_zCryoMin, private_zCryoMax, private_nBinsX, private_xCryoMin, private_xCryoMax);
    fBinHit_Light                 = calibDir.make<TH2D>("hist_BinHit_Light", "Each bin is filled only once per event when an OpHit is registered in that bin.",
        private_nBinsZ, private_zCryoMin, private_zCryoMax, private_nBinsX, private_xCryoMin, private_xCryoMax);
    fBinNHits_Charge              = calibDir.make<TH2D>("hist_nBinHit_Charge", "Counts the number of hits in each bin (no charge weighting).",
        private_nBinsZ, private_zCryoMin, private_zCryoMax, private_nBinsX, private_xCryoMin, private_xCryoMax);
    fBinNHits_Light               = calibDir.make<TH2D>("hist_nBinHit_Light", "Counts the number of OpHits in each bin (no PE weighting).",
        private_nBinsZ, private_zCryoMin, private_zCryoMax, private_nBinsX, private_xCryoMin, private_xCryoMax);

    //Buffers
    fBinHit_Charge_Buffer = new TH2D("hist_BinHit_Charge_Buffer", "Each bin is filled only once per event when a hit is registered in that bin.", 
        private_nBinsZ, private_zCryoMin, private_zCryoMax, private_nBinsX, private_xCryoMin, private_xCryoMax);
    fBinHit_Light_Buffer = new TH2D("hist_BinHit_Light_Buffer", "Each bin is filled only once per event when an OpHit is registered in that bin.", 
        private_nBinsZ, private_zCryoMin, private_zCryoMax, private_nBinsX, private_xCryoMin, private_xCryoMax);

  }
  

   //-----------------------------------------------------------------------
  void Nyarlathotep::endJob()
  {
    double scaleFrac = (1.0 / ( double(private_nEvt) ) );
    std::cout<<"Flow Check. End Job\n Scaling by "<<scaleFrac<<"\n" ; //DELETE THIS LINE

    fFlashPEsHist          -> Scale( scaleFrac );
    private_trackLength    -> Scale( scaleFrac );
    fOpHitPEs              -> Scale( scaleFrac );
    fHitCharge             -> Scale( scaleFrac );
    fPeakAmpVXZ            -> Scale( scaleFrac );
    fIntegralVXZ           -> Scale( scaleFrac );
    fHitPEVXZ              -> Scale( scaleFrac );
 
    //Unscalled Histograms
    fPeakAmpVXhist         -> Scale( scaleFrac );
    fIntegralVXhist        -> Scale( scaleFrac );
    fHitPEVXhist           -> Scale( scaleFrac );
    fFlashPEsHist          -> Scale( scaleFrac );

    //Scaled histograms
    fPeakAmpScaledVXhist     -> Scale( scaleFrac );
    fIntegralScaledVXhist    -> Scale( scaleFrac );
    fHitPEScaledVXhist       -> Scale( scaleFrac );

    //TDir Hist Pointer Containers

    for(std::map<std::string, TH1D*>::iterator itr=private_generatorHitIntegralsVX.begin(); itr!=private_generatorHitIntegralsVX.end(); ++itr){
      (itr->second)->Scale( scaleFrac );
    }
    for(std::map<std::string, TH1D*>::iterator itr=private_generatorHitIntegrals.begin(); itr!=private_generatorHitIntegrals.end(); ++itr){
      (itr->second)->Scale( scaleFrac );
    }
    for(std::map<std::string, TH1D*>::iterator itr=private_generatorHitsPerEvt.begin(); itr!=private_generatorHitsPerEvt.end(); ++itr){
      (itr->second)->Scale( scaleFrac );
    }
    for(std::map<std::string, TH1D*>::iterator itr=private_generatorOpHitPEs.begin(); itr!=private_generatorOpHitPEs.end(); ++itr){
      (itr->second)->Scale( scaleFrac );
    }
  }
  //-----------------------------------------------------------------------
  void Nyarlathotep::reconfigure(fhicl::ParameterSet const& parameterSet)
  {
    // Read parameters from the .fcl file. The names in the arguments
    std::cout<<"reconfigure\n";
    private_trackLabel           = parameterSet.get< std::string >("TrackLabel", "pmtracktc");
    private_hitLabel             = parameterSet.get< std::string >("HitLabel");
    private_opHitLabel           = parameterSet.get< std::string >("PhotLabel");
    private_opFlashLabel         = parameterSet.get< std::string >("FlashLabel", "opflash");
    private_BTLabel              = parameterSet.get< std::string >("BackTrackerLabel");
    private_PBTLabel             = parameterSet.get< std::string >("PhotonBackTrackerLabel");
     
     
     
    //Setting Global Boundaries
    private_xCryoMin=0.0;
    private_xCryoMax=0.0;
    private_yCryoMin=0.0;
    private_yCryoMax=0.0;
    private_zCryoMin=0.0;
    private_zCryoMax=0.0;
    private_xTpcMin=0.0;
    private_xTpcMax=0.0;
    private_yTpcMin=0.0;
    private_yTpcMax=0.0;
    private_zTpcMin=0.0;
    private_zTpcMax=0.0;
    private_edgeDelta = 5.0;
  }

  //-----------------------------------------------------------------------
  void Nyarlathotep::analyze(const art::Event& evt) 
  {
    private_nEvt++;
    global_event = evt.id().event();
    art::Handle< std::vector< recob::Hit > > hitHandle;
    std::vector< art::Ptr< recob::Hit > > hitList;
    if (evt.getByLabel(private_hitLabel, hitHandle) )
      art::fill_ptr_vector(hitList, hitHandle);

    art::Handle< std::vector< recob::OpHit > > opHitHandle;
    std::vector< art::Ptr< recob::OpHit > > opHitList;
    if (evt.getByLabel(private_opHitLabel, opHitHandle) )
      art::fill_ptr_vector(opHitList, opHitHandle);
    
    art::Handle< std::vector< recob::OpFlash > > opFlashHandle;
    std::vector< art::Ptr< recob::OpFlash > > opFlashList;
    if (evt.getByLabel(private_opFlashLabel, opFlashHandle) ){
      art::fill_ptr_vector(opFlashList, opFlashHandle);
    }

    art::Handle< std::vector< recob::Track> > trackHandle;
    std::vector< art::Ptr< recob::Track > >   trackList;
    if(evt.getByLabel(private_trackLabel, trackHandle) ){
      art::fill_ptr_vector(trackList, trackHandle);
    }

    std::vector< art::Handle< std::vector< simb::MCTruth > > > mcHandles;
    evt.getManyByType(mcHandles);
    std::map< long int, std::string> tIdToLabel;
    
    for( auto const& mcHandle : mcHandles ){
      const std::string& sModuleLabel = mcHandle.provenance()->moduleLabel();
      //I have the label here. I can use that to mimic DAQSimAna
      auto mcTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(sModuleLabel);
      art::FindManyP<simb::MCParticle> findMCParts(mcTrue, evt, "largeant" );
      std::vector<art::Ptr<simb::MCParticle > > mcParts = findMCParts.at(0);
      for(const art::Ptr<simb::MCParticle > ptr : mcParts){
        simb::MCParticle part = *ptr;
        int track = ptr->TrackId();
        //tIdToLabel[part] = sModuleLabel;
        tIdToLabel[track] = sModuleLabel;
      }
      std::string private_intPlotTitle = "Plot of Integrated Charge collected from generator ";
      std::string private_hitPlotTitle = "Plot of Number of Hits per Event collected from generator ";
      std::string private_pePlotTitle = "Plot of PEs collected from generator ";
      private_intPlotTitle.append(sModuleLabel);
      std::string intVxName = "IntVX"; std::string intName = "Int"; std::string opName = "PE"; std::string hitsName="Hits";
      if(private_generatorHitIntegralsVX.find(sModuleLabel) == private_generatorHitIntegralsVX.end()){
        TH1D* tmpHist = generatorPlotsPtr->make<TH1D>((intVxName.append(sModuleLabel)).c_str(),  (private_intPlotTitle.append(" Vs X")).c_str(), private_nBinsX, private_xCryoMin, private_xCryoMax);
        private_generatorHitIntegralsVX[sModuleLabel] = std::move(tmpHist);
      }
      if(private_generatorHitIntegrals.find(sModuleLabel) == private_generatorHitIntegrals.end()){
        TH1D* tmpHist = generatorPlotsPtr->make<TH1D>((intName.append(sModuleLabel)).c_str(),  private_intPlotTitle.c_str(), 5000, 0.0, 1000.0);
        private_generatorHitIntegrals[sModuleLabel] = std::move(tmpHist);
      }
      if(private_generatorHitsPerEvt.find(sModuleLabel) == private_generatorHitsPerEvt.end()){
        TH1D* tmpHist = generatorPlotsPtr->make<TH1D>((hitsName.append(sModuleLabel)).c_str(),  private_hitPlotTitle.c_str(), 100000, -0.5, 99999.5);
        private_generatorHitsPerEvt[sModuleLabel] = std::move(tmpHist);
      }
      if(private_generatorOpHitPEs.find(sModuleLabel) == private_generatorOpHitPEs.end()){
        TH1D* tmpHist = generatorPlotsPtr->make<TH1D>((opName.append(sModuleLabel)).c_str(),  private_pePlotTitle.c_str(), 100, 0.0, 50.0);
        private_generatorOpHitPEs[sModuleLabel] = std::move(tmpHist);
      }
    }

    //Call on Reco Tracks
    for( auto track : trackList ){
//      if(track->Length()!=0){std::cout<<"Track Length: "<<(track->Length())<<"\n";}
      private_trackLength-> Fill(track->Length());
    }


    //Call Charge events.
    std::map<std::string, int> nHitsPerGen;
    double chargeEventTotal =0.0;
    int nHitsEventTotal=0;
    for( auto ptrHit: hitList){
      try{
        recob::Hit hit = *ptrHit;
        double hitInt = hit.Integral();
        double hitAmp = hit.PeakAmplitude();
        double scaleTot = 0.;

        std::vector<double> xyzPos;
        xyzPos = (private_service_bt->HitToXYZ(ptrHit));
        //Fill Charge Buffer (Note we use set bin content instead of Fill, since we want each bin to be BOOLE.
        fBinHit_Charge_Buffer->SetBinContent(fBinHit_Charge_Buffer->FindBin(xyzPos[2],xyzPos[0]), 1.0);
        fBinNHits_Charge->Fill(xyzPos[2],xyzPos[0],1.0);
        fZVT_Charge->Fill(xyzPos[2],hit.PeakTime(),hitInt);

        chargeEventTotal+=hit.Integral();
        nHitsEventTotal++;

        std::vector<sim::TrackIDE>  eveIDEs = private_service_bt->HitToEveID(ptrHit);
        for(const sim::TrackIDE eveIDE : eveIDEs){
          auto tID = eveIDE.trackID;
          const simb::MCParticle* particle=private_service_bt->TrackIDToParticle(tID);
          if(particle){
            std::string generatorName=tIdToLabel[tID];
            if(nHitsPerGen.find(generatorName) == nHitsPerGen.end()){
              nHitsPerGen.emplace(generatorName, 1);
            }else{
              auto tmp = nHitsPerGen.find(generatorName);
              tmp->second = (tmp->second) + 1;
            }
            //Look up the Label of the generator of this particle
            //std::string generatorName = tIdToLabel[*particle];
            TH1D* genIntVxHist = private_generatorHitIntegralsVX[generatorName];
            TH1D* genIntHist   = private_generatorHitIntegrals[generatorName];
            double scaleFrac   = eveIDE.energyFrac;
            double scaleEn     = eveIDE.energy;
            scaleTot += scaleEn;
            std::vector<double> const evePos = { particle->Trajectory().X(0),  particle->Trajectory().Y(0), particle->Trajectory().Z(0) };
    

            genIntVxHist-> Fill(xyzPos.at(0), hitInt*scaleFrac);
            //Review. Do I want scale frac or not on genIntHist... 
            genIntHist-> Fill(hitInt*scaleFrac);
          }//End If Particle
        }
        if(scaleTot<1.0e-10){
          scaleTot=1.0;
        }
        //Write globals

        fHitCharge            -> Fill(hitInt);
        fPeakAmpVXhist        -> Fill(xyzPos.at(0), hitAmp);
        fIntegralVXhist       -> Fill(xyzPos.at(0), hitInt);
        fPeakAmpScaledVXhist  -> Fill(xyzPos.at(0), hitAmp/scaleTot);
        fIntegralScaledVXhist -> Fill(xyzPos.at(0), hitInt/scaleTot);

        fPeakAmpVXZ-> Fill(xyzPos.at(0), xyzPos.at(2), hitAmp);
        fIntegralVXZ-> Fill(xyzPos.at(0), xyzPos.at(2), hitInt);
      }catch(...){

      }//end trycatch
    }//end for(hits)
    //Loop over hits counted and add them to the correct histograms.
    for(std::map<std::string, int>::iterator itr=nHitsPerGen.begin(); itr!=nHitsPerGen.end(); ++itr){
      TH1D* genHitHist  = private_generatorHitsPerEvt[itr->first];
      genHitHist->Fill(itr->second);
    }

    //Add fBinHit_Charge_Buffer to fBinHitCharge.
//    if(fBinHit_Charge->GetSumw2N() == 0 ){fBinHit_Charge->Sumw2(1);}
//    if(fBinHit_Charge_Buffer->GetSumw2N() == 0 ){fBinHit_Charge_Buffer->Sumw2(1);}

    fBinHit_Charge->Add(fBinHit_Charge, fBinHit_Charge_Buffer, 1.0, 1.0);
    fBinHit_Charge_Buffer->Reset();
    //delete fBinHit_Charge_Buffer;

    fEventCharge->Fill(chargeEventTotal);
    fNHits_Charge->Fill(double(nHitsEventTotal));


    //Call Op Events

    int nOpHitsEventTotal=0;
    double lightEventTotal=0.0;

    for(const art::Ptr<recob::OpHit>& ptrOpHit: opHitList){
      try{
        recob::OpHit opHit = *ptrOpHit;
        double nPE=opHit.PE();
        double scaleTot=0.;
        lightEventTotal+=nPE;
        nOpHitsEventTotal++;
        //need position,
        std::vector<double> xyzPos          = private_service_pbt  -> OpHitToXYZ(ptrOpHit);
        std::vector<sim::TrackSDP>  eveSDPs = private_service_pbt  -> OpHitToEveSDPs(ptrOpHit);

        //Fill Charge Buffer (Note we use set bin content instead of Fill, since we want each bin to be BOOLE.
        fBinHit_Light_Buffer->SetBinContent(fBinHit_Light_Buffer->FindBin(xyzPos[2],xyzPos[0]), 1.0);
        fBinNHits_Light->Fill(xyzPos[2],xyzPos[0],1.0);

        for(const sim::TrackSDP eveSDP : eveSDPs){
          auto tID = eveSDP.trackID;
          const simb::MCParticle* particle = private_service_pbt->TrackIDToParticle(eveSDP.trackID);
          if(particle){
            std::string generatorName = tIdToLabel[ tID ];
            TH1D* genPEHist = private_generatorOpHitPEs[generatorName];

            double scaleFrac = eveSDP.energyFrac;
            double scaleEn   = eveSDP.energy;
            scaleTot += scaleEn;
            std::vector<double> const evePos = { particle->Trajectory().X(0),  particle->Trajectory().Y(0), particle->Trajectory().Z(0) };
  
            //Review. Do I want scale frac or not on genPEHist
            genPEHist      -> Fill(nPE*scaleFrac);
          }//End if particle
        }//End for (eveSDP)
        if(scaleTot<1.0e-10){
          scaleTot=1.0;
        }
        //Write globals
        fOpHitPEs          -> Fill(nPE);
        fHitPEScaledVXhist -> Fill(xyzPos.at(0), nPE/scaleTot);
        fHitPEVXhist       -> Fill(xyzPos.at(0), nPE);
        fHitPEVXZ          -> Fill(xyzPos.at(0), xyzPos.at(2), nPE);

      }//end try
      catch(...){//Lazy catch. Need to catch -No sim SDPs found- error.
   //     throw;
      }
    }//End for(OpHit)

    //if( fBinHit_Light->GetSumw2N() == 0 ){fBinHit_Light->Sumw2();}
    //if( fBinHit_Light_Buffer->GetSumw2N() == 0 ){fBinHit_Light_Buffer->Sumw2();}
    fBinHit_Light->Add(fBinHit_Light, fBinHit_Light_Buffer, 1.0, 1.0);
    fBinHit_Light_Buffer->Reset();
    fEventLight->Fill(lightEventTotal);
    fNHits_Light->Fill(double(nOpHitsEventTotal));


    fEventChargeVLight->Fill(chargeEventTotal,lightEventTotal,1.0);

    //Call OpFlashes
    for( auto ptrOpFlash : opFlashList ){
      double nPEcount = ptrOpFlash->TotalPE();
      std::vector<double> nPEs = ptrOpFlash->PEs();
    //  for( double nPE : nPEs ){
    //    if( abs(pid)==13 ){
    //      muon_FlashPEsHist           -> Fill(nPEcount);
    //    }else if( abs(pid)==11 ){
    //      beta_FlashPEsHist           -> Fill(nPEcount);
    //    }else if( abs(pid)==22){
    //      gamma_FlashPEsHist          -> Fill(nPEcount);
    //    }else if( pid==1000020040 ){
    //      alpha_FlashPEsHist          -> Fill(nPEcount);
    //    }else if ( abs(pid) == 2212 ){
    //      neutron_FlashPEsHist        -> Fill(nPEcount);
    //    }
    //  }
      fFlashPEsHist-> Fill( nPEcount );
    }

    

  } // Nyarlathotep::analyze()
  
  
  // This macro has to be defined for this module to be invoked from a
  // .fcl file; see Nyarlathotep.fcl for more information.
  DEFINE_ART_MODULE(Nyarlathotep)



} // namespace Nyarlathotep

#endif // Nyarlathotep_Module

