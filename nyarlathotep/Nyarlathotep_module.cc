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
//#include "messagefacility/Utilities/eception.h"

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

    TH2D* fAlpha_PeakAmpVXZ;
    TH2D* fAlpha_IntegralVXZ;
    TH2D* fAlpha_HitPEVXZ;

    TH2D* fBeta_PeakAmpVXZ;
    TH2D* fBeta_IntegralVXZ;
    TH2D* fBeta_HitPEVXZ;

    TH2D* fGamma_PeakAmpVXZ;
    TH2D* fGamma_IntegralVXZ;
    TH2D* fGamma_HitPEVXZ;

    TH2D* fNeutron_PeakAmpVXZ;
    TH2D* fNeutron_IntegralVXZ;
    TH2D* fNeutron_HitPEVXZ;

    //Unscalled Histograms
    TH1D* fPeakAmpVXhist;
    TH1D* fIntegralVXhist;
    TH1D* fHitPEVXhist;
    TH1D* fFlashPEsHist;
    TH1D* muon_PeakAmpVXhist;
    TH1D* muon_IntegralVXhist;
    TH1D* muon_HitPEVXhist;
//    TH1D* muon_FlashPEsHist;
    TH1D* alpha_PeakAmpVXhist;
    TH1D* alpha_IntegralVXhist;
    TH1D* alpha_HitPEVXhist;
//    TH1D* alpha_FlashPEsHist;
    TH1D* beta_PeakAmpVXhist;
    TH1D* beta_IntegralVXhist;
    TH1D* beta_HitPEVXhist;
//    TH1D* beta_FlashPEsHist;
    TH1D* gamma_PeakAmpVXhist;
    TH1D* gamma_IntegralVXhist;
    TH1D* gamma_HitPEVXhist;
//    TH1D* gamma_FlashPEsHist;
    TH1D* neutron_PeakAmpVXhist;
    TH1D* neutron_IntegralVXhist;
    TH1D* neutron_HitPEVXhist;
//    TH1D* neutron_FlashPEsHist;

    //Scaled histograms
    TH1D* fPeakAmpScaledVXhist;
    TH1D* fIntegralScaledVXhist;
    TH1D* fHitPEScaledVXhist;
    TH1D* muon_PeakAmpScaledVXhist;
    TH1D* muon_IntegralScaledVXhist;
    TH1D* muon_HitPEScaledVXhist;
    TH1D* alpha_PeakAmpScaledVXhist;
    TH1D* alpha_IntegralScaledVXhist;
    TH1D* alpha_HitPEScaledVXhist;
    TH1D* beta_PeakAmpScaledVXhist;
    TH1D* beta_IntegralScaledVXhist;
    TH1D* beta_HitPEScaledVXhist;
    TH1D* gamma_PeakAmpScaledVXhist;
    TH1D* gamma_IntegralScaledVXhist;
    TH1D* gamma_HitPEScaledVXhist;
    TH1D* neutron_PeakAmpScaledVXhist;
    TH1D* neutron_IntegralScaledVXhist;
    TH1D* neutron_HitPEScaledVXhist;

  



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
    art::TFileDirectory muonDir                 = private_service_tfs->mkdir("muons"       , "Contains histograms of muon hit information. ");
    art::TFileDirectory alphaDir                = private_service_tfs->mkdir("alphas"      , "Contains histograms with of alpha hit information. ");
    art::TFileDirectory betaDir                 = private_service_tfs->mkdir("betas"       , "Contains histograms of beta hit information. ");
    art::TFileDirectory gammaDir                = private_service_tfs->mkdir("gammas"      , "Contains histograms of gamma hit information. ");
    art::TFileDirectory neutronDir              = private_service_tfs->mkdir("neutrons"    , "Contains histograms of neutron hit information");

    generatorPlotsPtr = 
      std::make_unique<art::TFileDirectory>
      (allParticleDir.mkdir("generators", 
                            "Contains histograms for each " 
                            "generator process"));
    //generatorPlotsPtr = new art::TFileDirectory(allParticleDir.mkdir("generators", "Contains histograms for each generator process"));

    art::TFileDirectory scaledAllParticlesDir   = allParticleDir.mkdir("scale_AllParticles", "Contains scaled histograms with all hit information ");
    art::TFileDirectory scaledMuonDir           = muonDir.mkdir("scale_Muons",        "Contains scaled histograms of muon hit information. ");
    art::TFileDirectory scaledAlphaDir          = alphaDir.mkdir("scale_Alphas", "Contains histograms with of alpha hit information. ");
    art::TFileDirectory scaledBetaDir           = betaDir.mkdir("scale_Betas", "Contains histograms of beta hit information. ");
    art::TFileDirectory scaledGammaDir          = gammaDir.mkdir("scale_Gammas", "Contains histograms of gamma hit information. ");
    art::TFileDirectory scaledNeutronDir        = neutronDir.mkdir("scale_Neutrons", "Contains histograms of neutron hit information");

    art::TFileDirectory unscaledAllParticlesDir = allParticleDir.mkdir("unscaled_AllParticles", "Contains histograms with all hit information ");
    art::TFileDirectory unscaledMuonDir         = muonDir.mkdir("unscaled_Muons",        "Contains histograms of muon hit information. ");
    art::TFileDirectory unscaledAlphaDir        = alphaDir.mkdir("unscaled_Alphas", "Contains histograms with of alpha hit information. ");
    art::TFileDirectory unscaledBetaDir         = betaDir.mkdir("unscaled_Betas", "Contains histograms of beta hit information. ");
    art::TFileDirectory unscaledGammaDir        = gammaDir.mkdir("unscaled_Gammas", "Contains histograms of gamma hit information. ");
    art::TFileDirectory unscaledNeutronDir      = neutronDir.mkdir("unscaled_Neutrons", "Contains histograms of neutron hit information");


    /*
    art::TFileDirectory scaled   =       private_service_tfs->mkdir("scaledHistograms", "Contains the histograms rescaled by energy deposited.");
    art::TFileDirectory unscaled =       private_service_tfs->mkdir("unscaledHistograms", "Contains the histograms without rescaling by energy deposited.");
    art::TFileDirectory unscaled2D =     private_service_tfs->mkdir("unscaled2D", "Contains the 2D XZ histograms without rescaling by energy deposited.");
    */


    private_deltaTime     =           private_service_tfs->make<TH1D>    ("h_deltaTime", 
        "Plot of Times between hits. Restricted to times > 500ns", 4500, 500, 5000);
    private_trackLength   =           allParticleDir.make<TH1D>("h_trackLengths",
        "Histogram of Track Lengths",1000, 0.0, 100.0  );
    fOpHitPEs             =           allParticleDir.make<TH1D>("h_OpHitPEs",
        "Histogram of PEs from each OpHit.",21, -0.5, 20.5  );
    fHitCharge            =           allParticleDir.make<TH1D>("h_HitCharge",
        "Histogram of #int ADC Units from each Hit",5000, 0.0, 1000.0  );
    fPeakAmpVXZ           =           allParticleDir.make<TH2D>("h_PeakAmpVXZ",
        "Histogram of HitPeakAmplitude vs Z Position (cm) and X Position (cm)",
        private_nBinsX, private_xCryoMin, private_xCryoMax, private_nBinsZ,private_zCryoMin,private_zCryoMax );
    fPeakAmpVXZ -> SetYTitle("Z Position (cm)");
    fPeakAmpVXZ -> SetXTitle("X Position (cm)");

    fIntegralVXZ          =           allParticleDir.make<TH2D>("fIntegralVXZ",
        "Histogram of HitIntegral vs Z Position (cm) and X Position (cm)",
        private_nBinsXPE, private_xCryoMin, private_xCryoMax,        private_nBinsZ, private_zCryoMin, private_zCryoMax        );
    fIntegralVXZ-> SetYTitle("Z Position (cm)");
    fIntegralVXZ-> SetZTitle("X Position (cm)");

    fHitPEVXZ             =           allParticleDir.make<TH2D>("fHitPEVXZ",
        "Histogram of PhotoElectrons Detected vs Z Position (cm) and X Position (cm)",
        private_nBinsXPE, private_xCryoMin, private_xCryoMax,        private_nBinsZ, private_zCryoMin, private_zCryoMax        );
    fHitPEVXZ-> SetYTitle("Z Position (cm)");
    fHitPEVXZ-> SetXTitle("X Position (cm)");

    fAlpha_PeakAmpVXZ     =           alphaDir.make<TH2D>("fAlpha_PeakAmpVXZ",
        "Histogram of HitPeakAmplitude vs Z Position (cm) and X Position (cm) from #alpha's",
        private_nBinsX, private_xCryoMin, private_xCryoMax,        private_nBinsZ, private_zCryoMin, private_zCryoMax        );
    fAlpha_PeakAmpVXZ-> SetYTitle("Z Position (cm)");
    fAlpha_PeakAmpVXZ-> SetXTitle("X Position (cm)");

    fAlpha_IntegralVXZ    =           alphaDir.make<TH2D>("fAlpha_IntegralVXZ",
        "Histogram of HitIntegral vs Z Position (cm) and X Position (cm) from #alpha's",
        private_nBinsX, private_xCryoMin, private_xCryoMax,        private_nBinsZ, private_zCryoMin, private_zCryoMax         );
    fAlpha_IntegralVXZ-> SetYTitle("Z Position (cm)");
    fAlpha_IntegralVXZ-> SetXTitle("X Position (cm)");

    fAlpha_HitPEVXZ       =           alphaDir.make<TH2D>("fAlpha_HitPEVXZ",
        "Histogram of PhotoElectrons Detected vs Z Position (cm) and X Position (cm) from #alpha's",
        private_nBinsXPE, private_xCryoMin, private_xCryoMax,        private_nBinsZ, private_zCryoMin, private_zCryoMax        );
    fAlpha_HitPEVXZ-> SetYTitle("Z Position (cm)");
    fAlpha_HitPEVXZ-> SetXTitle("X Position (cm)");

    fBeta_PeakAmpVXZ      =           betaDir.make<TH2D>("fBeta_PeakAmpVXZ",
        "Histogram of HitPeakAmplitude vs Z Position (cm) and X Position (cm) from #beta's", 
        private_nBinsX, private_xCryoMin, private_xCryoMax,        private_nBinsZ, private_zCryoMin, private_zCryoMax        );
    fBeta_PeakAmpVXZ-> SetYTitle("Z Position (cm)");
    fBeta_PeakAmpVXZ-> SetXTitle("X Position (cm)");

    fBeta_IntegralVXZ     =           betaDir.make<TH2D>("fBeta_IntegralVXZ",
        "Histogram of HitIntegral vs Z Position (cm) and X Position (cm) from #beta's",
        private_nBinsX, private_xCryoMin, private_xCryoMax,        private_nBinsZ, private_zCryoMin, private_zCryoMax        );
    fBeta_IntegralVXZ-> SetYTitle("Z Position (cm)");
    fBeta_IntegralVXZ-> SetXTitle("X Position (cm)");

    fBeta_HitPEVXZ        =           betaDir.make<TH2D>("fBeta_HitPEVXZ",
        "Histogram of PhotoElectrons Detected vs Z Position (cm) and X Position (cm) from #beta's",
        private_nBinsXPE, private_xCryoMin, private_xCryoMax,        private_nBinsZ, private_zCryoMin, private_zCryoMax        );
    fBeta_HitPEVXZ-> SetYTitle("Z Position (cm)");
    fBeta_HitPEVXZ-> SetXTitle("X Position (cm)");

    fGamma_PeakAmpVXZ     =           gammaDir.make<TH2D>("fGamma_PeakAmpVXZ",
        "Histogram of HitPeakAmplitude vs Z Position (cm) and X Position (cm) from #gamma's",
        private_nBinsX, private_xCryoMin, private_xCryoMax,        private_nBinsZ, private_zCryoMin, private_zCryoMax        );
    fGamma_PeakAmpVXZ-> SetYTitle("Z Position (cm)");
    fGamma_PeakAmpVXZ-> SetXTitle("X Position (cm)");

    fGamma_IntegralVXZ    =           gammaDir.make<TH2D>("fGamma_IntegralVXZ",
        "Histogram of HitIntegral vs Z Position (cm) and X Position (cm) from #gamma's",
        private_nBinsX, private_xCryoMin, private_xCryoMax,        private_nBinsZ, private_zCryoMin, private_zCryoMax         );
    fGamma_IntegralVXZ-> SetYTitle("Z Position (cm)");
    fGamma_IntegralVXZ-> SetXTitle("X Position (cm)");

    fGamma_HitPEVXZ       =           gammaDir.make<TH2D>("fGamma_HitPEVXZ",
        "Histogram of PhotoElectrons Detected vs Z Position (cm) and X Position (cm) from #gamma's",
        private_nBinsXPE, private_xCryoMin, private_xCryoMax,        private_nBinsZ, private_zCryoMin, private_zCryoMax        );
    fGamma_HitPEVXZ-> SetYTitle("Z Position (cm)");
    fGamma_HitPEVXZ-> SetXTitle("X Position (cm)");

    fNeutron_PeakAmpVXZ   =           neutronDir.make<TH2D>("fNeutron_PeakAmpVXZ",
        "Histogram of HitPeakAmplitude vs Z Position (cm) and X Position (cm) from neutrons",
        private_nBinsX, private_xCryoMin, private_xCryoMax,        private_nBinsZ, private_zCryoMin, private_zCryoMax        );
    fNeutron_PeakAmpVXZ-> SetYTitle("Z Position (cm)");
    fNeutron_PeakAmpVXZ-> SetXTitle("X Position (cm)");

    fNeutron_IntegralVXZ  =           neutronDir.make<TH2D>("fNeutron_IntegralVXZ",
        "Histogram of HitIntegral vs Z Position (cm) and X Position (cm) from neutrons",
        private_nBinsX, private_xCryoMin, private_xCryoMax,        private_nBinsZ, private_zCryoMin, private_zCryoMax         );
    fNeutron_IntegralVXZ-> SetYTitle("Z Position (cm)");
    fNeutron_IntegralVXZ-> SetXTitle("X Position (cm)");

    fNeutron_HitPEVXZ     =           neutronDir.make<TH2D>("fNeutron_HitPEVXZ",
        "Histogram of PhotoElectrons Detected vs Z Position (cm) and X Position (cm) from neutrons",
        private_nBinsXPE, private_xCryoMin, private_xCryoMax,        private_nBinsZ, private_zCryoMin, private_zCryoMax        );
    fNeutron_HitPEVXZ-> SetYTitle("Z Position (cm)");
    fNeutron_HitPEVXZ-> SetXTitle("X Position (cm)");

    //make scaled histograms
    fPeakAmpScaledVXhist         =    scaledAllParticlesDir.make<TH1D>("fPeakAmpScaledVXhist", 
        "Plot of Summed Hit Peak Amplitudes vs X position (Distance from APA in cm)", private_nBinsX, private_xCryoMin, private_xCryoMax);
    fIntegralScaledVXhist        =    scaledAllParticlesDir.make<TH1D>("fIntegralScaledVXhist", 
        "Plot of Summed Hit Integral Values vs X position (Distance from APA in cm)", private_nBinsX, private_xCryoMin, private_xCryoMax);
    fHitPEScaledVXhist           =    scaledAllParticlesDir.make<TH1D>("fHitPEScaledVXhist", 
        "Plot of the number of Detected Photo Electrons vs X position (Distance from APA in cm)", private_nBinsXPE, private_xCryoMin, private_xCryoMax);
    muon_PeakAmpScaledVXhist     =    scaledMuonDir.make<TH1D>("muon_PeakAmpScaledVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Muons", private_nBinsX, private_xCryoMin, private_xCryoMax);
    muon_IntegralScaledVXhist    =    scaledMuonDir.make<TH1D>("muon_IntegralScaledVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Muons", private_nBinsX, private_xCryoMin, private_xCryoMax);
    muon_HitPEScaledVXhist       =    scaledMuonDir.make<TH1D>("muon_HitPEScaledVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Muons", private_nBinsXPE, private_xCryoMin, private_xCryoMax);
    alpha_PeakAmpScaledVXhist    =    scaledAlphaDir.make<TH1D>("alpha_PeakAmpScaledVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Alphas", private_nBinsX, private_xCryoMin, private_xCryoMax);
    alpha_IntegralScaledVXhist   =    scaledAlphaDir.make<TH1D>("alpha_IntegralScaledVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Alphas", private_nBinsX, private_xCryoMin, private_xCryoMax);
    alpha_HitPEScaledVXhist      =    scaledAlphaDir.make<TH1D>("alpha_HitPEScaledVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Alphas", private_nBinsXPE, private_xCryoMin, private_xCryoMax);
    beta_PeakAmpScaledVXhist     =    scaledBetaDir.make<TH1D>("beta_PeakAmpScaledVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Betas", private_nBinsX, private_xCryoMin, private_xCryoMax);
    beta_IntegralScaledVXhist    =    scaledBetaDir.make<TH1D>("beta_IntegralScaledVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Betas", private_nBinsX, private_xCryoMin, private_xCryoMax);
    beta_HitPEScaledVXhist       =    scaledBetaDir.make<TH1D>("beta_HitPEScaledVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Betas", private_nBinsXPE, private_xCryoMin, private_xCryoMax);
    gamma_PeakAmpScaledVXhist    =    scaledGammaDir.make<TH1D>("gamma_PeakAmpScaledVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for gammas", private_nBinsX, private_xCryoMin, private_xCryoMax);
    gamma_IntegralScaledVXhist   =    scaledGammaDir.make<TH1D>("gamma_IntegralScaledVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for gammas", private_nBinsX, private_xCryoMin, private_xCryoMax);
    gamma_HitPEScaledVXhist      =    scaledGammaDir.make<TH1D>("gamma_HitPEScaledVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for gammas", private_nBinsX, private_xCryoMin, private_xCryoMax);
    neutron_PeakAmpScaledVXhist  =    scaledNeutronDir.make<TH1D>("neutron_PeakAmpScaledVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for neutrons", private_nBinsX, private_xCryoMin, private_xCryoMax);
    neutron_IntegralScaledVXhist =    scaledNeutronDir.make<TH1D>("neutron_IntegralScaledVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for neutrons", private_nBinsX, private_xCryoMin, private_xCryoMax);
    neutron_HitPEScaledVXhist    =    scaledNeutronDir.make<TH1D>("neutron_HitPEScaledVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for neutrons", private_nBinsXPE, private_xCryoMin, private_xCryoMax);

    

    //Make unscaled Histograms
    fPeakAmpVXhist               =    unscaledAllParticlesDir.make<TH1D>("fPeakAmpVXhist", 
        "Plot of Summed Hit Peak Amplitudes vs X position (Distance from APA in cm)", private_nBinsX, private_xCryoMin, private_xCryoMax);
    fIntegralVXhist              =    unscaledAllParticlesDir.make<TH1D>("fIntegralVXhist", 
        "Plot of Summed Hit Integral Values vs X position (Distance from APA in cm)", private_nBinsX, private_xCryoMin, private_xCryoMax);
    fHitPEVXhist                 =    unscaledAllParticlesDir.make<TH1D>("fHitPEVXhist", 
        "Plot of the number of Detected Photo Electrons vs X position (Distance from APA in cm)", private_nBinsX, private_xCryoMin, private_xCryoMax);
    fFlashPEsHist                =    unscaledAllParticlesDir.make<TH1D>("fFlashPEsHist", 
        "Plot of the number of Detected Photo Electrons per Flash", private_nPEsBins, private_nPEsMin, private_nPEsMax);
    muon_PeakAmpVXhist           =    unscaledMuonDir.make<TH1D>("muon_PeakAmpVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Muons", private_nBinsX, private_xCryoMin, private_xCryoMax);
    muon_IntegralVXhist =     unscaledMuonDir.make<TH1D>("muon_IntegralVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Muons", private_nBinsX, private_xCryoMin, private_xCryoMax);
    muon_HitPEVXhist =           unscaledMuonDir.make<TH1D>("muon_HitPEVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Muons", private_nBinsXPE, private_xCryoMin, private_xCryoMax);
//    muon_FlashPEsHist =             unscaledMuonDir.make<TH1D>("muon_FlashPEsHist",
  //      "PhotoElectrons per Flash for Muons", private_nPEsBins, private_nPEsMin, private_xCryoMax);
    alpha_PeakAmpVXhist =     unscaledAlphaDir.make<TH1D>("alpha_PeakAmpVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Alphas", private_nBinsX, private_xCryoMin, private_xCryoMax);
    alpha_IntegralVXhist =    unscaledAlphaDir.make<TH1D>("alpha_IntegralVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Alphas", private_nBinsX, private_xCryoMin, private_xCryoMax);
    alpha_HitPEVXhist =          unscaledAlphaDir.make<TH1D>("alpha_HitPEVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Alphas", private_nBinsXPE, private_xCryoMin, private_xCryoMax);
//    alpha_FlashPEsHist =          unscaledAlphaDir.make<TH1D>("alpha_FlashPEsHist",
  //      "PhotoElectrons per Flash for Alphas", private_nPEsBins, private_nPEsMin, private_nPEsMax);
    beta_PeakAmpVXhist =      unscaledBetaDir.make<TH1D>("beta_PeakAmpVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Betas", private_nBinsX, private_xCryoMin, private_xCryoMax);
    beta_IntegralVXhist =     unscaledBetaDir.make<TH1D>("beta_IntegralVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Betas", private_nBinsX, private_xCryoMin, private_xCryoMax);
    beta_HitPEVXhist =           unscaledBetaDir.make<TH1D>("beta_HitPEVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Betas", private_nBinsXPE, private_xCryoMin, private_xCryoMax);
//    beta_FlashPEsHist =           unscaledBetaDir.make<TH1D>("beta_FlashPEsHist",
  //      "PhotoElectrons per OpFlash for Betas", private_nPEsBins, private_nPEsMin, private_nPEsMax);
    gamma_PeakAmpVXhist =     unscaledGammaDir.make<TH1D>("gamma_PeakAmpVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for gammas", private_nBinsX, private_xCryoMin, private_xCryoMax);
    gamma_IntegralVXhist =    unscaledGammaDir.make<TH1D>("gamma_IntegralVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for gammas", private_nBinsX, private_xCryoMin, private_xCryoMax);
    gamma_HitPEVXhist =          unscaledGammaDir.make<TH1D>("gamma_HitPEVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for gammas", private_nBinsXPE, private_xCryoMin, private_xCryoMax);
//    gamma_FlashPEsHist =          unscaledGammaDir.make<TH1D>("gamma_FlashPEsHist",
  //      "PhotoElectrons per OpFlash for gammas", private_nPEsBins, private_nPEsMin, private_nPEsMax);
    neutron_PeakAmpVXhist =     scaledNeutronDir.make<TH1D>("neutron_PeakAmpVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for neutrons", private_nBinsX, private_xCryoMin, private_xCryoMax);
    neutron_IntegralVXhist =    scaledNeutronDir.make<TH1D>("neutron_IntegralVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for neutrons", private_nBinsX, private_xCryoMin, private_xCryoMax);
    neutron_HitPEVXhist =          scaledNeutronDir.make<TH1D>("neutron_HitPEVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for neutrons", private_nBinsX, private_xCryoMin, private_xCryoMax);
//    neutron_FlashPEsHist =          scaledNeutronDir.make<TH1D>("neutron_FlashPEsHist",
  //      "PhotoElectrons per OpFlash for neutrons", private_nPEsBins, private_nPEsMin, private_nPEsMax);

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
    fAlpha_PeakAmpVXZ      -> Scale( scaleFrac );
    fAlpha_IntegralVXZ     -> Scale( scaleFrac );
    fAlpha_HitPEVXZ        -> Scale( scaleFrac );
    fBeta_PeakAmpVXZ       -> Scale( scaleFrac );
    fBeta_IntegralVXZ      -> Scale( scaleFrac );
    fBeta_HitPEVXZ         -> Scale( scaleFrac );
    fGamma_PeakAmpVXZ      -> Scale( scaleFrac );
    fGamma_IntegralVXZ     -> Scale( scaleFrac );
    fGamma_HitPEVXZ        -> Scale( scaleFrac );
 
    fNeutron_PeakAmpVXZ    -> Scale( scaleFrac );
    fNeutron_IntegralVXZ   -> Scale( scaleFrac );
    fNeutron_HitPEVXZ      -> Scale( scaleFrac );

    //Unscalled Histograms
    fPeakAmpVXhist         -> Scale( scaleFrac );
    fIntegralVXhist        -> Scale( scaleFrac );
    fHitPEVXhist           -> Scale( scaleFrac );
    fFlashPEsHist          -> Scale( scaleFrac );
    muon_PeakAmpVXhist     -> Scale( scaleFrac );
    muon_IntegralVXhist    -> Scale( scaleFrac );
    muon_HitPEVXhist       -> Scale( scaleFrac );
    alpha_PeakAmpVXhist    -> Scale( scaleFrac );
    alpha_IntegralVXhist   -> Scale( scaleFrac );
    alpha_HitPEVXhist      -> Scale( scaleFrac );
    beta_PeakAmpVXhist     -> Scale( scaleFrac );
    beta_IntegralVXhist    -> Scale( scaleFrac );
    beta_HitPEVXhist       -> Scale( scaleFrac );
    gamma_PeakAmpVXhist    -> Scale( scaleFrac );
    gamma_IntegralVXhist   -> Scale( scaleFrac );
    gamma_HitPEVXhist        -> Scale( scaleFrac );
    neutron_PeakAmpVXhist    -> Scale( scaleFrac );
    neutron_IntegralVXhist   -> Scale( scaleFrac );
    neutron_HitPEVXhist      -> Scale( scaleFrac );

    //Scaled histograms
    fPeakAmpScaledVXhist     -> Scale( scaleFrac );
    fIntegralScaledVXhist    -> Scale( scaleFrac );
    fHitPEScaledVXhist       -> Scale( scaleFrac );
    muon_PeakAmpScaledVXhist -> Scale( scaleFrac );
    muon_IntegralScaledVXhist -> Scale( scaleFrac );
    muon_HitPEScaledVXhist -> Scale( scaleFrac );
    alpha_PeakAmpScaledVXhist -> Scale( scaleFrac );
    alpha_IntegralScaledVXhist -> Scale( scaleFrac );
    alpha_HitPEScaledVXhist -> Scale( scaleFrac );
    beta_PeakAmpScaledVXhist -> Scale( scaleFrac );
    beta_IntegralScaledVXhist -> Scale( scaleFrac );
    beta_HitPEScaledVXhist -> Scale( scaleFrac );
    gamma_PeakAmpScaledVXhist -> Scale( scaleFrac );
    gamma_IntegralScaledVXhist -> Scale( scaleFrac );
    gamma_HitPEScaledVXhist -> Scale( scaleFrac );
    neutron_PeakAmpScaledVXhist -> Scale( scaleFrac );
    neutron_IntegralScaledVXhist -> Scale( scaleFrac );
    neutron_HitPEScaledVXhist -> Scale( scaleFrac );

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
    std::cout<<"Analyze\n";
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
    std::cout<<"TrackList length: "<<trackList.size()<<"\n";
    for( auto track : trackList ){
      std::cout<<"Track Length: "<<(track->Length())<<"\n";
      private_trackLength-> Fill(track->Length());
    }


    //Call Charge events.
    std::map<std::string, int> nHitsPerGen;
    for( auto ptrHit: hitList){
      try{
        recob::Hit hit = *ptrHit;
        double hitInt = hit.Integral();
        double hitAmp = hit.PeakAmplitude();
        double scaleTot = 0.;
        std::vector<double> xyzPos;
        xyzPos = (private_service_bt->HitToXYZ(ptrHit));
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
            int pid =  particle->PdgCode();
            double scaleFrac   = eveIDE.energyFrac;
            double scaleEn     = eveIDE.energy;
            scaleTot += scaleEn;
            std::vector<double> const evePos = { particle->Trajectory().X(0),  particle->Trajectory().Y(0), particle->Trajectory().Z(0) };
    

            if( abs(pid)==13 ){                                      //muon

              muon_PeakAmpScaledVXhist      -> Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              muon_IntegralScaledVXhist     -> Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);


              muon_PeakAmpVXhist            -> Fill(xyzPos.at(0), hitAmp*scaleFrac);
              muon_IntegralVXhist           -> Fill(xyzPos.at(0), hitInt*scaleFrac);

            }else if( abs(pid)==11 ){                                //beta
              fBeta_PeakAmpVXZ              -> Fill(xyzPos.at(0), xyzPos.at(2), hitAmp*scaleFrac);
              fBeta_IntegralVXZ             -> Fill(xyzPos.at(0), xyzPos.at(2), hitInt*scaleFrac);


              beta_PeakAmpScaledVXhist      -> Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              beta_IntegralScaledVXhist     -> Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);


              beta_PeakAmpVXhist            -> Fill(xyzPos.at(0), hitAmp*scaleFrac);
              beta_IntegralVXhist           -> Fill(xyzPos.at(0), hitInt*scaleFrac);

            }else if( abs(pid)==22){                                 //gamma
              fGamma_PeakAmpVXZ             -> Fill(xyzPos.at(0), xyzPos.at(2), hitAmp*scaleFrac);
              fGamma_IntegralVXZ            -> Fill(xyzPos.at(0), xyzPos.at(2), hitInt*scaleFrac);


              gamma_PeakAmpScaledVXhist     -> Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              gamma_IntegralScaledVXhist    -> Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);

              
              gamma_PeakAmpVXhist           -> Fill(xyzPos.at(0), hitAmp*scaleFrac);
              gamma_IntegralVXhist          -> Fill(xyzPos.at(0), hitInt*scaleFrac);

            }else if( pid==1000020040 ){                             //alpha
              fAlpha_PeakAmpVXZ             -> Fill(xyzPos.at(0), xyzPos.at(2), hitAmp*scaleFrac);
              fAlpha_IntegralVXZ            -> Fill(xyzPos.at(0), xyzPos.at(2), hitInt*scaleFrac);


              alpha_PeakAmpScaledVXhist     -> Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              alpha_IntegralScaledVXhist    -> Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);


              alpha_PeakAmpVXhist           -> Fill(xyzPos.at(0), hitAmp*scaleFrac);
              alpha_IntegralVXhist          -> Fill(xyzPos.at(0), hitInt*scaleFrac);

            }else if( abs(pid) == 2112 ){
              fNeutron_PeakAmpVXZ           -> Fill(xyzPos.at(0), xyzPos.at(2), hitAmp*scaleFrac);
              fNeutron_IntegralVXZ          -> Fill(xyzPos.at(0), xyzPos.at(2), hitInt*scaleFrac);


              neutron_PeakAmpScaledVXhist   -> Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              neutron_IntegralScaledVXhist  -> Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);


              neutron_PeakAmpVXhist         -> Fill(xyzPos.at(0), hitAmp*scaleFrac);
              neutron_IntegralVXhist        -> Fill(xyzPos.at(0), hitInt*scaleFrac);
            }
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


    //Call Op Events

    for(const art::Ptr<recob::OpHit>& ptrOpHit: opHitList){
      try{
        recob::OpHit opHit = *ptrOpHit;
        double nPE=opHit.PE();
        double scaleTot=0.;
        //need position,
        std::vector<double> xyzPos          = private_service_pbt  -> OpHitToXYZ(ptrOpHit);
        std::vector<sim::TrackSDP>  eveSDPs = private_service_pbt  -> OpHitToEveSDPs(ptrOpHit);
        for(const sim::TrackSDP eveSDP : eveSDPs){
          auto tID = eveSDP.trackID;
          const simb::MCParticle* particle = private_service_pbt->TrackIDToParticle(eveSDP.trackID);
          if(particle){
            std::string generatorName = tIdToLabel[ tID ];
            TH1D* genPEHist = private_generatorOpHitPEs[generatorName];

            int pid =  particle->PdgCode();
            double scaleFrac = eveSDP.energyFrac;
            double scaleEn   = eveSDP.energy;
            scaleTot += scaleEn;
            std::vector<double> const evePos = { particle->Trajectory().X(0),  particle->Trajectory().Y(0), particle->Trajectory().Z(0) };
  
            if( abs(pid)==13 ){
              muon_HitPEScaledVXhist   -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              muon_HitPEVXhist         -> Fill(xyzPos.at(0),nPE*scaleFrac);
            }else if( abs(pid)==11 ){
              fBeta_HitPEVXZ           -> Fill(xyzPos.at(0), xyzPos.at(2), nPE*scaleFrac);
              beta_HitPEScaledVXhist   -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              beta_HitPEVXhist         -> Fill(xyzPos.at(0),nPE*scaleFrac);
            }else if( abs(pid)==22){
              fGamma_HitPEVXZ          -> Fill(xyzPos.at(0), xyzPos.at(2), nPE*scaleFrac);
              gamma_HitPEScaledVXhist  -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              gamma_HitPEVXhist        -> Fill(xyzPos.at(0),nPE*scaleFrac);
            }else if( pid==1000020040 ){
              fAlpha_HitPEVXZ          -> Fill(xyzPos.at(0), xyzPos.at(2), nPE*scaleFrac);
              alpha_HitPEScaledVXhist  -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              alpha_HitPEVXhist        -> Fill(xyzPos.at(0),nPE*scaleFrac);
            }else if ( abs(pid) == 2212 ){
              fNeutron_HitPEVXZ          -> Fill(xyzPos.at(0), xyzPos.at(2), nPE*scaleFrac);
              neutron_HitPEScaledVXhist  -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              neutron_HitPEVXhist        -> Fill(xyzPos.at(0),nPE*scaleFrac);
            }
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

