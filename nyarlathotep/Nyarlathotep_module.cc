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

    int nEvt =0;
    int nPEsBins = 101;
    int stdBinWidth = 2; //Standard bin width in centimeters
    double stdBinWidthPE = 0.5; //Standard bin width for PE plots in centimeters
    double nPEsMin = -0.5;
    double nPEsMax = 100.5;

/*    double modulo (double a, double b);
    double modulo (double a, int b);
    int modulo (int a, int b);
    int signF (int a);
    int signF (double a);*/

    std::string fOpHitLabel;
    std::string fOpFlashLabel;
    std::string fHitLabel;
    std::string fBTRLabel;
    std::string fPBTRLabel;

    art::ServiceHandle<art::TFileService>        tfs;
    art::ServiceHandle<cheat::PhotonBackTracker> pbt;
    art::ServiceHandle<cheat::BackTracker>       bt;
    std::unique_ptr<art::TFileDirectory> generatorPlotsPtr;
    geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();
//    art::ServiceHandle<geo::Geometry> geom;

    Double_t xCryoMin, xCryoMax, yCryoMin, yCryoMax, zCryoMin, zCryoMax, xTpcMin, xTpcMax, yTpcMin, yTpcMax, zTpcMin, zTpcMax, edgeDelta;

    int nBinsX = 1;

    std::map<std::string, TH1D*> generatorIntegrals;
//    std::map<std::string, TH1D*> generatorPEs;

    TH1D* fDeltaTime;

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

    for(geo::CryostatID const& cID: geom->IterateCryostatIDs()){
      Double_t* bounds = new double[6];
      geom->CryostatBoundaries(bounds, cID);
      std::cout<< "Bounds xMax " << bounds[1] << "min " << bounds[0] <<"\n";
      xCryoMin = std::min(xCryoMin, bounds[0]);
      xCryoMax = std::max(xCryoMax, bounds[1]);
      yCryoMin = std::min(yCryoMin, bounds[2]);
      yCryoMax = std::max(yCryoMax, bounds[3]);
      zCryoMin = std::min(zCryoMin, bounds[4]);
      zCryoMax = std::max(zCryoMax, bounds[5]);
    }
    //Setting TPC boundaries
    for(geo::TPCGeo const& TPC: geom->IterateTPCs()){
      xTpcMin = std::min(xTpcMin, TPC.MinX());
      xTpcMax = std::max(xTpcMax, TPC.MaxX());
      yTpcMin = std::min(yTpcMin, TPC.MinY());
      yTpcMax = std::max(yTpcMax, TPC.MaxY());
      zTpcMin = std::min(zTpcMin, TPC.MinZ());
      zTpcMax = std::max(zTpcMax, TPC.MaxZ());
    }
    nBinsX = 1;
    int nBinsXPE = 1;
    //int nBinsX = int((xCryoMax - xCryoMin)/10);
    {
      int nStdBins = int((xCryoMax - xCryoMin)/stdBinWidth);
      if(nStdBins<100){ nBinsX=100; }
      else{ nBinsX=nStdBins;}
      int nStdBinsPE = int((xCryoMax-xCryoMin)/stdBinWidthPE);
      if(nStdBinsPE<100){ nBinsXPE=100; }
      else{ nBinsXPE=nStdBinsPE;}

    }
    int nBinsZ = int((zCryoMax - zCryoMin)/stdBinWidth);
    
    //TDirectories
    
    art::TFileDirectory allParticleDir      = tfs->mkdir("allParticles", "Contains histograms with all hit information ");
    art::TFileDirectory muonDir               = tfs->mkdir("muons",        "Contains histograms of muon hit information. ");
    art::TFileDirectory alphaDir              = tfs->mkdir("alphas", "Contains histograms with of alpha hit information. ");
    art::TFileDirectory betaDir               = tfs->mkdir("betas", "Contains histograms of beta hit information. ");
    art::TFileDirectory gammaDir              = tfs->mkdir("gammas", "Contains histograms of gamma hit information. ");
    art::TFileDirectory neutronDir            = tfs->mkdir("neutrons", "Contains histograms of neutron hit information");

    generatorPlotsPtr = std::make_unique<art::TFileDirectory>(allParticleDir.mkdir("generators", "Contains histograms for each generator process"));

    art::TFileDirectory scaledAllParticlesDir      = allParticleDir.mkdir("scale_AllParticles", "Contains scaled histograms with all hit information ");
    art::TFileDirectory scaledMuonDir             = muonDir.mkdir("scale_Muons",        "Contains scaled histograms of muon hit information. ");
    art::TFileDirectory scaledAlphaDir            = alphaDir.mkdir("scale_Alphas", "Contains histograms with of alpha hit information. ");
    art::TFileDirectory scaledBetaDir             = betaDir.mkdir("scale_Betas", "Contains histograms of beta hit information. ");
    art::TFileDirectory scaledGammaDir            = gammaDir.mkdir("scale_Gammas", "Contains histograms of gamma hit information. ");
    art::TFileDirectory scaledNeutronDir          = neutronDir.mkdir("scale_Neutrons", "Contains histograms of neutron hit information");

    art::TFileDirectory unscaledAllParticlesDir        = allParticleDir.mkdir("unscaled_AllParticles", "Contains histograms with all hit information ");
    art::TFileDirectory unscaledMuonDir               = muonDir.mkdir("unscaled_Muons",        "Contains histograms of muon hit information. ");
    art::TFileDirectory unscaledAlphaDir              = alphaDir.mkdir("unscaled_Alphas", "Contains histograms with of alpha hit information. ");
    art::TFileDirectory unscaledBetaDir               = betaDir.mkdir("unscaled_Betas", "Contains histograms of beta hit information. ");
    art::TFileDirectory unscaledGammaDir              = gammaDir.mkdir("unscaled_Gammas", "Contains histograms of gamma hit information. ");
    art::TFileDirectory unscaledNeutronDir            = neutronDir.mkdir("unscaled_Neutrons", "Contains histograms of neutron hit information");


    /*
    art::TFileDirectory scaled   =       tfs->mkdir("scaledHistograms", "Contains the histograms rescaled by energy deposited.");
    art::TFileDirectory unscaled =       tfs->mkdir("unscaledHistograms", "Contains the histograms without rescaling by energy deposited.");
    art::TFileDirectory unscaled2D =     tfs->mkdir("unscaled2D", "Contains the 2D XZ histograms without rescaling by energy deposited.");
    */


    fDeltaTime =              tfs->make<TH1D>    ("fDeltaTime", 
        "Plot of Times between hits. Restricted to times > 500ns", 4500, 500, 5000);
    fPeakAmpVXZ  =            allParticleDir.make<TH2D>("fPeakAmpVXZ",
        "Histogram of HitPeakAmplitude vs Z Position (cm) and X Position (cm)",
        nBinsX, xCryoMin, xCryoMax,        nBinsZ,zCryoMin,zCryoMax        );
    fPeakAmpVXZ->SetYTitle("Z Position (cm)");
    fPeakAmpVXZ->SetXTitle("X Position (cm)");

    fIntegralVXZ =            allParticleDir.make<TH2D>("fIntegralVXZ",
        "Histogram of HitIntegral vs Z Position (cm) and X Position (cm)",
        nBinsXPE, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fIntegralVXZ->SetYTitle("Z Position (cm)");
    fIntegralVXZ->SetZTitle("X Position (cm)");

    fHitPEVXZ       =            allParticleDir.make<TH2D>("fHitPEVXZ",
        "Histogram of PhotoElectrons Detected vs Z Position (cm) and X Position (cm)",
        nBinsXPE, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fHitPEVXZ->SetYTitle("Z Position (cm)");
    fHitPEVXZ->SetXTitle("X Position (cm)");

    fAlpha_PeakAmpVXZ  =            alphaDir.make<TH2D>("fAlpha_PeakAmpVXZ",
        "Histogram of HitPeakAmplitude vs Z Position (cm) and X Position (cm) from #alpha's",
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fAlpha_PeakAmpVXZ->SetYTitle("Z Position (cm)");
    fAlpha_PeakAmpVXZ->SetXTitle("X Position (cm)");

    fAlpha_IntegralVXZ =            alphaDir.make<TH2D>("fAlpha_IntegralVXZ",
        "Histogram of HitIntegral vs Z Position (cm) and X Position (cm) from #alpha's",
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax         );
    fAlpha_IntegralVXZ->SetYTitle("Z Position (cm)");
    fAlpha_IntegralVXZ->SetXTitle("X Position (cm)");

    fAlpha_HitPEVXZ       =            alphaDir.make<TH2D>("fAlpha_HitPEVXZ",
        "Histogram of PhotoElectrons Detected vs Z Position (cm) and X Position (cm) from #alpha's",
        nBinsXPE, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fAlpha_HitPEVXZ->SetYTitle("Z Position (cm)");
    fAlpha_HitPEVXZ->SetXTitle("X Position (cm)");

    fBeta_PeakAmpVXZ  =            betaDir.make<TH2D>("fBeta_PeakAmpVXZ",
        "Histogram of HitPeakAmplitude vs Z Position (cm) and X Position (cm) from #beta's", 
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fBeta_PeakAmpVXZ->SetYTitle("Z Position (cm)");
    fBeta_PeakAmpVXZ->SetXTitle("X Position (cm)");

    fBeta_IntegralVXZ =            betaDir.make<TH2D>("fBeta_IntegralVXZ",
        "Histogram of HitIntegral vs Z Position (cm) and X Position (cm) from #beta's",
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fBeta_IntegralVXZ->SetYTitle("Z Position (cm)");
    fBeta_IntegralVXZ->SetXTitle("X Position (cm)");

    fBeta_HitPEVXZ       =            betaDir.make<TH2D>("fBeta_HitPEVXZ",
        "Histogram of PhotoElectrons Detected vs Z Position (cm) and X Position (cm) from #beta's",
        nBinsXPE, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fBeta_HitPEVXZ->SetYTitle("Z Position (cm)");
    fBeta_HitPEVXZ->SetXTitle("X Position (cm)");

    fGamma_PeakAmpVXZ  =            gammaDir.make<TH2D>("fGamma_PeakAmpVXZ",
        "Histogram of HitPeakAmplitude vs Z Position (cm) and X Position (cm) from #gamma's",
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fGamma_PeakAmpVXZ->SetYTitle("Z Position (cm)");
    fGamma_PeakAmpVXZ->SetXTitle("X Position (cm)");

    fGamma_IntegralVXZ =            gammaDir.make<TH2D>("fGamma_IntegralVXZ",
        "Histogram of HitIntegral vs Z Position (cm) and X Position (cm) from #gamma's",
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax         );
    fGamma_IntegralVXZ->SetYTitle("Z Position (cm)");
    fGamma_IntegralVXZ->SetXTitle("X Position (cm)");

    fGamma_HitPEVXZ       =            gammaDir.make<TH2D>("fGamma_HitPEVXZ",
        "Histogram of PhotoElectrons Detected vs Z Position (cm) and X Position (cm) from #gamma's",
        nBinsXPE, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fGamma_HitPEVXZ->SetYTitle("Z Position (cm)");
    fGamma_HitPEVXZ->SetXTitle("X Position (cm)");

    fNeutron_PeakAmpVXZ  =            neutronDir.make<TH2D>("fNeutron_PeakAmpVXZ",
        "Histogram of HitPeakAmplitude vs Z Position (cm) and X Position (cm) from neutrons",
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fNeutron_PeakAmpVXZ->SetYTitle("Z Position (cm)");
    fNeutron_PeakAmpVXZ->SetXTitle("X Position (cm)");

    fNeutron_IntegralVXZ =            neutronDir.make<TH2D>("fNeutron_IntegralVXZ",
        "Histogram of HitIntegral vs Z Position (cm) and X Position (cm) from neutrons",
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax         );
    fNeutron_IntegralVXZ->SetYTitle("Z Position (cm)");
    fNeutron_IntegralVXZ->SetXTitle("X Position (cm)");

    fNeutron_HitPEVXZ       =            neutronDir.make<TH2D>("fNeutron_HitPEVXZ",
        "Histogram of PhotoElectrons Detected vs Z Position (cm) and X Position (cm) from neutrons",
        nBinsXPE, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fNeutron_HitPEVXZ->SetYTitle("Z Position (cm)");
    fNeutron_HitPEVXZ->SetXTitle("X Position (cm)");

    //make scaled histograms
    fPeakAmpScaledVXhist =          scaledAllParticlesDir.make<TH1D>("fPeakAmpScaledVXhist", 
        "Plot of Summed Hit Peak Amplitudes vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fIntegralScaledVXhist =         scaledAllParticlesDir.make<TH1D>("fIntegralScaledVXhist", 
        "Plot of Summed Hit Integral Values vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fHitPEScaledVXhist =               scaledAllParticlesDir.make<TH1D>("fHitPEScaledVXhist", 
        "Plot of the number of Detected Photo Electrons vs X position (Distance from APA in cm)", nBinsXPE, xCryoMin, xCryoMax);
    muon_PeakAmpScaledVXhist =      scaledMuonDir.make<TH1D>("muon_PeakAmpScaledVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_IntegralScaledVXhist =     scaledMuonDir.make<TH1D>("muon_IntegralScaledVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_HitPEScaledVXhist =           scaledMuonDir.make<TH1D>("muon_HitPEScaledVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Muons", nBinsXPE, xCryoMin, xCryoMax);
    alpha_PeakAmpScaledVXhist =     scaledAlphaDir.make<TH1D>("alpha_PeakAmpScaledVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_IntegralScaledVXhist =    scaledAlphaDir.make<TH1D>("alpha_IntegralScaledVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_HitPEScaledVXhist =          scaledAlphaDir.make<TH1D>("alpha_HitPEScaledVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Alphas", nBinsXPE, xCryoMin, xCryoMax);
    beta_PeakAmpScaledVXhist =      scaledBetaDir.make<TH1D>("beta_PeakAmpScaledVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_IntegralScaledVXhist =     scaledBetaDir.make<TH1D>("beta_IntegralScaledVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_HitPEScaledVXhist =           scaledBetaDir.make<TH1D>("beta_HitPEScaledVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Betas", nBinsXPE, xCryoMin, xCryoMax);
    gamma_PeakAmpScaledVXhist =     scaledGammaDir.make<TH1D>("gamma_PeakAmpScaledVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_IntegralScaledVXhist =    scaledGammaDir.make<TH1D>("gamma_IntegralScaledVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_HitPEScaledVXhist =          scaledGammaDir.make<TH1D>("gamma_HitPEScaledVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    neutron_PeakAmpScaledVXhist =   scaledNeutronDir.make<TH1D>("neutron_PeakAmpScaledVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for neutrons", nBinsX, xCryoMin, xCryoMax);
    neutron_IntegralScaledVXhist =  scaledNeutronDir.make<TH1D>("neutron_IntegralScaledVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for neutrons", nBinsX, xCryoMin, xCryoMax);
    neutron_HitPEScaledVXhist =        scaledNeutronDir.make<TH1D>("neutron_HitPEScaledVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for neutrons", nBinsXPE, xCryoMin, xCryoMax);

    

    //Make unscaled Histograms
    fPeakAmpVXhist =          unscaledAllParticlesDir.make<TH1D>("fPeakAmpVXhist", 
        "Plot of Summed Hit Peak Amplitudes vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fIntegralVXhist =         unscaledAllParticlesDir.make<TH1D>("fIntegralVXhist", 
        "Plot of Summed Hit Integral Values vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fHitPEVXhist =               unscaledAllParticlesDir.make<TH1D>("fHitPEVXhist", 
        "Plot of the number of Detected Photo Electrons vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fFlashPEsHist =                 unscaledAllParticlesDir.make<TH1D>("fFlashPEsHist", 
        "Plot of the number of Detected Photo Electrons per Flash", nPEsBins, nPEsMin, nPEsMax);
    muon_PeakAmpVXhist =      unscaledMuonDir.make<TH1D>("muon_PeakAmpVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_IntegralVXhist =     unscaledMuonDir.make<TH1D>("muon_IntegralVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_HitPEVXhist =           unscaledMuonDir.make<TH1D>("muon_HitPEVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Muons", nBinsXPE, xCryoMin, xCryoMax);
//    muon_FlashPEsHist =             unscaledMuonDir.make<TH1D>("muon_FlashPEsHist",
  //      "PhotoElectrons per Flash for Muons", nPEsBins, nPEsMin, xCryoMax);
    alpha_PeakAmpVXhist =     unscaledAlphaDir.make<TH1D>("alpha_PeakAmpVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_IntegralVXhist =    unscaledAlphaDir.make<TH1D>("alpha_IntegralVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_HitPEVXhist =          unscaledAlphaDir.make<TH1D>("alpha_HitPEVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Alphas", nBinsXPE, xCryoMin, xCryoMax);
//    alpha_FlashPEsHist =          unscaledAlphaDir.make<TH1D>("alpha_FlashPEsHist",
  //      "PhotoElectrons per Flash for Alphas", nPEsBins, nPEsMin, nPEsMax);
    beta_PeakAmpVXhist =      unscaledBetaDir.make<TH1D>("beta_PeakAmpVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_IntegralVXhist =     unscaledBetaDir.make<TH1D>("beta_IntegralVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_HitPEVXhist =           unscaledBetaDir.make<TH1D>("beta_HitPEVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Betas", nBinsXPE, xCryoMin, xCryoMax);
//    beta_FlashPEsHist =           unscaledBetaDir.make<TH1D>("beta_FlashPEsHist",
  //      "PhotoElectrons per OpFlash for Betas", nPEsBins, nPEsMin, nPEsMax);
    gamma_PeakAmpVXhist =     unscaledGammaDir.make<TH1D>("gamma_PeakAmpVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_IntegralVXhist =    unscaledGammaDir.make<TH1D>("gamma_IntegralVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_HitPEVXhist =          unscaledGammaDir.make<TH1D>("gamma_HitPEVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for gammas", nBinsXPE, xCryoMin, xCryoMax);
//    gamma_FlashPEsHist =          unscaledGammaDir.make<TH1D>("gamma_FlashPEsHist",
  //      "PhotoElectrons per OpFlash for gammas", nPEsBins, nPEsMin, nPEsMax);
    neutron_PeakAmpVXhist =     scaledNeutronDir.make<TH1D>("neutron_PeakAmpVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for neutrons", nBinsX, xCryoMin, xCryoMax);
    neutron_IntegralVXhist =    scaledNeutronDir.make<TH1D>("neutron_IntegralVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for neutrons", nBinsX, xCryoMin, xCryoMax);
    neutron_HitPEVXhist =          scaledNeutronDir.make<TH1D>("neutron_HitPEVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for neutrons", nBinsX, xCryoMin, xCryoMax);
//    neutron_FlashPEsHist =          scaledNeutronDir.make<TH1D>("neutron_FlashPEsHist",
  //      "PhotoElectrons per OpFlash for neutrons", nPEsBins, nPEsMin, nPEsMax);

  }
  

   //-----------------------------------------------------------------------
  void Nyarlathotep::endJob()
  {
    std::cout<<"Flow Check. End Job\n Scaling by"<<(1.0/( double (nEvt) ))<<"\n" ;
    fFlashPEsHist->Scale( (1.0/( double (nEvt) )) );
    fPeakAmpVXZ->Scale( (1.0/( double (nEvt) )) );
    fIntegralVXZ->Scale( (1.0/( double (nEvt) )) );
    fHitPEVXZ->Scale( (1.0/( double (nEvt) )) );
    fAlpha_PeakAmpVXZ->Scale( (1.0/( double (nEvt) )) );
    fAlpha_IntegralVXZ->Scale( (1.0/( double (nEvt) )) );
    fAlpha_HitPEVXZ->Scale( (1.0/( double (nEvt) )) );
    fBeta_PeakAmpVXZ->Scale( (1.0/( double (nEvt) )) );
    fBeta_IntegralVXZ->Scale( (1.0/( double (nEvt) )) );
    fBeta_HitPEVXZ->Scale( (1.0/( double (nEvt) )) );
    fGamma_PeakAmpVXZ->Scale( (1.0/( double (nEvt) )) );
    fGamma_IntegralVXZ->Scale( (1.0/( double (nEvt) )) );
    fGamma_HitPEVXZ->Scale( (1.0/( double (nEvt) )) );

    fNeutron_PeakAmpVXZ->Scale( (1.0/( double (nEvt) )) );
    fNeutron_IntegralVXZ->Scale( (1.0/( double (nEvt) )) );
    fNeutron_HitPEVXZ->Scale( (1.0/( double (nEvt) )) );

    //Unscalled Histograms
    fPeakAmpVXhist->Scale( (1.0/( double (nEvt) )) );
    fIntegralVXhist->Scale( (1.0/( double (nEvt) )) );
    fHitPEVXhist->Scale( (1.0/( double (nEvt) )) );
    fFlashPEsHist->Scale( (1.0/( double (nEvt) )) );
    muon_PeakAmpVXhist->Scale( (1.0/( double (nEvt) )) );
    muon_IntegralVXhist->Scale( (1.0/( double (nEvt) )) );
    muon_HitPEVXhist->Scale( (1.0/( double (nEvt) )) );
    alpha_PeakAmpVXhist->Scale( (1.0/( double (nEvt) )) );
    alpha_IntegralVXhist->Scale( (1.0/( double (nEvt) )) );
    alpha_HitPEVXhist->Scale( (1.0/( double (nEvt) )) );
    beta_PeakAmpVXhist->Scale( (1.0/( double (nEvt) )) );
    beta_IntegralVXhist->Scale( (1.0/( double (nEvt) )) );
    beta_HitPEVXhist->Scale( (1.0/( double (nEvt) )) );
    gamma_PeakAmpVXhist->Scale( (1.0/( double (nEvt) )) );
    gamma_IntegralVXhist->Scale( (1.0/( double (nEvt) )) );
    gamma_HitPEVXhist->Scale( (1.0/( double (nEvt) )) );
    neutron_PeakAmpVXhist->Scale( (1.0/( double (nEvt) )) );
    neutron_IntegralVXhist->Scale( (1.0/( double (nEvt) )) );
    neutron_HitPEVXhist->Scale( (1.0/( double (nEvt) )) );

    //Scaled histograms
    fPeakAmpScaledVXhist->Scale( (1.0/( double (nEvt) )) );
    fIntegralScaledVXhist->Scale( (1.0/( double (nEvt) )) );
    fHitPEScaledVXhist->Scale( (1.0/( double (nEvt) )) );
    muon_PeakAmpScaledVXhist->Scale( (1.0/( double (nEvt) )) );
    muon_IntegralScaledVXhist->Scale( (1.0/( double (nEvt) )) );
    muon_HitPEScaledVXhist->Scale( (1.0/( double (nEvt) )) );
    alpha_PeakAmpScaledVXhist->Scale( (1.0/( double (nEvt) )) );
    alpha_IntegralScaledVXhist->Scale( (1.0/( double (nEvt) )) );
    alpha_HitPEScaledVXhist->Scale( (1.0/( double (nEvt) )) );
    beta_PeakAmpScaledVXhist->Scale( (1.0/( double (nEvt) )) );
    beta_IntegralScaledVXhist->Scale( (1.0/( double (nEvt) )) );
    beta_HitPEScaledVXhist->Scale( (1.0/( double (nEvt) )) );
    gamma_PeakAmpScaledVXhist->Scale( (1.0/( double (nEvt) )) );
    gamma_IntegralScaledVXhist->Scale( (1.0/( double (nEvt) )) );
    gamma_HitPEScaledVXhist->Scale( (1.0/( double (nEvt) )) );
    neutron_PeakAmpScaledVXhist->Scale( (1.0/( double (nEvt) )) );
    neutron_IntegralScaledVXhist->Scale( (1.0/( double (nEvt) )) );
    neutron_HitPEScaledVXhist->Scale( (1.0/( double (nEvt) )) );
  }
  //-----------------------------------------------------------------------
  void Nyarlathotep::reconfigure(fhicl::ParameterSet const& parameterSet)
  {
    // Read parameters from the .fcl file. The names in the arguments
    std::cout<<"reconfigure\n";
    fHitLabel             = parameterSet.get< std::string >("HitLabel");
    fOpHitLabel           = parameterSet.get< std::string >("PhotLabel");
    fOpFlashLabel         = parameterSet.get< std::string >("FlashLabel", "opflash");
    fBTRLabel             = parameterSet.get< std::string >("BackTrackerLabel");
    fPBTRLabel            = parameterSet.get< std::string >("PhotonBackTrackerLabel");
     
     
     
    //Setting Global Boundaries
    xCryoMin=0.0;
    xCryoMax=0.0;
    yCryoMin=0.0;
    yCryoMax=0.0;
    zCryoMin=0.0;
    zCryoMax=0.0;
    xTpcMin=0.0;
    xTpcMax=0.0;
    yTpcMin=0.0;
    yTpcMax=0.0;
    zTpcMin=0.0;
    zTpcMax=0.0;
    edgeDelta = 5.0;
  }

  //-----------------------------------------------------------------------
  void Nyarlathotep::analyze(const art::Event& evt) 
  {
    std::cout<<"Analyze\n";
    nEvt++;
    global_event = evt.id().event();
    art::Handle< std::vector< recob::Hit > > hitHandle;
    std::vector< art::Ptr< recob::Hit > > hitList;
    if (evt.getByLabel(fHitLabel, hitHandle) )
      art::fill_ptr_vector(hitList, hitHandle);

    art::Handle< std::vector< recob::OpHit > > opHitHandle;
    std::vector< art::Ptr< recob::OpHit > > opHitList;
    if (evt.getByLabel(fOpHitLabel, opHitHandle) )
      art::fill_ptr_vector(opHitList, opHitHandle);
    
    art::Handle< std::vector< recob::OpFlash > > opFlashHandle;
    std::vector< art::Ptr< recob::OpFlash > > opFlashList;
    if (evt.getByLabel(fOpFlashLabel, opFlashHandle) ){
      art::fill_ptr_vector(opFlashList, opFlashHandle);
    }


    std::vector< art::Handle< std::vector< simb::MCParticle > > > mcHandles;
    evt.getManyByType(mcHandles);
    std::map< simb::MCParticle, std::string> mPartToLabel;
    
    for( auto const& mcHandle : mcHandles ){
      const std::string& sModuleLabel = mcHandle.provenance()->moduleLabel();
      //I have the label here. I can use that to mimic DAQSimAna
      auto mcTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(sModuleLabel);
      art::FindManyP<simb::MCParticle> findMCParts(mcTrue, evt, "largeant" );
      std::vector<art::Ptr<simb::MCParticle > > mcParts = findMCParts.at(1);
      for(const art::Ptr<simb::MCParticle > ptr : mcParts){
        simb::MCParticle part = *ptr;
        mPartToLabel[part] = sModuleLabel;
      }
      std::string sPlotTitle = "Plot of Integrated Charge collected from generator ";
      sPlotTitle.append(sModuleLabel);
      if(generatorIntegrals.find(sModuleLabel) == generatorIntegrals.end()){
        TH1D* tmpHist = generatorPlotsPtr->make<TH1D>(sModuleLabel.c_str(),  sPlotTitle.c_str(), nBinsX, xCryoMin, xCryoMax);
        generatorIntegrals[sModuleLabel] = std::move(tmpHist);
      }
    }


    //Call Charge events.

    for( auto ptrHit: hitList){
      try{
        recob::Hit hit = *ptrHit;
        double hitInt = hit.Integral();
        double hitAmp = hit.PeakAmplitude();
        double scaleTot = 0.;
        std::vector<double> xyzPos;
        xyzPos = (bt->HitToXYZ(ptrHit));
        std::vector<sim::TrackIDE>  eveIDEs = bt->HitToEveID(ptrHit);
        for(const sim::TrackIDE eveIDE : eveIDEs){
          auto tID = eveIDE.trackID;
          const simb::MCParticle* particle=bt->TrackIDToParticle(tID);
          if(particle){
            //Look up the Label of the generator of this particle
            std::string processName = mPartToLabel[*particle];
            TH1D* genIntHist = generatorIntegrals[processName];
            int pid =  particle->PdgCode();
            double scaleFrac = eveIDE.energyFrac;
            double scaleEn   = eveIDE.energy;
            scaleTot += scaleEn;
            std::vector<double> const evePos = { particle->Trajectory().X(0),  particle->Trajectory().Y(0), particle->Trajectory().Z(0) };
    

            if( (evePos.at(2)>= zTpcMin-edgeDelta && evePos.at(2)<= zTpcMin + edgeDelta ) || (evePos.at(2)>= zTpcMax - edgeDelta && evePos.at(2)<= zTpcMax+edgeDelta)){//side
            }else if( (evePos.at(1)>= yTpcMax - edgeDelta && evePos.at(1) <= yTpcMax + edgeDelta) ){  //Top
            }else if( (evePos.at(0) >= 0-edgeDelta && evePos.at(0) <= 0+edgeDelta ) ){                //APA
            }else if( (evePos.at(0) >= xTpcMax-edgeDelta && evePos.at(0) <= xTpcMax+edgeDelta ) ){    //CPA
            }
            if( abs(pid)==13 ){                                      //muon

              muon_PeakAmpScaledVXhist  ->Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              muon_IntegralScaledVXhist ->Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);


              muon_PeakAmpVXhist  ->Fill(xyzPos.at(0), hitAmp*scaleFrac);
              muon_IntegralVXhist ->Fill(xyzPos.at(0), hitInt*scaleFrac);

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
              fNeutron_PeakAmpVXZ             -> Fill(xyzPos.at(0), xyzPos.at(2), hitAmp*scaleFrac);
              fNeutron_IntegralVXZ            -> Fill(xyzPos.at(0), xyzPos.at(2), hitInt*scaleFrac);


              neutron_PeakAmpScaledVXhist     -> Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              neutron_IntegralScaledVXhist    -> Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);


              neutron_PeakAmpVXhist           -> Fill(xyzPos.at(0), hitAmp*scaleFrac);
              neutron_IntegralVXhist          -> Fill(xyzPos.at(0), hitInt*scaleFrac);
            }
            genIntHist->Fill(xyzPos.at(0), hitInt*scaleFrac);
          }//End If Particle
        }
        if(scaleTot<1.0e-10){
          scaleTot=1.0;
        }
        //Write globals
        //write globals

        fPeakAmpVXhist ->Fill(xyzPos.at(0), hitAmp);
        fIntegralVXhist->Fill(xyzPos.at(0), hitInt);
        fPeakAmpScaledVXhist ->Fill(xyzPos.at(0), hitAmp/scaleTot);
        fIntegralScaledVXhist->Fill(xyzPos.at(0), hitInt/scaleTot);

        fPeakAmpVXZ->Fill(xyzPos.at(0), xyzPos.at(2), hitAmp);
        fIntegralVXZ->Fill(xyzPos.at(0), xyzPos.at(2), hitInt);
      }catch(...){

      }//end trycatch
    }//end for(hits)


    //Call Op Events

    for(const art::Ptr<recob::OpHit>& ptrOpHit: opHitList){
      try{
        recob::OpHit opHit = *ptrOpHit;
        //need position,
        std::vector<double> xyzPos = pbt->OpHitToXYZ(ptrOpHit);
        std::vector<sim::TrackSDP>  eveSDPs=pbt->OpHitToEveSDPs(ptrOpHit);
        double nPE=opHit.PE();
        double scaleTot=0.;
        for(const sim::TrackSDP eveSDP : eveSDPs){
          const simb::MCParticle* particle = pbt->TrackIDToParticle(eveSDP.trackID);
          if(particle){
            //Calculate Process Name Goes HERE!!!!!!!!!!!!
            /*
            std::string processName = "";
            std::string processPlotDescription = "Plot of PEs generated by ";
            processPlotDescription.append(processName);
            TH1D* tmpHist;
            if(generatorPEs.find(processName)==generatorPEs.end() ){
              generatorPEs[processName] = generatorPlots..make<TH1D>(processName, processPlotDescription, 
                  nBinsX, xCryoMin, xCryoMax);
            }
            */
            /////////////////////////////////////////////////Maintinence HERE!!!!!!!!!?///////////////////////

            int pid =  particle->PdgCode();
            double scaleFrac = eveSDP.energyFrac;
            double scaleEn   = eveSDP.energy;
            scaleTot += scaleEn;
            std::vector<double> const evePos = { particle->Trajectory().X(0),  particle->Trajectory().Y(0), particle->Trajectory().Z(0) };
            if( (evePos.at(2)>= zTpcMin-edgeDelta && evePos.at(2)<= zTpcMin + edgeDelta ) || (evePos.at(2)>= zTpcMax - edgeDelta && evePos.at(2)<= zTpcMax+edgeDelta)){ //side
            }else if( (evePos.at(1)>= yTpcMax - edgeDelta && evePos.at(1) <= yTpcMax + edgeDelta) ){ //Top
            }else if( (evePos.at(0) >= 0-edgeDelta && evePos.at(0) <= 0+edgeDelta ) ){ //APA
            }else if( (evePos.at(0) >= xTpcMax-edgeDelta && evePos.at(0) <= xTpcMax+edgeDelta ) ){ //CPA
            }
  
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
          }//End for(eveSDP)
        }
        if(scaleTot<1.0e-10){
          scaleTot=1.0;
        }
        //Write globals
        fHitPEScaledVXhist -> Fill(xyzPos.at(0), nPE/scaleTot);
        fHitPEVXhist -> Fill(xyzPos.at(0), nPE);

        fHitPEVXZ->Fill(xyzPos.at(0), xyzPos.at(2), nPE);

      }//end try
      catch(...){//Lazy catch. Need to catch -No sim SDPs found- error.
   //     throw;
      }
    }//End for(OpHit)

    //Call OpFlashes
    for( auto ptrOpFlash : opFlashList ){
      double nPEcount = ptrOpFlash->TotalPE();
     /*
      *std::vector<double> nPEs = ptrOpFlash->PEs();
      *for( double nPE : nPEs ){
      *  if( abs(pid)==13 ){
      *    muon_FlashPEsHist           -> Fill(nPEcount);
      *  }else if( abs(pid)==11 ){
      *    beta_FlashPEsHist           -> Fill(nPEcount);
      *  }else if( abs(pid)==22){
      *    gamma_FlashPEsHist          -> Fill(nPEcount);
      *  }else if( pid==1000020040 ){
      *    alpha_FlashPEsHist          -> Fill(nPEcount);
      *  }else if ( abs(pid) == 2212 ){
      *    neutron_FlashPEsHist        -> Fill(nPEcount);
      *  }
      *}
      */
      fFlashPEsHist->Fill( nPEcount );
    }

    

  } // Nyarlathotep::analyze()
  
  
  // This macro has to be defined for this module to be invoked from a
  // .fcl file; see Nyarlathotep.fcl for more information.
  DEFINE_ART_MODULE(Nyarlathotep)



} // namespace Nyarlathotep

#endif // Nyarlathotep_Module

