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
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
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
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/Exception.h"
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
#include "TProfile.h"
#include "TProfile3D.h"
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

    int nPEsBins = 51;
    double nPEsMin = -0.5;
    double nPEsMax = 50.5;

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
    geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();
//    art::ServiceHandle<geo::Geometry> geom;

    Double_t xCryoMin, xCryoMax, yCryoMin, yCryoMax, zCryoMin, zCryoMax, xTpcMin, xTpcMax, yTpcMin, yTpcMax, zTpcMin, zTpcMax, edgeDelta;


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
    TProfile* fPeakAmpVX;
    TProfile* fIntegralVX;
    TProfile* fHitPEVX;
    TProfile* localTop_PeakAmpVX;
    TProfile* localTop_IntegralVX;
    TProfile* localTop_HitPEVX;
    TProfile* localSides_PeakAmpVX;
    TProfile* localSides_IntegralVX;
    TProfile* localSides_HitPEVX;
    TProfile* localAPA_PeakAmpVX;
    TProfile* localAPA_IntegralVX;
    TProfile* localAPA_HitPEVX;
    TProfile* localCPA_PeakAmpVX;
    TProfile* localCPA_IntegralVX;
    TProfile* localCPA_HitPEVX;
    TProfile* muon_PeakAmpVX;
    TProfile* muon_IntegralVX;
    TProfile* muon_HitPEVX;
    TProfile* alpha_PeakAmpVX;
    TH1D* alpha_PeakAmpVXhist;
    TProfile* alpha_IntegralVX;
    TH1D* alpha_IntegralVXhist;
    TProfile* alpha_HitPEVX;
    TH1D* alpha_HitPEVXhist;
//    TH1D* alpha_FlashPEsHist;
    TProfile* beta_PeakAmpVX;
    TH1D* beta_PeakAmpVXhist;
    TProfile* beta_IntegralVX;
    TH1D* beta_IntegralVXhist;
    TProfile* beta_HitPEVX;
    TH1D* beta_HitPEVXhist;
//    TH1D* beta_FlashPEsHist;
    TProfile* gamma_PeakAmpVX;
    TH1D* gamma_PeakAmpVXhist;
    TProfile* gamma_IntegralVX;
    TH1D* gamma_IntegralVXhist;
    TProfile* gamma_HitPEVX;
    TH1D* gamma_HitPEVXhist;
//    TH1D* gamma_FlashPEsHist;
    TProfile* neutron_PeakAmpVX;
    TH1D* neutron_PeakAmpVXhist;
    TProfile* neutron_IntegralVX;
    TH1D* neutron_IntegralVXhist;
    TProfile* neutron_HitPEVX;
    TH1D* neutron_HitPEVXhist;
//    TH1D* neutron_FlashPEsHist;

    //Scaled histograms
    TProfile* fPeakAmpScaledVX;
    TH1D* fPeakAmpScaledVXhist;
    TProfile* fIntegralScaledVX;
    TH1D* fIntegralScaledVXhist;
    TProfile* fHitPEScaledVX;
    TH1D* fHitPEScaledVXhist;
    TProfile* localTop_PeakAmpScaledVX;
    TProfile* localTop_IntegralScaledVX;
    TProfile* localTop_HitPEScaledVX;
    TProfile* localSides_PeakAmpScaledVX;
    TProfile* localSides_IntegralScaledVX;
    TProfile* localSides_HitPEScaledVX;
    TProfile* localAPA_PeakAmpScaledVX;
    TProfile* localAPA_IntegralScaledVX;
    TProfile* localAPA_HitPEScaledVX;
    TProfile* localCPA_PeakAmpScaledVX;
    TProfile* localCPA_IntegralScaledVX;
    TProfile* localCPA_HitPEScaledVX;
    TProfile* muon_PeakAmpScaledVX;
    TH1D* muon_PeakAmpScaledVXhist;
    TProfile* muon_IntegralScaledVX;
    TH1D* muon_IntegralScaledVXhist;
    TProfile* muon_HitPEScaledVX;
    TH1D* muon_HitPEScaledVXhist;
    TProfile* alpha_PeakAmpScaledVX;
    TH1D* alpha_PeakAmpScaledVXhist;
    TProfile* alpha_IntegralScaledVX;
    TH1D* alpha_IntegralScaledVXhist;
    TProfile* alpha_HitPEScaledVX;
    TH1D* alpha_HitPEScaledVXhist;
    TProfile* beta_PeakAmpScaledVX;
    TH1D* beta_PeakAmpScaledVXhist;
    TProfile* beta_IntegralScaledVX;
    TH1D* beta_IntegralScaledVXhist;
    TProfile* beta_HitPEScaledVX;
    TH1D* beta_HitPEScaledVXhist;
    TProfile* gamma_PeakAmpScaledVX;
    TH1D* gamma_PeakAmpScaledVXhist;
    TProfile* gamma_IntegralScaledVX;
    TH1D* gamma_IntegralScaledVXhist;
    TProfile* gamma_HitPEScaledVX;
    TH1D* gamma_HitPEScaledVXhist;
    TProfile* neutron_PeakAmpScaledVX;
    TH1D* neutron_PeakAmpScaledVXhist;
    TProfile* neutron_IntegralScaledVX;
    TH1D* neutron_IntegralScaledVXhist;
    TProfile* neutron_HitPEScaledVX;
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
    int nBinsX = 1;
    //int nBinsX = int((xCryoMax - xCryoMin)/10);
    {
      int n1cmBins = int((xCryoMax - xCryoMin)/10);
      if(n1cmBins<100) nBinsX=100;
      else nBinsX=n1cmBins;

    }
    int nBinsZ = int((zCryoMax - zCryoMin)/10);
    
    //TDirectories
    
    art::TFileDirectory allParticleDir      = tfs->mkdir("allParticles", "Contains histograms with all hit information ");
    art::TFileDirectory muonDir               = tfs->mkdir("muons",        "Contains histograms of muon hit information. ");
    art::TFileDirectory alphaDir              = tfs->mkdir("alphas", "Contains histograms with of alpha hit information. ");
    art::TFileDirectory betaDir               = tfs->mkdir("betas", "Contains histograms of beta hit information. ");
    art::TFileDirectory gammaDir              = tfs->mkdir("gammas", "Contains histograms of gamma hit information. ");
    art::TFileDirectory neutronDir            = tfs->mkdir("neutrons", "Contains histograms of neutron hit information");

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
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fIntegralVXZ->SetYTitle("Z Position (cm)");
    fIntegralVXZ->SetZTitle("X Position (cm)");

    fHitPEVXZ       =            allParticleDir.make<TH2D>("fHitPEVXZ",
        "Histogram of PhotoElectrons Detected vs Z Position (cm) and X Position (cm)",
        nBinsX*20, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
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
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
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
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
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
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
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
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fNeutron_HitPEVXZ->SetYTitle("Z Position (cm)");
    fNeutron_HitPEVXZ->SetXTitle("X Position (cm)");

    //make scaled histograms
    fPeakAmpScaledVX =              scaledAllParticlesDir.make<TProfile>("fPeakAmpScaledVX", 
        "Plot of Hit Peak Amplitudes scaled by Track Deposited Energy vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fPeakAmpScaledVXhist =          scaledAllParticlesDir.make<TH1D>("fPeakAmpScaledVXhist", 
        "Plot of Summed Hit Peak Amplitudes vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fIntegralScaledVX =             scaledAllParticlesDir.make<TProfile>("fIntegralScaledVX", 
        "Plot of the Hit Integral Values vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fIntegralScaledVXhist =         scaledAllParticlesDir.make<TH1D>("fIntegralScaledVXhist", 
        "Plot of Summed Hit Integral Values vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fHitPEScaledVX =                   scaledAllParticlesDir.make<TProfile>("fHitPEScaledVX", 
        "Plot of the number of Detected Photo Electrons per OpHit vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fHitPEScaledVXhist =               scaledAllParticlesDir.make<TH1D>("fHitPEScaledVXhist", 
        "Plot of the number of Detected Photo Electrons vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    localTop_PeakAmpScaledVX =      scaledAllParticlesDir.make<TProfile>("fromTop_PeakAmpScaledVX", 
        "Peak Amplitude per hit vs X (distance from APA in cm) for particles from the top of the TPC", nBinsX, xCryoMin, xCryoMax);
    localTop_IntegralScaledVX =     scaledAllParticlesDir.make<TProfile>("fromTop_IntegralScaledVX", 
        "Integral per hit vs X (distance from APA in cm) for particles from the top of the TPC", nBinsX, xCryoMin, xCryoMax);
    localTop_HitPEScaledVX =           scaledAllParticlesDir.make<TProfile>("fromTop_PEsVsXpos", 
        "Photo Electrons per ophit vs X (distance from APA) for particles from the top of the TPC", nBinsX, xCryoMin, xCryoMax);
    localSides_PeakAmpScaledVX =    scaledAllParticlesDir.make<TProfile>("fromSides_PeakAmpScaledVX",
        "Plot of the Hit Peak Amplitudes vs X ositions (distance from APA in cm) for particles from the sides of the TPC", nBinsX, xCryoMin, xCryoMax);
    localSides_IntegralScaledVX =   scaledAllParticlesDir.make<TProfile>("fromSides_IntegralScaledVX",
        "Integral values per hit vs X position (distance from APA in cm) for particles from the TPC sides" , nBinsX, xCryoMin, xCryoMax);
    localSides_HitPEScaledVX =         scaledAllParticlesDir.make<TProfile>("fromSides_HitPEScaledVX", 
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for particles from the TPC sides" , nBinsX, xCryoMin, xCryoMax);
    localAPA_PeakAmpScaledVX =      scaledAllParticlesDir.make<TProfile>("fromAPA_PeakAmpScaledVX", 
        "Peak Amp per hit vs X position (distance from APA in cm) for particles from the APA", nBinsX, xCryoMin, xCryoMax);
    localAPA_IntegralScaledVX =     scaledAllParticlesDir.make<TProfile>("fromAPA_IntegralScaledVX",
        "Integral per hit Vs X position (distance from APA in cm) for particles from the APA", nBinsX, xCryoMin, xCryoMax);
    localAPA_HitPEScaledVX =           scaledAllParticlesDir.make<TProfile>("fromAPA_HitPEScaledVX",
        "PhotoElectrons per OpHit Vs X position (distance from APA in cm) for particles from the APA", nBinsX, xCryoMin, xCryoMax);
    muon_PeakAmpScaledVX =          scaledMuonDir.make<TProfile>("muon_PeakAmpScaledVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_PeakAmpScaledVXhist =      scaledMuonDir.make<TH1D>("muon_PeakAmpScaledVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_IntegralScaledVX =         scaledMuonDir.make<TProfile>("muon_IntegralScaledVX",
        "Integral per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_IntegralScaledVXhist =     scaledMuonDir.make<TH1D>("muon_IntegralScaledVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_HitPEScaledVX =               scaledMuonDir.make<TProfile>("muon_HitPEScaledVX",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_HitPEScaledVXhist =           scaledMuonDir.make<TH1D>("muon_HitPEScaledVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    alpha_PeakAmpScaledVX =         scaledAlphaDir.make<TProfile>("alpha_PeakAmpScaledVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_PeakAmpScaledVXhist =     scaledAlphaDir.make<TH1D>("alpha_PeakAmpScaledVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_IntegralScaledVX =        scaledAlphaDir.make<TProfile>("alpha_IntegralScaledVX",
        "Integral per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_IntegralScaledVXhist =    scaledAlphaDir.make<TH1D>("alpha_IntegralScaledVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_HitPEScaledVX =              scaledAlphaDir.make<TProfile>("alpha_HitPEScaledVX",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_HitPEScaledVXhist =          scaledAlphaDir.make<TH1D>("alpha_HitPEScaledVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    beta_PeakAmpScaledVX =          scaledBetaDir.make<TProfile>("beta_PeakAmpScaledVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_PeakAmpScaledVXhist =      scaledBetaDir.make<TH1D>("beta_PeakAmpScaledVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_IntegralScaledVX =         scaledBetaDir.make<TProfile>("beta_IntegralScaledVX",
        "Integral per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_IntegralScaledVXhist =     scaledBetaDir.make<TH1D>("beta_IntegralScaledVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_HitPEScaledVX =               scaledBetaDir.make<TProfile>("beta_HitPEScaledVX",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_HitPEScaledVXhist =           scaledBetaDir.make<TH1D>("beta_HitPEScaledVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    gamma_PeakAmpScaledVX =         scaledGammaDir.make<TProfile>("gamma_PeakAmpScaledVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_PeakAmpScaledVXhist =     scaledGammaDir.make<TH1D>("gamma_PeakAmpScaledVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_IntegralScaledVX =        scaledGammaDir.make<TProfile>("gamma_IntegralScaledVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_IntegralScaledVXhist =    scaledGammaDir.make<TH1D>("gamma_IntegralScaledVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_HitPEScaledVX =              scaledGammaDir.make<TProfile>("gamma_HitPEScaledVX",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_HitPEScaledVXhist =          scaledGammaDir.make<TH1D>("gamma_HitPEScaledVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    neutron_PeakAmpScaledVX =       scaledNeutronDir.make<TProfile>("neutron_PeakAmpScaledVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for neutrons", nBinsX, xCryoMin, xCryoMax);
    neutron_PeakAmpScaledVXhist =   scaledNeutronDir.make<TH1D>("neutron_PeakAmpScaledVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for neutrons", nBinsX, xCryoMin, xCryoMax);
    neutron_IntegralScaledVX =      scaledNeutronDir.make<TProfile>("neutron_IntegralScaledVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for neutrons", nBinsX, xCryoMin, xCryoMax);
    neutron_IntegralScaledVXhist =  scaledNeutronDir.make<TH1D>("neutron_IntegralScaledVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for neutrons", nBinsX, xCryoMin, xCryoMax);
    neutron_HitPEScaledVX =            scaledNeutronDir.make<TProfile>("neutron_HitPEScaledVX",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for neutrons", nBinsX, xCryoMin, xCryoMax);
    neutron_HitPEScaledVXhist =        scaledNeutronDir.make<TH1D>("neutron_HitPEScaledVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for neutrons", nBinsX, xCryoMin, xCryoMax);

    

    //Make unscaled Histograms
    fPeakAmpVX =              unscaledAllParticlesDir.make<TProfile>("fPeakAmpVX", 
        "Plot of Hit Peak Amplitudes vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fPeakAmpVXhist =          unscaledAllParticlesDir.make<TH1D>("fPeakAmpVXhist", 
        "Plot of Summed Hit Peak Amplitudes vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fIntegralVX =             unscaledAllParticlesDir.make<TProfile>("fIntegralVX", 
        "Plot of the Hit Integral Values vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fIntegralVXhist =         unscaledAllParticlesDir.make<TH1D>("fIntegralVXhist", 
        "Plot of Summed Hit Integral Values vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fHitPEVX =                   unscaledAllParticlesDir.make<TProfile>("fHitPEVX", 
        "Plot of the number of Detected Photo Electrons per OpHit vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fHitPEVXhist =               unscaledAllParticlesDir.make<TH1D>("fHitPEVXhist", 
        "Plot of the number of Detected Photo Electrons vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fFlashPEsHist =                 unscaledAllParticlesDir.make<TH1D>("fFlashPEsHist", 
        "Plot of the number of Detected Photo Electrons per Flash", nPEsBins, nPEsMin, nPEsMax);
    localTop_PeakAmpVX =      unscaledAllParticlesDir.make<TProfile>("fromTop_PeakAmpVX", 
        "Peak Amplitude per hit vs X (distance from APA in cm) for particles from the top of the TPC", nBinsX, xCryoMin, xCryoMax);
    localTop_IntegralVX =     unscaledAllParticlesDir.make<TProfile>("fromTop_IntegralVX", 
        "Integral per hit vs X (distance from APA in cm) for particles from the top of the TPC", nBinsX, xCryoMin, xCryoMax);
    localTop_HitPEVX =           unscaledAllParticlesDir.make<TProfile>("fromTop_PEsVsXpos", 
        "Photo Electrons per ophit vs X (distance from APA) for particles from the top of the TPC", nBinsX, xCryoMin, xCryoMax);
    localSides_PeakAmpVX =    unscaledAllParticlesDir.make<TProfile>("fromSides_PeakAmpVX",
        "Plot of the Hit Peak Amplitudes vs X ositions (distance from APA in cm) for particles from the sides of the TPC", nBinsX, xCryoMin, xCryoMax);
    localSides_IntegralVX =   unscaledAllParticlesDir.make<TProfile>("fromSides_IntegralVX",
        "Integral values per hit vs X position (distance from APA in cm) for particles from the TPC sides" , nBinsX, xCryoMin, xCryoMax);
    localSides_HitPEVX =         unscaledAllParticlesDir.make<TProfile>("fromSides_HitPEVX", 
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for particles from the TPC sides" , nBinsX, xCryoMin, xCryoMax);
    localAPA_PeakAmpVX =      unscaledAllParticlesDir.make<TProfile>("fromAPA_PeakAmpVX", 
        "Peak Amp per hit vs X position (distance from APA in cm) for particles from the APA", nBinsX, xCryoMin, xCryoMax);
    localAPA_IntegralVX =     unscaledAllParticlesDir.make<TProfile>("fromAPA_IntegralVX",
        "Integral per hit Vs X position (distance from APA in cm) for particles from the APA", nBinsX, xCryoMin, xCryoMax);
    localAPA_HitPEVX =           unscaledAllParticlesDir.make<TProfile>("fromAPA_HitPEVX",
        "PhotoElectrons per OpHit Vs X position (distance from APA in cm) for particles from the APA", nBinsX, xCryoMin, xCryoMax);
    muon_PeakAmpVX =          unscaledMuonDir.make<TProfile>("muon_PeakAmpVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_PeakAmpVXhist =      unscaledMuonDir.make<TH1D>("muon_PeakAmpVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_IntegralVX =         unscaledMuonDir.make<TProfile>("muon_IntegralVX",
        "Integral per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_IntegralVXhist =     unscaledMuonDir.make<TH1D>("muon_IntegralVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_HitPEVX =               unscaledMuonDir.make<TProfile>("muon_HitPEVX",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_HitPEVXhist =           unscaledMuonDir.make<TH1D>("muon_HitPEVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
//    muon_FlashPEsHist =             unscaledMuonDir.make<TH1D>("muon_FlashPEsHist",
  //      "PhotoElectrons per Flash for Muons", nPEsBins, nPEsMin, xCryoMax);
    alpha_PeakAmpVX =         unscaledAlphaDir.make<TProfile>("alpha_PeakAmpVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_PeakAmpVXhist =     unscaledAlphaDir.make<TH1D>("alpha_PeakAmpVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_IntegralVX =        unscaledAlphaDir.make<TProfile>("alpha_IntegralVX",
        "Integral per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_IntegralVXhist =    unscaledAlphaDir.make<TH1D>("alpha_IntegralVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_HitPEVX =              unscaledAlphaDir.make<TProfile>("alpha_HitPEVX",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_HitPEVXhist =          unscaledAlphaDir.make<TH1D>("alpha_HitPEVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
//    alpha_FlashPEsHist =          unscaledAlphaDir.make<TH1D>("alpha_FlashPEsHist",
  //      "PhotoElectrons per Flash for Alphas", nPEsBins, nPEsMin, nPEsMax);
    beta_PeakAmpVX =          unscaledBetaDir.make<TProfile>("beta_PeakAmpVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_PeakAmpVXhist =      unscaledBetaDir.make<TH1D>("beta_PeakAmpVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_IntegralVX =         unscaledBetaDir.make<TProfile>("beta_IntegralVX",
        "Integral per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_IntegralVXhist =     unscaledBetaDir.make<TH1D>("beta_IntegralVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_HitPEVX =               unscaledBetaDir.make<TProfile>("beta_HitPEVX",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_HitPEVXhist =           unscaledBetaDir.make<TH1D>("beta_HitPEVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
//    beta_FlashPEsHist =           unscaledBetaDir.make<TH1D>("beta_FlashPEsHist",
  //      "PhotoElectrons per OpFlash for Betas", nPEsBins, nPEsMin, nPEsMax);
    gamma_PeakAmpVX =         unscaledGammaDir.make<TProfile>("gamma_PeakAmpVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_PeakAmpVXhist =     unscaledGammaDir.make<TH1D>("gamma_PeakAmpVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_IntegralVX =        unscaledGammaDir.make<TProfile>("gamma_IntegralVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_IntegralVXhist =    unscaledGammaDir.make<TH1D>("gamma_IntegralVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_HitPEVX =              unscaledGammaDir.make<TProfile>("gamma_HitPEVX",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_HitPEVXhist =          unscaledGammaDir.make<TH1D>("gamma_HitPEVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
//    gamma_FlashPEsHist =          unscaledGammaDir.make<TH1D>("gamma_FlashPEsHist",
  //      "PhotoElectrons per OpFlash for gammas", nPEsBins, nPEsMin, nPEsMax);
    neutron_PeakAmpVX =         scaledNeutronDir.make<TProfile>("neutron_PeakAmpVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for neutrons", nBinsX, xCryoMin, xCryoMax);
    neutron_PeakAmpVXhist =     scaledNeutronDir.make<TH1D>("neutron_PeakAmpVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for neutrons", nBinsX, xCryoMin, xCryoMax);
    neutron_IntegralVX =        scaledNeutronDir.make<TProfile>("neutron_IntegralVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for neutrons", nBinsX, xCryoMin, xCryoMax);
    neutron_IntegralVXhist =    scaledNeutronDir.make<TH1D>("neutron_IntegralVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for neutrons", nBinsX, xCryoMin, xCryoMax);
    neutron_HitPEVX =              scaledNeutronDir.make<TProfile>("neutron_HitPEVX",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for neutrons", nBinsX, xCryoMin, xCryoMax);
    neutron_HitPEVXhist =          scaledNeutronDir.make<TH1D>("neutron_HitPEVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for neutrons", nBinsX, xCryoMin, xCryoMax);
//    neutron_FlashPEsHist =          scaledNeutronDir.make<TH1D>("neutron_FlashPEsHist",
  //      "PhotoElectrons per OpFlash for neutrons", nPEsBins, nPEsMin, nPEsMax);

  }
  

   //-----------------------------------------------------------------------
  void Nyarlathotep::endJob()
  {

  }
  //-----------------------------------------------------------------------
  void Nyarlathotep::reconfigure(fhicl::ParameterSet const& parameterSet)
  {
    // Read parameters from the .fcl file. The names in the arguments
    std::cout<<"reconfigure\n";
    fHitLabel             = parameterSet.get< std::string >("HitLabel");
    fOpHitLabel           = parameterSet.get< std::string >("PhotLabel");
    fOpFlashLabel           = parameterSet.get< std::string >("PhotLabel", "opflash");
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
    if (evt.getByLabel(fOpFlashLabel, opFlashHandle) )
      art::fill_ptr_vector(opFlashList, opFlashHandle);

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
            int pid =  particle->PdgCode();
            double scaleFrac = eveIDE.energyFrac;
            double scaleEn   = eveIDE.energy;
            scaleTot += scaleEn;
            std::vector<double> const evePos = { particle->Trajectory().X(0),  particle->Trajectory().Y(0), particle->Trajectory().Z(0) };
    

            if( (evePos.at(2)>= zTpcMin-edgeDelta && evePos.at(2)<= zTpcMin + edgeDelta ) || (evePos.at(2)>= zTpcMax - edgeDelta && evePos.at(2)<= zTpcMax+edgeDelta)){//side
              localSides_PeakAmpScaledVX ->Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              localSides_IntegralScaledVX->Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);
              localSides_PeakAmpVX ->Fill(xyzPos.at(0), hitAmp*scaleFrac);
              localSides_IntegralVX->Fill(xyzPos.at(0), hitInt*scaleFrac);
            }else if( (evePos.at(1)>= yTpcMax - edgeDelta && evePos.at(1) <= yTpcMax + edgeDelta) ){  //Top
              localTop_PeakAmpScaledVX   ->Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              localTop_IntegralScaledVX  ->Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);
              localTop_PeakAmpVX   ->Fill(xyzPos.at(0), hitAmp*scaleFrac);
              localTop_IntegralVX  ->Fill(xyzPos.at(0), hitInt*scaleFrac);
            }else if( (evePos.at(0) >= 0-edgeDelta && evePos.at(0) <= 0+edgeDelta ) ){                //APA
              localAPA_PeakAmpScaledVX   ->Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              localAPA_IntegralScaledVX  ->Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);
              localAPA_PeakAmpVX   ->Fill(xyzPos.at(0), hitAmp*scaleFrac);
              localAPA_IntegralVX  ->Fill(xyzPos.at(0), hitInt*scaleFrac);
            }else if( (evePos.at(0) >= xTpcMax-edgeDelta && evePos.at(0) <= xTpcMax+edgeDelta ) ){    //CPA
              localCPA_PeakAmpScaledVX   ->Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              localCPA_IntegralScaledVX  ->Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);
              localCPA_PeakAmpVX   ->Fill(xyzPos.at(0), hitAmp*scaleFrac);
              localCPA_IntegralVX  ->Fill(xyzPos.at(0), hitInt*scaleFrac);
            }
            if( abs(pid)==13 ){                                      //muon
              muon_PeakAmpScaledVX  ->Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              muon_IntegralScaledVX ->Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);

              muon_PeakAmpScaledVXhist  ->Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              muon_IntegralScaledVXhist ->Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);

              muon_PeakAmpVX  ->Fill(xyzPos.at(0), hitAmp*scaleFrac);
              muon_IntegralVX ->Fill(xyzPos.at(0), hitInt*scaleFrac);

              muon_PeakAmpVXhist  ->Fill(xyzPos.at(0), hitAmp*scaleFrac);
              muon_IntegralVXhist ->Fill(xyzPos.at(0), hitInt*scaleFrac);

            }else if( abs(pid)==11 ){                                //beta
              fBeta_PeakAmpVXZ              -> Fill(xyzPos.at(0), xyzPos.at(2), hitAmp*scaleFrac);
              fBeta_IntegralVXZ             -> Fill(xyzPos.at(0), xyzPos.at(2), hitInt*scaleFrac);

              beta_PeakAmpScaledVX          -> Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              beta_IntegralScaledVX         -> Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);

              beta_PeakAmpScaledVXhist      -> Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              beta_IntegralScaledVXhist     -> Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);

              beta_PeakAmpVX                -> Fill(xyzPos.at(0), hitAmp*scaleFrac);
              beta_IntegralVX               -> Fill(xyzPos.at(0), hitInt*scaleFrac);

              beta_PeakAmpVXhist            -> Fill(xyzPos.at(0), hitAmp*scaleFrac);
              beta_IntegralVXhist           -> Fill(xyzPos.at(0), hitInt*scaleFrac);

            }else if( abs(pid)==22){                                 //gamma
              fGamma_PeakAmpVXZ             -> Fill(xyzPos.at(0), xyzPos.at(2), hitAmp*scaleFrac);
              fGamma_IntegralVXZ            -> Fill(xyzPos.at(0), xyzPos.at(2), hitInt*scaleFrac);

              gamma_PeakAmpScaledVX         -> Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              gamma_IntegralScaledVX        -> Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);

              gamma_PeakAmpScaledVXhist     -> Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              gamma_IntegralScaledVXhist    -> Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);

              gamma_PeakAmpVX               -> Fill(xyzPos.at(0), hitAmp*scaleFrac);
              gamma_IntegralVX              -> Fill(xyzPos.at(0), hitInt*scaleFrac);
              
              gamma_PeakAmpVXhist           -> Fill(xyzPos.at(0), hitAmp*scaleFrac);
              gamma_IntegralVXhist          -> Fill(xyzPos.at(0), hitInt*scaleFrac);

            }else if( pid==1000020040 ){                             //alpha
              fAlpha_PeakAmpVXZ             -> Fill(xyzPos.at(0), xyzPos.at(2), hitAmp*scaleFrac);
              fAlpha_IntegralVXZ            -> Fill(xyzPos.at(0), xyzPos.at(2), hitInt*scaleFrac);

              alpha_PeakAmpScaledVX         -> Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              alpha_IntegralScaledVX        -> Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);

              alpha_PeakAmpScaledVXhist     -> Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              alpha_IntegralScaledVXhist    -> Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);

              alpha_PeakAmpVX               -> Fill(xyzPos.at(0), hitAmp*scaleFrac);
              alpha_IntegralVX              -> Fill(xyzPos.at(0), hitInt*scaleFrac);

              alpha_PeakAmpVXhist           -> Fill(xyzPos.at(0), hitAmp*scaleFrac);
              alpha_IntegralVXhist          -> Fill(xyzPos.at(0), hitInt*scaleFrac);

            }else if( abs(pid) == 2112 ){
              fNeutron_PeakAmpVXZ             -> Fill(xyzPos.at(0), xyzPos.at(2), hitAmp*scaleFrac);
              fNeutron_IntegralVXZ            -> Fill(xyzPos.at(0), xyzPos.at(2), hitInt*scaleFrac);

              neutron_PeakAmpScaledVX         -> Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              neutron_IntegralScaledVX        -> Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);

              neutron_PeakAmpScaledVXhist     -> Fill(xyzPos.at(0), hitAmp*scaleFrac/scaleEn);
              neutron_IntegralScaledVXhist    -> Fill(xyzPos.at(0), hitInt*scaleFrac/scaleEn);

              neutron_PeakAmpVX               -> Fill(xyzPos.at(0), hitAmp*scaleFrac);
              neutron_IntegralVX              -> Fill(xyzPos.at(0), hitInt*scaleFrac);

              neutron_PeakAmpVXhist           -> Fill(xyzPos.at(0), hitAmp*scaleFrac);
              neutron_IntegralVXhist          -> Fill(xyzPos.at(0), hitInt*scaleFrac);
            }
          }
        }
        if(scaleTot<1.0e-10){
          scaleTot=1.0;
        }
        //Write globals
        //write globals
        fPeakAmpVX ->Fill(xyzPos.at(0), hitAmp);
        fIntegralVX->Fill(xyzPos.at(0), hitInt);
        fPeakAmpVXhist ->Fill(xyzPos.at(0), hitAmp);
        fIntegralVXhist->Fill(xyzPos.at(0), hitInt);
        fPeakAmpScaledVX ->Fill(xyzPos.at(0), hitAmp/scaleTot);
        fIntegralScaledVX->Fill(xyzPos.at(0), hitInt/scaleTot);
        fPeakAmpScaledVXhist ->Fill(xyzPos.at(0), hitAmp/scaleTot);
        fIntegralScaledVXhist->Fill(xyzPos.at(0), hitInt/scaleTot);

        fPeakAmpVXZ->Fill(xyzPos.at(0), xyzPos.at(2), hitAmp);
        fIntegralVXZ->Fill(xyzPos.at(0), xyzPos.at(2), hitInt);
      }catch(...){

      }//end trycatch
    }


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
            int pid =  particle->PdgCode();
            double scaleFrac = eveSDP.energyFrac;
            double scaleEn   = eveSDP.energy;
            scaleTot += scaleEn;
            std::vector<double> const evePos = { particle->Trajectory().X(0),  particle->Trajectory().Y(0), particle->Trajectory().Z(0) };
            if( (evePos.at(2)>= zTpcMin-edgeDelta && evePos.at(2)<= zTpcMin + edgeDelta ) || (evePos.at(2)>= zTpcMax - edgeDelta && evePos.at(2)<= zTpcMax+edgeDelta)){ //side
              localSides_HitPEScaledVX  ->Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              localSides_HitPEVX  ->Fill(xyzPos.at(0),nPE*scaleFrac);
            }else if( (evePos.at(1)>= yTpcMax - edgeDelta && evePos.at(1) <= yTpcMax + edgeDelta) ){ //Top
              localTop_HitPEScaledVX    ->Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              localTop_HitPEVX    ->Fill(xyzPos.at(0),nPE*scaleFrac);
            }else if( (evePos.at(0) >= 0-edgeDelta && evePos.at(0) <= 0+edgeDelta ) ){ //APA
              localAPA_HitPEScaledVX    ->Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              localAPA_HitPEVX    ->Fill(xyzPos.at(0),nPE*scaleFrac);
            }else if( (evePos.at(0) >= xTpcMax-edgeDelta && evePos.at(0) <= xTpcMax+edgeDelta ) ){ //CPA
              localCPA_HitPEScaledVX    ->Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              localCPA_HitPEVX    ->Fill(xyzPos.at(0),nPE*scaleFrac);
            }
  
            if( abs(pid)==13 ){
              muon_HitPEScaledVX       -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              muon_HitPEScaledVXhist   -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              muon_HitPEVX             -> Fill(xyzPos.at(0),nPE*scaleFrac);
              muon_HitPEVXhist         -> Fill(xyzPos.at(0),nPE*scaleFrac);
            }else if( abs(pid)==11 ){
              fBeta_HitPEVXZ           -> Fill(xyzPos.at(0), xyzPos.at(2), nPE*scaleFrac);
              beta_HitPEScaledVX       -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              beta_HitPEScaledVXhist   -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              beta_HitPEVX             -> Fill(xyzPos.at(0),nPE*scaleFrac);
              beta_HitPEVXhist         -> Fill(xyzPos.at(0),nPE*scaleFrac);
            }else if( abs(pid)==22){
              fGamma_HitPEVXZ          -> Fill(xyzPos.at(0), xyzPos.at(2), nPE*scaleFrac);
              gamma_HitPEScaledVX      -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              gamma_HitPEScaledVXhist  -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              gamma_HitPEVX            -> Fill(xyzPos.at(0),nPE*scaleFrac);
              gamma_HitPEVXhist        -> Fill(xyzPos.at(0),nPE*scaleFrac);
            }else if( pid==1000020040 ){
              fAlpha_HitPEVXZ          -> Fill(xyzPos.at(0), xyzPos.at(2), nPE*scaleFrac);
              alpha_HitPEScaledVX      -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              alpha_HitPEScaledVXhist  -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              alpha_HitPEVX            -> Fill(xyzPos.at(0),nPE*scaleFrac);
              alpha_HitPEVXhist        -> Fill(xyzPos.at(0),nPE*scaleFrac);
            }else if ( abs(pid) == 2212 ){
              fNeutron_HitPEVXZ          -> Fill(xyzPos.at(0), xyzPos.at(2), nPE*scaleFrac);
              neutron_HitPEScaledVX      -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              neutron_HitPEScaledVXhist  -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              neutron_HitPEVX            -> Fill(xyzPos.at(0),nPE*scaleFrac);
              neutron_HitPEVXhist        -> Fill(xyzPos.at(0),nPE*scaleFrac);
            }
          }//End for(eveSDP)
        }
        if(scaleTot<1.0e-10){
          scaleTot=1.0;
        }
        //Write globals
        fHitPEScaledVX -> Fill(xyzPos.at(0), nPE/scaleTot);
        fHitPEScaledVXhist -> Fill(xyzPos.at(0), nPE/scaleTot);
        fHitPEVX -> Fill(xyzPos.at(0), nPE);
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

