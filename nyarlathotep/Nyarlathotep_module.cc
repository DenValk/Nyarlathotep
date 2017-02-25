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

/*    double modulo (double a, double b);
    double modulo (double a, int b);
    int modulo (int a, int b);
    int signF (int a);
    int signF (double a);*/

    std::string fOpHitLabel;
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
    TH2D* fPEVXZ;

    TH2D* fAlpha_PeakAmpVXZ;
    TH2D* fAlpha_IntegralVXZ;
    TH2D* fAlpha_PEVXZ;

    TH2D* fBeta_PeakAmpVXZ;
    TH2D* fBeta_IntegralVXZ;
    TH2D* fBeta_PEVXZ;

    TH2D* fGamma_PeakAmpVXZ;
    TH2D* fGamma_IntegralVXZ;
    TH2D* fGamma_PEVXZ;

    //Unscalled Histograms
    TH1D* fPeakAmpVXhist;
    TH1D* fIntegralVXhist;
    TH1D* fPEVXhist;
    TH1D* muon_PeakAmpVXhist;
    TH1D* muon_IntegralVXhist;
    TH1D* muon_PEVXhist;
    TProfile* fPeakAmpVX;
    TProfile* fIntegralVX;
    TProfile* fPEVX;
    TProfile* localTop_PeakAmpVX;
    TProfile* localTop_IntegralVX;
    TProfile* localTop_PEVX;
    TProfile* localSides_PeakAmpVX;
    TProfile* localSides_IntegralVX;
    TProfile* localSides_PEVX;
    TProfile* localAPA_PeakAmpVX;
    TProfile* localAPA_IntegralVX;
    TProfile* localAPA_PEVX;
    TProfile* localCPA_PeakAmpVX;
    TProfile* localCPA_IntegralVX;
    TProfile* localCPA_PEVX;
    TProfile* muon_PeakAmpVX;
    TProfile* muon_IntegralVX;
    TProfile* muon_PEVX;
    TProfile* alpha_PeakAmpVX;
    TH1D* alpha_PeakAmpVXhist;
    TProfile* alpha_IntegralVX;
    TH1D* alpha_IntegralVXhist;
    TProfile* alpha_PEVX;
    TH1D* alpha_PEVXhist;
    TProfile* beta_PeakAmpVX;
    TH1D* beta_PeakAmpVXhist;
    TProfile* beta_IntegralVX;
    TH1D* beta_IntegralVXhist;
    TProfile* beta_PEVX;
    TH1D* beta_PEVXhist;
    TProfile* gamma_PeakAmpVX;
    TH1D* gamma_PeakAmpVXhist;
    TProfile* gamma_IntegralVX;
    TH1D* gamma_IntegralVXhist;
    TProfile* gamma_PEVX;
    TH1D* gamma_PEVXhist;

    //Scaled histograms
    TProfile* fPeakAmpScaledVX;
    TH1D* fPeakAmpScaledVXhist;
    TProfile* fIntegralScaledVX;
    TH1D* fIntegralScaledVXhist;
    TProfile* fPEScaledVX;
    TH1D* fPEScaledVXhist;
    TProfile* localTop_PeakAmpScaledVX;
    TProfile* localTop_IntegralScaledVX;
    TProfile* localTop_PEScaledVX;
    TProfile* localSides_PeakAmpScaledVX;
    TProfile* localSides_IntegralScaledVX;
    TProfile* localSides_PEScaledVX;
    TProfile* localAPA_PeakAmpScaledVX;
    TProfile* localAPA_IntegralScaledVX;
    TProfile* localAPA_PEScaledVX;
    TProfile* localCPA_PeakAmpScaledVX;
    TProfile* localCPA_IntegralScaledVX;
    TProfile* localCPA_PEScaledVX;
    TProfile* muon_PeakAmpScaledVX;
    TH1D* muon_PeakAmpScaledVXhist;
    TProfile* muon_IntegralScaledVX;
    TH1D* muon_IntegralScaledVXhist;
    TProfile* muon_PEScaledVX;
    TH1D* muon_PEScaledVXhist;
    TProfile* alpha_PeakAmpScaledVX;
    TH1D* alpha_PeakAmpScaledVXhist;
    TProfile* alpha_IntegralScaledVX;
    TH1D* alpha_IntegralScaledVXhist;
    TProfile* alpha_PEScaledVX;
    TH1D* alpha_PEScaledVXhist;
    TProfile* beta_PeakAmpScaledVX;
    TH1D* beta_PeakAmpScaledVXhist;
    TProfile* beta_IntegralScaledVX;
    TH1D* beta_IntegralScaledVXhist;
    TProfile* beta_PEScaledVX;
    TH1D* beta_PEScaledVXhist;
    TProfile* gamma_PeakAmpScaledVX;
    TH1D* gamma_PeakAmpScaledVXhist;
    TProfile* gamma_IntegralScaledVX;
    TH1D* gamma_IntegralScaledVXhist;
    TProfile* gamma_PEScaledVX;
    TH1D* gamma_PEScaledVXhist;

  



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
    int nBinsX = int((xCryoMax - xCryoMin)/10);
    int nBinsZ = int((zCryoMax - zCryoMin)/10);
    
    //TDirectories
    art::TFileDirectory scaled   =       tfs->mkdir("scaledHistograms", "Contains the histograms rescaled by energy deposited.");
    art::TFileDirectory unscaled =       tfs->mkdir("unscaledHistograms", "Contains the histograms without rescaling by energy deposited.");
    art::TFileDirectory unscaled2D =     tfs->mkdir("unscaled2D", "Contains the 2D XZ histograms without rescaling by energy deposited.");


    fDeltaTime =              tfs->make<TH1D>    ("fDeltaTime", 
        "Plot of Times between hits. Restricted to times > 500ns", 4500, 500, 5000);
    fPeakAmpVXZ  =            unscaled2D.make<TH2D>("fPeakAmpVXZ",
        "Histogram of HitPeakAmplitude vs Z Position (cm) and X Position (cm)",
        nBinsX, xCryoMin, xCryoMax,        nBinsZ,zCryoMin,zCryoMax        );
    fPeakAmpVXZ->SetYTitle("Z Position (cm)");
    fPeakAmpVXZ->SetXTitle("X Position (cm)");

    fIntegralVXZ =            unscaled2D.make<TH2D>("fIntegralVXZ",
        "Histogram of HitIntegral vs Z Position (cm) and X Position (cm)",
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fIntegralVXZ->SetYTitle("Z Position (cm)");
    fIntegralVXZ->SetZTitle("X Position (cm)");

    fPEVXZ       =            unscaled2D.make<TH2D>("fPEVXZ",
        "Histogram of PhotoElectrons Detected vs Z Position (cm) and X Position (cm)",
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fPEVXZ->SetYTitle("Z Position (cm)");
    fPEVXZ->SetXTitle("X Position (cm)");

    fAlpha_PeakAmpVXZ  =            unscaled2D.make<TH2D>("fAlpha_PeakAmpVXZ",
        "Histogram of HitPeakAmplitude vs Z Position (cm) and X Position (cm) from #alpha's",
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fAlpha_PeakAmpVXZ->SetYTitle("Z Position (cm)");
    fAlpha_PeakAmpVXZ->SetXTitle("X Position (cm)");

    fAlpha_IntegralVXZ =            unscaled2D.make<TH2D>("fAlpha_IntegralVXZ",
        "Histogram of HitIntegral vs Z Position (cm) and X Position (cm) from #alpha's",
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax         );
    fAlpha_IntegralVXZ->SetYTitle("Z Position (cm)");
    fAlpha_IntegralVXZ->SetXTitle("X Position (cm)");

    fAlpha_PEVXZ       =            unscaled2D.make<TH2D>("fAlpha_PEVXZ",
        "Histogram of PhotoElectrons Detected vs Z Position (cm) and X Position (cm) from #alpha's",
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fAlpha_PEVXZ->SetYTitle("Z Position (cm)");
    fAlpha_PEVXZ->SetXTitle("X Position (cm)");

    fBeta_PeakAmpVXZ  =            unscaled2D.make<TH2D>("fBeta_PeakAmpVXZ",
        "Histogram of HitPeakAmplitude vs Z Position (cm) and X Position (cm) from #beta's", 
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fBeta_PeakAmpVXZ->SetYTitle("Z Position (cm)");
    fBeta_PeakAmpVXZ->SetXTitle("X Position (cm)");

    fBeta_IntegralVXZ =            unscaled2D.make<TH2D>("fBeta_IntegralVXZ",
        "Histogram of HitIntegral vs Z Position (cm) and X Position (cm) from #beta's",
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fBeta_IntegralVXZ->SetYTitle("Z Position (cm)");
    fBeta_IntegralVXZ->SetXTitle("X Position (cm)");

    fBeta_PEVXZ       =            unscaled2D.make<TH2D>("fBeta_PEVXZ",
        "Histogram of PhotoElectrons Detected vs Z Position (cm) and X Position (cm) from #beta's",
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fBeta_PEVXZ->SetYTitle("Z Position (cm)");
    fBeta_PEVXZ->SetXTitle("X Position (cm)");

    fGamma_PeakAmpVXZ  =            unscaled2D.make<TH2D>("fGamma_PeakAmpVXZ",
        "Histogram of HitPeakAmplitude vs Z Position (cm) and X Position (cm) from #gamma's",
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fGamma_PeakAmpVXZ->SetYTitle("Z Position (cm)");
    fGamma_PeakAmpVXZ->SetXTitle("X Position (cm)");

    fGamma_IntegralVXZ =            unscaled2D.make<TH2D>("fGamma_IntegralVXZ",
        "Histogram of HitIntegral vs Z Position (cm) and X Position (cm) from #gamma's",
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax         );
    fGamma_IntegralVXZ->SetYTitle("Z Position (cm)");
    fGamma_IntegralVXZ->SetXTitle("X Position (cm)");

    fGamma_PEVXZ       =            unscaled2D.make<TH2D>("fGamma_PEVXZ",
        "Histogram of PhotoElectrons Detected vs Z Position (cm) and X Position (cm) from #gamma's",
        nBinsX, xCryoMin, xCryoMax,        nBinsZ, zCryoMin, zCryoMax        );
    fGamma_PEVXZ->SetYTitle("Z Position (cm)");
    fGamma_PEVXZ->SetXTitle("X Position (cm)");

    //make scaled histograms
    fPeakAmpScaledVX =              scaled.make<TProfile>("fPeakAmpScaledVX", 
        "Plot of Hit Peak Amplitudes scaled by Track Deposited Energy vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fPeakAmpScaledVXhist =          scaled.make<TH1D>("fPeakAmpScaledVXhist", 
        "Plot of Summed Hit Peak Amplitudes vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fIntegralScaledVX =             scaled.make<TProfile>("fIntegralScaledVX", 
        "Plot of the Hit Integral Values vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fIntegralScaledVXhist =         scaled.make<TH1D>("fIntegralScaledVXhist", 
        "Plot of Summed Hit Integral Values vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fPEScaledVX =                   scaled.make<TProfile>("fPEScaledVX", 
        "Plot of the number of Detected Photo Electrons per OpHit vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fPEScaledVXhist =               scaled.make<TH1D>("fPEScaledVXhist", 
        "Plot of the number of Detected Photo Electrons per OpHit vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    localTop_PeakAmpScaledVX =      scaled.make<TProfile>("fromTop_PeakAmpScaledVX", 
        "Peak Amplitude per hit vs X (distance from APA in cm) for particles from the top of the TPC", nBinsX, xCryoMin, xCryoMax);
    localTop_IntegralScaledVX =     scaled.make<TProfile>("fromTop_IntegralScaledVX", 
        "Integral per hit vs X (distance from APA in cm) for particles from the top of the TPC", nBinsX, xCryoMin, xCryoMax);
    localTop_PEScaledVX =           scaled.make<TProfile>("fromTop_PEsVsXpos", 
        "Photo Electrons per ophit vs X (distance from APA) for particles from the top of the TPC", nBinsX, xCryoMin, xCryoMax);
    localSides_PeakAmpScaledVX =    scaled.make<TProfile>("fromSides_PeakAmpScaledVX",
        "Plot of the Hit Peak Amplitudes vs X ositions (distance from APA in cm) for particles from the sides of the TPC", nBinsX, xCryoMin, xCryoMax);
    localSides_IntegralScaledVX =   scaled.make<TProfile>("fromSides_IntegralScaledVX",
        "Integral values per hit vs X position (distance from APA in cm) for particles from the TPC sides" , nBinsX, xCryoMin, xCryoMax);
    localSides_PEScaledVX =         scaled.make<TProfile>("fromSides_PEScaledVX", 
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for particles from the TPC sides" , nBinsX, xCryoMin, xCryoMax);
    localAPA_PeakAmpScaledVX =      scaled.make<TProfile>("fromAPA_PeakAmpScaledVX", 
        "Peak Amp per hit vs X position (distance from APA in cm) for particles from the APA", nBinsX, xCryoMin, xCryoMax);
    localAPA_IntegralScaledVX =     scaled.make<TProfile>("fromAPA_IntegralScaledVX",
        "Integral per hit Vs X position (distance from APA in cm) for particles from the APA", nBinsX, xCryoMin, xCryoMax);
    localAPA_PEScaledVX =           scaled.make<TProfile>("fromAPA_PEScaledVX",
        "PhotoElectrons per OpHit Vs X position (distance from APA in cm) for particles from the APA", nBinsX, xCryoMin, xCryoMax);
    muon_PeakAmpScaledVX =          scaled.make<TProfile>("muon_PeakAmpScaledVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_PeakAmpScaledVXhist =      scaled.make<TH1D>("muon_PeakAmpScaledVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_IntegralScaledVX =         scaled.make<TProfile>("muon_IntegralScaledVX",
        "Integral per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_IntegralScaledVXhist =     scaled.make<TH1D>("muon_IntegralScaledVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_PEScaledVX =               scaled.make<TProfile>("muon_PEScaledVX",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_PEScaledVXhist =           scaled.make<TH1D>("muon_PEScaledVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    alpha_PeakAmpScaledVX =         scaled.make<TProfile>("alpha_PeakAmpScaledVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_PeakAmpScaledVXhist =     scaled.make<TH1D>("alpha_PeakAmpScaledVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_IntegralScaledVX =        scaled.make<TProfile>("alpha_IntegralScaledVX",
        "Integral per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_IntegralScaledVXhist =    scaled.make<TH1D>("alpha_IntegralScaledVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_PEScaledVX =              scaled.make<TProfile>("alpha_PEScaledVX",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_PEScaledVXhist =          scaled.make<TH1D>("alpha_PEScaledVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    beta_PeakAmpScaledVX =          scaled.make<TProfile>("beta_PeakAmpScaledVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_PeakAmpScaledVXhist =      scaled.make<TH1D>("beta_PeakAmpScaledVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_IntegralScaledVX =         scaled.make<TProfile>("beta_IntegralScaledVX",
        "Integral per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_IntegralScaledVXhist =     scaled.make<TH1D>("beta_IntegralScaledVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_PEScaledVX =               scaled.make<TProfile>("beta_PEScaledVX",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_PEScaledVXhist =           scaled.make<TH1D>("beta_PEScaledVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    gamma_PeakAmpScaledVX =         scaled.make<TProfile>("gamma_PeakAmpScaledVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_PeakAmpScaledVXhist =     scaled.make<TH1D>("gamma_PeakAmpScaledVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_IntegralScaledVX =        scaled.make<TProfile>("gamma_IntegralScaledVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_IntegralScaledVXhist =    scaled.make<TH1D>("gamma_IntegralScaledVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_PEScaledVX =              scaled.make<TProfile>("gamma_PEScaledVX",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_PEScaledVXhist =          scaled.make<TH1D>("gamma_PEScaledVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);

    

    //Make unscaled Histograms
    fPeakAmpVX =              unscaled.make<TProfile>("fPeakAmpVX", 
        "Plot of Hit Peak Amplitudes vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fPeakAmpVXhist =          unscaled.make<TH1D>("fPeakAmpVXhist", 
        "Plot of Summed Hit Peak Amplitudes vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fIntegralVX =             unscaled.make<TProfile>("fIntegralVX", 
        "Plot of the Hit Integral Values vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fIntegralVXhist =         unscaled.make<TH1D>("fIntegralVXhist", 
        "Plot of Summed Hit Integral Values vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fPEVX =                   unscaled.make<TProfile>("fPEVX", 
        "Plot of the number of Detected Photo Electrons per OpHit vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    fPEVXhist =               unscaled.make<TH1D>("fPEVXhist", 
        "Plot of the number of Detected Photo Electrons per OpHit vs X position (Distance from APA in cm)", nBinsX, xCryoMin, xCryoMax);
    localTop_PeakAmpVX =      unscaled.make<TProfile>("fromTop_PeakAmpVX", 
        "Peak Amplitude per hit vs X (distance from APA in cm) for particles from the top of the TPC", nBinsX, xCryoMin, xCryoMax);
    localTop_IntegralVX =     unscaled.make<TProfile>("fromTop_IntegralVX", 
        "Integral per hit vs X (distance from APA in cm) for particles from the top of the TPC", nBinsX, xCryoMin, xCryoMax);
    localTop_PEVX =           unscaled.make<TProfile>("fromTop_PEsVsXpos", 
        "Photo Electrons per ophit vs X (distance from APA) for particles from the top of the TPC", nBinsX, xCryoMin, xCryoMax);
    localSides_PeakAmpVX =    unscaled.make<TProfile>("fromSides_PeakAmpVX",
        "Plot of the Hit Peak Amplitudes vs X ositions (distance from APA in cm) for particles from the sides of the TPC", nBinsX, xCryoMin, xCryoMax);
    localSides_IntegralVX =   unscaled.make<TProfile>("fromSides_IntegralVX",
        "Integral values per hit vs X position (distance from APA in cm) for particles from the TPC sides" , nBinsX, xCryoMin, xCryoMax);
    localSides_PEVX =         unscaled.make<TProfile>("fromSides_PEVX", 
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for particles from the TPC sides" , nBinsX, xCryoMin, xCryoMax);
    localAPA_PeakAmpVX =      unscaled.make<TProfile>("fromAPA_PeakAmpVX", 
        "Peak Amp per hit vs X position (distance from APA in cm) for particles from the APA", nBinsX, xCryoMin, xCryoMax);
    localAPA_IntegralVX =     unscaled.make<TProfile>("fromAPA_IntegralVX",
        "Integral per hit Vs X position (distance from APA in cm) for particles from the APA", nBinsX, xCryoMin, xCryoMax);
    localAPA_PEVX =           unscaled.make<TProfile>("fromAPA_PEVX",
        "PhotoElectrons per OpHit Vs X position (distance from APA in cm) for particles from the APA", nBinsX, xCryoMin, xCryoMax);
    muon_PeakAmpVX =          unscaled.make<TProfile>("muon_PeakAmpVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_PeakAmpVXhist =      unscaled.make<TH1D>("muon_PeakAmpVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_IntegralVX =         unscaled.make<TProfile>("muon_IntegralVX",
        "Integral per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_IntegralVXhist =     unscaled.make<TH1D>("muon_IntegralVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_PEVX =               unscaled.make<TProfile>("muon_PEVX",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    muon_PEVXhist =           unscaled.make<TH1D>("muon_PEVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Muons", nBinsX, xCryoMin, xCryoMax);
    alpha_PeakAmpVX =         unscaled.make<TProfile>("alpha_PeakAmpVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_PeakAmpVXhist =     unscaled.make<TH1D>("alpha_PeakAmpVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_IntegralVX =        unscaled.make<TProfile>("alpha_IntegralVX",
        "Integral per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_IntegralVXhist =    unscaled.make<TH1D>("alpha_IntegralVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_PEVX =              unscaled.make<TProfile>("alpha_PEVX",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    alpha_PEVXhist =          unscaled.make<TH1D>("alpha_PEVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Alphas", nBinsX, xCryoMin, xCryoMax);
    beta_PeakAmpVX =          unscaled.make<TProfile>("beta_PeakAmpVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_PeakAmpVXhist =      unscaled.make<TH1D>("beta_PeakAmpVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_IntegralVX =         unscaled.make<TProfile>("beta_IntegralVX",
        "Integral per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_IntegralVXhist =     unscaled.make<TH1D>("beta_IntegralVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_PEVX =               unscaled.make<TProfile>("beta_PEVX",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    beta_PEVXhist =           unscaled.make<TH1D>("beta_PEVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for Betas", nBinsX, xCryoMin, xCryoMax);
    gamma_PeakAmpVX =         unscaled.make<TProfile>("gamma_PeakAmpVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_PeakAmpVXhist =     unscaled.make<TH1D>("gamma_PeakAmpVXhist",
        "Peak Amp per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_IntegralVX =        unscaled.make<TProfile>("gamma_IntegralVX",
        "Peak Amp per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_IntegralVXhist =    unscaled.make<TH1D>("gamma_IntegralVXhist",
        "Integral per Hit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_PEVX =              unscaled.make<TProfile>("gamma_PEVX",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);
    gamma_PEVXhist =          unscaled.make<TH1D>("gamma_PEVXhist",
        "PhotoElectrons per OpHit vs X position (distance from APA in cm) for gammas", nBinsX, xCryoMin, xCryoMax);

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
    std::vector< art::Ptr< recob::Hit > > hitlist;
    if (evt.getByLabel(fHitLabel, hitHandle) )
      art::fill_ptr_vector(hitlist, hitHandle);

    art::Handle< std::vector< recob::OpHit > > ophitHandle;
    std::vector< art::Ptr< recob::OpHit > > ophitlist;
    if (evt.getByLabel(fOpHitLabel, ophitHandle) )
      art::fill_ptr_vector(ophitlist, ophitHandle);

    //Call Charge events.
//    for(art::Ptr<recob::Hit> const& ptrHit: hitlist){
    for( auto ptrHit: hitlist){
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

    for(const art::Ptr<recob::OpHit>& ptrOpHit: ophitlist){
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
              localSides_PEScaledVX  ->Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              localSides_PEVX  ->Fill(xyzPos.at(0),nPE*scaleFrac);
            }else if( (evePos.at(1)>= yTpcMax - edgeDelta && evePos.at(1) <= yTpcMax + edgeDelta) ){ //Top
              localTop_PEScaledVX    ->Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              localTop_PEVX    ->Fill(xyzPos.at(0),nPE*scaleFrac);
            }else if( (evePos.at(0) >= 0-edgeDelta && evePos.at(0) <= 0+edgeDelta ) ){ //APA
              localAPA_PEScaledVX    ->Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              localAPA_PEVX    ->Fill(xyzPos.at(0),nPE*scaleFrac);
            }else if( (evePos.at(0) >= xTpcMax-edgeDelta && evePos.at(0) <= xTpcMax+edgeDelta ) ){ //CPA
              localCPA_PEScaledVX    ->Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              localCPA_PEVX    ->Fill(xyzPos.at(0),nPE*scaleFrac);
            }
  
            if( abs(pid)==13 ){
              muon_PEScaledVX       -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              muon_PEScaledVXhist   -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              muon_PEVX             -> Fill(xyzPos.at(0),nPE*scaleFrac);
              muon_PEVXhist         -> Fill(xyzPos.at(0),nPE*scaleFrac);
            }else if( abs(pid)==11 ){
              fBeta_PEVXZ           -> Fill(xyzPos.at(0), xyzPos.at(2), nPE*scaleFrac);
              beta_PEScaledVX       -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              beta_PEScaledVXhist   -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              beta_PEVX             -> Fill(xyzPos.at(0),nPE*scaleFrac);
              beta_PEVXhist         -> Fill(xyzPos.at(0),nPE*scaleFrac);
            }else if( abs(pid)==22){
              fGamma_PEVXZ          -> Fill(xyzPos.at(0), xyzPos.at(2), nPE*scaleFrac);
              gamma_PEScaledVX      -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              gamma_PEScaledVXhist  -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              gamma_PEVX            -> Fill(xyzPos.at(0),nPE*scaleFrac);
              gamma_PEVXhist        -> Fill(xyzPos.at(0),nPE*scaleFrac);
            }else if( pid==1000020040 ){
              fAlpha_PEVXZ          -> Fill(xyzPos.at(0), xyzPos.at(2), nPE*scaleFrac);
              alpha_PEScaledVX      -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              alpha_PEScaledVXhist  -> Fill(xyzPos.at(0),nPE*scaleFrac/scaleEn);
              alpha_PEVX            -> Fill(xyzPos.at(0),nPE*scaleFrac);
              alpha_PEVXhist        -> Fill(xyzPos.at(0),nPE*scaleFrac);
            }
          }//End for(eveSDP)
        }
        if(scaleTot<1.0e-10){
          scaleTot=1.0;
        }
        //Write globals
        fPEScaledVX -> Fill(xyzPos.at(0), nPE/scaleTot);
        fPEScaledVXhist -> Fill(xyzPos.at(0), nPE/scaleTot);
        fPEVX -> Fill(xyzPos.at(0), nPE);
        fPEVXhist -> Fill(xyzPos.at(0), nPE);

        fPEVXZ->Fill(xyzPos.at(0), xyzPos.at(2), nPE);

      }//end try
      catch(...){//Lazy catch. Need to catch -No sim SDPs found- error.
   //     throw;
      }
    }//End for(OpHit)
    

  } // Nyarlathotep::analyze()
  
  
  // This macro has to be defined for this module to be invoked from a
  // .fcl file; see Nyarlathotep.fcl for more information.
  DEFINE_ART_MODULE(Nyarlathotep)



} // namespace Nyarlathotep

#endif // Nyarlathotep_Module

