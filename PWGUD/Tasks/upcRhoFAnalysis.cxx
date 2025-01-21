// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief  task for analysis of rho in UPCs using UD tables (from SG producer)
///         includes event tagging based on ZN information, track selection, reconstruction,
///         and also some basic stuff for decay phi anisotropy studies
/// \author Cesar, cesar.ramirez@cern.ch

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h" // similiar to TLorentzVector (which is now legacy apparently)
#include "random"

#include "Common/DataModel/PIDResponse.h"

#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using FullUDSgCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced, aod::SGCollisions>::iterator;
using FullUDTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksFlags>;



namespace o2::aod
{
namespace fourpi
{



//for event
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);

//for Rhos
DECLARE_SOA_COLUMN(M, m, double);
DECLARE_SOA_COLUMN(Pt, pt, double);
DECLARE_SOA_COLUMN(Eta, eta, double);
DECLARE_SOA_COLUMN(Phi, phi, double);

//for vertex
DECLARE_SOA_COLUMN(PosX, posX, double);
DECLARE_SOA_COLUMN(PosY, posY, double);
DECLARE_SOA_COLUMN(PosZ, posZ, double);

//for other
//DECLARE_SOA_COLUMN(TotalCharge, totalCharge, int);

//for pion tracks
//DECLARE_SOA_COLUMN(TrackSign, trackSign, std::vector<int>);
//DECLARE_SOA_COLUMN(TrackM, trackM, std::vector<double>);
//DECLARE_SOA_COLUMN(TrackPt, trackPt, std::vector<double>);
//DECLARE_SOA_COLUMN(TrackEta, trackEta, std::vector<double>);
//DECLARE_SOA_COLUMN(TrackPhi, trackPhi, std::vector<double>);



} // namespace fourpi
DECLARE_SOA_TABLE(SYSTEMTREE, "AOD", "SystemTree", fourpi::RunNumber, fourpi::M, fourpi::Pt, fourpi::Eta, fourpi::Phi)
                  fourpi::PosX, fourpi::PosY, fourpi::PosZ); //, fourpi::TotalCharge); //fourpi::TrackSign,
                  //fourpi::TrackM, fourpi::TrackPt, fourpi::TrackEta, fourpi::TrackPhi);
} // namespace o2::aod

struct upcRhoFAnalysis {
  Produces<aod::SYSTEMTREE>systemTree;

  double PcEtaCut = 0.9; // physics coordination recommendation

  Configurable<bool> specifyGapSide{"specifyGapSide", true, "specify gap side for SG/DG produced data"};
  Configurable<int> gapSide{"gapSide", 2, "gap side for SG produced data"};
  Configurable<bool> requireTof{"requireTof", false, "require TOF signal"};

  Configurable<double> collisionsPosZMaxCut{"collisionsPosZMaxCut", 10.0, "max Z position cut on collisions"};
  Configurable<double> ZNcommonEnergyCut{"ZNcommonEnergyCut", 0.0, "ZN common energy cut"};
  Configurable<double> ZNtimeCut{"ZNtimeCut", 2.0, "ZN time cut"};

  Configurable<double> tracksTpcNSigmaPiCut{"tracksTpcNSigmaPiCut", 3.0, "TPC nSigma pion cut"};
  Configurable<double> tracksDcaMaxCut{"tracksDcaMaxCut", 1.0, "max DCA cut on tracks"};

  Configurable<double> systemMassMinCut{"systemMassMinCut", 0.5, "min M cut for reco system"};
  Configurable<double> systemMassMaxCut{"systemMassMaxCut", 1.2, "max M cut for reco system"};
  Configurable<double> systemPtCut{"systemPtMaxCut", 0.1, "max pT cut for reco system"};
  Configurable<double> systemYCut{"systemYCut", 0.9, "rapiditiy cut for reco system"};

  ConfigurableAxis mAxis{"mAxis", {1000, 0.0, 10.0}, "m (GeV/#it{c}^{2})"};
  ConfigurableAxis mCutAxis{"mCutAxis", {70, 0.5, 1.2}, "m (GeV/#it{c}^{2})"};
  ConfigurableAxis ptAxis{"ptAxis", {1000, 0.0, 10.0}, "p_{T} (GeV/#it{c})"};
  ConfigurableAxis ptCutAxis{"ptCutAxis", {300, 0.0, 0.3}, "p_{T} (GeV/#it{c})"};
  ConfigurableAxis pt2Axis{"pt2Axis", {300, 0.0, 0.09}, "p_{T}^{2} (GeV^{2}/#it{c}^{2})"};
  ConfigurableAxis etaAxis{"etaAxis", {180, -0.9, 0.9}, "#eta"};
  ConfigurableAxis yAxis{"yAxis", {180, -0.9, 0.9}, "y"};
  ConfigurableAxis phiAxis{"phiAxis", {180, 0.0, o2::constants::math::TwoPI}, "#phi"};
  ConfigurableAxis phiAsymmAxis{"phiAsymmAxis", {182, -o2::constants::math::PI, o2::constants::math::PI}, "#phi"};
  ConfigurableAxis momentumFromPhiAxis{"momentumFromPhiAxis", {400, -0.1, 0.1}, "p (GeV/#it{c})"};
  ConfigurableAxis ptQuantileAxis{"ptQuantileAxis", {0, 0.0181689, 0.0263408, 0.0330488, 0.0390369, 0.045058, 0.0512604, 0.0582598, 0.066986, 0.0788085, 0.1}, "p_{T} (GeV/#it{c})"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    // selection counter
    std::vector<std::string> selectionCounterLabels = {"all tracks", "PV contributor", "ITS + TPC hit", "TOF requirement", "DCA cut", "#eta cut", "2D TPC n#sigma_{#pi} cut"};

    // 4PI SYSTEM
    registry.add("4pi/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("4pi/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    //registry.add("4pi/hEta", ";Eta (1);counts", kTH1D, {etaAxis});
    //registry.add("4pi/hPhi", ";Phi ();counts", kTH1D, {phiAxis});
  }

  template <typename C>
  bool collisionPassesCuts(const C& collision) // collision cuts
  {
    if (std::abs(collision.posZ()) > collisionsPosZMaxCut)
      return false;
    if (specifyGapSide && collision.gapSide() != gapSide)
      return false;
    return true;
  }

  template <typename T>
  bool trackPassesCuts(const T& track) // track cuts (PID done separately)
  {
    if (!track.isPVContributor())
      return false;
    //registry.fill(HIST("QC/tracks/hSelectionCounter"), 1);
    if (!track.hasITS() || !track.hasTPC())
      return false;
    //registry.fill(HIST("QC/tracks/hSelectionCounter"), 2);
    if (requireTof && !track.hasTOF())
      return false;
    //registry.fill(HIST("QC/tracks/hSelectionCounter"), 3);
    if (std::abs(track.dcaZ()) > tracksDcaMaxCut || std::abs(track.dcaXY()) > (0.0182 + 0.0350 / std::pow(track.pt(), 1.01))) // Run 2 dynamic DCA cut
      return false;
    //registry.fill(HIST("QC/tracks/hSelectionCounter"), 4);
    if (std::abs(eta(track.px(), track.py(), track.pz())) > PcEtaCut)
      return false;
    //registry.fill(HIST("QC/tracks/hSelectionCounter"), 5);
    return true;
  }

  template <typename T>
  bool tracksPassPiPID(const T& cutTracks) // n-dimensional PID cut
  {
    double radius = 0.0;
    for (const auto& track : cutTracks)
      radius += std::pow(track.tpcNSigmaPi(), 2);
    return radius < std::pow(tracksTpcNSigmaPiCut, 2);
  }

  template <typename T>
  double tracksTotalCharge(const T& cutTracks) // total charge of selected tracks
  {
    double charge = 0.0;
    for (const auto& track : cutTracks)
      charge += track.sign();
    return charge;
  }

  bool systemPassCuts(const ROOT::Math::PxPyPzMVector& system) // system cuts
  {
    if (system.M() < systemMassMinCut || system.M() > systemMassMaxCut)
      return false;
    if (system.Pt() > systemPtCut)
      return false;
    if (std::abs(system.Rapidity()) > systemYCut)
      return false;
    return true;
  }

  ROOT::Math::PxPyPzMVector reconstructSystem(const std::vector<ROOT::Math::PxPyPzMVector>& cutTracks4Vecs) // reconstruct system from 4-vectors
  {
    ROOT::Math::PxPyPzMVector system;
    for (const auto& track4Vec : cutTracks4Vecs)
      system += track4Vec;
    return system;
  }

  double deltaPhi(const ROOT::Math::PxPyPzMVector& p1, const ROOT::Math::PxPyPzMVector& p2)
  {
    double dPhi = p1.Phi() - p2.Phi();
    if (dPhi > o2::constants::math::PI)
      dPhi -= o2::constants::math::TwoPI;
    else if (dPhi < -o2::constants::math::PI)
      dPhi += o2::constants::math::TwoPI;
    return dPhi; // calculate delta phi in (-pi, pi)
  }

  double getPhiRandom(const std::vector<ROOT::Math::PxPyPzMVector>& cutTracks4Vecs) // decay phi anisotropy
  {                                                                                 // two possible definitions of phi: randomize the tracks
    std::vector<int> indices = {0, 1};
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();    // get time-based seed
    std::shuffle(indices.begin(), indices.end(), std::default_random_engine(seed)); // shuffle indices
    // calculate phi
    ROOT::Math::PxPyPzMVector pOne = cutTracks4Vecs[indices[0]];
    ROOT::Math::PxPyPzMVector pTwo = cutTracks4Vecs[indices[1]];
    ROOT::Math::PxPyPzMVector pPlus = pOne + pTwo;
    ROOT::Math::PxPyPzMVector pMinus = pOne - pTwo;
    return deltaPhi(pPlus, pMinus);
  }

  template <typename T>
  double getPhiCharge(const T& cutTracks, const std::vector<ROOT::Math::PxPyPzMVector>& cutTracks4Vecs)
  { // two possible definitions of phi: charge-based assignment
    ROOT::Math::PxPyPzMVector pOne, pTwo;
    if (cutTracks[0].sign() > 0) {
      pOne = cutTracks4Vecs[0];
      pTwo = cutTracks4Vecs[1];
    } else {
      pOne = cutTracks4Vecs[1];
      pTwo = cutTracks4Vecs[0];
    }
    ROOT::Math::PxPyPzMVector pPlus = pOne + pTwo;
    ROOT::Math::PxPyPzMVector pMinus = pOne - pTwo;
    return deltaPhi(pPlus, pMinus);
  }

  void processReco(FullUDSgCollision const& collision, FullUDTracks const& tracks)
  {
    

    if (!collisionPassesCuts(collision))
      return;

    // event tagging
    bool XnXn = false, OnOn = false, XnOn = false, OnXn = false; // note: On == 0n...
    if (collision.energyCommonZNA() < ZNcommonEnergyCut && collision.energyCommonZNC() < ZNcommonEnergyCut)
      OnOn = true;
    if (collision.energyCommonZNA() > ZNcommonEnergyCut && std::abs(collision.timeZNA()) < ZNtimeCut &&
        collision.energyCommonZNC() > ZNcommonEnergyCut && std::abs(collision.timeZNC()) < ZNtimeCut)
      XnXn = true;
    if (collision.energyCommonZNA() > ZNcommonEnergyCut && std::abs(collision.timeZNA()) < ZNtimeCut && collision.energyCommonZNC() < ZNcommonEnergyCut)
      XnOn = true;
    if (collision.energyCommonZNA() < ZNcommonEnergyCut && collision.energyCommonZNC() > ZNcommonEnergyCut && std::abs(collision.timeZNC()) < ZNtimeCut)
      OnXn = true;
    // vectors for storing selected tracks and their 4-vectors
    std::vector<decltype(tracks.begin())> cutTracks;
    std::vector<ROOT::Math::PxPyPzMVector> cutTracks4Vecs;

    int trackCounter = 0;
    for (const auto& track : tracks) {
   

      if (!trackPassesCuts(track))
        continue;
      trackCounter++;
      cutTracks.push_back(track);
      cutTracks4Vecs.push_back(ROOT::Math::PxPyPzMVector(track.px(), track.py(), track.pz(), o2::constants::physics::MassPionCharged)); // apriori assume pion mass
      
    }


    if (!tracksPassPiPID(cutTracks))
      return;
    // reonstruct system and calculate total charge, save commonly used values into variables
    ROOT::Math::PxPyPzMVector system = reconstructSystem(cutTracks4Vecs);
    int totalCharge = tracksTotalCharge(cutTracks);
    int nTracks = cutTracks.size();
    double mass = system.M();
    double pT = system.Pt();
    double pTsquare = pT * pT;
    double rapidity = system.Rapidity();
    double systemPhi = system.Phi() + o2::constants::math::PI;

    if (nTracks == 4 && tracksTotalCharge(cutTracks) == 0) { // 4pi system
      systemTree(collision.runNumber(), mass, pT, rapidity, systemPhi, collision.posX(), collision.posY(), collision.posZ());
      registry.fill(HIST("4pi/hM"), mass);
      registry.fill(HIST("4pi/hPt"), pT);
      //registry.fill(HIST("4pi/hEta"), rapiditiy);
      //registry.fill(HIST("4pi/hPhi"), system);
    } 
    //std::cout<<"Hola Mundo"<<std::endl;
  }
  PROCESS_SWITCH(upcRhoFAnalysis, processReco, "analyse reco tracks", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    o2::framework::adaptAnalysisTask<upcRhoFAnalysis>(cfgc)};
}