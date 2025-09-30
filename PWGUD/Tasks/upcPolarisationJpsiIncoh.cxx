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

/// \file upcPolarisationJpsiIncoh.cxx
/// \brief Workflow to analyse UPC forward events and perform J/psi polarization selections
/// \author Niveditha Ram, IP2I <niv.ram@cern.ch>
/// \ingroup PWGUD
/// executable name: o2-analysis-ud-upc-polarisation-jpsiincoh

#include "PWGUD/DataModel/UDTables.h"

#include "Common/Core/RecoDecay.h"


#include "CommonConstants/PhysicsConstants.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TSystem.h"

#include <unordered_map>
#include <vector>

using namespace ROOT::Math;

// table for saving tree with info on data
namespace dimu
{
// dimuon
  DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
  DECLARE_SOA_COLUMN(M, m, float);
  DECLARE_SOA_COLUMN(Energy, energy, float);
  DECLARE_SOA_COLUMN(Px, px, float);
  DECLARE_SOA_COLUMN(Py, py, float);
  DECLARE_SOA_COLUMN(Pz, pz, float);
  DECLARE_SOA_COLUMN(Pt, pt, float);
  DECLARE_SOA_COLUMN(Rap, rap, float);
  DECLARE_SOA_COLUMN(Phi, phi, float);

  //tracks positive (p) and negative (n) stored
  DECLARE_SOA_COLUMN(EnergyP, energyP, float);
  DECLARE_SOA_COLUMN(Pxp, pxp, float);
  DECLARE_SOA_COLUMN(Pyp, pyp, float);
  DECLARE_SOA_COLUMN(Pzp, pzp, float);
  DECLARE_SOA_COLUMN(Ptp, ptp, float);
  DECLARE_SOA_COLUMN(Etap, etap, float);
  DECLARE_SOA_COLUMN(Phip, phip, float);

  DECLARE_SOA_COLUMN(EnergyN, energyN, float);
  DECLARE_SOA_COLUMN(Pxn, pxn, float);
  DECLARE_SOA_COLUMN(Pyn, pyn, float);
  DECLARE_SOA_COLUMN(Pzn, pzn, float);
  DECLARE_SOA_COLUMN(Ptn, ptn, float);
  DECLARE_SOA_COLUMN(Etan, etan, float);
  DECLARE_SOA_COLUMN(Phin, phin, float);
  
} // namespace dimu

namespace o2::aod
{
  DECLARE_SOA_TABLE(DiMu, "AOD", "DIMU",
                  dimu::RunNumber,
                  dimu::M, dimu::Energy, dimu::Px, dimu::Py, dimu::Pz, dimu::Pt, dimu::Rap, dimu::Phi,
		    dimu::EnergyP, dimu::Pxp, dimu::Pyp, dimu::Pzp, dimu::Ptp, dimu::Etap, dimu::Phip,
		    dimu::EnergyN, dimu::Pxn, dimu::Pyn, dimu::Pzn, dimu::Ptn, dimu::Etan, dimu::Phin);


} // namespace o2::aod
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// constants used in the track selection
const float kRAbsMin = 17.6;
const float kRAbsMax = 89.5;
const float kPDca = 200.;
float kEtaMin = -4.0;
float kEtaMax = -2.5;
const float kPtMin = 0.;
const float kMaxAmpV0A = 100.;
const int kReqMatchMIDTracks = 2;
const int kReqMatchMFTTracks = 2;
const int kMaxChi2MFTMatch = 30;
const float kMaxZDCTime = 2.;
const float kMaxZDCTimeHisto = 10.;



struct UpcPolarisationJpsiIncoh {

  using CandidatesFwd = soa::Join<o2::aod::UDCollisions, o2::aod::UDCollisionsSelsFwd>;
  using ForwardTracks = soa::Join<o2::aod::UDFwdTracks, o2::aod::UDFwdTracksExtra>;
  using CompleteFwdTracks = soa::Join<ForwardTracks, o2::aod::UDMcFwdTrackLabels>;

  Produces<o2::aod::DiMu> dimuSel;
  // defining histograms using histogram registry: different histos for the different process functions
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  // CONFIGURABLES
  static constexpr double Pi = o2::constants::math::PI;
  // pT of muon pairs
  Configurable<int> nBinsPt{"nBinsPt", 250, "N bins in pT histo"};
  Configurable<float> lowPt{"lowPt", 0., "lower limit in pT histo"};
  Configurable<float> highPt{"highPt", 2, "upper limit in pT histo"};
  // mass of muon pairs
  Configurable<int> nBinsMass{"nBinsMass", 500, "N bins in mass histo"};
  Configurable<float> lowMass{"lowMass", 0., "lower limit in mass histo"};
  Configurable<float> highMass{"highMass", 10., "upper limit in mass histo"};
  // eta of muon pairs
  Configurable<int> nBinsEta{"nBinsEta", 600, "N bins in eta histo"};
  Configurable<float> lowEta{"lowEta", -10., "lower limit in eta histo"};
  Configurable<float> highEta{"highEta", -2., "upper limit in eta histo"};
  // rapidity of muon pairs
  Configurable<int> nBinsRapidity{"nBinsRapidity", 250, "N bins in rapidity histo"};
  Configurable<float> lowRapidity{"lowRapidity", -4.5, "lower limit in rapidity histo"};
  Configurable<float> highRapidity{"highRapidity", -2., "upper limit in rapidity histo"};
  // phi of muon pairs
  Configurable<int> nBinsPhi{"nBinsPhi", 600, "N bins in phi histo"};
  Configurable<float> lowPhi{"lowPhi", -Pi, "lower limit in phi histo"};
  Configurable<float> highPhi{"highPhi", Pi, "upper limit in phi histo"};
 // pT of single muons                                                                                                    
  Configurable<int> nBinsPtSingle{"nBinsPtSingle", 500, "N bins in pT histo single muon"};
  Configurable<float> lowPtSingle{"lowPtSingle", 0., "lower limit in pT histo single muon"};
  Configurable<float> highPtSingle{"highPtSingle", 2., "upper limit in pT histo single muon"};
  // eta of single muons                                                                                                   
  Configurable<int> nBinsEtaSingle{"nBinsEtaSingle", 250, "N bins in eta histo single muon"};
  Configurable<float> lowEtaSingle{"lowEtaSingle", -4.5, "lower limit in eta histo single muon"};
  Configurable<float> highEtaSingle{"highEtaSingle", -2., "upper limit in eta histo single muon"};
  // phi of single muons                                                                                                   
  Configurable<int> nBinsPhiSingle{"nBinsPhiSingle", 600, "N bins in phi histo single muon"};
  Configurable<float> lowPhiSingle{"lowPhiSingle", -Pi, "lower limit in phi histo single muon"};
  Configurable<float> highPhiSingle{"highPhiSingle", Pi, "upper limit in phi histo single muon"};
  // ZDC                                                                                                                   
  Configurable<int> nBinsZDCen{"nBinsZDCen", 200, "N bins in ZN energy"};
  Configurable<float> lowEnZN{"lowEnZN", -50., "lower limit in ZN energy histo"};
  Configurable<float> highEnZN{"highEnZN", 250., "upper limit in ZN energy histo"};
  // Analysis cuts
  Configurable<float> minJpsiMass{"minJpsiMass", 3.0, "Minimum of the jpsi peak for peak cut"};
  Configurable<float> maxJpsiMass{"maxJpsiMass", 3.18, "Maximum of the jpsi peak for peak cut"};

  // my track type
  // 0 = MCH-MID-MFT
  // 1 = MCH-MID
  Configurable<int> myTrackType{"myTrackType", 1, "My track type"};

  void init(InitContext&)
  {
    // axis
    const AxisSpec axisPt{nBinsPt, lowPt, highPt, "#it{p}_{T} GeV/#it{c}"};
    const AxisSpec axisMass{nBinsMass, lowMass, highMass, "m_{#mu#mu} GeV/#it{c}^{2}"};
    const AxisSpec axisEta{nBinsEta, lowEta, highEta, "#eta"};
    const AxisSpec axisRapidity{nBinsRapidity, lowRapidity, highRapidity, "Rapidity"};
    const AxisSpec axisPhi{nBinsPhi, lowPhi, highPhi, "#varphi"};
    const AxisSpec axisPtSingle{nBinsPtSingle, lowPtSingle, highPtSingle, "#it{p}_{T}_{ trk} GeV/#it{c}"};
    const AxisSpec axisEtaSingle{nBinsEtaSingle, lowEtaSingle, highEtaSingle, "#eta_{trk}"};
    const AxisSpec axisPhiSingle{nBinsPhiSingle, lowPhiSingle, highPhiSingle, "#varphi_{trk}"};

    // histos
    // data and reco MC
    registry.add("hMass", "Invariant mass of muon pairs;;#counts", kTH1D, {axisMass});
    registry.add("hPt", "Transverse momentum of muon pairs;;#counts", kTH1D, {axisPt});
    registry.add("hEta", "Pseudorapidty of muon pairs;;#counts", kTH1D, {axisEta});
    registry.add("hRapidity", "Rapidty of muon pairs;;#counts", kTH1D, {axisRapidity});
    registry.add("hPhi", "#varphi of muon pairs;;#counts", kTH1D, {axisPhi});
    registry.add("hPtTrkPos", "Pt of positive muons;;#counts", kTH1D, {axisPtSingle});
    registry.add("hPtTrkNeg", "Pt of negative muons;;#counts", kTH1D, {axisPtSingle});
    registry.add("hEtaTrkPos", "#eta of positive muons;;#counts", kTH1D, {axisEtaSingle});
    registry.add("hEtaTrkNeg", "#eta of negative muons;;#counts", kTH1D, {axisEtaSingle});
    registry.add("hPhiTrkPos", "#varphi of positive muons;;#counts", kTH1D, {axisPhiSingle});
    registry.add("hPhiTrkNeg", "#varphi of negative muons;;#counts", kTH1D, {axisPhiSingle});


  }

  // template function that fills a map with the collision id of each udcollision as key
  // and a vector with the tracks
  // map == (key, element) == (udCollisionId, vector of trks)
  template <typename TTracks>
  void collectCandIDs(std::unordered_map<int32_t, std::vector<int32_t>>& tracksPerCand, TTracks& tracks)
  {
    for (const auto& tr : tracks) {
      int32_t candId = tr.udCollisionId();
      if (candId < 0) {
        continue;
      }
      tracksPerCand[candId].push_back(tr.globalIndex());
    }
  }

  // template function that fills a map with the collision id of each udmccollision as key
  // and a vector with the tracks
  // map == (key, element) == (udMcCollisionId, vector of mc particles)
  template <typename TTracks>
  void collectMcCandIDs(std::unordered_map<int32_t, std::vector<int32_t>>& tracksPerCand, TTracks& tracks)
  {
    for (const auto& tr : tracks) {
      int32_t candId = tr.udMcCollisionId();
      if (candId < 0) {
        continue;
      }
      tracksPerCand[candId].push_back(tr.globalIndex());
    }
  }

  // struct used to store the ZDC info in a map
  struct ZDCinfo {
    float timeA;
    float timeC;
    float enA;
    float enC;
    int32_t id;
  };

  // function that fills a map with the collision id of each udcollision as key
  // and a ZDCinfo struct with the ZDC information
  void collectCandZDCInfo(std::unordered_map<int32_t, ZDCinfo>& zdcPerCand, o2::aod::UDZdcsReduced const& ZDCs)
  {

    for (const auto& zdc : ZDCs) {
      int32_t candId = zdc.udCollisionId();
      if (candId < 0) {
        continue;
      }

      zdcPerCand[candId].timeA = zdc.timeZNA();
      zdcPerCand[candId].timeC = zdc.timeZNC();
      zdcPerCand[candId].enA = zdc.energyCommonZNA();
      zdcPerCand[candId].enC = zdc.energyCommonZNC();

      // take care of the infinity
      if (std::isinf(zdcPerCand[candId].timeA))
        zdcPerCand[candId].timeA = -999;
      if (std::isinf(zdcPerCand[candId].timeC))
        zdcPerCand[candId].timeC = -999;
      if (std::isinf(zdcPerCand[candId].enA))
        zdcPerCand[candId].enA = -999;
      if (std::isinf(zdcPerCand[candId].enC))
        zdcPerCand[candId].enC = -999;
    }
  }

  // function to select muon tracks
  template <typename TTracks>
  bool isMuonSelected(const TTracks& fwdTrack)
  {
    float rAbs = fwdTrack.rAtAbsorberEnd();
    float pDca = fwdTrack.pDca();
    float pt = RecoDecay::pt(fwdTrack.px(), fwdTrack.py());
    float eta = RecoDecay::eta(std::array{fwdTrack.px(), fwdTrack.py(), fwdTrack.pz()});
    if (eta < kEtaMin || eta > kEtaMax)
      return false;
    if (pt < kPtMin)
      return false;
    if (rAbs < kRAbsMin || rAbs > kRAbsMax)
      return false;
    if (pDca > kPDca)
      return false;
    return true;
  }

  // function that processes the candidates:
  // it applies V0 selection, trk selection, kine selection, and fills the histograms
  // it also divides the data in neutron classes
  // used for real data
  void processCand(CandidatesFwd::iterator const& cand,
                   ForwardTracks::iterator const& tr1, ForwardTracks::iterator const& tr2)
  {
    // V0 selection
    const auto& ampsV0A = cand.amplitudesV0A();
    const auto& ampsRelBCsV0A = cand.ampRelBCsV0A();
    for (unsigned int i = 0; i < ampsV0A.size(); ++i) {
      if (std::abs(ampsRelBCsV0A[i]) <= 1) {
        if (ampsV0A[i] > kMaxAmpV0A)
          return;
      }
    }

    // select opposite charge events only                                                                                  
    if (cand.netCharge() != 0) {
      return;
    }
    
    // MCH-MID match selection
    int nMIDs = 0;
    if (tr1.chi2MatchMCHMID() > 0)
      nMIDs++;
    if (tr2.chi2MatchMCHMID() > 0)
      nMIDs++;
    if (nMIDs != kReqMatchMIDTracks)
      return;
    // MFT-MID match selection (if MFT is requested by the trackType)
    if (myTrackType == 0) {
      // if MFT is requested check that the tracks is inside the MFT acceptance
      kEtaMin = -3.6;
      kEtaMax = -2.5;

      int nMFT = 0;
      if (tr1.chi2MatchMCHMFT() > 0 && tr1.chi2MatchMCHMFT() < kMaxChi2MFTMatch)
        nMFT++;
      if (tr2.chi2MatchMCHMFT() > 0 && tr2.chi2MatchMCHMFT() < kMaxChi2MFTMatch)
        nMFT++;
      if (nMFT != kReqMatchMFTTracks)
        return;
    }
    // track selection
    if (!isMuonSelected(*tr1))
      return;
    if (!isMuonSelected(*tr2))
      return;

    // form Lorentz vectors
    auto mMu = o2::constants::physics::MassMuonMinus;
    LorentzVector<PxPyPzM4D<float>> p1(tr1.px(), tr1.py(), tr1.pz(), mMu);
    LorentzVector<PxPyPzM4D<float>> p2(tr2.px(), tr2.py(), tr2.pz(), mMu);
    LorentzVector p = p1 + p2;

    // cut on pair kinematics
    // select mass
    if (p.M() < lowMass)
      return;
    if (p.M() > highMass)
      return;
    // select pt
    if (p.Pt() < lowPt)
      return;
    if (p.Pt() > highPt)
      return;
    // select rapidity
    if (p.Rapidity() < lowRapidity)
      return;
    if (p.Rapidity() > highRapidity)
      return;
    // fill the histos without looking at neutron emission
    registry.fill(HIST("hMass"), p.M());
    registry.fill(HIST("hPt"), p.Pt());
    registry.fill(HIST("hEta"), p.Eta());
    registry.fill(HIST("hRapidity"), p.Rapidity());
    registry.fill(HIST("hPhi"), p.Phi());
    registry.fill(HIST("hPtTrkPos"), p1.Pt());
    registry.fill(HIST("hPtTrkNeg"), p2.Pt());
    registry.fill(HIST("hEtaTrkPos"), p1.Eta());
    registry.fill(HIST("hEtaTrkNeg"), p2.Eta());
    registry.fill(HIST("hPhiTrkPos"), p1.Phi());
    registry.fill(HIST("hPhiTrkNeg"), p2.Phi());

    
 // store the event to save it into a tree                                                                                 
    if (tr1.sign() > 0) {
      dimuSel(cand.runNumber(),
              p.M(), p.E(), p.Px(), p.Py(), p.Pz(), p.Pt(), p.Rapidity(), p.Phi(),
	      p1.E(), p1.Px(), p1.Py(), p1.Pz(), p1.Pt(), p1.Eta(), p1.Phi(),
              p2.E(), p2.Px(), p2.Py(), p2.Pz(), p2.Pt(), p2.Eta(), p2.Phi());
    } else {
      dimuSel(cand.runNumber(),
              p.M(), p.E(), p.Px(), p.Py(), p.Pz(), p.Pt(), p.Rapidity(), p.Phi(),
              p2.E(), p2.Px(), p2.Py(), p2.Pz(), p2.Pt(), p2.Eta(), p2.Phi(),
              p1.E(), p1.Px(), p1.Py(), p1.Pz(), p1.Pt(), p1.Eta(), p1.Phi());
    }
    
  }
  // PROCESS FUNCTION
  void processData(CandidatesFwd const& eventCandidates,
                   o2::aod::UDZdcsReduced const& ZDCs,
                   ForwardTracks const& fwdTracks)
  {

    // map with the tracks
    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    // takes a tracks table with a coloumn of collision ID and makes it into a map of collision ID to each track.
    collectCandIDs(tracksPerCand, fwdTracks);

    // map with the ZDC info
    std::unordered_map<int32_t, ZDCinfo> zdcPerCand;
    collectCandZDCInfo(zdcPerCand, ZDCs);

    // loop over the candidates
    for (const auto& item : tracksPerCand) {
      int32_t trId1 = item.second[0];
      int32_t trId2 = item.second[1];
      int32_t candID = item.first;
      auto cand = eventCandidates.iteratorAt(candID);
      auto tr1 = fwdTracks.iteratorAt(trId1);
      auto tr2 = fwdTracks.iteratorAt(trId2);
      processCand(cand, tr1, tr2);
    }
  }

  PROCESS_SWITCH(UpcPolarisationJpsiIncoh, processData, "", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcPolarisationJpsiIncoh>(cfgc),
  };
}
