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

/// \file doublephitable.cxx
/// \brief Selection of events with triplets and pairs for femtoscopic studies
///
/// \author Sourav Kundu, sourav.kundu@cern.ch

#include "PWGLF/DataModel/ReducedDoublePhiTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/PID.h>
#include <ReconstructionDataFormats/Track.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TMath.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct doublephitable {
  // Produce derived tables
  Produces<aod::RedPhiEvents> redPhiEvents;
  Produces<aod::PhiTracks> phiTrack;
  Produces<aod::PhiPhiPairs> phiPhiPair;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  /// Event selection
  struct : ConfigurableGroup {
    std::string prefix = std::string("EventSel");
    Configurable<bool> ConfEvtSelectZvtx{"ConfEvtSelectZvtx", true, "Event selection includes max. z-Vertex"};
    Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
    Configurable<bool> ConfEvtSel8{"ConfEvtSel8", true, "Event selection sel8"};
  } EventSel;

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  // Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 0.0f, "Accepted maximum Centrality"};
  // Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 100.0f, "Accepted minimum Centrality"};
  // track
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "ccdb-url"};
  Configurable<bool> useTrigger{"useTrigger", true, "use Trigger"};
  Configurable<bool> useGlobalTrack{"useGlobalTrack", true, "use Global track"};
  Configurable<float> cfgCutTOFBeta{"cfgCutTOFBeta", 0.0, "cut TOF beta"};
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "cut on Charge"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", -2.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutTOF{"nsigmaCutTOF", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<bool> isDeepAngle{"isDeepAngle", true, "Deep Angle cut"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {120, 0.98, 1.1}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.0, 10.}, "#it{p}_{T} (GeV/#it{c})"};
  Configurable<float> minPhiMass{"minPhiMass", 1.01, "Minimum phi mass"};
  Configurable<float> maxPhiMass{"maxPhiMass", 1.03, "Maximum phi mass"};
  Configurable<float> nsigmaCutTPCPreSel{"nsigmacutTPCPreSel", 3.0, "Value of the TPC Nsigma cut Pre selection"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path for Run 3 magnetic field"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "Propagate tracks to the fitted common vertex"};
  Configurable<bool> useAbsDCAFit{"useAbsDCAFit", false, "Use absolute-distance minimisation"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Use covariance-weighted final PCA"};
  Configurable<float> maxRFit{"maxRFit", 200.f, "Maximum fitted-vertex radius"};
  Configurable<float> maxDZIniFit{"maxDZIniFit", 4.f, "Maximum initial longitudinal separation"};
  Configurable<float> maxDXYIniFit{"maxDXYIniFit", 4.f, "Maximum initial transverse separation"};
  Configurable<float> maxChi2Fit{"maxChi2Fit", 1.e9f, "Maximum four-track fit chi2"};
  Configurable<float> minParamChangeFit{"minParamChangeFit", 1.e-3f, "Fitter convergence threshold"};
  Configurable<float> minRelChi2ChangeFit{"minRelChi2ChangeFit", 0.9f, "Relative chi2 convergence threshold"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  // Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  Filter PIDcutFilter = nabs(aod::pidtpc::tpcNSigmaKa) < nsigmaCutTPCPreSel;

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>>;
  using TrackCandidates = soa::Filtered<
    soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTOFbeta, aod::pidTPCFullKa, aod::pidTOFFullKa>>;

  SliceCache cache;
  Partition<TrackCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  Partition<TrackCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  o2::vertexing::DCAFitterN<4> df4;
  int fieldRunNumber = -999;
  float bz = 0.f;

  OutputObj<TH1D> hProcessedEvents{TH1D("hProcessedEvents", ";; Number of events", 4, 0.0f, 4.0f)};
  // Histogram
  HistogramRegistry qaRegistry{"QAHistos",
                               {
                                 {"hEventstat", "hEventstat", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
                                 {"hInvMassPhi", "hInvMassPhi", {HistType::kTH2F, {{40, 1.0f, 1.04f}, {100, 0.0f, 10.0f}}}},
                                 {"hNsigmaPtkaonTPCITS", "hNsigmaPtkaonTPCITS", {HistType::kTH3F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                                 {"hNsigmaPtkaonTPC", "hNsigmaPtkaonTPC", {HistType::kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                                 {"hNsigmaPtkaonTOF", "hNsigmaPtkaonTOF", {HistType::kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                               },
                               OutputObjHandlingPolicy::AnalysisObject};

  HistogramRegistry pairQaRegistry{"PairVertexQAHistos",
                                   {
                                     {"hPairMassPtFit",
                                      "Refitted pair;M_{#phi#phi} (GeV/#it{c}^{2});"
                                      "#it{p}_{T}^{#phi#phi} (GeV/#it{c})",
                                      {HistType::kTH2F, {{400, 2.2f, 3.4f}, {200, 0.0f, 20.0f}}}},
                                     {"hFitChi2Ndf", "Four-kaon common-vertex fit;#chi^{2}/NDF;counts", {HistType::kTH1F, {{500, 0.f, 100.f}}}},
                                     {"hSumDcaChi2", "Four-kaon DCA compatibility;#Sigma_{i=1}^{4} #chi^{2}_{DCA,i};counts", {HistType::kTH1F, {{500, 0.f, 200.f}}}},
                                     {"hVertexL3DSig", "Fitted vertex relative to PV;L_{3D}/#sigma_{L_{3D}};counts", {HistType::kTH1F, {{500, 0.f, 100.f}}}},
                                     {"hRmsDcaSig", "Four-kaon DCA compatibility;#sqrt{#Sigma#chi^{2}_{DCA}/8};counts", {HistType::kTH1F, {{300, 0.0f, 15.0f}}}},
                                   },
                                   OutputObjHandlingPolicy::AnalysisObject};

  double massKa = o2::constants::physics::MassKPlus;

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(url.value);
    ccdbApi.init(url);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    zorroSummary.setObject(zorro.getZorroSummary());
    hProcessedEvents->GetXaxis()->SetBinLabel(1, "All Trigger events");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "Events with Double Phi Trigger");
    hProcessedEvents->GetXaxis()->SetBinLabel(3, "Events eith trigger and Evsel");
    hProcessedEvents->GetXaxis()->SetBinLabel(4, "Final Event");
    df4.setPropagateToPCA(propagateToPCA);
    df4.setMaxR(maxRFit);
    df4.setMaxDZIni(maxDZIniFit);
    df4.setMaxDXYIni(maxDXYIniFit);
    df4.setMaxChi2(maxChi2Fit);
    df4.setMinParamChange(minParamChangeFit);
    df4.setMinRelChi2Change(minRelChi2ChangeFit);
    df4.setUseAbsDCA(useAbsDCAFit);
    df4.setWeightedFinalPCA(useWeightedFinalPCA);
    df4.setMatCorrType(o2::base::Propagator::MatCorrType::USEMatCorrNONE);
  }

  template <typename T>
  bool isSelectedEvent(T const& col)
  {
    if (EventSel.ConfEvtSelectZvtx && std::abs(col.posZ()) > EventSel.ConfEvtZvtx) {
      return false;
    }
    if (EventSel.ConfEvtSel8 && !col.sel8()) {
      return false;
    }
    return true;
  }
  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (useGlobalTrack && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsCrossedRows() > cfgTPCcluster)) {
      return false;
    }
    return true;
  }
  template <typename T>
  bool selectionPID(const T& candidate)
  {
    if (!candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
      return true;
    }
    if (candidate.hasTOF() && candidate.beta() > cfgCutTOFBeta && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && TMath::Abs(candidate.tofNSigmaKa()) < nsigmaCutTOF) {
      return true;
    }
    return false;
  }

  int currentRunNumber = -999;
  int lastRunNumber = -999;
  ROOT::Math::PxPyPzMVector KaonPlus, KaonMinus, PhiMesonMother, PhiVectorDummy, Phid1dummy, Phid2dummy;
  void processPhiReducedTable(EventCandidates::iterator const& collision, TrackCandidates const&, aod::BCsWithTimestamps const&)
  {
    o2::aod::ITSResponse itsResponse;
    bool keepEventDoublePhi = false;
    int numberPhi = 0;
    std::vector<int64_t> Phid1Index = {};
    std::vector<int64_t> Phid2Index = {};
    std::vector<float> Phid1Charge = {};
    std::vector<float> Phid2Charge = {};
    std::vector<float> Phid1TPC = {};
    std::vector<float> Phid2TPC = {};
    std::vector<float> Phid1TOF = {};
    std::vector<float> Phid2TOF = {};
    std::vector<int> Phid1TOFHit = {};
    std::vector<int> Phid2TOFHit = {};
    std::vector<ROOT::Math::PtEtaPhiMVector> phiresonance, phiresonanced1, phiresonanced2;
    int Npostrack = 0;
    int Nnegtrack = 0;
    float centrality = collision.centFT0M();
    currentRunNumber = collision.bc_as<aod::BCsWithTimestamps>().runNumber();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    hProcessedEvents->Fill(0.5);
    bool zorroSelected = false;
    if (currentRunNumber != lastRunNumber) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), "fTriggerEventDoublePhi");
      zorro.populateHistRegistry(qaRegistry, bc.runNumber());
      lastRunNumber = currentRunNumber;
    }
    if (useTrigger) {
      zorroSelected = zorro.isSelected(collision.template bc_as<aod::BCsWithTimestamps>().globalBC());
    } else {
      zorroSelected = true;
    }
    if (zorroSelected) {
      hProcessedEvents->Fill(1.5);
    }
    if (zorroSelected && isSelectedEvent(collision)) {
      hProcessedEvents->Fill(2.5);
      auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      for (auto track1 : posThisColl) {
        // track selection
        if (!selectionTrack(track1)) {
          continue;
        }
        qaRegistry.fill(HIST("hNsigmaPtkaonTPCITS"), itsResponse.nSigmaITS<o2::track::PID::Kaon>(track1), track1.tpcNSigmaKa(), track1.pt());
        // PID check
        if (!selectionPID(track1)) {
          continue;
        }
        if (!(itsResponse.nSigmaITS<o2::track::PID::Kaon>(track1) > -2.0 && itsResponse.nSigmaITS<o2::track::PID::Kaon>(track1) < 3.0)) {
          continue;
        }
        Npostrack = Npostrack + 1;
        qaRegistry.fill(HIST("hNsigmaPtkaonTPC"), track1.tpcNSigmaKa(), track1.pt());
        if (track1.hasTOF()) {
          qaRegistry.fill(HIST("hNsigmaPtkaonTOF"), track1.tofNSigmaKa(), track1.pt());
        }
        auto track1ID = track1.globalIndex();
        for (auto track2 : negThisColl) {
          auto track2ID = track2.globalIndex();
          if (track2ID == track1ID) {
            continue;
          }
          // track selection
          if (!selectionTrack(track2)) {
            continue;
          }
          // PID check
          if (!selectionPID(track2)) {
            continue;
          }
          if (!(itsResponse.nSigmaITS<o2::track::PID::Kaon>(track2) > -2.0 && itsResponse.nSigmaITS<o2::track::PID::Kaon>(track2) < 3.0)) {
            continue;
          }
          if (Npostrack == 1) {
            Nnegtrack = Nnegtrack + 1;
          }
          KaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
          KaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
          PhiMesonMother = KaonPlus + KaonMinus;
          if (PhiMesonMother.M() > minPhiMass && PhiMesonMother.M() < maxPhiMass) {
            numberPhi = numberPhi + 1;
            ROOT::Math::PtEtaPhiMVector temp1(track1.pt(), track1.eta(), track1.phi(), massKa);
            ROOT::Math::PtEtaPhiMVector temp2(track2.pt(), track2.eta(), track2.phi(), massKa);
            ROOT::Math::PtEtaPhiMVector temp3(PhiMesonMother.pt(), PhiMesonMother.eta(), PhiMesonMother.phi(), PhiMesonMother.M());
            phiresonanced1.push_back(temp1);
            phiresonanced2.push_back(temp2);
            phiresonance.push_back(temp3);
            Phid1Index.push_back(track1.globalIndex());
            Phid2Index.push_back(track2.globalIndex());
            Phid1Charge.push_back(track1.sign());
            Phid2Charge.push_back(track2.sign());
            Phid1TPC.push_back(track1.tpcNSigmaKa());
            Phid2TPC.push_back(track2.tpcNSigmaKa());
            auto d1TOFHit = -1;
            auto d2TOFHit = -1;
            auto d1TOF = -999.0;
            auto d2TOF = -999.0;
            if (track1.hasTOF()) {
              d1TOFHit = 1;
              d1TOF = track1.tofNSigmaKa();
            }
            if (track2.hasTOF()) {
              d2TOFHit = 1;
              d2TOF = track2.tofNSigmaKa();
            }
            Phid1TOF.push_back(d1TOF);
            Phid2TOF.push_back(d2TOF);
            Phid1TOFHit.push_back(d1TOFHit);
            Phid2TOFHit.push_back(d2TOFHit);
            qaRegistry.fill(HIST("hInvMassPhi"), PhiMesonMother.M(), PhiMesonMother.Pt());
          }
        }
      }
    } // select collision
    if (numberPhi > 1 && zorroSelected && Nnegtrack > 1 && Npostrack > 1) {
      keepEventDoublePhi = true;
      hProcessedEvents->Fill(3.5);
    }
    qaRegistry.fill(HIST("hEventstat"), 0.5);
    if (keepEventDoublePhi && numberPhi > 1 && (phiresonance.size() == phiresonanced1.size()) && (phiresonance.size() == phiresonanced2.size())) {
      qaRegistry.fill(HIST("hEventstat"), 1.5);
      /////////// Fill collision table///////////////
      redPhiEvents(bc.globalBC(), currentRunNumber, bc.timestamp(), collision.posZ(), collision.numContrib(), Npostrack, Nnegtrack, centrality);
      auto indexEvent = redPhiEvents.lastIndex();
      //// Fill track table for Phi//////////////////
      for (auto if1 = phiresonance.begin(); if1 != phiresonance.end(); ++if1) {
        auto i5 = std::distance(phiresonance.begin(), if1);
        PhiVectorDummy = phiresonance.at(i5);
        Phid1dummy = phiresonanced1.at(i5);
        Phid2dummy = phiresonanced2.at(i5);
        phiTrack(indexEvent, PhiVectorDummy.Px(), PhiVectorDummy.Py(), PhiVectorDummy.Pz(), Phid1dummy.Px(), Phid1dummy.Py(), Phid1dummy.Pz(), Phid2dummy.Px(), Phid2dummy.Py(),
                 Phid2dummy.Pz(), PhiVectorDummy.M(), Phid1Index.at(i5), Phid2Index.at(i5), Phid1Charge.at(i5), Phid2Charge.at(i5), Phid1TPC.at(i5), Phid2TPC.at(i5), Phid1TOFHit.at(i5),
                 Phid2TOFHit.at(i5), Phid1TOF.at(i5), Phid2TOF.at(i5));
      }
    }
  } // process
  PROCESS_SWITCH(doublephitable, processPhiReducedTable, "Process table creation for double phi", false);

  struct PhiCandidateVtx {
    TrackCandidates::iterator kPlus;
    TrackCandidates::iterator kMinus;
    ROOT::Math::PxPyPzMVector phiPreFit;
    float itsPlus = -999.f;
    float itsMinus = -999.f;
  };

  struct KaonPairPayload {
    int64_t index = -1;
    int8_t charge = 0;
    float dcaXY = -999.f;
    float dcaZ = -999.f;
    float dcaXYSig = -999.f;
    float dcaZSig = -999.f;
    float dcaChi2 = -1.f;
    float px = 0.f;
    float py = 0.f;
    float pz = 0.f;
    int8_t tofHit = -1;
    float tpc = -999.f;
    float tof = -999.f;
    float its = -999.f;
  };

  struct PhiPhiPairPayload {
    size_t phi1Index = 0;
    size_t phi2Index = 0;
    float pairMass = 0.f;
    float pairPx = 0.f;
    float pairPy = 0.f;
    float pairPz = 0.f;
    float phi1Mass = 0.f;
    float phi1Px = 0.f;
    float phi1Py = 0.f;
    float phi1Pz = 0.f;
    float phi2Mass = 0.f;
    float phi2Px = 0.f;
    float phi2Py = 0.f;
    float phi2Pz = 0.f;
    KaonPairPayload k1;
    KaonPairPayload k2;
    KaonPairPayload k3;
    KaonPairPayload k4;
    int8_t fitStatus = -1;
    float fitChi2 = -1.f;

    // Nominal geometric NDF for four tracks fitted to one 3D point:
    // 2 * Ntracks - 3 = 5.
    int8_t fitNdf = 5;
    float fitChi2Ndf = -1.f;
    float pvX = 0.f;
    float pvY = 0.f;
    float pvZ = 0.f;
    float vtxX = 0.f;
    float vtxY = 0.f;
    float vtxZ = 0.f;
    float deltaVtxX = 0.f;
    float deltaVtxY = 0.f;
    float deltaVtxZ = 0.f;
    float vertexLxy = 0.f;
    float vertexL3D = 0.f;
    float vertexLxyErr = 0.f;
    float vertexL3DErr = 0.f;
    float vertexLxySig = 0.f;
    float vertexL3DSig = 0.f;
    int8_t nValidDca = 0;
    float sumDcaXYSig2 = 0.f;
    float sumDcaZSig2 = 0.f;
    float sumDcaChi2 = 0.f;
    float rmsDcaSig = -1.f;
    float maxDcaChi2 = -1.f;
    float maxAbsDcaXYSig = -1.f;
    float maxAbsDcaZSig = -1.f;
  };

  template <typename T>
  bool selectionITSKaon(const T& candidate, float& nSigmaITS)
  {
    nSigmaITS = o2::aod::ITSResponse::nSigmaITS<o2::track::PID::Kaon>(candidate);
    if (!(nSigmaITS > -2.0f && nSigmaITS < 3.0f)) {
      return false;
    }

    return true;
  }

  void updateMagneticField(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (fieldRunNumber == bc.runNumber()) {
      return;
    }

    auto* grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrpMag.value, bc.timestamp());
    if (grpmag == nullptr) {
      LOGF(fatal, "GRPMagField unavailable for run %d, timestamp %llu", bc.runNumber(), static_cast<unsigned long long>(bc.timestamp()));
    }

    o2::base::Propagator::initFieldFromGRP(grpmag);
    bz = o2::base::Propagator::Instance()->getNominalBz();
    df4.setBz(bz);
    fieldRunNumber = bc.runNumber();
  }

  template <typename T, typename V>
  KaonPairPayload makeKaonPayload(const T& track, const std::array<float, 3>& momentumAtVertex, float itsNSigma, const V& primaryVertex)
  {
    KaonPairPayload out;
    out.index = track.globalIndex();
    out.charge = static_cast<int8_t>(track.sign());
    out.dcaXY = track.dcaXY();
    out.dcaZ = track.dcaZ();
    out.px = momentumAtVertex[0];
    out.py = momentumAtVertex[1];
    out.pz = momentumAtVertex[2];
    out.tofHit = track.hasTOF() ? 1 : -1;
    out.tpc = track.tpcNSigmaKa();
    out.tof = track.hasTOF() ? track.tofNSigmaKa() : -999.f;
    out.its = itsNSigma;
    auto trackParCov = getTrackParCov(track);
    o2::dataformats::DCA dca;
    if (trackParCov.propagateToDCA(primaryVertex, bz, &dca)) {
      out.dcaXY = dca.getY();
      out.dcaZ = dca.getZ();
      if (dca.getSigmaY2() > 0.f) {
        out.dcaXYSig = out.dcaXY / std::sqrt(dca.getSigmaY2());
      }

      if (dca.getSigmaZ2() > 0.f) {
        out.dcaZSig = out.dcaZ / std::sqrt(dca.getSigmaZ2());
      }

      const float chi2 = dca.calcChi2();
      if (std::isfinite(chi2) && chi2 >= 0.f) {
        out.dcaChi2 = chi2;
      }
    }

    return out;
  }

  template <typename CovPV, typename CovVtx>
  static void calculateVertexQuantities(PhiPhiPairPayload& pair, const CovPV& covPV, const CovVtx& covVtx)
  {
    const float cxx = covPV[0] + covVtx[0];
    const float cxy = covPV[1] + covVtx[1];
    const float cyy = covPV[2] + covVtx[2];
    const float cxz = covPV[3] + covVtx[3];
    const float cyz = covPV[4] + covVtx[4];
    const float czz = covPV[5] + covVtx[5];
    pair.vertexLxy = std::hypot(pair.deltaVtxX, pair.deltaVtxY);
    pair.vertexL3D = std::sqrt(pair.deltaVtxX * pair.deltaVtxX + pair.deltaVtxY * pair.deltaVtxY + pair.deltaVtxZ * pair.deltaVtxZ);
    float varLxy = 0.5f * (cxx + cyy);
    if (pair.vertexLxy > 1.e-12f) {
      varLxy = (pair.deltaVtxX * pair.deltaVtxX * cxx + 2.f * pair.deltaVtxX * pair.deltaVtxY * cxy + pair.deltaVtxY * pair.deltaVtxY * cyy) / (pair.vertexLxy * pair.vertexLxy);
    }

    float varL3D = (cxx + cyy + czz) / 3.f;
    if (pair.vertexL3D > 1.e-12f) {
      varL3D = (pair.deltaVtxX * pair.deltaVtxX * cxx + pair.deltaVtxY * pair.deltaVtxY * cyy + pair.deltaVtxZ * pair.deltaVtxZ * czz +
                2.f * pair.deltaVtxX * pair.deltaVtxY * cxy + 2.f * pair.deltaVtxX * pair.deltaVtxZ * cxz + 2.f * pair.deltaVtxY * pair.deltaVtxZ * cyz) /
               (pair.vertexL3D * pair.vertexL3D);
    }

    pair.vertexLxyErr = varLxy > 0.f ? std::sqrt(varLxy) : 0.f;
    pair.vertexL3DErr = varL3D > 0.f ? std::sqrt(varL3D) : 0.f;
    pair.vertexLxySig = pair.vertexLxyErr > 0.f ? pair.vertexLxy / pair.vertexLxyErr : 0.f;
    pair.vertexL3DSig = pair.vertexL3DErr > 0.f ? pair.vertexL3D / pair.vertexL3DErr : 0.f;
  }

  static void calculateDcaSummary(PhiPhiPairPayload& pair)
  {
    const std::array<const KaonPairPayload*, 4> kaons{&pair.k1, &pair.k2, &pair.k3, &pair.k4};
    pair.nValidDca = 0;
    pair.sumDcaXYSig2 = 0.f;
    pair.sumDcaZSig2 = 0.f;
    pair.sumDcaChi2 = 0.f;
    pair.maxDcaChi2 = -1.f;
    pair.maxAbsDcaXYSig = -1.f;
    pair.maxAbsDcaZSig = -1.f;
    for (const auto* kaon : kaons) {
      if (std::isfinite(kaon->dcaXYSig) && kaon->dcaXYSig > -900.f) {
        pair.sumDcaXYSig2 += kaon->dcaXYSig * kaon->dcaXYSig;
        pair.maxAbsDcaXYSig = std::max(pair.maxAbsDcaXYSig, std::abs(kaon->dcaXYSig));
      }

      if (std::isfinite(kaon->dcaZSig) && kaon->dcaZSig > -900.f) {
        pair.sumDcaZSig2 += kaon->dcaZSig * kaon->dcaZSig;
        pair.maxAbsDcaZSig = std::max(pair.maxAbsDcaZSig, std::abs(kaon->dcaZSig));
      }

      if (std::isfinite(kaon->dcaChi2) && kaon->dcaChi2 >= 0.f) {
        ++pair.nValidDca;
        pair.sumDcaChi2 += kaon->dcaChi2;
        pair.maxDcaChi2 = std::max(pair.maxDcaChi2, kaon->dcaChi2);
      }
    }

    if (pair.nValidDca == 4) {
      pair.rmsDcaSig = std::sqrt(pair.sumDcaChi2 / 8.f);
    }
  }
  void processPhiPairVertexTable(EventCandidates::iterator const& collision, TrackCandidates const&, aod::BCsWithTimestamps const&)
  {
    qaRegistry.fill(HIST("hEventstat"), 0.5);
    hProcessedEvents->Fill(0.5);
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    currentRunNumber = bc.runNumber();
    if (currentRunNumber != lastRunNumber) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), "fTriggerEventDoublePhi");
      zorro.populateHistRegistry(qaRegistry, bc.runNumber());
      lastRunNumber = currentRunNumber;
    }

    const bool zorroSelected = useTrigger ? zorro.isSelected(bc.globalBC()) : true;
    if (!zorroSelected) {
      return;
    }

    hProcessedEvents->Fill(1.5);
    if (!isSelectedEvent(collision)) {
      return;
    }

    hProcessedEvents->Fill(2.5);
    updateMagneticField(bc);
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    std::vector<TrackCandidates::iterator> selectedPos;
    std::vector<TrackCandidates::iterator> selectedNeg;
    std::vector<float> selectedPosITS;
    std::vector<float> selectedNegITS;
    selectedPos.reserve(posThisColl.size());
    selectedNeg.reserve(negThisColl.size());
    selectedPosITS.reserve(posThisColl.size());
    selectedNegITS.reserve(negThisColl.size());
    for (auto track : posThisColl) {
      if (!selectionTrack(track)) {
        continue;
      }

      if (!selectionPID(track)) {
        continue;
      }

      float nSigmaITS = -999.f;
      if (!selectionITSKaon(track, nSigmaITS)) {
        continue;
      }

      selectedPos.push_back(track);
      selectedPosITS.push_back(nSigmaITS);
      qaRegistry.fill(HIST("hNsigmaPtkaonTPC"), track.tpcNSigmaKa(), track.pt());
      if (track.hasTOF()) {
        qaRegistry.fill(HIST("hNsigmaPtkaonTOF"), track.tofNSigmaKa(), track.pt());
      }
    }

    for (auto track : negThisColl) {
      if (!selectionTrack(track)) {
        continue;
      }

      if (!selectionPID(track)) {
        continue;
      }

      float nSigmaITS = -999.f;
      if (!selectionITSKaon(track, nSigmaITS)) {
        continue;
      }

      selectedNeg.push_back(track);
      selectedNegITS.push_back(nSigmaITS);
      qaRegistry.fill(HIST("hNsigmaPtkaonTPC"), track.tpcNSigmaKa(), track.pt());
      if (track.hasTOF()) {
        qaRegistry.fill(HIST("hNsigmaPtkaonTOF"), track.tofNSigmaKa(), track.pt());
      }
    }

    const int nPosTrack = static_cast<int>(selectedPos.size());
    const int nNegTrack = static_cast<int>(selectedNeg.size());
    if (nPosTrack < 2 || nNegTrack < 2) {
      return;
    }

    std::vector<PhiCandidateVtx> phiCandidates;
    phiCandidates.reserve(selectedPos.size() * selectedNeg.size());
    for (size_t iPos = 0; iPos < selectedPos.size(); ++iPos) {
      for (size_t iNeg = 0; iNeg < selectedNeg.size(); ++iNeg) {
        const auto& trackPlus = selectedPos[iPos];
        const auto& trackMinus = selectedNeg[iNeg];
        ROOT::Math::PxPyPzMVector kaonPlus(trackPlus.px(), trackPlus.py(), trackPlus.pz(), massKa);
        ROOT::Math::PxPyPzMVector kaonMinus(trackMinus.px(), trackMinus.py(), trackMinus.pz(), massKa);
        const auto phi = kaonPlus + kaonMinus;
        if (phi.M() <= minPhiMass || phi.M() >= maxPhiMass) {
          continue;
        }

        phiCandidates.push_back({trackPlus, trackMinus, phi, selectedPosITS[iPos], selectedNegITS[iNeg]});
        qaRegistry.fill(HIST("hInvMassPhi"), phi.M(), phi.Pt());
      }
    }

    if (phiCandidates.size() < 2) {
      return;
    }

    const auto primaryVertex = getPrimaryVertex(collision);
    const auto covPV = primaryVertex.getCov();
    std::vector<PhiPhiPairPayload> acceptedPairs;
    for (size_t iPhi = 0; iPhi < phiCandidates.size(); ++iPhi) {
      for (size_t jPhi = iPhi + 1; jPhi < phiCandidates.size(); ++jPhi) {
        const auto& phi1 = phiCandidates[iPhi];
        const auto& phi2 = phiCandidates[jPhi];
        const int64_t id1 = phi1.kPlus.globalIndex();
        const int64_t id2 = phi1.kMinus.globalIndex();
        const int64_t id3 = phi2.kPlus.globalIndex();
        const int64_t id4 = phi2.kMinus.globalIndex();
        if (id1 == id2 || id1 == id3 || id1 == id4 || id2 == id3 || id2 == id4 || id3 == id4) {
          continue;
        }
        auto trackPar1 = getTrackParCov(phi1.kPlus);
        auto trackPar2 = getTrackParCov(phi1.kMinus);
        auto trackPar3 = getTrackParCov(phi2.kPlus);
        auto trackPar4 = getTrackParCov(phi2.kMinus);
        int nFitCandidates = 0;
        try {
          nFitCandidates = df4.process(trackPar1, trackPar2, trackPar3, trackPar4);
        } catch (const std::runtime_error& error) {
          LOGF(debug, "Four-kaon DCAFitterN failed: %s", error.what());
          continue;
        } catch (...) {
          LOG(debug) << "Four-kaon DCAFitterN failed with an unknown exception";
          continue;
        }

        if (nFitCandidates == 0) {
          continue;
        }

        std::array<float, 3> pK1{0.f, 0.f, 0.f};
        std::array<float, 3> pK2{0.f, 0.f, 0.f};
        std::array<float, 3> pK3{0.f, 0.f, 0.f};
        std::array<float, 3> pK4{0.f, 0.f, 0.f};
        df4.getTrack(0).getPxPyPzGlo(pK1);
        df4.getTrack(1).getPxPyPzGlo(pK2);
        df4.getTrack(2).getPxPyPzGlo(pK3);
        df4.getTrack(3).getPxPyPzGlo(pK4);
        const ROOT::Math::PxPyPzMVector k1Fit(pK1[0], pK1[1], pK1[2], massKa);
        const ROOT::Math::PxPyPzMVector k2Fit(pK2[0], pK2[1], pK2[2], massKa);
        const ROOT::Math::PxPyPzMVector k3Fit(pK3[0], pK3[1], pK3[2], massKa);
        const ROOT::Math::PxPyPzMVector k4Fit(pK4[0], pK4[1], pK4[2], massKa);
        const auto phi1Fit = k1Fit + k2Fit;
        const auto phi2Fit = k3Fit + k4Fit;
        const auto pairFit = phi1Fit + phi2Fit;

        // No pair-pT or pair-mass cut here.

        PhiPhiPairPayload pair;
        pair.phi1Index = iPhi;
        pair.phi2Index = jPhi;
        pair.pairMass = pairFit.M();
        pair.pairPx = pairFit.Px();
        pair.pairPy = pairFit.Py();
        pair.pairPz = pairFit.Pz();
        pair.phi1Mass = phi1Fit.M();
        pair.phi1Px = phi1Fit.Px();
        pair.phi1Py = phi1Fit.Py();
        pair.phi1Pz = phi1Fit.Pz();
        pair.phi2Mass = phi2Fit.M();
        pair.phi2Px = phi2Fit.Px();
        pair.phi2Py = phi2Fit.Py();
        pair.phi2Pz = phi2Fit.Pz();
        pair.k1 = makeKaonPayload(phi1.kPlus, pK1, phi1.itsPlus, primaryVertex);
        pair.k2 = makeKaonPayload(phi1.kMinus, pK2, phi1.itsMinus, primaryVertex);
        pair.k3 = makeKaonPayload(phi2.kPlus, pK3, phi2.itsPlus, primaryVertex);
        pair.k4 = makeKaonPayload(phi2.kMinus, pK4, phi2.itsMinus, primaryVertex);
        pair.fitStatus = static_cast<int8_t>(df4.getFitStatus());
        pair.fitChi2 = df4.getChi2AtPCACandidate();
        pair.fitChi2Ndf = pair.fitChi2 / static_cast<float>(pair.fitNdf);
        const auto& fittedVertex = df4.getPCACandidate();
        const auto covVtx = df4.calcPCACovMatrixFlat();
        pair.pvX = collision.posX();
        pair.pvY = collision.posY();
        pair.pvZ = collision.posZ();
        pair.vtxX = fittedVertex[0];
        pair.vtxY = fittedVertex[1];
        pair.vtxZ = fittedVertex[2];
        pair.deltaVtxX = pair.vtxX - pair.pvX;
        pair.deltaVtxY = pair.vtxY - pair.pvY;
        pair.deltaVtxZ = pair.vtxZ - pair.pvZ;
        calculateVertexQuantities(pair, covPV, covVtx);
        calculateDcaSummary(pair);
        acceptedPairs.push_back(pair);
      }
    }

    if (acceptedPairs.empty()) {
      return;
    }

    redPhiEvents(bc.globalBC(), currentRunNumber, bc.timestamp(), collision.posZ(), collision.numContrib(), nPosTrack, nNegTrack, collision.centFT0M());
    const auto indexEvent = redPhiEvents.lastIndex();
    std::vector<int64_t> phiRows(phiCandidates.size(), -1);
    for (size_t iPhi = 0; iPhi < phiCandidates.size(); ++iPhi) {
      const auto& phi = phiCandidates[iPhi];
      const int d1TOFHit = phi.kPlus.hasTOF() ? 1 : -1;
      const int d2TOFHit = phi.kMinus.hasTOF() ? 1 : -1;
      const float d1TOF = phi.kPlus.hasTOF() ? phi.kPlus.tofNSigmaKa() : -999.f;
      const float d2TOF = phi.kMinus.hasTOF() ? phi.kMinus.tofNSigmaKa() : -999.f;
      phiTrack(indexEvent, phi.phiPreFit.Px(), phi.phiPreFit.Py(), phi.phiPreFit.Pz(), phi.kPlus.px(), phi.kPlus.py(), phi.kPlus.pz(), phi.kMinus.px(), phi.kMinus.py(),
               phi.kMinus.pz(), phi.phiPreFit.M(), phi.kPlus.globalIndex(), phi.kMinus.globalIndex(), phi.kPlus.sign(), phi.kMinus.sign(), phi.kPlus.tpcNSigmaKa(),
               phi.kMinus.tpcNSigmaKa(), d1TOFHit, d2TOFHit, d1TOF, d2TOF);

      phiRows[iPhi] = phiTrack.lastIndex();
    }

    for (const auto& pair : acceptedPairs) {
      const auto& k1 = pair.k1;
      const auto& k2 = pair.k2;
      const auto& k3 = pair.k3;
      const auto& k4 = pair.k4;
      pairQaRegistry.fill(HIST("hPairMassPtFit"), pair.pairMass, std::hypot(pair.pairPx, pair.pairPy));
      pairQaRegistry.fill(HIST("hFitChi2Ndf"), pair.fitChi2Ndf);
      pairQaRegistry.fill(HIST("hSumDcaChi2"), pair.sumDcaChi2);
      pairQaRegistry.fill(HIST("hRmsDcaSig"), pair.rmsDcaSig);
      pairQaRegistry.fill(HIST("hVertexL3DSig"), pair.vertexL3DSig);
      phiPhiPair(indexEvent, phiRows[pair.phi1Index], phiRows[pair.phi2Index], pair.pairMass, pair.pairPx, pair.pairPy, pair.pairPz, pair.phi1Mass, pair.phi1Px, pair.phi1Py,
                 pair.phi1Pz, pair.phi2Mass, pair.phi2Px, pair.phi2Py, pair.phi2Pz, k1.index, k1.charge, k1.dcaXY, k1.dcaZ, k1.dcaXYSig, k1.dcaZSig, k1.dcaChi2, k1.px, k1.py, k1.pz,
                 k1.tofHit, k1.tpc, k1.tof, k1.its, k2.index, k2.charge, k2.dcaXY, k2.dcaZ, k2.dcaXYSig, k2.dcaZSig, k2.dcaChi2, k2.px, k2.py, k2.pz, k2.tofHit, k2.tpc, k2.tof, k2.its,
                 k3.index, k3.charge, k3.dcaXY, k3.dcaZ, k3.dcaXYSig, k3.dcaZSig, k3.dcaChi2, k3.px, k3.py, k3.pz, k3.tofHit, k3.tpc, k3.tof, k3.its, k4.index, k4.charge, k4.dcaXY,
                 k4.dcaZ, k4.dcaXYSig, k4.dcaZSig, k4.dcaChi2, k4.px, k4.py, k4.pz, k4.tofHit, k4.tpc, k4.tof, k4.its, pair.fitStatus, pair.fitChi2, pair.fitNdf, pair.fitChi2Ndf,
                 pair.pvX, pair.pvY, pair.pvZ, pair.vtxX, pair.vtxY, pair.vtxZ, pair.deltaVtxX, pair.deltaVtxY, pair.deltaVtxZ, pair.vertexLxy, pair.vertexL3D, pair.vertexLxyErr,
                 pair.vertexL3DErr, pair.vertexLxySig, pair.vertexL3DSig, pair.nValidDca, pair.sumDcaXYSig2, pair.sumDcaZSig2, pair.sumDcaChi2, pair.rmsDcaSig, pair.maxDcaChi2,
                 pair.maxAbsDcaXYSig, pair.maxAbsDcaZSig);
    }

    hProcessedEvents->Fill(3.5);
    qaRegistry.fill(HIST("hEventstat"), 1.5);
  }
  PROCESS_SWITCH(doublephitable, processPhiPairVertexTable, "Process collision, phi and four-kaon pair tables", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<doublephitable>(cfg)};
}
