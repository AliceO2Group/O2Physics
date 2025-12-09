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

#include "Common/Core/TrackSelection.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/MathConstants.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include <Framework/Configurable.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h>
#include <TDatabasePDG.h> // FIXME
#include <TMath.h>
#include <TPDGCode.h> // FIXME

#include <fairlogger/Logger.h>

#include <iostream>
#include <iterator>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct doublephitable {

  // Produce derived tables
  Produces<aod::RedPhiEvents> redPhiEvents;
  Produces<aod::PhiTracks> phiTrack;

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  // Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 0.0f, "Accepted maximum Centrality"};
  // Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 100.0f, "Accepted minimum Centrality"};
  // track
  Configurable<bool> useGlobalTrack{"useGlobalTrack", true, "use Global track"};
  Configurable<float> cfgCutTOFBeta{"cfgCutTOFBeta", 0.0, "cut TOF beta"};
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "cut on Charge"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutTOF{"nsigmaCutTOF", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<bool> isDeepAngle{"isDeepAngle", true, "Deep Angle cut"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {120, 0.98, 1.1}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.0, 10.}, "#it{p}_{T} (GeV/#it{c})"};
  Configurable<float> minPhiMass{"minPhiMass", 1.01, "Minimum phi mass"};
  Configurable<float> maxPhiMass{"maxPhiMass", 1.03, "Maximum phi mass"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  // Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  Filter PIDcutFilter = nabs(aod::pidtpc::tpcNSigmaKa) < nsigmaCutTPC;

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTOFbeta, aod::pidTPCFullKa, aod::pidTOFFullKa>>;

  SliceCache cache;
  Partition<TrackCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  Partition<TrackCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  // Histogram
  HistogramRegistry qaRegistry{"QAHistos", {
                                             {"hEventstat", "hEventstat", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
                                             {"hInvMassPhi", "hInvMassPhi", {HistType::kTH2F, {{40, 1.0f, 1.04f}, {100, 0.0f, 10.0f}}}},
                                             {"hNsigmaPtkaonTPC", "hNsigmaPtkaonTPC", {HistType::kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                                             {"hNsigmaPtkaonTOF", "hNsigmaPtkaonTOF", {HistType::kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                                           },
                               OutputObjHandlingPolicy::AnalysisObject};

  double massKa = o2::constants::physics::MassKPlus;

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (useGlobalTrack && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster)) {
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
  // deep angle cut on pair to remove photon conversion
  template <typename T1, typename T2>
  bool selectionPair(const T1& candidate1, const T2& candidate2)
  {
    double pt1, pt2, pz1, pz2, p1, p2, angle;
    pt1 = candidate1.pt();
    pt2 = candidate2.pt();
    pz1 = candidate1.pz();
    pz2 = candidate2.pz();
    p1 = candidate1.p();
    p2 = candidate2.p();
    angle = TMath::ACos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    if (isDeepAngle && angle < cfgDeepAngle) {
      return false;
    }
    return true;
  }

  ROOT::Math::PxPyPzMVector KaonPlus, KaonMinus, PhiMesonMother, PhiVectorDummy, Phid1dummy, Phid2dummy;
  void processPhiReducedTable(EventCandidates::iterator const& collision, TrackCandidates const&, aod::BCsWithTimestamps const&)
  {
    o2::aod::ITSResponse itsResponse;
    bool keepEventDoublePhi = false;
    int numberPhi = 0;
    auto currentRunNumber = collision.bc_as<aod::BCsWithTimestamps>().runNumber();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
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

    if (collision.sel8() && collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(aod::evsel::kNoITSROFrameBorder) && collision.selection_bit(aod::evsel::kNoSameBunchPileup) && collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      for (auto track1 : posThisColl) {
        // track selection
        if (!selectionTrack(track1)) {
          continue;
        }
        // PID check
        if (!selectionPID(track1)) {
          continue;
        }
        if (!(itsResponse.nSigmaITS<o2::track::PID::Kaon>(track1) > -3.0 && itsResponse.nSigmaITS<o2::track::PID::Kaon>(track1) < 3.0)) {
          continue;
        }
        Npostrack = Npostrack + 1;
        qaRegistry.fill(HIST("hNsigmaPtkaonTPC"), track1.tpcNSigmaKa(), track1.pt());
        if (track1.hasTOF()) {
          qaRegistry.fill(HIST("hNsigmaPtkaonTOF"), track1.tofNSigmaKa(), track1.pt());
        }
        auto track1ID = track1.globalIndex();
        for (auto track2 : negThisColl) {
          // track selection
          if (!selectionTrack(track2)) {
            continue;
          }
          // PID check
          if (!selectionPID(track2)) {
            continue;
          }
          if (Npostrack == 1) {
            Nnegtrack = Nnegtrack + 1;
          }
          auto track2ID = track2.globalIndex();
          if (track2ID == track1ID) {
            continue;
          }
          if (!(itsResponse.nSigmaITS<o2::track::PID::Kaon>(track2) > -3.0 && itsResponse.nSigmaITS<o2::track::PID::Kaon>(track2) < 3.0)) {
            continue;
          }
          if (!selectionPair(track1, track2)) {
            continue;
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
    if (numberPhi > 1 && Npostrack > 1 && Nnegtrack > 1) {
      keepEventDoublePhi = true;
    }
    qaRegistry.fill(HIST("hEventstat"), 0.5);
    if (keepEventDoublePhi && numberPhi > 1 && Npostrack > 1 && Nnegtrack > 1 && (phiresonance.size() == phiresonanced1.size()) && (phiresonance.size() == phiresonanced2.size())) {
      qaRegistry.fill(HIST("hEventstat"), 1.5);
      /////////// Fill collision table///////////////
      redPhiEvents(bc.globalBC(), currentRunNumber, bc.timestamp(), collision.posZ(), collision.numContrib(), Npostrack, Nnegtrack);
      auto indexEvent = redPhiEvents.lastIndex();
      //// Fill track table for Phi//////////////////
      for (auto if1 = phiresonance.begin(); if1 != phiresonance.end(); ++if1) {
        auto i5 = std::distance(phiresonance.begin(), if1);
        PhiVectorDummy = phiresonance.at(i5);
        Phid1dummy = phiresonanced1.at(i5);
        Phid2dummy = phiresonanced2.at(i5);
        phiTrack(indexEvent, PhiVectorDummy.Px(), PhiVectorDummy.Py(), PhiVectorDummy.Pz(), Phid1dummy.Px(), Phid1dummy.Py(), Phid1dummy.Pz(), Phid2dummy.Px(), Phid2dummy.Py(), Phid2dummy.Pz(),
                 PhiVectorDummy.M(),
                 Phid1Index.at(i5), Phid2Index.at(i5),
                 Phid1Charge.at(i5), Phid2Charge.at(i5),
                 Phid1TPC.at(i5), Phid2TPC.at(i5), Phid1TOFHit.at(i5), Phid2TOFHit.at(i5), Phid1TOF.at(i5), Phid2TOF.at(i5));
      }
    }
  } // process
  PROCESS_SWITCH(doublephitable, processPhiReducedTable, "Process table creation for double phi", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<doublephitable>(cfg)};
}
