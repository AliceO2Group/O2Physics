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

/// \file kstarzeroanalysis.cxx
/// \recostruction f K*(892) resonance
///
///
/// \author Marta Urioni <marta.urioni@cern.ch>

#include <TLorentzVector.h>
#include "TF1.h"

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "DataFormatsParameters/GRPObject.h"
#include "CommonConstants/PhysicsConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

struct kstarzeroanalysis {
  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  // bins and axis
  ConfigurableAxis CfgPtAxis{"CfgPtAxis", {50, 0.0f, 10.f}, "Pt axis"};
  ConfigurableAxis CfgMinvAxis{"CfgMinvAxis", {300, 0.6f, 1.5f}, "Invariant mass axis"};
  ConfigurableAxis CfgCentrAxis{"CfgCentrAxis", {110, 0.0f, 110.0f}, "centrality axis"};

  // mixed event
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.f, 1.f, 5.f, 10.f, 20.f, 30.f, 50.f, 70.f, 100.f, 110.f}, "multiplicity bins"};
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};

  // event selections
  Configurable<int> cMaxCentr{"cMaxCentr", 80, "Maximum centrality of event"};

  // track cuts
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.5, "Track maximum DCAr to PV cut"};
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track maximum DCAz to PV cut"};
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};           // PV Contributor
  Configurable<float> cfgITScluster{"cfgITScluster", 0, "ITS min number of clusters"};
  Configurable<float> cfgTPCcluster{"cfgTPCcluster", 70, "TPC min number of clusters"};
  Configurable<float> cfgRCRFC{"cfgRCRFC", 0.8f, "TPC minimum number of crossed raws over findable"};

  // PID selection
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutTOF{"nsigmaCutTOF", 3.0, "Value of the TOF Nsigma cut"};

  // cuts on mother
  Configurable<float> cRapidityMotherCut{"cRapidityMotherCut", 0.5, "Maximum rapidity of the mother"};
  Configurable<double> cMaxPtMotherCut{"cMaxPtMotherCut", 15.0, "Maximum pt of mother cut"};
  Configurable<double> cMaxMinvMotherCut{"cMaxMinvMotherCut", 1.5, "Maximum Minv of mother cut"};

  // constants
  double massKa = MassKPlus;
  double massPi = MassPiPlus;

  void init(o2::framework::InitContext&)
  {
    // Test on Mixed event
    histos.add("TestME/hCollisionIndexSameE", "coll index sameE", HistType::kTH1F, {{500, 0.0f, 500.0f}});
    histos.add("TestME/hCollisionIndexMixedE", "coll index mixedE", HistType::kTH1F, {{500, 0.0f, 500.0f}});
    histos.add("TestME/hnTrksSameE", "n tracks per event SameE", HistType::kTH1F, {{1000, 0.0f, 1000.0f}});
    histos.add("TestME/hnTrksMixedE", "n tracks per event MixedE", HistType::kTH1F, {{1000, 0.0f, 1000.0f}});
    histos.add("TestME/hPairsCounterSameE", "tot n pairs sameE", HistType::kTH1F, {{1, 0.5f, 1.5f}});
    histos.add("TestME/hPairsCounterMixedE", "tot n pairs mixedE", HistType::kTH1F, {{1, 0.5f, 1.5f}});

    // event histograms
    histos.add("QAevent/hEvtCounterSameE", "Number of analyzed Same Events", HistType::kTH1F, {{1, 0.5, 1.5}});
    histos.add("QAevent/hVertexZSameE", "Collision Vertex Z position", HistType::kTH1F, {{100, -15., 15.}});
    histos.add("QAevent/hMultiplicityPercentSameE", "Multiplicity percentile of collision", HistType::kTH1F, {{120, 0.0f, 120.0f}});

    histos.add("QAevent/hEvtCounterMixedE", "Number of analyzed Mixed Events", HistType::kTH1F, {{1, 0.5, 1.5}});
    histos.add("QAevent/hVertexZMixedE", "Collision Vertex Z position", HistType::kTH1F, {{100, -15., 15.}});
    histos.add("QAevent/hMultiplicityPercentMixedE", "Multiplicity percentile of collision", HistType::kTH1F, {{120, 0.0f, 120.0f}});

    // track QA before cuts
    histos.add("QAbeforeCuts/hEta_BC", "Eta distribution before cuts", kTH1F, {{250, -1.0f, 1.0f}});
    histos.add("QAbeforeCuts/hDcaxy_BC", "DCAxy distribution before cuts", kTH1F, {{250, -1.0f, 1.0f}});
    histos.add("QAbeforeCuts/hDcaz_BC", "DCAz distribution before  cuts", kTH1F, {{250, -1.0f, 1.0f}});

    // kaon QA histograms
    histos.add("QAafterCuts/hEtaKaon_AC", "Eta distribution Kaon after cuts", kTH1F, {{250, -1.0f, 1.0f}});
    histos.add("QAafterCuts/hDcaxyKaon_AC", "DCAxy Kaon distribution after cuts", kTH1F, {{250, -1.0f, 1.0f}});
    histos.add("QAafterCuts/hDcazKaon_AC", "DCAz Kaon distribution after  cuts", kTH1F, {{250, -1.0f, 1.0f}});

    histos.add("QAbeforeCuts/hNsigmaKaonTPC_BC", "NsigmaTPC Kaon before cuts", kTH2F, {CfgPtAxis, {250, -6.0f, 6.0f}});
    histos.add("QAbeforeCuts/hNsigmaKaonTOF_BC", "NsigmaTOF Kaon before cuts", kTH2F, {CfgPtAxis, {250, -6.0f, 6.0f}});
    histos.add("QAafterCuts/hNsigmaKaonTPC_AC", "NsigmaTPC Kaon after cuts", kTH2F, {CfgPtAxis, {250, -6.0f, 6.0f}});
    histos.add("QAafterCuts/hNsigmaKaonTOF_AC", "NsigmaTOF Kaon after cuts", kTH2F, {CfgPtAxis, {250, -6.0f, 6.0f}});

    // Pion QA histograms
    histos.add("QAafterCuts/hEtaPion_AC", "Eta distribution Pion after cuts", kTH1F, {{250, -1.0f, 1.0f}});
    histos.add("QAafterCuts/hDcaxyPion_AC", "DCAxy Pion distribution after cuts", kTH1F, {{250, -1.0f, 1.0f}});
    histos.add("QAafterCuts/hDcazPion_AC", "DCAz Pion distribution after  cuts", kTH1F, {{250, -1.0f, 1.0f}});

    histos.add("QAbeforeCuts/hNsigmaPionTPC_BC", "NsigmaTPC Pion before cuts", kTH2F, {CfgPtAxis, {250, -6.0f, 6.0f}});
    histos.add("QAbeforeCuts/hNsigmaPionTOF_BC", "NsigmaTOF Pion before cuts", kTH2F, {CfgPtAxis, {250, -6.0f, 6.0f}});
    histos.add("QAafterCuts/hNsigmaPionTPC_AC", "NsigmaTPC Pion after cuts", kTH2F, {CfgPtAxis, {250, -6.0f, 6.0f}});
    histos.add("QAafterCuts/hNsigmaPionTOF_AC", "NsigmaTOF Pion after cuts", kTH2F, {CfgPtAxis, {250, -6.0f, 6.0f}});

    // invariant mass istograms
    histos.add("h1Kstar0InvMassUnlikeSignPM", "Invariant mass of Kstar0 meson Unlike Sign Kaon+ Pion-", kTH1F, {CfgMinvAxis});
    histos.add("h1Kstar0InvMassUnlikeSignMP", "Invariant mass of Kstar0 meson Unlike Sign Kaon- Pion+", kTH1F, {CfgMinvAxis});
    histos.add("h1Kstar0InvMassLikeSignPP", "Invariant mass of Kstar0 meson Like Sign positive", kTH1F, {CfgMinvAxis});
    histos.add("h1Kstar0InvMassLikeSignMM", "Invariant mass of Kstar0 meson Like Sign negative", kTH1F, {CfgMinvAxis});
    histos.add("h1Kstar0InvMassMixedUnlikeSignPM", "Invariant mass of Kstar0 meson Unlike Sign Kaon+ Pion- for mixedevent", kTH1F, {CfgMinvAxis});
    histos.add("h1Kstar0InvMassMixedUnlikeSignMP", "Invariant mass of Kstar0 meson Unlike Sign Kaon- Pion+ for mixedevent", kTH1F, {CfgMinvAxis});

    histos.add("h3Kstar0InvMassUnlikeSignPM", "Invariant mass of Kstar0 meson Unlike Sign Kaon+ Pion-", kTH3F, {CfgCentrAxis, CfgPtAxis, CfgMinvAxis});
    histos.add("h3Kstar0InvMassUnlikeSignMP", "Invariant mass of Kstar0 meson Unlike Sign Kaon- Pion+", kTH3F, {CfgCentrAxis, CfgPtAxis, CfgMinvAxis});
    histos.add("h3Kstar0InvMassLikeSignPP", "Invariant mass of Kstar0 meson Like Sign positive", kTH3F, {CfgCentrAxis, CfgPtAxis, CfgMinvAxis});
    histos.add("h3Kstar0InvMassLikeSignMM", "Invariant mass of Kstar0 meson Like Sign negative", kTH3F, {CfgCentrAxis, CfgPtAxis, CfgMinvAxis});
    histos.add("h3Kstar0InvMassMixedUnlikeSignPM", "Invariant mass of Kstar0 meson Unlike Sign Kaon+ Pion- for mixedevent", kTH3F, {CfgCentrAxis, CfgPtAxis, CfgMinvAxis});
    histos.add("h3Kstar0InvMassMixedUnlikeSignMP", "Invariant mass of Kstar0 meson Unlike Sign Kaon- Pion+ for mixedevent", kTH3F, {CfgCentrAxis, CfgPtAxis, CfgMinvAxis});
  }

  template <typename TCollision>
  bool eventSelection(TCollision collision, const float& centrality)
  {
    if (centrality > cMaxCentr)
      return false;

    return true;
  }

  template <typename TrackType>
  bool trackCut(const TrackType track)
  {
    if (std::abs(track.pt()) < cMinPtcut)
      return false;
    if (std::abs(track.dcaXY()) > cMaxDCArToPVcut)
      return false;
    if (std::abs(track.dcaZ()) > cMaxDCAzToPVcut)
      return false;
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgGlobalWoDCATrack && !(track.isGlobalTrackWoDCA() && track.itsNCls() > cfgITScluster && track.tpcNClsFound() > cfgTPCcluster && track.tpcCrossedRowsOverFindableCls() > cfgRCRFC))
      return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;

    return true;
  }

  template <typename T>
  bool selectionPIDKaon(const T& candidate)
  {
    bool tpcPass = std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC;
    bool tofPass = (candidate.hasTOF()) ? std::abs(candidate.tofNSigmaKa()) < nsigmaCutTOF : true;
    if (tpcPass && tofPass) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selectionPIDPion(const T& candidate)
  {
    bool tpcPass = std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC;
    bool tofPass = (candidate.hasTOF()) ? std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOF : true;
    if (tpcPass && tofPass) {
      return true;
    }
    return false;
  }

  double rapidity, mass, pT, paircharge;
  TLorentzVector daughter1, daughter2, mother;
  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2)
  {
    auto multiplicity = collision.cent();
    if (!eventSelection(collision, multiplicity)) // selection on event
      return;

    if constexpr (!IsMix) {
      histos.fill(HIST("TestME/hCollisionIndexSameE"), collision.globalIndex());
      histos.fill(HIST("TestME/hnTrksSameE"), dTracks1.size());
    } else {
      histos.fill(HIST("TestME/hCollisionIndexMixedE"), collision.globalIndex());
      histos.fill(HIST("TestME/hnTrksMixedE"), dTracks1.size());
    }

    for (auto const& [trkKa, trkPi] : soa::combinations(soa::CombinationsFullIndexPolicy(dTracks1, dTracks2))) {
      if constexpr (!IsMix) {
        histos.fill(HIST("TestME/hPairsCounterSameE"), 1.0);
      } else {
        histos.fill(HIST("TestME/hPairsCounterMixedE"), 1.0);
      }

      if (trkKa.index() == trkPi.index()) // Exclude same index tracks.
        continue;

      if (!trackCut(trkKa) || !trackCut(trkPi)) // track selections
        continue;

      if constexpr (!IsMix) { // Same event
        // QA plots before cuts
        histos.fill(HIST("QAbeforeCuts/hEta_BC"), trkKa.eta());
        histos.fill(HIST("QAbeforeCuts/hDcaxy_BC"), trkKa.dcaXY());
        histos.fill(HIST("QAbeforeCuts/hDcaz_BC"), trkKa.dcaZ());
        histos.fill(HIST("QAbeforeCuts/hNsigmaKaonTPC_BC"), trkKa.pt(), trkKa.tpcNSigmaKa());
        histos.fill(HIST("QAbeforeCuts/hNsigmaPionTPC_BC"), trkKa.pt(), trkKa.tpcNSigmaPi());
        if (trkKa.hasTOF()) {
          histos.fill(HIST("QAbeforeCuts/hNsigmaKaonTOF_BC"), trkKa.pt(), trkKa.tofNSigmaKa());
          histos.fill(HIST("QAbeforeCuts/hNsigmaPionTOF_BC"), trkKa.pt(), trkKa.tofNSigmaPi());
        }
      }

      if (!selectionPIDKaon(trkKa) || !selectionPIDPion(trkPi)) // PID selections
        continue;

      if constexpr (!IsMix) { // Same event
        // QA plots Kaon
        histos.fill(HIST("QAafterCuts/hEtaKaon_AC"), trkKa.eta());
        histos.fill(HIST("QAafterCuts/hDcaxyKaon_AC"), trkKa.dcaXY());
        histos.fill(HIST("QAafterCuts/hDcazKaon_AC"), trkKa.dcaZ());
        histos.fill(HIST("QAafterCuts/hNsigmaKaonTPC_AC"), trkKa.pt(), trkKa.tpcNSigmaKa());
        if (trkKa.hasTOF()) {
          histos.fill(HIST("QAafterCuts/hNsigmaKaonTOF_AC"), trkKa.pt(), trkKa.tofNSigmaKa());
        }
        // QA plots Pion
        histos.fill(HIST("QAafterCuts/hEtaPion_AC"), trkPi.eta());
        histos.fill(HIST("QAafterCuts/hDcaxyPion_AC"), trkPi.dcaXY());
        histos.fill(HIST("QAafterCuts/hDcazPion_AC"), trkPi.dcaZ());
        histos.fill(HIST("QAafterCuts/hNsigmaPionTPC_AC"), trkPi.pt(), trkPi.tpcNSigmaPi());
        if (trkPi.hasTOF()) {
          histos.fill(HIST("QAafterCuts/hNsigmaPionTOF_AC"), trkPi.pt(), trkPi.tofNSigmaPi());
        }
      }

      daughter1.SetXYZM(trkKa.px(), trkKa.py(), trkKa.pz(), massKa); // set the daughter1 4-momentum
      daughter2.SetXYZM(trkPi.px(), trkPi.py(), trkPi.pz(), massPi); // set the daughter2 4-momentum
      mother = daughter1 + daughter2;                                // calculate the mother 4-momentum
      mass = mother.M();
      pT = mother.Pt();
      rapidity = mother.Rapidity();

      if (std::abs(rapidity) > cRapidityMotherCut)
        continue;                // rapidity cut
      if (pT >= cMaxPtMotherCut) // excluding candidates in overflow
        continue;
      if (mass >= cMaxMinvMotherCut) // excluding candidates in overflow
        continue;

      if constexpr (!IsMix) {                     // Same event
        if (trkKa.sign() > 0 && trkPi.sign() < 0) // unlike sign PM
        {
          histos.fill(HIST("h3Kstar0InvMassUnlikeSignPM"), multiplicity, pT, mass);
          histos.fill(HIST("h1Kstar0InvMassUnlikeSignPM"), mass);
        } else if (trkKa.sign() < 0 && trkPi.sign() > 0) // unlike sign MP
        {
          histos.fill(HIST("h3Kstar0InvMassUnlikeSignMP"), multiplicity, pT, mass);
          histos.fill(HIST("h1Kstar0InvMassUnlikeSignMP"), mass);
        } else if (trkKa.sign() > 0 && trkPi.sign() > 0) // like sign PP
        {
          histos.fill(HIST("h3Kstar0InvMassLikeSignPP"), multiplicity, pT, mass);
          histos.fill(HIST("h1Kstar0InvMassLikeSignPP"), mass);
        } else // like sign MM
        {
          histos.fill(HIST("h3Kstar0InvMassLikeSignMM"), multiplicity, pT, mass);
          histos.fill(HIST("h1Kstar0InvMassLikeSignMM"), mass);
        }
      } else {                                    // Mixed event
        if (trkKa.sign() > 0 && trkPi.sign() < 0) // unlike sign PM
        {
          histos.fill(HIST("h3Kstar0InvMassMixedUnlikeSignPM"), multiplicity, pT, mass);
          histos.fill(HIST("h1Kstar0InvMassMixedUnlikeSignPM"), mass);
        } else if (trkKa.sign() < 0 && trkPi.sign() > 0) // unlike sign MP
        {
          histos.fill(HIST("h3Kstar0InvMassMixedUnlikeSignMP"), multiplicity, pT, mass);
          histos.fill(HIST("h1Kstar0InvMassMixedUnlikeSignMP"), mass);
        }
      }
    } // end loop on tracks
  }   // end fillHistograms

  // Process the data
  void processSE(aod::ResoCollision& collision, aod::ResoTracks const& resotracks)
  {
    // Fill the event counter
    histos.fill(HIST("QAevent/hEvtCounterSameE"), 1.0);
    histos.fill(HIST("QAevent/hVertexZSameE"), collision.posZ());
    histos.fill(HIST("QAevent/hMultiplicityPercentSameE"), collision.cent());
    fillHistograms<false, false>(collision, resotracks, resotracks); // Fill histograms, no MC, no mixing
  }
  PROCESS_SWITCH(kstarzeroanalysis, processSE, "Process same event", false);

  // Processing Event Mixing
  using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processME(o2::aod::ResoCollisions& collisions, aod::ResoTracks const& resotracks)
  {
    auto tracksTuple = std::make_tuple(resotracks);
    BinningTypeVtxZT0M colBinning{{CfgVtxBins, CfgMultBins}, true};
    SameKindPair<aod::ResoCollisions, aod::ResoTracks, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (auto& [collision1, tracks1, collision2, tracks2] : pairs) { // loop over all pairs
      // Fill the event counter
      histos.fill(HIST("QAevent/hEvtCounterMixedE"), 1.0);
      histos.fill(HIST("QAevent/hVertexZMixedE"), collision1.posZ());
      histos.fill(HIST("QAevent/hMultiplicityPercentMixedE"), collision1.cent());
      fillHistograms<false, true>(collision1, tracks1, tracks2); // Fill histograms, no MC, mixing
    }
  }
  PROCESS_SWITCH(kstarzeroanalysis, processME, "Process EventMixing for combinatorial background", false); // Event Mixing
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<kstarzeroanalysis>(cfgc, TaskName{"lf-kstarzeroanalysis"})};
}
