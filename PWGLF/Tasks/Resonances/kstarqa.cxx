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
/// \brief this is a code for the kstarqa resonance
/// \author prottay das, sawan
/// \since 13/03/2024

#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstdlib>

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct kstarqa {

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{
    "ccdb-no-later-than",
    std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::system_clock::now().time_since_epoch())
      .count(),
    "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080",
                                "url of the ccdb repository"};

  SliceCache cache;

  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection",
                                    {},
                                    OutputObjHandlingPolicy::AnalysisObject,
                                    true,
                                    true};
  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Confugrable for QA histograms
  Configurable<bool> QAbefore{"QAbefore", false, "QAbefore"};
  Configurable<bool> QAafter{"QAafter", false, "QAafter"};
  Configurable<bool> QAv0{"QAv0", false, "QAv0"};
  Configurable<bool> onlyTPC{"onlyTPC", false, "only TPC tracks"};

  // Configurable for event selection
  Configurable<float>
    cutzvertex{"cutzvertex", 10.0f,
               "Accepted z-vertex range (cm)"};

  // Configurables for track selections
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2f, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f,
                                  "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0,
                                   "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0,
                                        "Value of the Combined Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5,
                                     "Number of mixed events per event"};
  Configurable<bool> cfgMultFT0{"cfgMultFT0", false, "cfgMultFT0"};
  Configurable<bool> cfgCentFT0C{"cfgCentFT0C", true, "cfgCentFT0C"};
  Configurable<bool> iscustomDCAcut{"iscustomDCAcut", false, "iscustomDCAcut"};
  Configurable<bool> ismanualDCAcut{"ismanualDCAcut", true, "ismanualDCAcut"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  ConfigurableAxis cMixMultBins{"cMixMultBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f}, "Mixing bins - multiplicity"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};

  void init(InitContext const&)
  {
    // Axes
    AxisSpec vertexZAxis = {nBins, -10., 10., "vrtx_{Z} [cm]"};
    AxisSpec ptAxis = {200, 0.0f, 20.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec multAxis = {100, 0.0f, 100.0f, "Multiplicity"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec",
                        {HistType::kTH1F, {vertexZAxis}});
    rEventSelection.add("hmult", "Centrality distribution", kTH1F,
                        {{200, 0.0f, 200.0f}});

    // for primary tracks
    if (QAbefore && QAafter) {
      histos.add("hNsigmaPionTPC_before", "NsigmaPion TPC distribution before",
                 kTH1F, {{200, -10.0f, 10.0f}});
      histos.add("hNsigmaPionTOF_before", "NsigmaPion TOF distribution before",
                 kTH1F, {{200, -10.0f, 10.0f}});

      histos.add("hEta_after", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
      histos.add("hDcaxy_after", "Dcaxy distribution", kTH1F,
                 {{200, -10.0f, 10.0f}});
      histos.add("hDcaz_after", "Dcaz distribution", kTH1F,
                 {{200, -10.0f, 10.0f}});
      histos.add("hNsigmaPionTPC_after", "NsigmaPion TPC distribution", kTH1F,
                 {{200, -10.0f, 10.0f}});
      histos.add("hNsigmaPionTOF_after", "NsigmaPion TOF distribution", kTH1F,
                 {{200, -10.0f, 10.0f}});
    }

    // CKStar histograms
    histos.add("h3CKSInvMassUnlikeSign",
               "Invariant mass of CKS meson Unlike Sign", kTHnSparseF,
               {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {90, 0.6, 1.5}}, true);
    histos.add("h3CKSInvMasslikeSign",
               "Invariant mass of CKS meson like Sign", kTHnSparseF,
               {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {90, 0.6, 1.5}}, true);
    histos.add("h3CKSInvMassMixed", "Invariant mass of CKS meson Mixed",
               kTHnSparseF,
               {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {90, 0.6, 1.5}}, true);
  }

  double massPi = TDatabasePDG::Instance()
                    ->GetParticle(kPiPlus)
                    ->Mass(); // FIXME: Get from the common header
  double massKa = o2::constants::physics::MassKPlus;
  ROOT::Math::PtEtaPhiMVector CKSVector;

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (iscustomDCAcut &&
        (!candidate.isGlobalTrack() || !candidate.isPVContributor() ||
         candidate.itsNCls() < cfgITScluster)) {
      return false;
    }
    if (ismanualDCAcut &&
        !(candidate.isGlobalTrackWoDCA() && candidate.isPVContributor() &&
          std::abs(candidate.dcaXY()) < cfgCutDCAxy &&
          std::abs(candidate.dcaZ()) < cfgCutDCAz &&
          candidate.itsNCls() > cfgITScluster)) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool selectionPID(const T& candidate, int PID)
  {
    if (PID == 0) {
      if (candidate.hasTOF() &&
          (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() +
           candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) <
            (nsigmaCutCombined * nsigmaCutCombined)) {
        return true;
      }
      if (onlyTPC) {
        if (!candidate.hasTOF() &&
            std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
          return true;
        }
      }
    } else if (PID == 1) {
      if (candidate.hasTOF() &&
          (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() +
           candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) <
            (nsigmaCutCombined * nsigmaCutCombined)) {
        return true;
      }
      if (onlyTPC) {
        if (!candidate.hasTOF() &&
            std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
      }
    }
    return false;
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection
  // requirements
  // Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);

  Filter acceptanceFilter =
    (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) &&
                        (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<
    soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs,
              aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;
  using TrackCandidates = soa::Filtered<
    soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
              aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa>>;
  using V0TrackCandidate = aod::V0Datas;

  ConfigurableAxis axisVertex{
    "axisVertex",
    {20, -10, 10},
    "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{
    "axisMultiplicityClass",
    {10, 0, 100},
    "multiplicity percentile for bin"};
  ConfigurableAxis axisMultiplicity{
    "axisMultiplicity",
    {2000, 0, 10000},
    "TPC multiplicity  for bin"};

  using BinningTypeTPCMultiplicity =
    ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  // using BinningTypeVertexContributor =
  // ColumnBinningPolicy<aod::collision::PosZ, aod::collision::NumContrib>;
  using BinningTypeCentralityM =
    ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypeVertexContributor =
    ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;

  BinningTypeVertexContributor binningOnPositions{
    {axisVertex, axisMultiplicity},
    true};

  SameKindPair<EventCandidates, TrackCandidates,
               BinningTypeVertexContributor>
    pair{binningOnPositions, cfgNoMixedEvents, -1, &cache};

  void processSE(EventCandidates::iterator const& collision,
                 TrackCandidates const& tracks,
                 aod::BCs const&)

  {

    if (!collision.sel8()) {
      return;
    }

    if (timFrameEvsel && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return;
    }

    std::vector<ROOT::Math::PtEtaPhiMVector> pions, kaons;
    std::vector<int64_t> PionIndex = {};
    std::vector<int64_t> KaonIndex = {};
    std::vector<int64_t> PioncollIndex = {};
    std::vector<int64_t> KaoncollIndex = {};
    std::vector<int64_t> PionSign = {};
    std::vector<int64_t> KaonSign = {};

    float multiplicity = 0.0f;
    /*  if (cfgMultFT0)
      multiplicity = collision.multZeqFT0A() + collision.multZeqFT0C();
    if (cfgMultFT0 == 0 && cfgCentFT0C == 1)
      multiplicity = collision.centFT0C();
      if (cfgMultFT0 == 0 && cfgCentFT0C == 0)*/
    multiplicity = collision.centFT0M();

    // Fill the event counter
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    rEventSelection.fill(HIST("hmult"), multiplicity);

    for (auto track1 : tracks) {
      /*
      if (QAbefore) {
        histos.fill(HIST("hNsigmaPionTPC_before"), track1.tpcNSigmaPi());
        histos.fill(HIST("hNsigmaPionTOF_before"), track1.tofNSigmaPi());
      }
      */

      if (!selectionPID(track1, 1))
        continue; // for primary particle PID

      if (!selectionTrack(track1)) {
        continue;
      }

      /*
      if (QAafter) {
        histos.fill(HIST("hEta_after"), track1.eta());
        histos.fill(HIST("hDcaxy_after"), track1.dcaXY());
        histos.fill(HIST("hDcaz_after"), track1.dcaZ());
        histos.fill(HIST("hNsigmaPionTPC_after"), track1.tpcNSigmaPi());
        histos.fill(HIST("hNsigmaPionTOF_after"), track1.tofNSigmaPi());
      }
      */
      ROOT::Math::PtEtaPhiMVector temp1(track1.pt(), track1.eta(), track1.phi(),
                                        massKa);
      kaons.push_back(temp1);
      KaonIndex.push_back(track1.globalIndex());
      KaoncollIndex.push_back(track1.collisionId());
      KaonSign.push_back(track1.sign());

    } // track loop ends

    for (auto track2 : tracks) {
      /*
      if (QAbefore) {
        histos.fill(HIST("hNsigmaPionTPC_before"), track1.tpcNSigmaPi());
        histos.fill(HIST("hNsigmaPionTOF_before"), track1.tofNSigmaPi());
      }
      */

      if (!selectionPID(track2, 0))
        continue; // for primary particle PID

      if (!selectionTrack(track2)) {
        continue;
      }

      /*
      if (QAafter) {
        histos.fill(HIST("hEta_after"), track1.eta());
        histos.fill(HIST("hDcaxy_after"), track1.dcaXY());
        histos.fill(HIST("hDcaz_after"), track1.dcaZ());
        histos.fill(HIST("hNsigmaPionTPC_after"), track1.tpcNSigmaPi());
        histos.fill(HIST("hNsigmaPionTOF_after"), track1.tofNSigmaPi());
      }
      */
      ROOT::Math::PtEtaPhiMVector temp2(track2.pt(), track2.eta(), track2.phi(),
                                        massPi);
      pions.push_back(temp2);
      PionIndex.push_back(track2.globalIndex());
      PioncollIndex.push_back(track2.collisionId());
      PionSign.push_back(track2.sign());

    } // track loop ends

    if (pions.size() != 0 && kaons.size() != 0) {
      for (auto ikaon = kaons.begin(); ikaon != kaons.end(); ++ikaon) {
        auto i1 = std::distance(kaons.begin(), ikaon);
        for (auto ipion = pions.begin(); ipion != pions.end();
             ++ipion) {
          auto i3 = std::distance(pions.begin(), ipion);

          if (PionIndex.at(i3) <= KaonIndex.at(i1))
            continue;
          CKSVector = kaons.at(i1) + pions.at(i3);

          if (TMath::Abs(CKSVector.Rapidity()) < 0.5) {
            if (PionSign.at(i3) * KaonSign.at(i1) < 0)
              histos.fill(HIST("h3CKSInvMassUnlikeSign"), multiplicity,
                          CKSVector.Pt(), CKSVector.M());
            else if (PionSign.at(i3) * KaonSign.at(i1) > 0)
              histos.fill(HIST("h3CKSInvMasslikeSign"), multiplicity,
                          CKSVector.Pt(), CKSVector.M());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(kstarqa, processSE, "Process Same event", true);

  void processME(EventCandidates const& collisions,
                 TrackCandidates const& tracks)

  {

    for (auto& [c1, tracks1, c2, tracks2] : pair) {

      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }

      if (timFrameEvsel && (!c1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c2.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
        return;
      }

      float multiplicity = 0.0f;
      /*      if (cfgMultFT0)
        multiplicity = c1.multZeqFT0A() + c1.multZeqFT0C();
      if (cfgMultFT0 == 0 && cfgCentFT0C == 1)
        multiplicity = c1.centFT0C();
  if (cfgMultFT0 == 0 && cfgCentFT0C == 0)*/
      multiplicity = c1.centFT0M();

      for (auto& [t1, t2] : o2::soa::combinations(
             o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (!selectionTrack(t1))
          continue;
        if (!selectionTrack(t2))
          continue;
        if (!selectionPID(t1, 1))
          continue;
        if (!selectionPID(t2, 0))
          continue;

        TLorentzVector KAON;
        KAON.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), massKa);
        TLorentzVector PION;
        PION.SetPtEtaPhiM(t2.pt(), t2.eta(), t2.phi(), massPi);

        TLorentzVector CKSmix = KAON + PION;

        if (TMath::Abs(CKSmix.Rapidity()) < 0.5) {
          if (t1.sign() * t2.sign() < 0)
            histos.fill(HIST("h3CKSInvMassMixed"), multiplicity, CKSmix.Pt(),
                        CKSmix.M());
        }
      }
    }
  }

  PROCESS_SWITCH(kstarqa, processME, "Process Mixed event", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<kstarqa>(cfgc)};
}
