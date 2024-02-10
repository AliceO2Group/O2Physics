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

/// \file skimCreatorMini.cxx
/// \brief Mini version of the track index skim creator
///
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

// O2
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

// O2Physics
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"

// PWGHF
#include "Tutorials/PWGHF/DataModelMini.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Track selection =====================================================================

/// Track selection
struct HfTrackIndexSkimCreatorTagSelTracks {
  Produces<aod::HfSelTrack> rowSelectedTrack;

  // 2-prong cuts
  Configurable<double> ptTrackMin{"ptTrackMin", 0.3, "min. track pT for 2 prong candidate"};
  Configurable<double> etaTrackMax{"etaTrackMax", 0.8, "max. pseudorapidity for 2 prong candidate"};
  Configurable<double> dcaTrackMin{"dcaTrackMin", 0.0025, "min. DCA for 2 prong candidate"};

  using TracksWDcaSel = soa::Join<aod::Tracks, aod::TracksDCA, aod::TrackSelection>;

  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext&)
  {
    const TString strTitle = "D^{0} candidates";
    const TString strPt = "#it{p}_{T}^{track} (GeV/#it{c})";
    const TString strEntries = "entries";
    registry.add("hPtNoCuts", "all tracks;" + strPt + ";" + strEntries, {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("hPtCuts2Prong", "tracks selected for 2-prong vertexing;" + strPt + ";" + strEntries, {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("hPtVsDcaXYToPvCuts2Prong", "tracks selected for 2-prong vertexing;" + strPt + ";" + "DCAxy to prim. vtx. (cm)" + ";" + strEntries, {HistType::kTH2F, {{100, 0., 10.}, {400, -2., 2.}}});
    registry.add("hEtaCuts2Prong", "tracks selected for 2-prong vertexing;#it{#eta};" + strEntries, {HistType::kTH1F, {{static_cast<int>(1.2 * etaTrackMax * 100), -1.2 * etaTrackMax, 1.2 * etaTrackMax}}});
  }

  void process(TracksWDcaSel const& tracks)
  {
    for (const auto& track : tracks) {
      bool statusProng = true;

      auto ptTrack = track.pt();
      registry.fill(HIST("hPtNoCuts"), ptTrack);

      // pT cut
      if (ptTrack < ptTrackMin) {
        statusProng = false;
      }

      // eta cut
      auto etaTrack = track.eta();
      if (statusProng && std::abs(etaTrack) > etaTrackMax) {
        statusProng = false;
      }

      // quality cut
      if (!track.isGlobalTrackWoDCA()) {
        statusProng = false;
      }

      // DCA cut
      auto dcaXY = track.dcaXY();
      if (statusProng && std::abs(dcaXY) < dcaTrackMin) {
        statusProng = false;
      }

      // fill histograms
      if (statusProng) {
        registry.fill(HIST("hPtCuts2Prong"), ptTrack);
        registry.fill(HIST("hEtaCuts2Prong"), etaTrack);
        registry.fill(HIST("hPtVsDcaXYToPvCuts2Prong"), ptTrack, dcaXY);
      }

      // fill table row
      rowSelectedTrack(statusProng);
    }
  }
};

// Track index skimming =====================================================================

/// Track index skim creator
/// Pre-selection of 2-prong secondary vertices
struct HfTrackIndexSkimCreator {
  Produces<aod::HfTrackIndexProng2> rowTrackIndexProng2;

  // vertexing parameters
  Configurable<double> magneticField{"magneticField", 5., "magnetic field [kG]"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};

  o2::vertexing::DCAFitterN<2> fitter; // 2-prong vertex fitter

  using SelectedTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksCov, aod::HfSelTrack>>;

  Filter filterSelectTracks = aod::hf_seltrack::isSelProng == true;

  HistogramRegistry registry{
    "registry",
    {// 2-prong histograms
     {"hVtx2ProngX", "2-prong candidates;#it{x}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -2., 2.}}}},
     {"hVtx2ProngY", "2-prong candidates;#it{y}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -2., 2.}}}},
     {"hVtx2ProngZ", "2-prong candidates;#it{z}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -20., 20.}}}},
     {"hMassD0ToPiK", "D^{0} candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}}}}};

  void init(InitContext&)
  {
    // Configure the vertexer
    fitter.setBz(magneticField);
    fitter.setPropagateToPCA(propagateToPCA);
    fitter.setMaxR(maxR);
    fitter.setMaxDZIni(maxDZIni);
    fitter.setMinParamChange(minParamChange);
    fitter.setMinRelChi2Change(minRelChi2Change);
    fitter.setUseAbsDCA(useAbsDCA);
  }

  void process(aod::Collision const&,
               SelectedTracks const& tracks)
  {
    // loop over positive tracks
    for (const auto& trackPos1 : tracks) {
      if (trackPos1.signed1Pt() < 0) {
        continue;
      }
      auto trackParVarPos1 = getTrackParCov(trackPos1);

      // loop over negative tracks
      for (const auto& trackNeg1 : tracks) {
        if (trackNeg1.signed1Pt() > 0) {
          continue;
        }
        auto trackParVarNeg1 = getTrackParCov(trackNeg1);

        // secondary vertex reconstruction and further 2-prong selections
        if (fitter.process(trackParVarPos1, trackParVarNeg1) == 0) {
          continue;
        }
        //  get secondary vertex
        const auto& secondaryVertex = fitter.getPCACandidate();
        // get track momenta
        std::array<float, 3> pVec0;
        std::array<float, 3> pVec1;
        fitter.getTrack(0).getPxPyPzGlo(pVec0);
        fitter.getTrack(1).getPxPyPzGlo(pVec1);

        // fill table row
        rowTrackIndexProng2(trackPos1.globalIndex(),
                            trackNeg1.globalIndex());

        // fill histograms
        registry.fill(HIST("hVtx2ProngX"), secondaryVertex[0]);
        registry.fill(HIST("hVtx2ProngY"), secondaryVertex[1]);
        registry.fill(HIST("hVtx2ProngZ"), secondaryVertex[2]);
        std::array<std::array<float, 3>, 2> arrMom = {pVec0, pVec1};
        auto mass2Prong = RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus});
        registry.fill(HIST("hMassD0ToPiK"), mass2Prong);
      }
    }
  }
};

// Add all tasks in the workflow specification.
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTrackIndexSkimCreatorTagSelTracks>(cfgc),
    adaptAnalysisTask<HfTrackIndexSkimCreator>(cfgc)};
}
