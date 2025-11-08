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

#include "Tutorials/PWGHF/DataModelMini.h"
//
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TString.h>

#include <array>
#include <cstdlib>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

// Track selection =====================================================================

/// Track selection
struct HfSkimCreatorMiniTagSelTracks {
  Produces<aod::HfTSelTrack> rowSelectedTrack;

  // 2-prong cuts
  Configurable<float> ptTrackMin{"ptTrackMin", 0.3f, "Min. track pT for 2-prong candidate [GeV/c]"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8f, "Max. pseudorapidity for 2-prong candidate"};
  Configurable<float> dcaTrackMin{"dcaTrackMin", 0.0025f, "Min. DCA for 2-prong candidate [cm]"};

  using TracksWDcaSel = soa::Join<aod::Tracks, aod::TracksDCA, aod::TrackSelection>;

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    const TString strPt = "#it{p}_{T}^{track} (GeV/#it{c})";
    const TString strEntries = "entries";
    registry.add("hPtNoCuts", "all tracks;" + strPt + ";" + strEntries, {HistType::kTH1D, {{100, 0., 10.}}});
    registry.add("hPtCuts2Prong", "tracks selected for 2-prong vertexing;" + strPt + ";" + strEntries, {HistType::kTH1D, {{100, 0., 10.}}});
    registry.add("hPtVsDcaXYToPvCuts2Prong", "tracks selected for 2-prong vertexing;" + strPt + ";" + "DCAxy to prim. vtx. (cm)" + ";" + strEntries, {HistType::kTH2F, {{100, 0., 10.}, {400, -2., 2.}}});
    registry.add("hEtaCuts2Prong", "tracks selected for 2-prong vertexing;#it{#eta};" + strEntries, {HistType::kTH1D, {{static_cast<int>(1.2 * etaTrackMax * 100), -1.2 * etaTrackMax, 1.2 * etaTrackMax}}});
  }

  void process(TracksWDcaSel::iterator const& track)
  {
    const auto ptTrack{track.pt()};
    registry.fill(HIST("hPtNoCuts"), ptTrack);

    bool statusProng{true};

    // pT cut
    if (ptTrack < ptTrackMin) {
      statusProng = false;
    }

    // eta cut
    const auto etaTrack{track.eta()};
    if (statusProng && std::abs(etaTrack) > etaTrackMax) {
      statusProng = false;
    }

    // quality cut
    if (!track.isGlobalTrackWoDCA()) {
      statusProng = false;
    }

    // DCA cut
    const auto dcaXY{track.dcaXY()};
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
};

// Track index skimming =====================================================================

/// Track index skim creator
/// Pre-selection of 2-prong secondary vertices
struct HfSkimCreatorMini {
  Produces<aod::HfT2Prongs> rowTrackIndexProng2;

  // vertexing parameters
  Configurable<float> magneticField{"magneticField", 5.f, "magnetic field [kG]"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "Create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<float> maxR{"maxR", 200.f, "Reject PCA's above this radius"};
  Configurable<float> maxDZIni{"maxDZIni", 4.f, "Reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<float> minParamChange{"minParamChange", 1.e-3f, "Stop iterations if largest change of any X is smaller than this"};
  Configurable<float> minRelChi2Change{"minRelChi2Change", 0.9f, "Stop iterations if chi2/chi2old > this"};

  o2::vertexing::DCAFitterN<2> fitter{}; // 2-prong vertex fitter

  using SelectedTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksCov, aod::HfTSelTrack>>;

  Filter filterSelectTracks = aod::hf_seltrack::isSelProng == true;

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    // 2-prong histograms
    registry.add("hVtx2ProngX", "2-prong candidates;#it{x}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -2., 2.}}});
    registry.add("hVtx2ProngY", "2-prong candidates;#it{y}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -2., 2.}}});
    registry.add("hVtx2ProngZ", "2-prong candidates;#it{z}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -20., 20.}}});
    registry.add("hMassD0ToPiK", "D^{0} candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}});

    // Configure the vertexer.
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
    for (const auto& trackPos : tracks) {
      if (trackPos.signed1Pt() < 0) {
        continue;
      }
      const auto trackParVarPos{getTrackParCov(trackPos)};

      // loop over negative tracks
      for (const auto& trackNeg : tracks) {
        if (trackNeg.signed1Pt() > 0) {
          continue;
        }
        const auto trackParVarNeg{getTrackParCov(trackNeg)};

        // secondary-vertex reconstruction and further 2-prong selections
        int nVtxFromFitter = 0;
        try {
          nVtxFromFitter = fitter.process(trackParVarPos, trackParVarNeg);
        } catch (...) {
        }
        if (nVtxFromFitter == 0) {
          continue;
        }
        // get secondary vertex
        const auto& secondaryVertex = fitter.getPCACandidate();
        // get track momenta
        std::array<float, 3> pVec0{};
        std::array<float, 3> pVec1{};
        fitter.getTrack(0).getPxPyPzGlo(pVec0);
        fitter.getTrack(1).getPxPyPzGlo(pVec1);

        // fill table row
        rowTrackIndexProng2(trackPos.globalIndex(),
                            trackNeg.globalIndex());

        // fill histograms
        registry.fill(HIST("hVtx2ProngX"), secondaryVertex[0]);
        registry.fill(HIST("hVtx2ProngY"), secondaryVertex[1]);
        registry.fill(HIST("hVtx2ProngZ"), secondaryVertex[2]);
        const std::array arrayMomenta{pVec0, pVec1};
        const auto massPiK{RecoDecay::m(arrayMomenta, std::array{MassPiPlus, MassKPlus})};
        registry.fill(HIST("hMassD0ToPiK"), massPiK);
        // const auto massKPi{RecoDecay::m(arrayMomenta, std::array{MassKPlus, MassPiPlus})};
        // registry.fill(HIST("hMassD0ToPiK"), massKPi);
      }
    }
  }
};

// Add all tasks in the workflow specification.
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfSkimCreatorMiniTagSelTracks>(cfgc),
    adaptAnalysisTask<HfSkimCreatorMini>(cfgc)};
}
