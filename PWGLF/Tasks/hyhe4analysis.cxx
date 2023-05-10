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
//
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//  3-body Hyperhelium 4 analysis
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//
// with contributions from:
// --- navneet.kumar@cern.ch
// --- natasha.sharma@cern.ch
// --- kumar.lokesh@cern.ch
// --- david.dobrigkeit.chinellato@cern.ch

#include <cmath>
#include <array>
#include <cstdlib>
#include <map>
#include <iterator>
#include <utility>

#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DCAFitter/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "PWGLF/DataModel/LFHyperhelium4Tables.h"

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// use parameters + cov mat non-propagated, aux info + (extension propagated)
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
using TracksWithExtra = soa::Join<aod::Tracks, aod::TracksExtra>;

// For MC and dE/dx association
using TracksExtraWithPIDandLabels = soa::Join<aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullHe, aod::McTrackLabels>;

// For MC association in pre-selection
using LabeledTracksExtra = soa::Join<aod::TracksExtra, aod::McTrackLabels>;

struct hyhefouranalysis {
  // Basic selection criteria
  Configurable<float> selHyHe4daudca{"selHyHe4daudca", 1, "DCA between HyHe4 daughters"};

  // storing output
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext& context)
  {
    const AxisSpec axisMassHyHe4{(int)400, 3.7f, 4.1f, "Hyperhelium4 Mass Distribution (GeV/c^{2})"};
    const AxisSpec axisNCandidates{(int)100, -0.5f, 99.5f, "Number of 3-body candidates"};

    const AxisSpec axisEventCounter{(int)1, 0.0f, 1.0f, "Number of events"};
    const AxisSpec axisPt{(int)1000, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisDCADaughters{(int)1000, 0.0f, 10.0f, "DCA Daughters (cm)"};
    const AxisSpec axisDCAtoPV{(int)1000, -10.0f, 10.0f, "DCA to PV (cm)"};

    // Base counters
    histos.add("hNEvents", "hNEvents", kTH1F, {axisEventCounter});
    histos.add("hNCandidates", "hNCandidates", kTH1F, {axisNCandidates});

    // Kinematics
    histos.add("hPtHyHe4_Helium3", "hPtHyHe4_Helium3", kTH1F, {axisPt});
    histos.add("hPtHyHe4_Proton", "hPtHyHe4_Proton", kTH1F, {axisPt});
    histos.add("hPtHyHe4_Pion", "hPtHyHe4_Pion", kTH1F, {axisPt});
    histos.add("hPtAntiHyHe4_Helium3", "hPtAntiHyHe4_Helium3", kTH1F, {axisPt});
    histos.add("hPtAntiHyHe4_Proton", "hPtAntiHyHe4_Proton", kTH1F, {axisPt});
    histos.add("hPtAntiHyHe4_Pion", "hPtAntiHyHe4_Pion", kTH1F, {axisPt});

    // Topological variables
    histos.add("hPtHyHe4_DCADaughters", "hPtHyHe4_DCADaughters", kTH1F, {axisDCADaughters});
    histos.add("hPtHyHe4_DCAxy", "hPtHyHe4_DCAxy", kTH1F, {axisDCAtoPV});
    histos.add("hPtHyHe4_DCAz", "hPtHyHe4_DCAz", kTH1F, {axisDCAtoPV});
    histos.add("hPtHyHe4_Helium3_DCA", "hPtHyHe4_Helium3_DCA", kTH1F, {axisDCAtoPV});
    histos.add("hPtHyHe4_Proton_DCA", "hPtHyHe4_Proton_DCA", kTH1F, {axisDCAtoPV});
    histos.add("hPtHyHe4_Pion_DCA", "hPtHyHe4_Pion_DCA", kTH1F, {axisDCAtoPV});
    histos.add("hPtAntiHyHe4_DCADaughters", "hPtAntiHyHe4_DCADaughters", kTH1F, {axisDCADaughters});
    histos.add("hPtAntiHyHe4_DCAxy", "hPtAntiHyHe4_DCAxy", kTH1F, {axisDCAtoPV});
    histos.add("hPtAntiHyHe4_DCAz", "hPtAntiHyHe4_DCAz", kTH1F, {axisDCAtoPV});
    histos.add("hPtAntiHyHe4_Helium3_DCA", "hPtAntiHyHe4_Helium3_DCA", kTH1F, {axisDCAtoPV});
    histos.add("hPtAntiHyHe4_Proton_DCA", "hPtAntiHyHe4_Proton_DCA", kTH1F, {axisDCAtoPV});
    histos.add("hPtAntiHyHe4_Pion_DCA", "hPtAntiHyHe4_Pion_DCA", kTH1F, {axisDCAtoPV});

    // Base mass
    histos.add("hMassVsPtHyHe4", "hMassVsPtHyHe4", kTH2F, {axisPt, axisMassHyHe4});
    histos.add("hMassVsPtAntiHyHe4", "hMassVsPtAntiHyHe4", kTH2F, {axisPt, axisMassHyHe4});
  }

  void process(aod::Collision const& collision, aod::HyHe4Datas const& hyhe4candidates, FullTracksExtIU const&, aod::BCsWithTimestamps const&)
  {
    histos.fill(HIST("hNEvents"), 0.5);
    histos.fill(HIST("hNCandidates"), hyhe4candidates.size());
    for (auto const& hyhe4cand : hyhe4candidates) {
      // process this particular candidate with existing data model
      if (fabs(hyhe4cand.yHyHe4()) > 0.5) {
        continue;
      }
      if (hyhe4cand.sign() > 0) {
        histos.fill(HIST("hMassVsPtHyHe4"), hyhe4cand.pt(), hyhe4cand.m());
        histos.fill(HIST("hPtHyHe4_Helium3"), hyhe4cand.ptHelium3());
        histos.fill(HIST("hPtHyHe4_Proton"), hyhe4cand.ptProton());
        histos.fill(HIST("hPtHyHe4_Pion"), hyhe4cand.ptPion());

        histos.fill(HIST("hPtHyHe4_DCADaughters"), hyhe4cand.dcaDaughters());
        histos.fill(HIST("hPtHyHe4_DCAxy"), hyhe4cand.dcaxyHyHe4ToPV());
        histos.fill(HIST("hPtHyHe4_DCAz"), hyhe4cand.dcazHyHe4ToPV());
        histos.fill(HIST("hPtHyHe4_Helium3_DCA"), hyhe4cand.dcaHelium3ToPV());
        histos.fill(HIST("hPtHyHe4_Proton_DCA"), hyhe4cand.dcaProtonToPV());
        histos.fill(HIST("hPtHyHe4_Pion_DCA"), hyhe4cand.dcaPionToPV());
      }

      if (hyhe4cand.sign() < 0) {
        histos.fill(HIST("hMassVsPtAntiHyHe4"), hyhe4cand.pt(), hyhe4cand.m());
        histos.fill(HIST("hPtAntiHyHe4_Helium3"), hyhe4cand.ptHelium3());
        histos.fill(HIST("hPtAntiHyHe4_Proton"), hyhe4cand.ptProton());
        histos.fill(HIST("hPtAntiHyHe4_Pion"), hyhe4cand.ptPion());

        histos.fill(HIST("hPtAntiHyHe4_DCADaughters"), hyhe4cand.dcaDaughters());
        histos.fill(HIST("hPtAntiHyHe4_DCAxy"), hyhe4cand.dcaxyHyHe4ToPV());
        histos.fill(HIST("hPtAntiHyHe4_DCAz"), hyhe4cand.dcazHyHe4ToPV());
        histos.fill(HIST("hPtAntiHyHe4_Helium3_DCA"), hyhe4cand.dcaHelium3ToPV());
        histos.fill(HIST("hPtAntiHyHe4_Proton_DCA"), hyhe4cand.dcaProtonToPV());
        histos.fill(HIST("hPtAntiHyHe4_Pion_DCA"), hyhe4cand.dcaPionToPV());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hyhefouranalysis>(cfgc)};
}
