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
/// \file   qaCascades.cxx
/// \author Maria Barlou m.barlou@cern.ch
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  QA task for basic quantities on cascades
///

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

using PIDTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullPi>;
using PIDTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullPi>;
using CollisionCandidates = soa::Join<aod::Collisions, aod::EvSels>;

struct qaCascades {
  // Configurables
  ConfigurableAxis invMassXi{"invMassXi", {1000, 0.f, 5.f}, "Binning for the Xi mass histograms"};
  ConfigurableAxis invMassOmega{"invMassOmega", {1000, 0.f, 5.f}, "Binning for the Omega mass histograms"};
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply event selection cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  HistogramRegistry registry{"qaCascades"};

#define fillHistogram(name, ...) registry.fill(HIST(name), __VA_ARGS__)

  void init(InitContext const&)
  {
    // Axis definition
    const AxisSpec axisMassXi{invMassXi, "invMassXi"};
    const AxisSpec axisMassOmega{invMassOmega, "invMassOmega"};

    // Histogram definition
    registry.add("massXi", "massXi", kTH1D, {axisMassXi});
    registry.add("massOmega", "massOmega", kTH1D, {axisMassOmega});

    LOG(info) << "qaCascades initialized";
    registry.print();
  }

  Filter eventFilter = (applyEvSel.node() == 0) ||
                       ((applyEvSel.node() == 1) && (o2::aod::evsel::sel7 == true)) ||
                       ((applyEvSel.node() == 2) && (o2::aod::evsel::sel8 == true));

  void process(soa::Filtered<CollisionCandidates>::iterator const&, aod::CascDataExt const& cascades)
  {
    // Monitor average number of cascades per event (for trending vs run number)
    for (auto& casc : cascades) {
      registry.fill(HIST("massXi"), casc.mXi());
      registry.fill(HIST("massOmega"), casc.mOmega());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<qaCascades>(cfgc)};
}
