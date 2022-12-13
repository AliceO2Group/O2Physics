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
// Task producing QA histograms to study track (pre-)propagation
//

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/trackUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct propagatorQa {
  Configurable<float> windowDCA{"windowDCA", 50, "windowDCA"};
  Configurable<int> Nbins{"Nbins", 10000, "Nbins"};

  // Momentum distribution
  OutputObj<TH1F> hPtBeforeXcut{TH1F("hPtBeforeXcut", "hPtBeforeXcut", Nbins, 0, 10)};
  OutputObj<TH1F> hPt{TH1F("hPt", "hPt", Nbins, 0, 10)};
  OutputObj<TH1F> hPtusedInSVertexer{TH1F("hPtusedInSVertexer", "hPtusedInSVertexer", Nbins, 0, 10)};

  // IU radii, also from svertexer
  OutputObj<TH1F> hUpdateRadii{TH1F("hUpdateRadii", "hUpdateRadii", 5000, 0, 100)};
  OutputObj<TH1F> hUpdateRadiiusedInSVertexer{TH1F("hUpdateRadiiusedInSVertexer", "hUpdateRadii", 5000, 0, 100)};

  // DCA
  OutputObj<TH1F> hdcaXYall{TH1F("hdcaXYall", "hdcaXYall", Nbins, -windowDCA, windowDCA)};
  OutputObj<TH1F> hdcaXYusedInSVertexer{TH1F("hdcaXYusedInSVertexer", "hdcaXYusedInSVertexer", Nbins, -windowDCA, windowDCA)};

  o2::track::TrackPar lTrackParametrization;

  void process(aod::Collision const& collision, aod::V0s const& V0s, aod::Cascades const& cascades, soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA> const& tracks)
  {
    for (auto& track : tracks) {
      hPtBeforeXcut->Fill(track.pt());

      float maxXtoConsider = o2::constants::geom::XTPCInnerRef + 0.1;
      if (track.trackType() != aod::track::TrackIU && track.x() > maxXtoConsider)
        continue;

      std::array<float, 3> pos;
      lTrackParametrization = getTrackPar(track);
      lTrackParametrization.getXYZGlo(pos);
      float lRadiusOfLastUpdate = TMath::Sqrt(pos[0] * pos[0] + pos[1] * pos[1]);

      hUpdateRadii->Fill(lRadiusOfLastUpdate);

      // fill kinematic variables
      hPt->Fill(track.pt());

      float lDCA = track.dcaXY();

      hdcaXYall->Fill(lDCA);

      // determine if track was used in svertexer
      bool usedInSVertexer = false;
      bool lUsedByV0 = false, lUsedByCascade = false;
      for (auto& V0 : V0s) {
        if (V0.posTrackId() == track.globalIndex()) {
          lUsedByV0 = true;
          break;
        }
        if (V0.negTrackId() == track.globalIndex()) {
          lUsedByV0 = true;
          break;
        }
      }
      for (auto& cascade : cascades) {
        if (cascade.bachelorId() == track.globalIndex()) {
          lUsedByCascade = true;
          break;
        }
      }
      if (lUsedByV0 || lUsedByCascade)
        usedInSVertexer = true;

      if (usedInSVertexer)
        hUpdateRadiiusedInSVertexer->Fill(lRadiusOfLastUpdate);
      if (usedInSVertexer)
        hdcaXYusedInSVertexer->Fill(lDCA);
      if (usedInSVertexer)
        hPtusedInSVertexer->Fill(track.pt());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<propagatorQa>(cfgc)};
}
