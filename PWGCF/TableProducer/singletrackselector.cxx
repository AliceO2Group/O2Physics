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
/// \brief create a table applying some basic cuts on the ITS and DCA.
/// \author Sofia Tomassini
/// \since 31 May 2023

#include <Framework/AnalysisDataModel.h>
#include <fairlogger/Logger.h>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"

#include "PWGCF/DataModel/singletrackselector.h"

#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;
using namespace o2::aod;

struct singleTrackSelector {

  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply rapidity cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<float> cutDcaXy{"cutDcaXy", 0.12, ""};
  Configurable<float> cutPtMin{"cutPtMin", 0.4, "Minimum cut in pT"};
  Configurable<float> cutTPCNSigmaPr{"cutTPCNSigmaPr", 5.f, "Cut on the TPC nsigma for protons"};
  Configurable<float> cutTPCNSigmaDe{"cutTPCNSigmaDe", 5.f, "Cut on the TPC nsigma for deuteron"};

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidEvTimeFlags, aod::TracksDCA,
                         aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                         aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr,
                         aod::pidTPCFullDe, aod::pidTOFFullDe, aod::pidTOFbeta,
                         aod::TrackSelection>;
  using Coll = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::FT0sCorrected>;

  Produces<o2::aod::SingleTrackSel> tableRow;
  Produces<o2::aod::SingleCollSel> tableRowColl;

  Filter eventFilter = (applyEvSel.node() == 0) ||
                       ((applyEvSel.node() == 1) && (aod::evsel::sel7 == true)) ||
                       ((applyEvSel.node() == 2) && (aod::evsel::sel8 == true));
  Filter vertexFilter = ((o2::aod::collision::posZ < 15.f) && (o2::aod::collision::posZ > -15.f));
  Filter trackFilter = ((o2::aod::track::itsChi2NCl <= 36.f) && (o2::aod::track::itsChi2NCl >= 0.f) && (o2::aod::track::tpcChi2NCl >= 0.f) && (o2::aod::track::tpcChi2NCl <= 4.f));

  void process(soa::Filtered<Coll>::iterator const& collision, soa::Filtered<Trks> const& tracks)
  {
    tableRow.reserve(tracks.size());
    tableRowColl(collision.globalIndex(),
                 collision.multTPC(),
                 collision.posZ());

    for (auto& track : tracks) {
      if (track.pt() < cutPtMin) {
        continue;
      }
      if (abs(track.dcaXY()) > cutDcaXy) {
        continue;
      }
      if (abs(track.tpcNSigmaPr()) < cutTPCNSigmaPr || abs(track.tpcNSigmaDe()) < cutTPCNSigmaDe) {

        tableRow(tableRowColl.lastIndex(),
                 track.hasITS(),
                 track.hasTOF(),
                 track.px(),
                 track.py(),
                 track.pz(),
                 track.tpcInnerParam(),
                 track.tpcSignal(),
                 track.beta(),
                 track.dcaXY(),
                 track.dcaZ(),
                 track.tpcNClsFound(),
                 track.tpcFoundOverFindableCls(),
                 track.tpcChi2NCl(),
                 track.itsNCls(),
                 track.itsChi2NCl(),
                 track.sign(),
                 track.eta(),
                 track.phi(),
                 singletrackselector::packInTableOffset<singletrackselector::storedcrossedrows::binning>(track.tpcNClsCrossedRows()),
                 singletrackselector::packInTable<singletrackselector::nsigma::binning>(track.tofNSigmaPr()),
                 singletrackselector::packInTable<singletrackselector::nsigma::binning>(track.tpcNSigmaPr()),
                 singletrackselector::packInTable<singletrackselector::nsigma::binning>(track.tofNSigmaDe()),
                 singletrackselector::packInTable<singletrackselector::nsigma::binning>(track.tpcNSigmaDe()));
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<singleTrackSelector>(cfgc)};
}
