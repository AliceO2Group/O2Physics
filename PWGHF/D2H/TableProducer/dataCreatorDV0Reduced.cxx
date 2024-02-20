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

/// \file dataCreatorDV0Reduced.cxx
/// \brief Creation of D-V0 pairs
///
/// \author Luca Aglietta <luca.aglietta@cern.ch>, UniTO Turin

#include <cmath>
#include <map>

#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

// event types
enum Event : uint8_t {
  Processed = 0,
  NoDV0Selected,
  DV0Selected,
  kNEvent
};

/// Creation of D-V0 pairs
struct HfDataCreatorDV0Reduced {
  // Produces AOD tables to store track information
  Produces<aod::StraCollisions> lfReducedCollision; //Defined in PWGLF/DataModel/LFStrangenessTables.h
  Produces<aod::HfOrigColCounts> hfCollisionCounter; // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  Produces<aod::HfRedVzeros> hfCandV0; // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  Produces<aod::HfRedD> hfCandD; // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h

  // selection  
  Configurable<int> selectionFlagD{"selectionFlagD", 7, "Selection Flag for D"};

  HfHelper hfHelper;

  // O2DatabasePDG service
  Service<o2::framework::O2DatabasePDG> pdg;

  bool isHfCandResoConfigFilled = false;

  using TracksPidAll = soa::Join<aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using TracksPIDWithSel = soa::Join<aod::TracksWCovDcaExtra, TracksPidAll, aod::TrackSelection>;
  using TracksPIDWithSelAndMc = soa::Join<TracksPIDWithSel, aod::McTrackLabels>;
  using CandsDFiltered = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagD);

  Preslice<CandsDFiltered> candsDPerCollision = aod::track_association::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<aod::V0Datas> candsV0PerCollision = aod::track_association::collisionId;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    // histograms
    constexpr int kNBinsEvents = kNEvent;
    std::string labels[kNBinsEvents];
    labels[Event::Processed] = "processed";
    labels[Event::NoDV0Selected] = "without DV0 pairs";
    labels[Event::DV0Selected] = "with DV0 pairs";
    static const AxisSpec axisEvents = {kNBinsEvents, 0.5, kNBinsEvents + 0.5, ""};
    registry.add("hEvents", "Events;;entries", HistType::kTH1F, {axisEvents});
    for (int iBin = 0; iBin < kNBinsEvents; iBin++) {
      registry.get<TH1>(HIST("hEvents"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }
    registry.add("hMassDToPiKPi", "D^{#minus} candidates;inv. mass (p^{#minus} K^{#plus} #pi^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, 1.7, 2}}});
    registry.add("hMassV0", "K0^{s} candidates;inv. mass (#pi^{#plus}#pi^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, 0.35, 0.65}}});
    registry.add("hPtD", "D^{#minus} candidates;D^{#minus} candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("hPtV0", "V0 candidates;V0 candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("hMassDs1", "Ds1 candidates;Ds1 candidate #inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, 2., 3.}}});
  }

  void runDataCreation(aod::Collision const& collision,
                       CandsDFiltered const& candsD,
                       aod::V0Datas const& V0s)
  {
    // helpers for ReducedTables filling
    int indexHfReducedCollision = lfReducedCollision.lastIndex() + 1;
    // std::map where the key is the V0.globalIndex() and
    // the value is the V0 index in the table of the selected v0s
    std::map<int64_t, int64_t> selectedV0s;
    bool fillHfReducedCollision = false;

    for (const auto& candD : candsD) {
      bool fillHfCandD = false;
      float invMassD = hfHelper.invMassDplusToPiKPi(candD);
      std::array<float, 3> pVecD = candD.pVector();

      for (const auto& v0 : V0s) {
        std::array<float, 3> pVecV0 = {v0.pxpos() + v0.pxneg(), v0.pypos() + v0.pyneg(), v0.pzpos() + v0.pzneg()};
        float ptV0 = sqrt(pVecV0[0]*pVecV0[0] + pVecV0[1]*pVecV0[1]);
        auto invMass2DV0 = RecoDecay::m2(std::array{pVecD, pVecV0}, std::array{MassDPlus, MassK0Short});
        registry.fill(HIST("hMassDs1"), sqrt(invMass2DV0));
     
        // fill V0 table
        // if information on V0 already stored, go to next V0
        if (!selectedV0s.count(v0.globalIndex())) {
          hfCandV0(v0.posTrackId(), v0.negTrackId(), 
                   indexHfReducedCollision,
                   v0.x(), v0.y(), v0.z(),
                   v0.mK0Short(),
                   pVecV0[0], pVecV0[1], pVecV0[2]);
          selectedV0s[v0.globalIndex()] = hfCandV0.lastIndex();
          registry.fill(HIST("hPtV0"), ptV0);
          registry.fill(HIST("hMassV0"), v0.mK0Short());
        }
        fillHfCandD = true;
      } // V0 loop
      
      if (fillHfCandD) { // fill candDplus table only once per D candidate, only if at least one V0 is found
        hfCandD(candD.prong0Id(), candD.prong1Id(), candD.prong2Id(),
                indexHfReducedCollision,
                candD.xSecondaryVertex(), candD.ySecondaryVertex(), candD.zSecondaryVertex(),
                invMassD,
                pVecD[0],pVecD[1],pVecD[2]);
        fillHfReducedCollision = true;
        registry.fill(HIST("hMassDToPiKPi"), invMassD);
        registry.fill(HIST("hPtD"), candD.pt());
      }
    } // candsD loop
    registry.fill(HIST("hEvents"), 1 + Event::Processed);
    if (!fillHfReducedCollision) {
      registry.fill(HIST("hEvents"), 1 + Event::NoDV0Selected);
      return;
    }
    registry.fill(HIST("hEvents"), 1 + Event::DV0Selected);
    // fill collision table if it contains a DPi pair a minima
    lfReducedCollision(collision.posX(), collision.posY(), collision.posZ());
  }

  void processData(aod::Collisions const& collisions,
                   CandsDFiltered const& candsD,
                   aod::TrackAssoc const& trackIndices,
                   aod::V0Datas const& V0s)
  {
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize());

    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);
      auto V0sThisColl = V0s.sliceBy(candsV0PerCollision, thisCollId);
      runDataCreation(collision, candsDThisColl, V0sThisColl);
    }
  }
  PROCESS_SWITCH(HfDataCreatorDV0Reduced, processData, "Process without MC info and without ML info", true);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDataCreatorDV0Reduced>(cfgc)};
}
