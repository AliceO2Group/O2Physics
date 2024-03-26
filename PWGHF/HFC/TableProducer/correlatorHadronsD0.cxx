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

/// \author Xu wang <wangxuwx@mails.ccnu.edu.cn>.
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
/// _________________________________________
/// Table for storing trigger track indices
namespace triggerTracks
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                       //!
DECLARE_SOA_INDEX_COLUMN_FULL(Track, track, int, Tracks, "_Trigger"); //!
} // namespace triggerTracks
DECLARE_SOA_TABLE(TriggerTracks, "AOD", "TRIGGERTRACKS", o2::soa::Index<>, triggerTracks::CollisionId, triggerTracks::TrackId);
} // namespace o2::aod
/// _________________________________________
/// Table for storing assoc track indices

// histogram binning definition
const int massAxisNBins = 200;
const double massAxisMin = 1.3848;
const double massAxisMax = 2.3848;
AxisSpec axisDeltaEta = {100, -2., 2., ""};
AxisSpec axisDeltaPhi = {64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf, ""};
AxisSpec axisPtD = {10, 0., 10., ""};
AxisSpec axisPtHadron = {11, 0., 11., ""};
AxisSpec axisPoolBin = {9, 0., 9., ""};
std::vector<AxisSpec> axes = {axisDeltaPhi, axisDeltaEta, axisPtD, axisPtHadron, axisPoolBin};

// definition of ME variables and new types
std::vector<double> zBins{VARIABLE_WIDTH, -10.0, -2.5, 2.5, 10.0};
std::vector<double> multBins{VARIABLE_WIDTH, 0., 200., 500.0, 5000.};
std::vector<double> multBinsMcGen{VARIABLE_WIDTH, 0., 20., 50.0, 500.}; // In MCGen multiplicity is defined by counting primaries
using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFV0M<aod::mult::MultFV0A, aod::mult::MultFV0C>>;
BinningType corrBinning{{zBins, multBins}, true};

double getDeltaPhi(double phiHadron, double phiD)
{
  return RecoDecay::constrainAngle(phiD - phiHadron, -o2::constants::math::PIHalf);
}

struct HfCorrelatorHadronsD0Selection {
  Produces<aod::DmesonSelection> d0Sel;
  Produces<aod::TriggerTracks> triggerTrack;

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<float> yCandMax{"yCandMax", 4.0, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", -1., "min. cand. pT"};
  Configurable<float> triggerEtaMin{"triggerEtaCutMin", -0.8, "triggeretamin"};
  Configurable<float> triggerEtaMax{"triggerEtaCutMax", 0.8, "triggeretamax"};
  Configurable<float> triggerPtCutMin{"triggerPtCutMin", 3, "triggerptmin"};
  Configurable<float> triggerPtCutMax{"triggerPtCutMax", 20, "triggerptmax"};

  HfHelper hfHelper;
  SliceCache cache;

  Preslice<aod::HfCand2Prong> perCol = aod::hf_cand::collisionId;
  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0>> selectedD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;

  void init(InitContext&)
  {
  }

  void processD0SelectionData(aod::Collision const& collision, soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& candidates)
  {
    bool isD0Found = 0;
    if (selectedD0Candidates.size() > 0) {
      auto selectedD0CandidatesGrouped = selectedD0Candidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);

      for (const auto& candidate1 : selectedD0CandidatesGrouped) {
        // check decay channel flag for candidate1
        if (!TESTBIT(candidate1.hfflag(), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          continue;
        }
        if (yCandMax >= 0. && std::abs(hfHelper.yD0(candidate1)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate1.pt() < ptCandMin) {
          continue;
        }
        isD0Found = 1;
      }
    }
    d0Sel(isD0Found);
  }
  PROCESS_SWITCH(HfCorrelatorHadronsD0Selection, processD0SelectionData, "Process D0 Selection Data", false);

  void processTriggerHardons(aod::Collision const& collision, aod::TracksWDca const& tracks)
  {
    for (auto const& track : tracks) {
      if (track.eta() > triggerEtaMax || track.eta() < triggerEtaMin) {
        continue;
      }
      if (track.pt() > triggerPtCutMax || track.pt() < triggerPtCutMin) {
        continue;
      }
      triggerTrack(track.collisionId(), track.globalIndex());
    }
  }
  PROCESS_SWITCH(HfCorrelatorHadronsD0Selection, processTriggerHardons, "Process Hadrons Trigger Selection Data", false);
};

struct HfCorrelatorHadronsD0 {

  SliceCache cache;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};

  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0>> selectedD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Filter collisionFilter = aod::hf_selection_dmeson_collision::dmesonSel == true;
  Filter d0Filter = (aod::hf_sel_candidate_d0::isSelD0 >= 1) || (aod::hf_sel_candidate_d0::isSelD0bar >= 1);

  void init(InitContext&)
  {
    histos.add("CorrelatorHadronsD0", "", kTHnF, axes);
    histos.add("hMass1D", "D0,D0bar candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisNBins, massAxisMin, massAxisMax}}});
    histos.add("hMassD01D", "D0,D0bar candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisNBins, massAxisMin, massAxisMax}}});
    histos.add("hMassD0bar1D", "D0,D0bar candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisNBins, massAxisMin, massAxisMax}}});
  };

  void processCorrelatorHadronD0(soa::Join<aod::Collisions, aod::Mults>::iterator const& collision,
                                 aod::TracksWDca const& tracks,
                                 soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& candidates,
                                 aod::TriggerTracks const& triggers)
  {
    int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multFV0M()));
    for (auto& triggerTrack : triggers) {
      auto trigg = triggerTrack.track_as<aod::TracksWCovDca>();
      auto selectedD0CandidatesGrouped = selectedD0Candidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      for (const auto& candidate1 : selectedD0CandidatesGrouped) {
        // Remove D0 daughters by checking track indices
        if ((candidate1.prong0Id() == trigg.globalIndex()) || (candidate1.prong1Id() == trigg.globalIndex())) {
          continue;
        }
        float deltaPhi = getDeltaPhi(trigg.phi(), candidate1.phi());
        float deltaEta = trigg.eta() - candidate1.eta();
        float ptCandidate = candidate1.pt();
        float ptTrigger = trigg.pt();
        histos.fill(HIST("CorrelatorHadronsD0"), deltaPhi, deltaEta, ptCandidate, ptTrigger, poolBin);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorHadronsD0, processCorrelatorHadronD0, "HfCorrelatorHadronsD0", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCorrelatorHadronsD0Selection>(cfgc),
    adaptAnalysisTask<HfCorrelatorHadronsD0>(cfgc)};
}
