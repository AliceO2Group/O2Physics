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

/// \file correlatorDsHadronsReduced.cxx
/// \brief Ds-Hadrons correlator task for ME offline
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>

#include <vector>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/HFC/DataModel/DerivedDataCorrelationTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
double getDeltaPhi(double phiHadron, double phiD)
{
  return RecoDecay::constrainAngle(phiHadron - phiD, -PIHalf);
}

// binning type
using BinningTypeDerived = ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Multiplicity>;

/// Ds-Hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfCorrelatorDsHadronsReduced {
  Produces<aod::DsHadronPair> entryDsHadronPair;
  Produces<aod::DsHadronRecoInfo> entryDsHadronRecoInfo;
  // Produces<aod::DsHadronGenInfo> entryDsHadronGenInfo;

  Configurable<bool> fillHistoData{"fillHistoData", true, "Flag for filling histograms in data processes"};
  Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "Number of events mixed in ME process"};

  SliceCache cache;

  // Preslice<aod::AssocTrackReds> tracksPerCol = aod::hf_assoc_track_reduced::hfcRedCollisionId;
  Preslice<aod::AssocTrackReds> tracksPerCol = aod::hf_candidate_reduced::hfcRedCollisionId;
  Preslice<aod::DsCandReduceds> candPerCol = aod::hf_candidate_reduced::hfcRedCollisionId;

  ConfigurableAxis zPoolBins{"zPoolBins", {VARIABLE_WIDTH, -10.0, -2.5, 2.5, 10.0}, "z vertex position pools"};
  ConfigurableAxis multPoolBins{"multPoolBins", {VARIABLE_WIDTH, 0., 900., 1800., 6000.}, "event multiplicity pools (FT0M)"};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {600, 0., 6000.}, "Multiplicity as FT0M signal amplitude"};
  ConfigurableAxis binsPosZ{"binsPosZ", {100, -10., 10.}, "primary vertex z coordinate"};
  ConfigurableAxis binsPoolBin{"binsPoolBin", {9, 0., 9.}, "PoolBin"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    AxisSpec axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    AxisSpec axisPosZ = {binsPosZ, "PosZ"};
    AxisSpec axisPoolBin = {binsPoolBin, "PoolBin"};

    // Histograms for data analysis
    if (fillHistoData) {
      registry.add("hMultFT0M", "Multiplicity FT0M", {HistType::kTH1F, {axisMultFT0M}});
      registry.add("hZVtx", "z vertex", {HistType::kTH1F, {axisPosZ}});
      registry.add("hDsPoolBin", "Ds candidates pool bin", {HistType::kTH1F, {axisPoolBin}});
      registry.add("hTracksPoolBin", "Particles associated pool bin", {HistType::kTH1F, {axisPoolBin}});
    }
  }

  void processDerivedDataME(aod::HfcRedCollisions const& collisions,
                            aod::DsCandReduceds const& candidates,
                            aod::AssocTrackReds const& tracks)
  {

    BinningTypeDerived corrBinning{{zPoolBins, multPoolBins}, true};

    auto tracksTuple = std::make_tuple(candidates, tracks);

    Pair<aod::HfcRedCollisions, aod::DsCandReduceds, aod::AssocTrackReds, BinningTypeDerived> pairData{corrBinning, numberEventsMixed, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairData) {
      if (tracks1.size() == 0) {
        continue;
      }

      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multiplicity()));
      int poolBinDs = corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multiplicity()));
      registry.fill(HIST("hMultFT0M"), c1.multiplicity());
      registry.fill(HIST("hZVtx"), c1.posZ());
      registry.fill(HIST("hTracksPoolBin"), poolBin);
      registry.fill(HIST("hDsPoolBin"), poolBinDs);

      for (const auto& [cand, pAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d)", cand.index(), pAssoc.index(), c1.index(), c2.index(), cand.hfcRedCollisionId(), pAssoc.hfcRedCollisionId());

        entryDsHadronPair(getDeltaPhi(pAssoc.phiAssocTrack(), cand.phiCand()),
                          pAssoc.etaAssocTrack() - cand.etaCand(),
                          cand.ptCand(),
                          pAssoc.ptAssocTrack(),
                          poolBin);
        entryDsHadronRecoInfo(cand.invMassDs(), false, false);
        // entryDsHadronGenInfo(false, false, 0);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadronsReduced, processDerivedDataME, "Process Mixed Event Derived Data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorDsHadronsReduced>(cfgc)};
}
