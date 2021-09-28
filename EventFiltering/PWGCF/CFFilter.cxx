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

/// \file CFFilter.cxx
/// \brief Selection of events with triplets for femtoscopic studies
///
/// \author Laura Serksnyte, TU München, laura.serksnyte@cern.ch

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

#include "../filterTables.h"

#include "../../PWGCF/DataModel/FemtoDerived.h"
#include "../../PWGCF/FemtoDream/FemtoDreamParticleHisto.h"
#include "../../PWGCF/FemtoDream/FemtoDreamPairCleaner.h"
#include "../../PWGCF/FemtoDream/FemtoDreamContainer.h"
#include "../../PWGCF/FemtoDream/FemtoDreamMath.h"
#include "../../PWGCF/FemtoDream/FemtoDreamPairCleaner.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include <cmath>
#include <string>

namespace
{

static constexpr int nTriplets{4};

enum CFTriggers {
  kPPP = 0,
  kPPL,
  kPLL,
  kLLL
};

static const std::vector<std::string> CfTriggerNames{"ppp", "ppL", "pLL", "LLL"};

} // namespace

namespace o2::aod
{
using FullCollision = soa::Join<aod::Collisions,
                                aod::EvSels,
                                aod::Mults>::iterator;
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoDream;

struct CFFilter {

  Produces<aod::CFFilters> tags;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  //FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCleanerTT; Currently not used, will be needed later
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> pairCleanerTV;

  void init(o2::framework::InitContext&)
  {
    registry.add("fProcessedEvents", "CF - event filtered;;events", HistType::kTH1F, {{6, -0.5, 4.5}});
    std::array<std::string, 6> eventTitles = {"all", "rejected", "p-p-p", "p-p-L", "p-L-L", "L-L-L"};
    for (size_t iBin = 0; iBin < eventTitles.size(); iBin++) {
      registry.get<TH1>(HIST("fProcessedEvents"))->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
    }
    registry.add("fSameEvent", "CF - same event distribution;;events", HistType::kTH1F, {{8000, 0, 8}});
  }

  float mMassProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  float mMassLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

  void process(aod::FullCollision const& col, o2::aod::FemtoDreamCollisions const& colsFemto,
               o2::aod::FemtoDreamParticles& partsFemto)
  {
    registry.get<TH1>(HIST("fProcessedEvents"))->Fill(0);
    bool keepEvent[nTriplets]{false};

    // check if current collision was selected by femto skimmer
    auto collsInFemto = colsFemto.sliceBy(o2::aod::femtodreamcollision::collisionId, col.globalIndex());
    if (collsInFemto.size() == 0) {
      keepEvent[kPPL] = false;
    } else {
      auto collInFemto = collsInFemto.begin();

      //uint8_t trackTypeSel = o2::aod::femtodreamparticle::ParticleType::kTrack; Fix this to work instead of below hardcoded lines
      //uint V0TypeSel = o2::aod::femtodreamparticle::ParticleType::kV0; Fix this to work instead of below hardcoded lines
      uint8_t Track = 0;      // Track
      uint8_t V0 = 1;         // V0
      uint8_t V0Daughter = 2; // V0
      uint32_t trackMask = 524288;
      uint32_t V0Mask = 64;
      uint32_t testedBit0 = 0;

      //Obtain particle and antiparticle candidates of protons and lambda hyperons for current femto collision
      Partition<o2::aod::FemtoDreamParticles> partsProton1 = (o2::aod::femtodreamparticle::femtoDreamCollisionId == collInFemto.globalIndex()) &&
                                                             (o2::aod::femtodreamparticle::partType == Track) && ((o2::aod::femtodreamparticle::cut & trackMask) == trackMask);
      Partition<o2::aod::FemtoDreamParticles> partsLambda1 = (o2::aod::femtodreamparticle::femtoDreamCollisionId == collInFemto.globalIndex()) &&
                                                             (o2::aod::femtodreamparticle::partType == V0) && ((o2::aod::femtodreamparticle::cut & V0Mask) == V0Mask);
      Partition<o2::aod::FemtoDreamParticles> partsProton0 = (o2::aod::femtodreamparticle::femtoDreamCollisionId == collInFemto.globalIndex()) &&
                                                             (o2::aod::femtodreamparticle::partType == Track) && ((o2::aod::femtodreamparticle::cut & trackMask) == testedBit0);
      Partition<o2::aod::FemtoDreamParticles> partsLambda0 = (o2::aod::femtodreamparticle::femtoDreamCollisionId == collInFemto.globalIndex()) &&
                                                             (o2::aod::femtodreamparticle::partType == V0) && ((o2::aod::femtodreamparticle::cut & V0Mask) == testedBit0);

      partsProton1.bindTable(partsFemto);
      partsLambda1.bindTable(partsFemto);
      partsProton0.bindTable(partsFemto);
      partsLambda0.bindTable(partsFemto);

      // This is the main trigger part for proton-proton-Lambda
      // pairCleanerTV -> Test if lambda hyperons don't have a daughter which is as well used as primary proton
      // Calculate Q3 and check if it is smaller than 0.6
      // If at collision has at least one triplet with Q3<0.6, the trigger value is set to true!
      // IMPORTANT: Include close pair rejection here
      int lowQ3Triplets = 0;
      if (partsLambda0.size() >= 1 && partsProton0.size() >= 2) {
        for (auto& partLambda : partsLambda0) {
          if (!pairCleanerTV.isCleanPair(partLambda, partLambda, partsFemto)) {
            continue;
          }
          for (auto& [p1, p2] : combinations(partsProton0, partsProton0)) {
            auto Q3 = FemtoDreamMath::getQ3(p1, mMassProton, p2, mMassProton, partLambda, mMassLambda);
            registry.get<TH1>(HIST("fSameEvent"))->Fill(Q3);
            if (Q3 < 1.5)
              lowQ3Triplets++; // real value 0.6. We use 1.5 here for testing locally because of statistics
          }
        }
      }
      if (partsLambda1.size() >= 1 && partsProton1.size() >= 2) {
        for (auto& partLambda : partsLambda1) {
          if (!pairCleanerTV.isCleanPair(partLambda, partLambda, partsFemto)) {
            continue;
          }
          for (auto& [p1, p2] : combinations(partsProton1, partsProton1)) {
            auto Q3 = FemtoDreamMath::getQ3(p1, mMassProton, p2, mMassProton, partLambda, mMassLambda);
            registry.get<TH1>(HIST("fSameEvent"))->Fill(Q3);
            if (Q3 < 1.5)
              lowQ3Triplets++; // real value 0.6. We use 1.5 here for testing locally because of statistics
          }
        }
      }
      if (lowQ3Triplets > 0)
        keepEvent[kPPL] = true;
    }

    tags(keepEvent[kPPP], keepEvent[kPPL], keepEvent[kPLL], keepEvent[kLLL]);

    if (!keepEvent[kPPP] && !keepEvent[kPPL] && !keepEvent[kPLL] && !keepEvent[kLLL]) {
      registry.get<TH1>(HIST("fProcessedEvents"))->Fill(1);
    } else {
      for (int iTrigger{0}; iTrigger < nTriplets; iTrigger++) {
        if (keepEvent[iTrigger]) {
          registry.get<TH1>(HIST("fProcessedEvents"))->Fill(iTrigger + 1);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<CFFilter>(cfg)};
}
