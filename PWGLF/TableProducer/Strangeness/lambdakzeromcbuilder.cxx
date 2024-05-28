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
//  Lambdakzero label builder task
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    david.dobrigkeit.chinellato@cern.ch
//

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
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

//*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
struct lambdakzeromcbuilder {
  Produces<aod::McV0Labels> v0labels; // MC labels for V0s
  Produces<aod::V0MCCores> v0mccores; // optionally aggregate information from MC side for posterior analysis (derived data)
  Produces<aod::V0CoreMCLabels> v0CoreMCLabels; // interlink V0Cores -> V0MCCores in asymmetric mode

  Configurable<bool> populateV0MCCoresSymmetric{"populateV0MCCoresSymmetric", false, "populate V0MCCores table for derived data analysis, keep V0MCCores joinable with V0Cores"};
  Configurable<bool> populateV0MCCoresAsymmetric{"populateV0MCCoresAsymmetric", false, "populate V0MCCores table for derived data analysis, create V0Cores -> V0MCCores interlink. Saves only labeled V0s."};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    if (populateV0MCCoresAsymmetric) {
      LOGF(info, "Asymmetric V0MCCores filling enabled!");
    }
    if (populateV0MCCoresSymmetric) {
      LOGF(info, "Symmetric V0MCCores filling enabled!");
    }
    if (populateV0MCCoresAsymmetric && populateV0MCCoresSymmetric) {
      LOGF(fatal, "Error in configuration: please select only one out of populateV0MCCoresAsymmetric and populateV0MCCoresSymmetric! Crashing!");
    }

    // for storing basic statistics
    auto h = histos.add<TH1>("hBuildingStatistics", "hBuildingStatistics", kTH1F, {{4, -0.5, 3.5f}});
    h->GetXaxis()->SetBinLabel(1, "V0Cores population");
    h->GetXaxis()->SetBinLabel(2, "V0MCCores population");
    h->GetXaxis()->SetBinLabel(3, "x check: duplicates");
    h->GetXaxis()->SetBinLabel(4, "x check: unique");
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // Helper struct to contain V0MCCore information prior to filling
  struct mcV0info {
    int label;
    int motherLabel;
    int pdgCode;
    int pdgCodeMother;
    int pdgCodePositive;
    int pdgCodeNegative;
    bool isPhysicalPrimary;
    std::array<float, 3> xyz;
    std::array<float, 3> posP;
    std::array<float, 3> negP;
    uint64_t packedMcParticleIndices;
  };
  mcV0info thisInfo;
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

  // prong index combiner
  uint64_t combineProngIndices(uint32_t low, uint32_t high)
  {
    return (((uint64_t)high) << 32) | ((uint64_t)low);
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // build V0 labels
  void process(aod::V0Datas const& v0table, aod::McTrackLabels const&, aod::McParticles const& /*particlesMC*/)
  {
    // to be used if using the populateV0MCCoresAsymmetric mode, kept empty otherwise
    std::vector<mcV0info> mcV0infos; // V0MCCore information

    for (auto& v0 : v0table) {
      thisInfo.packedMcParticleIndices = 0; // not de-referenced properly yet
      thisInfo.label = -1;
      thisInfo.motherLabel = -1;
      thisInfo.pdgCode = 0;
      thisInfo.pdgCodeMother = 0;
      thisInfo.pdgCodePositive = 0;
      thisInfo.pdgCodeNegative = 0;
      auto lNegTrack = v0.negTrack_as<aod::McTrackLabels>();
      auto lPosTrack = v0.posTrack_as<aod::McTrackLabels>();

      // Association check
      // There might be smarter ways of doing this in the future
      if (lNegTrack.has_mcParticle() && lPosTrack.has_mcParticle()) {
        thisInfo.packedMcParticleIndices = combineProngIndices(lPosTrack.mcParticleId(), lNegTrack.mcParticleId());
        auto lMCNegTrack = lNegTrack.mcParticle_as<aod::McParticles>();
        auto lMCPosTrack = lPosTrack.mcParticle_as<aod::McParticles>();
        thisInfo.pdgCodePositive = lMCPosTrack.pdgCode();
        thisInfo.pdgCodeNegative = lMCNegTrack.pdgCode();
        thisInfo.posP[0] = lMCPosTrack.px();
        thisInfo.posP[1] = lMCPosTrack.py();
        thisInfo.posP[2] = lMCPosTrack.pz();
        thisInfo.negP[0] = lMCNegTrack.px();
        thisInfo.negP[1] = lMCNegTrack.py();
        thisInfo.negP[2] = lMCNegTrack.pz();
        if (lMCNegTrack.has_mothers() && lMCPosTrack.has_mothers()) {
          for (auto& lNegMother : lMCNegTrack.mothers_as<aod::McParticles>()) {
            for (auto& lPosMother : lMCPosTrack.mothers_as<aod::McParticles>()) {
              if (lNegMother.globalIndex() == lPosMother.globalIndex()) {
                thisInfo.label = lNegMother.globalIndex();
                // acquire information
                thisInfo.xyz[0] = lMCPosTrack.vx();
                thisInfo.xyz[1] = lMCPosTrack.vy();
                thisInfo.xyz[2] = lMCPosTrack.vz();
                thisInfo.pdgCode = lNegMother.pdgCode();
                thisInfo.isPhysicalPrimary = lNegMother.isPhysicalPrimary();
                if (lNegMother.has_mothers()) {
                  for (auto& lNegGrandMother : lNegMother.mothers_as<aod::McParticles>()) {
                    thisInfo.pdgCodeMother = lNegGrandMother.pdgCode();
                    thisInfo.motherLabel = lNegGrandMother.globalIndex();
                  }
                }
              }
            }
          }
        }
      } // end association check
      // Construct label table (note: this will be joinable with V0Datas!)
      v0labels(
        thisInfo.label, thisInfo.motherLabel);

      // ---] Symmetric populate [---
      // in this approach, V0Cores will be joinable with V0MCCores.
      // this is the most pedagogical approach, but it is also more limited
      // and it might use more disk space unnecessarily.
      if (populateV0MCCoresSymmetric) {
        v0mccores(
          thisInfo.label, thisInfo.pdgCode,
          thisInfo.pdgCodeMother, thisInfo.pdgCodePositive, thisInfo.pdgCodeNegative,
          thisInfo.isPhysicalPrimary, thisInfo.xyz[0], thisInfo.xyz[1], thisInfo.xyz[2],
          thisInfo.posP[0], thisInfo.posP[1], thisInfo.posP[2],
          thisInfo.negP[0], thisInfo.negP[1], thisInfo.negP[2]);

        // n.b. placing the interlink index here allows for the writing of
        //      code that is agnostic with respect to the joinability of
        //      V0Cores and V0MCCores (always dereference -> safe)
        v0CoreMCLabels(v0.globalIndex()); // interlink index
      }
      // ---] Asymmetric populate [---
      // in this approach, V0Cores will NOT be joinable with V0MCCores.
      // an additional reference to V0MCCore that IS joinable with V0Cores
      // will be provided to the user.
      if (populateV0MCCoresAsymmetric) {
        int thisV0MCCoreIndex = -1;
        // step 1: check if this element is already provided in the table
        //         using the packedIndices variable calculated above
        for (uint32_t ii = 0; ii < mcV0infos.size(); ii++) {
          if (thisInfo.packedMcParticleIndices == mcV0infos[ii].packedMcParticleIndices && mcV0infos[ii].packedMcParticleIndices > 0) {
            thisV0MCCoreIndex = ii;
            histos.fill(HIST("hBuildingStatistics"), 2.0f); // found
            break;                                          // this exists already in list
          }
        }
        if (thisV0MCCoreIndex < 0) {
          // this V0MCCore does not exist yet. Create it and reference it
          histos.fill(HIST("hBuildingStatistics"), 3.0f); // new
          thisV0MCCoreIndex = mcV0infos.size();
          mcV0infos.push_back(thisInfo);
        }
        v0CoreMCLabels(thisV0MCCoreIndex); // interlink index
      }
    }

    // now populate V0MCCores if in asymmetric mode
    if (populateV0MCCoresAsymmetric) {
      for (auto info : mcV0infos) {
        v0mccores(
          info.label, info.pdgCode,
          info.pdgCodeMother, info.pdgCodePositive, info.pdgCodeNegative,
          info.isPhysicalPrimary, info.xyz[0], info.xyz[1], info.xyz[2],
          info.posP[0], info.posP[1], info.posP[2],
          info.negP[0], info.negP[1], info.negP[2]);
      }
    }

    // collect operating parameters
    histos.fill(HIST("hBuildingStatistics"), 0.0f, v0table.size());
    histos.fill(HIST("hBuildingStatistics"), 1.0f, mcV0infos.size());
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdakzeromcbuilder>(cfgc)};
}
