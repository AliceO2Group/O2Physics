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

  Configurable<bool> populateV0MCCores{"populateV0MCCores", false, "populate V0MCCores table for derived data analysis"};

  void init(InitContext const&) {}

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // build V0 labels
  void process(aod::V0Datas const& v0table, aod::McTrackLabels const&, aod::McParticles const& particlesMC)
  {
    for (auto& v0 : v0table) {
      int lLabel = -1;
      int pdgCode = -1, pdgCodeMother = -1, pdgCodePositive = -1, pdgCodeNegative = -1;
      bool isPhysicalPrimary = false;
      float xmc = -999.0f, ymc = -999.0f, zmc = -999.0f;
      float pxposmc = -999.0f, pyposmc = -999.0f, pzposmc = -999.0f;
      float pxnegmc = -999.0f, pynegmc = -999.0f, pznegmc = -999.0f;

      auto lNegTrack = v0.negTrack_as<aod::McTrackLabels>();
      auto lPosTrack = v0.posTrack_as<aod::McTrackLabels>();

      // Association check
      // There might be smarter ways of doing this in the future
      if (lNegTrack.has_mcParticle() && lPosTrack.has_mcParticle()) {
        auto lMCNegTrack = lNegTrack.mcParticle_as<aod::McParticles>();
        auto lMCPosTrack = lPosTrack.mcParticle_as<aod::McParticles>();
        pdgCodePositive = lMCPosTrack.pdgCode();
        pdgCodeNegative = lMCNegTrack.pdgCode();
        pxposmc = lMCPosTrack.px();
        pyposmc = lMCPosTrack.py();
        pzposmc = lMCPosTrack.pz();
        pxnegmc = lMCNegTrack.px();
        pynegmc = lMCNegTrack.py();
        pznegmc = lMCNegTrack.pz();
        if (lMCNegTrack.has_mothers() && lMCPosTrack.has_mothers()) {
          for (auto& lNegMother : lMCNegTrack.mothers_as<aod::McParticles>()) {
            for (auto& lPosMother : lMCPosTrack.mothers_as<aod::McParticles>()) {
              if (lNegMother.globalIndex() == lPosMother.globalIndex()) {
                lLabel = lNegMother.globalIndex();
                // acquire information
                xmc = lMCPosTrack.vx();
                ymc = lMCPosTrack.vy();
                zmc = lMCPosTrack.vz();
                pdgCode = lNegMother.pdgCode();
                isPhysicalPrimary = lNegMother.isPhysicalPrimary();
                if (lNegMother.has_mothers()) {
                  for (auto& lNegGrandMother : lNegMother.mothers_as<aod::McParticles>()) {
                    pdgCodeMother = lNegGrandMother.pdgCode();
                  }
                }
              }
            }
          }
        }
      } // end association check
      // Construct label table (note: this will be joinable with V0Datas!)
      v0labels(
        lLabel);
      if (populateV0MCCores) {
        v0mccores(
          pdgCode, pdgCodeMother, pdgCodePositive, pdgCodeNegative,
          isPhysicalPrimary, xmc, ymc, zmc,
          pxposmc, pyposmc, pzposmc,
          pxnegmc, pynegmc, pznegmc);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdakzeromcbuilder>(cfgc)};
}
