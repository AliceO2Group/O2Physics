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
// Cascade labeler task
// ====================
//
// This code loops over a CascData table and produces a
// standard table of Cascade -> McParticle indices. This
// is meant to simplify the acquisiton of the correct MC
// particle when associating at analysis level.
//
// The existence of these labels will also mean that
// users can easily select Cascades that actually have
// MC associations to them.
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    david.dobrigkeit.chinellato@cern.ch
//

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/StrangenessTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include <CCDB/BasicCCDBManager.h>

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;

/// Cascade builder task: rebuilds cascades
struct cascadeLabeler {
  Produces<aod::McCascLabels> casclabels;

  //for bookkeeping purposes: how many V0s come from same mother etc
  HistogramRegistry registry{
    "registry",
    {
      {"hLabelCounter", "hLabelCounter", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
    },
  };

  void process(aod::Collisions::iterator const& collision, aod::CascDataExt const& casctable, aod::V0sLinked const&, aod::V0Datas const& v0table, aod::Tracks const& tracks, aod::McParticles const& particlesMC)
  {
    for (auto& casc : casctable) {

      //Loop over those that actually have the corresponding V0 associated to them
      auto v0index = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0index.has_v0Data())) {
        continue; // skip those cascades for which V0 doesn't exist
      }
      auto v0 = v0index.v0Data(); // de-reference index to correct v0data in case it exists

      int lLabel = -1;
      float lFillVal = 0.5f; //all considered V0s

      //Acquire all three daughter tracks, please
      auto lBachTrack = casc.bachelor_as<LabeledTracks>();
      auto lNegTrack = v0.negTrack_as<LabeledTracks>();
      auto lPosTrack = v0.posTrack_as<LabeledTracks>();

      //Association check
      //There might be smarter ways of doing this in the future
      if (lNegTrack.has_mcParticle() && lPosTrack.has_mcParticle() && lBachTrack.has_mcParticle()) {
        auto lMCBachTrack = lNegTrack.mcParticle_as<aod::McParticles>();
        auto lMCNegTrack = lNegTrack.mcParticle_as<aod::McParticles>();
        auto lMCPosTrack = lPosTrack.mcParticle_as<aod::McParticles>();

        //Step 1: check if the mother is the same, go up a level
        if (lMCNegTrack.has_mothers() && lMCPosTrack.has_mothers()) {

          for (auto& lNegMother : lMCNegTrack.mothers_as<aod::McParticles>()) {
            for (auto& lPosMother : lMCPosTrack.mothers_as<aod::McParticles>()) {
              if (lNegMother == lPosMother) {
                //if we got to this level, it means the mother particle exists and is the same
                //now we have to go one level up and compare to the bachelor mother too
                for (auto& lV0Mother : lNegMother.mothers_as<aod::McParticles>()) {
                  for (auto& lBachMother : lMCBachTrack.mothers_as<aod::McParticles>()) {
                    if (lV0Mother == lBachMother) {
                      lLabel = lV0Mother.globalIndex();
                      lFillVal = 1.5f; //v0s with the same mother
                    }
                  }
                } //end conditional V0-bach pair
              }   //end neg = pos mother conditional
            }
          } //end loop neg/pos mothers
        }   //end conditional of mothers existing
      }     //end association check

      registry.fill(HIST("hLabelCounter"), lFillVal);

      //Construct label table (note: this will be joinable with V0Datas)
      casclabels(
        lLabel);
    } //end casctable loop
  }   //end process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascadeLabeler>(cfgc)};
}
