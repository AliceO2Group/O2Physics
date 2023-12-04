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
//  *+-+*+-+*+-+*+-+*+-+*+-+*
//  Cascade label builder
//  *+-+*+-+*+-+*+-+*+-+*+-+*
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

// For MC association in pre-selection
using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;

//*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
struct cascademcbuilder {
  Produces<aod::McCascLabels> casclabels;       // MC labels for cascades
  Produces<aod::McKFCascLabels> kfcasclabels;   // MC labels for tracked cascades
  Produces<aod::McTraCascLabels> tracasclabels; // MC labels for tracked cascades
  Produces<aod::McCascBBTags> bbtags;           // bb tags (inv structure tagging)
  Produces<aod::CascMCCores> cascmccores;       // optionally aggregate information from MC side for posterior analysis (derived data)

  Configurable<bool> populateCascMCCores{"populateCascMCCores", false, "populate CascMCCores table for derived data analysis"};

  void init(InitContext const&) {}

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // build cascade labels
  void processCascades(aod::CascDatas const& casctable, aod::V0sLinked const&, aod::V0Datas const& v0table, aod::McTrackLabels const&, aod::McParticles const&)
  {
    for (auto& casc : casctable) {
      int pdgCode = -1, pdgCodeMother = -1;
      int pdgCodePositive = -1, pdgCodeNegative = -1, pdgCodeBachelor = -1, pdgCodeV0 = -1;
      bool isPhysicalPrimary = false;
      float xmc = -999.0f, ymc = -999.0f, zmc = -999.0f;
      float xlmc = -999.0f, ylmc = -999.0f, zlmc = -999.0f;
      float pxposmc = -999.0f, pyposmc = -999.0f, pzposmc = -999.0f;
      float pxnegmc = -999.0f, pynegmc = -999.0f, pznegmc = -999.0f;
      float pxbachmc = -999.0f, pybachmc = -999.0f, pzbachmc = -999.0f;
      float px = -999.0f, py = -999.0f, pz = -999.0f;

      // Loop over those that actually have the corresponding V0 associated to them
      auto v0 = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0.has_v0Data())) {
        casclabels(-1);
        continue; // skip those cascades for which V0 doesn't exist (but: should never happen)
      }
      auto v0data = v0.v0Data(); // de-reference index to correct v0data in case it exists
      int lLabel = -1;

      // Acquire all three daughter tracks, please
      auto lBachTrack = casc.bachelor_as<aod::McTrackLabels>();
      auto lNegTrack = v0data.negTrack_as<aod::McTrackLabels>();
      auto lPosTrack = v0data.posTrack_as<aod::McTrackLabels>();

      // Association check
      // There might be smarter ways of doing this in the future
      if (lNegTrack.has_mcParticle() && lPosTrack.has_mcParticle() && lBachTrack.has_mcParticle()) {
        auto lMCBachTrack = lBachTrack.mcParticle_as<aod::McParticles>();
        auto lMCNegTrack = lNegTrack.mcParticle_as<aod::McParticles>();
        auto lMCPosTrack = lPosTrack.mcParticle_as<aod::McParticles>();

        pdgCodePositive = lMCPosTrack.pdgCode();
        pdgCodeNegative = lMCNegTrack.pdgCode();
        pdgCodeBachelor = lMCBachTrack.pdgCode();
        pxposmc = lMCPosTrack.px();
        pyposmc = lMCPosTrack.py();
        pzposmc = lMCPosTrack.pz();
        pxnegmc = lMCNegTrack.px();
        pynegmc = lMCNegTrack.py();
        pznegmc = lMCNegTrack.pz();
        pxbachmc = lMCBachTrack.px();
        pybachmc = lMCBachTrack.py();
        pzbachmc = lMCBachTrack.pz();

        // Step 1: check if the mother is the same, go up a level
        if (lMCNegTrack.has_mothers() && lMCPosTrack.has_mothers()) {
          for (auto& lNegMother : lMCNegTrack.mothers_as<aod::McParticles>()) {
            for (auto& lPosMother : lMCPosTrack.mothers_as<aod::McParticles>()) {
              if (lNegMother == lPosMother) {
                // acquire information
                xlmc = lMCPosTrack.vx();
                ylmc = lMCPosTrack.vy();
                zlmc = lMCPosTrack.vz();
                pdgCodeV0 = lNegMother.pdgCode();

                // if we got to this level, it means the mother particle exists and is the same
                // now we have to go one level up and compare to the bachelor mother too
                for (auto& lV0Mother : lNegMother.mothers_as<aod::McParticles>()) {
                  for (auto& lBachMother : lMCBachTrack.mothers_as<aod::McParticles>()) {
                    if (lV0Mother == lBachMother) {
                      lLabel = lV0Mother.globalIndex();
                      pdgCode = lV0Mother.pdgCode();
                      isPhysicalPrimary = lV0Mother.isPhysicalPrimary();
                      xmc = lMCBachTrack.vx();
                      ymc = lMCBachTrack.vy();
                      zmc = lMCBachTrack.vz();
                      px = lV0Mother.px();
                      py = lV0Mother.py();
                      pz = lV0Mother.pz();
                      if (lV0Mother.has_mothers()) {
                        for (auto& lV0GrandMother : lV0Mother.mothers_as<aod::McParticles>()) {
                          pdgCodeMother = lV0GrandMother.pdgCode();
                        }
                      }
                    }
                  }
                } // end conditional V0-bach pair
              }   // end neg = pos mother conditional
            }
          } // end loop neg/pos mothers
        }   // end conditional of mothers existing
      }     // end association check
      // Construct label table (note: this will be joinable with CascDatas)
      casclabels(
        lLabel);
      if (populateCascMCCores) {
        cascmccores(
          pdgCode, pdgCodeMother, pdgCodeV0, isPhysicalPrimary,
          pdgCodePositive, pdgCodeNegative, pdgCodeBachelor,
          xmc, ymc, zmc, xlmc, ylmc, zlmc,
          pxposmc, pyposmc, pzposmc,
          pxnegmc, pynegmc, pznegmc,
          pxbachmc, pybachmc, pzbachmc,
          px, py, pz);
      }
    } // end casctable loop
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // build kf cascade labels
  void processKFCascades(aod::KFCascDatas const& casctable, aod::V0sLinked const&, aod::V0Datas const& v0table, aod::McTrackLabels const&, aod::McParticles const&)
  {
    for (auto& casc : casctable) {
      // Loop over those that actually have the corresponding V0 associated to them
      auto v0 = casc.v0();
      int lLabel = -1;

      // Acquire all three daughter tracks, please
      auto lBachTrack = casc.bachelor_as<aod::McTrackLabels>();
      auto lNegTrack = v0.negTrack_as<aod::McTrackLabels>();
      auto lPosTrack = v0.posTrack_as<aod::McTrackLabels>();

      // Association check
      // There might be smarter ways of doing this in the future
      if (lNegTrack.has_mcParticle() && lPosTrack.has_mcParticle() && lBachTrack.has_mcParticle()) {
        auto lMCBachTrack = lBachTrack.mcParticle_as<aod::McParticles>();
        auto lMCNegTrack = lNegTrack.mcParticle_as<aod::McParticles>();
        auto lMCPosTrack = lPosTrack.mcParticle_as<aod::McParticles>();
        // Step 1: check if the mother is the same, go up a level
        if (lMCNegTrack.has_mothers() && lMCPosTrack.has_mothers()) {
          for (auto& lNegMother : lMCNegTrack.mothers_as<aod::McParticles>()) {
            for (auto& lPosMother : lMCPosTrack.mothers_as<aod::McParticles>()) {
              if (lNegMother == lPosMother) {
                // if we got to this level, it means the mother particle exists and is the same
                // now we have to go one level up and compare to the bachelor mother too
                for (auto& lV0Mother : lNegMother.mothers_as<aod::McParticles>()) {
                  for (auto& lBachMother : lMCBachTrack.mothers_as<aod::McParticles>()) {
                    if (lV0Mother == lBachMother) {
                      lLabel = lV0Mother.globalIndex();
                    }
                  }
                } // end conditional V0-bach pair
              }   // end neg = pos mother conditional
            }
          } // end loop neg/pos mothers
        }   // end conditional of mothers existing
      }     // end association check
      // Construct label table (note: this will be joinable with CascDatas)
      kfcasclabels(
        lLabel);
    } // end casctable loop
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // build tracked cascade labels
  void processTrackedCascades(aod::TraCascDatas const& casctable, aod::V0sLinked const&, aod::V0Datas const& v0table, aod::McTrackLabels const&, aod::McParticles const&)
  {
    for (auto& casc : casctable) {
      // Loop over those that actually have the corresponding V0 associated to them
      auto v0 = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0.has_v0Data())) {
        tracasclabels(-1);
        continue; // skip those cascades for which V0 doesn't exist
      }
      auto v0data = v0.v0Data(); // de-reference index to correct v0data in case it exists
      int lLabel = -1;

      // Acquire all three daughter tracks, please
      auto lBachTrack = casc.bachelor_as<aod::McTrackLabels>();
      auto lNegTrack = v0data.negTrack_as<aod::McTrackLabels>();
      auto lPosTrack = v0data.posTrack_as<aod::McTrackLabels>();

      // Association check
      // There might be smarter ways of doing this in the future
      if (lNegTrack.has_mcParticle() && lPosTrack.has_mcParticle() && lBachTrack.has_mcParticle()) {
        auto lMCBachTrack = lBachTrack.mcParticle_as<aod::McParticles>();
        auto lMCNegTrack = lNegTrack.mcParticle_as<aod::McParticles>();
        auto lMCPosTrack = lPosTrack.mcParticle_as<aod::McParticles>();
        // Step 1: check if the mother is the same, go up a level
        if (lMCNegTrack.has_mothers() && lMCPosTrack.has_mothers()) {
          for (auto& lNegMother : lMCNegTrack.mothers_as<aod::McParticles>()) {
            for (auto& lPosMother : lMCPosTrack.mothers_as<aod::McParticles>()) {
              if (lNegMother == lPosMother) {
                // if we got to this level, it means the mother particle exists and is the same
                // now we have to go one level up and compare to the bachelor mother too
                for (auto& lV0Mother : lNegMother.mothers_as<aod::McParticles>()) {
                  for (auto& lBachMother : lMCBachTrack.mothers_as<aod::McParticles>()) {
                    if (lV0Mother == lBachMother) {
                      lLabel = lV0Mother.globalIndex();
                    }
                  }
                } // end conditional V0-bach pair
              }   // end neg = pos mother conditional
            }
          } // end loop neg/pos mothers
        }   // end conditional of mothers existing
      }     // end association check
      // Construct label table (note: this will be joinable with CascDatas)
      tracasclabels(
        lLabel);
    } // end casctable loop
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // build cascade labels
  void processBBTags(aod::CascDatas const& casctable, aod::V0sLinked const&, aod::V0Datas const& v0table, aod::McTrackLabels const&, aod::McParticles const&)
  {
    for (auto& casc : casctable) {
      bool bbTag = false; // bachelor-baryon correlation tag to pass

      // Loop over those that actually have the corresponding V0 associated to them
      auto v0 = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0.has_v0Data())) {
        bbtags(bbTag);
        continue; // skip those cascades for which V0 doesn't exist
      }
      auto v0data = v0.v0Data(); // de-reference index to correct v0data in case it exists

      // Acquire all three daughter tracks, please
      auto lBachTrack = casc.bachelor_as<aod::McTrackLabels>();
      auto lNegTrack = v0data.negTrack_as<aod::McTrackLabels>();
      auto lPosTrack = v0data.posTrack_as<aod::McTrackLabels>();

      // Bachelor-baryon association checker
      // this will allow for analyses to pinpoint the effect of spurious, unwanted correlations!

      if (lBachTrack.has_mcParticle()) {
        auto bachelorParticle = lBachTrack.mcParticle_as<aod::McParticles>();
        if (bachelorParticle.pdgCode() == 211) { // pi+, look for antiproton in negative prong
          if (lNegTrack.has_mcParticle()) {
            auto baryonParticle = lNegTrack.mcParticle_as<aod::McParticles>();
            if (baryonParticle.has_mothers() && bachelorParticle.has_mothers() && baryonParticle.pdgCode() == -2212) {
              for (auto& baryonMother : baryonParticle.mothers_as<aod::McParticles>()) {
                for (auto& pionMother : bachelorParticle.mothers_as<aod::McParticles>()) {
                  if (baryonMother.globalIndex() == pionMother.globalIndex() && baryonMother.pdgCode() == -3122) {
                    bbTag = true;
                  }
                }
              }
            }
          }
        }                                         // end if-pion
        if (bachelorParticle.pdgCode() == -211) { // pi-, look for proton in positive prong
          if (lNegTrack.has_mcParticle()) {
            auto baryonParticle = lPosTrack.mcParticle_as<aod::McParticles>();
            if (baryonParticle.has_mothers() && bachelorParticle.has_mothers() && baryonParticle.pdgCode() == 2212) {
              for (auto& baryonMother : baryonParticle.mothers_as<aod::McParticles>()) {
                for (auto& pionMother : bachelorParticle.mothers_as<aod::McParticles>()) {
                  if (baryonMother.globalIndex() == pionMother.globalIndex() && baryonMother.pdgCode() == 3122) {
                    bbTag = true;
                  }
                }
              }
            }
          }
        } // end if-pion
      }   // end bachelor has mcparticle
      // Construct label table (note: this will be joinable with CascDatas)
      bbtags(bbTag);
    } // end casctable loop
  }

  PROCESS_SWITCH(cascademcbuilder, processCascades, "Produce regular cascade label tables", true);
  PROCESS_SWITCH(cascademcbuilder, processKFCascades, "Produce KF cascade label tables", false);
  PROCESS_SWITCH(cascademcbuilder, processTrackedCascades, "Produce tracked cascade label tables", false);
  PROCESS_SWITCH(cascademcbuilder, processBBTags, "Produce cascade bach-baryon correlation tags", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascademcbuilder>(cfgc)};
}
