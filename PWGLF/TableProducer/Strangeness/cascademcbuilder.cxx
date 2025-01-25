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
  Produces<aod::CascCoreMCLabels> cascCoreMClabels; // optionally aggregate information from MC side for posterior analysis (derived data)
  Produces<aod::CascMCCollRefs> cascmccollrefs;     // references MC collisions from MC cascades

  Configurable<bool> populateCascMCCoresSymmetric{"populateCascMCCoresSymmetric", false, "populate CascMCCores table for derived data analysis, keep CascMCCores joinable with CascCores"};
  Configurable<bool> populateCascMCCoresAsymmetric{"populateCascMCCoresAsymmetric", false, "populate CascMCCores table for derived data analysis, create CascCores -> CascMCCores interlink. Saves only labeled Cascades."};

  Configurable<bool> addGeneratedXiMinus{"addGeneratedXiMinus", false, "add CascMCCore entry for generated, not-recoed XiMinus"};
  Configurable<bool> addGeneratedXiPlus{"addGeneratedXiPlus", false, "add CascMCCore entry for generated, not-recoed XiPlus"};
  Configurable<bool> addGeneratedOmegaMinus{"addGeneratedOmegaMinus", false, "add CascMCCore entry for generated, not-recoed OmegaMinus"};
  Configurable<bool> addGeneratedOmegaPlus{"addGeneratedOmegaPlus", false, "add CascMCCore entry for generated, not-recoed OmegaPlus"};

  Configurable<bool> treatPiToMuDecays{"treatPiToMuDecays", true, "if true, will correctly capture pi -> mu and V0 label will still point to originating V0 decay in those cases. Nota bene: prong info will still be for the muon!"};

  Configurable<float> rapidityWindow{"rapidityWindow", 0.5, "rapidity window to save non-recoed candidates"};

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // Helper struct to contain CascMCCore information prior to filling
  struct mcCascinfo {
    int label;
    int motherLabel;
    int mcCollision;
    int pdgCode;
    int pdgCodeMother;
    int pdgCodeV0;
    int pdgCodePositive;
    int pdgCodeNegative;
    int pdgCodeBachelor;
    bool isPhysicalPrimary;
    int processPositive = -1;
    int processNegative = -1;
    int processBachelor = -1;
    std::array<float, 3> xyz;
    std::array<float, 3> lxyz;
    std::array<float, 3> posP;
    std::array<float, 3> negP;
    std::array<float, 3> bachP;
    std::array<float, 3> momentum;
    int mcParticlePositive;
    int mcParticleNegative;
    int mcParticleBachelor;
  };
  mcCascinfo thisInfo;
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

  // kink handling
  template <typename mcpart>
  int getOriginatingParticle(mcpart const& part, int& indexForPositionOfDecay)
  {
    int returnValue = -1;
    if (part.has_mothers()) {
      auto const& motherList = part.template mothers_as<aod::McParticles>();
      if (motherList.size() == 1) {
        for (const auto& mother : motherList) {
          if (std::abs(part.pdgCode()) == 13 && treatPiToMuDecays) {
            // muon decay, de-ref mother twice
            if (mother.has_mothers()) {
              auto grandMotherList = mother.template mothers_as<aod::McParticles>();
              if (grandMotherList.size() == 1) {
                for (const auto& grandMother : grandMotherList) {
                  returnValue = grandMother.globalIndex();
                  indexForPositionOfDecay = mother.globalIndex(); // for V0 decay position: grab muon
                }
              }
            }
          } else {
            returnValue = mother.globalIndex();
            indexForPositionOfDecay = part.globalIndex();
          }
        }
      }
    }
    return returnValue;
  }

  template <typename TCascadeTable, typename TMCParticleTable>
  void generateCascadeMCinfo(TCascadeTable cascTable, TMCParticleTable mcParticles)
  {

    // to be used if using the asymmetric mode, kept empty otherwise
    std::vector<mcCascinfo> mcCascinfos;                           // V0MCCore information
    std::vector<bool> mcParticleIsReco(mcParticles.size(), false); // mc Particle not recoed by V0s

    for (auto& casc : cascTable) {
      thisInfo.pdgCode = -1, thisInfo.pdgCodeMother = -1;
      thisInfo.pdgCodePositive = -1, thisInfo.pdgCodeNegative = -1;
      thisInfo.pdgCodeBachelor = -1, thisInfo.pdgCodeV0 = -1;
      thisInfo.isPhysicalPrimary = false;
      thisInfo.xyz[0] = -999.0f, thisInfo.xyz[1] = -999.0f, thisInfo.xyz[2] = -999.0f;
      thisInfo.lxyz[0] = -999.0f, thisInfo.lxyz[1] = -999.0f, thisInfo.lxyz[2] = -999.0f;
      thisInfo.posP[0] = -999.0f, thisInfo.posP[1] = -999.0f, thisInfo.posP[2] = -999.0f;
      thisInfo.negP[0] = -999.0f, thisInfo.negP[1] = -999.0f, thisInfo.negP[2] = -999.0f;
      thisInfo.bachP[0] = -999.0f, thisInfo.bachP[1] = -999.0f, thisInfo.bachP[2] = -999.0f;
      thisInfo.momentum[0] = -999.0f, thisInfo.momentum[1] = -999.0f, thisInfo.momentum[2] = -999.0f;
      thisInfo.label = -1, thisInfo.motherLabel = -1;
      thisInfo.mcParticlePositive = -1;
      thisInfo.mcParticleNegative = -1;
      thisInfo.mcParticleBachelor = -1;

      // Acquire all three daughter tracks, please
      auto lBachTrack = casc.template bachelor_as<aod::McTrackLabels>();
      auto lNegTrack = casc.template negTrack_as<aod::McTrackLabels>();
      auto lPosTrack = casc.template posTrack_as<aod::McTrackLabels>();

      // Association check
      // There might be smarter ways of doing this in the future
      if (lNegTrack.has_mcParticle() && lPosTrack.has_mcParticle() && lBachTrack.has_mcParticle()) {
        auto lMCBachTrack = lBachTrack.template mcParticle_as<aod::McParticles>();
        auto lMCNegTrack = lNegTrack.template mcParticle_as<aod::McParticles>();
        auto lMCPosTrack = lPosTrack.template mcParticle_as<aod::McParticles>();

        thisInfo.mcParticlePositive = lMCPosTrack.globalIndex();
        thisInfo.mcParticleNegative = lMCNegTrack.globalIndex();
        thisInfo.mcParticleBachelor = lMCBachTrack.globalIndex();
        thisInfo.pdgCodePositive = lMCPosTrack.pdgCode();
        thisInfo.pdgCodeNegative = lMCNegTrack.pdgCode();
        thisInfo.pdgCodeBachelor = lMCBachTrack.pdgCode();
        thisInfo.posP[0] = lMCPosTrack.px();
        thisInfo.posP[1] = lMCPosTrack.py();
        thisInfo.posP[2] = lMCPosTrack.pz();
        thisInfo.negP[0] = lMCNegTrack.px();
        thisInfo.negP[1] = lMCNegTrack.py();
        thisInfo.negP[2] = lMCNegTrack.pz();
        thisInfo.bachP[0] = lMCBachTrack.px();
        thisInfo.bachP[1] = lMCBachTrack.py();
        thisInfo.bachP[2] = lMCBachTrack.pz();
        thisInfo.processPositive = lMCPosTrack.getProcess();
        thisInfo.processNegative = lMCNegTrack.getProcess();
        thisInfo.processBachelor = lMCBachTrack.getProcess();

        // Step 0: treat pi -> mu + antineutrino
        // if present, de-reference original V0 correctly and provide label to original object
        // NOTA BENE: the prong info will still correspond to a muon, treat carefully!
        int negOriginating = -1, posOriginating = -1, bachOriginating = -1;
        int particleForLambdaDecayPositionIdx = -1, particleForCascadeDecayPositionIdx = -1;
        negOriginating = getOriginatingParticle(lMCNegTrack, particleForLambdaDecayPositionIdx);
        posOriginating = getOriginatingParticle(lMCPosTrack, particleForLambdaDecayPositionIdx);
        bachOriginating = getOriginatingParticle(lMCBachTrack, particleForCascadeDecayPositionIdx);

        if (negOriginating > -1 && negOriginating == posOriginating) {
          auto originatingV0 = mcParticles.rawIteratorAt(negOriginating);
          auto particleForLambdaDecayPosition = mcParticles.rawIteratorAt(particleForLambdaDecayPositionIdx);

          thisInfo.label = originatingV0.globalIndex();
          thisInfo.lxyz[0] = particleForLambdaDecayPosition.vx();
          thisInfo.lxyz[1] = particleForLambdaDecayPosition.vy();
          thisInfo.lxyz[2] = particleForLambdaDecayPosition.vz();
          thisInfo.pdgCodeV0 = originatingV0.pdgCode();

          if (originatingV0.has_mothers()) {
            for (auto& lV0Mother : originatingV0.template mothers_as<aod::McParticles>()) {
              if (lV0Mother.globalIndex() == bachOriginating) { // found mother particle
                thisInfo.label = lV0Mother.globalIndex();

                if (lV0Mother.has_mcCollision()) {
                  thisInfo.mcCollision = lV0Mother.mcCollisionId(); // save this reference, please
                }

                thisInfo.pdgCode = lV0Mother.pdgCode();
                thisInfo.isPhysicalPrimary = lV0Mother.isPhysicalPrimary();
                thisInfo.xyz[0] = originatingV0.vx();
                thisInfo.xyz[1] = originatingV0.vy();
                thisInfo.xyz[2] = originatingV0.vz();
                thisInfo.momentum[0] = lV0Mother.px();
                thisInfo.momentum[1] = lV0Mother.py();
                thisInfo.momentum[2] = lV0Mother.pz();
                if (lV0Mother.has_mothers()) {
                  for (auto& lV0GrandMother : lV0Mother.template mothers_as<aod::McParticles>()) {
                    thisInfo.pdgCodeMother = lV0GrandMother.pdgCode();
                    thisInfo.motherLabel = lV0GrandMother.globalIndex();
                  }
                }
              }
            } // end v0 mother loop
          } // end has_mothers check for V0
        } // end conditional of pos/neg originating being the same
      }     // end association check
      // Construct label table (note: this will be joinable with CascDatas)
      casclabels(
        thisInfo.label, thisInfo.motherLabel);

      // Mark mcParticle as recoed (no searching necessary afterwards)
      if (thisInfo.label > -1) {
        mcParticleIsReco[thisInfo.label] = true;
      }

      if (populateCascMCCoresSymmetric) {
        cascmccores(
          thisInfo.pdgCode, thisInfo.pdgCodeMother, thisInfo.pdgCodeV0, thisInfo.isPhysicalPrimary,
          thisInfo.pdgCodePositive, thisInfo.pdgCodeNegative, thisInfo.pdgCodeBachelor,
          thisInfo.xyz[0], thisInfo.xyz[1], thisInfo.xyz[2],
          thisInfo.lxyz[0], thisInfo.lxyz[1], thisInfo.lxyz[2],
          thisInfo.posP[0], thisInfo.posP[1], thisInfo.posP[2],
          thisInfo.negP[0], thisInfo.negP[1], thisInfo.negP[2],
          thisInfo.bachP[0], thisInfo.bachP[1], thisInfo.bachP[2],
          thisInfo.momentum[0], thisInfo.momentum[1], thisInfo.momentum[2]);
        cascmccollrefs(thisInfo.mcCollision);
      }

      if (populateCascMCCoresAsymmetric) {
        int thisCascMCCoreIndex = -1;
        // step 1: check if this element is already provided in the table
        //         using the packedIndices variable calculated above
        for (uint32_t ii = 0; ii < mcCascinfos.size(); ii++) {
          if (thisInfo.label == mcCascinfos[ii].label && mcCascinfos[ii].label > -1) {
            thisCascMCCoreIndex = ii;
            break; // this exists already in list
          }
        }
        if (thisCascMCCoreIndex < 0) {
          // this CascMCCore does not exist yet. Create it and reference it
          thisCascMCCoreIndex = mcCascinfos.size();
          mcCascinfos.push_back(thisInfo);
        }
        cascCoreMClabels(thisCascMCCoreIndex); // interlink: reconstructed -> MC index
      }
    } // end casctable loop

    // now populate V0MCCores if in asymmetric mode
    if (populateCascMCCoresAsymmetric) {
      // first step: add any un-recoed v0mmcores that were requested
      for (auto& mcParticle : mcParticles) {
        thisInfo.pdgCode = -1, thisInfo.pdgCodeMother = -1;
        thisInfo.pdgCodePositive = -1, thisInfo.pdgCodeNegative = -1;
        thisInfo.pdgCodeBachelor = -1, thisInfo.pdgCodeV0 = -1;
        thisInfo.isPhysicalPrimary = false;
        thisInfo.xyz[0] = 0.0f, thisInfo.xyz[1] = 0.0f, thisInfo.xyz[2] = 0.0f;
        thisInfo.lxyz[0] = 0.0f, thisInfo.lxyz[1] = 0.0f, thisInfo.lxyz[2] = 0.0f;
        thisInfo.posP[0] = 0.0f, thisInfo.posP[1] = 0.0f, thisInfo.posP[2] = 0.0f;
        thisInfo.negP[0] = 0.0f, thisInfo.negP[1] = 0.0f, thisInfo.negP[2] = 0.0f;
        thisInfo.bachP[0] = 0.0f, thisInfo.bachP[1] = 0.0f, thisInfo.bachP[2] = 0.0f;
        thisInfo.momentum[0] = 0.0f, thisInfo.momentum[1] = 0.0f, thisInfo.momentum[2] = 0.0f;
        thisInfo.label = -1, thisInfo.motherLabel = -1;
        thisInfo.mcParticlePositive = -1;
        thisInfo.mcParticleNegative = -1;
        thisInfo.mcParticleBachelor = -1;

        if (mcParticleIsReco[mcParticle.globalIndex()] == true)
          continue; // skip if already created in list

        if (TMath::Abs(mcParticle.y()) > rapidityWindow)
          continue; // skip outside midrapidity

        if (
          (addGeneratedXiMinus && mcParticle.pdgCode() == 3312) ||
          (addGeneratedXiPlus && mcParticle.pdgCode() == -3312) ||
          (addGeneratedOmegaMinus && mcParticle.pdgCode() == 3334) ||
          (addGeneratedOmegaPlus && mcParticle.pdgCode() == -3334)) {
          thisInfo.pdgCode = mcParticle.pdgCode();
          thisInfo.isPhysicalPrimary = mcParticle.isPhysicalPrimary();

          if (mcParticle.has_mcCollision()) {
            thisInfo.mcCollision = mcParticle.mcCollisionId(); // save this reference, please
          }
          thisInfo.momentum[0] = mcParticle.px();
          thisInfo.momentum[1] = mcParticle.py();
          thisInfo.momentum[2] = mcParticle.pz();
          thisInfo.label = mcParticle.globalIndex();

          if (mcParticle.has_daughters()) {
            auto const& daughters = mcParticle.template daughters_as<aod::McParticles>();
            for (auto& dau : daughters) {
              if (dau.getProcess() != 4) // check whether the daughter comes from a decay
                continue;

              if (TMath::Abs(dau.pdgCode()) == 211 || TMath::Abs(dau.pdgCode()) == 321) {
                thisInfo.pdgCodeBachelor = dau.pdgCode();
                thisInfo.bachP[0] = dau.px();
                thisInfo.bachP[1] = dau.py();
                thisInfo.bachP[2] = dau.pz();
                thisInfo.xyz[0] = dau.vx();
                thisInfo.xyz[1] = dau.vy();
                thisInfo.xyz[2] = dau.vz();
                thisInfo.mcParticleBachelor = dau.globalIndex();
              }
              if (TMath::Abs(dau.pdgCode()) == 2212) {
                thisInfo.pdgCodeV0 = dau.pdgCode();

                for (auto& v0Dau : dau.template daughters_as<aod::McParticles>()) {
                  if (v0Dau.getProcess() != 4)
                    continue;

                  if (v0Dau.pdgCode() > 0) {
                    thisInfo.pdgCodePositive = v0Dau.pdgCode();
                    thisInfo.processPositive = v0Dau.getProcess();
                    thisInfo.posP[0] = v0Dau.px();
                    thisInfo.posP[1] = v0Dau.py();
                    thisInfo.posP[2] = v0Dau.pz();
                    thisInfo.lxyz[0] = v0Dau.vx();
                    thisInfo.lxyz[1] = v0Dau.vy();
                    thisInfo.lxyz[2] = v0Dau.vz();
                    thisInfo.mcParticlePositive = v0Dau.globalIndex();
                  }
                  if (v0Dau.pdgCode() < 0) {
                    thisInfo.pdgCodeNegative = v0Dau.pdgCode();
                    thisInfo.processNegative = v0Dau.getProcess();
                    thisInfo.negP[0] = v0Dau.px();
                    thisInfo.negP[1] = v0Dau.py();
                    thisInfo.negP[2] = v0Dau.pz();
                    thisInfo.mcParticleNegative = v0Dau.globalIndex();
                  }
                }
              }
            }
          }

          // if I got here, it means this MC particle was not recoed and is of interest. Add it please
          mcCascinfos.push_back(thisInfo);
        }
      }

      for (auto thisInfo : mcCascinfos) {
        cascmccores( // a lot of the info below will be compressed in case of not-recoed MC (good!)
          thisInfo.pdgCode, thisInfo.pdgCodeMother, thisInfo.pdgCodeV0, thisInfo.isPhysicalPrimary,
          thisInfo.pdgCodePositive, thisInfo.pdgCodeNegative, thisInfo.pdgCodeBachelor,
          thisInfo.xyz[0], thisInfo.xyz[1], thisInfo.xyz[2],
          thisInfo.lxyz[0], thisInfo.lxyz[1], thisInfo.lxyz[2],
          thisInfo.posP[0], thisInfo.posP[1], thisInfo.posP[2],
          thisInfo.negP[0], thisInfo.negP[1], thisInfo.negP[2],
          thisInfo.bachP[0], thisInfo.bachP[1], thisInfo.bachP[2],
          thisInfo.momentum[0], thisInfo.momentum[1], thisInfo.momentum[2]);
        cascmccollrefs(thisInfo.mcCollision);
      }
    }
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // build cascade labels
  void processCascades(aod::CascDatas const& casctable, aod::V0sLinked const&, aod::V0Datas const& /*v0table*/, aod::McTrackLabels const&, aod::McParticles const& mcParticles)
  {
    generateCascadeMCinfo(casctable, mcParticles);
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // build findable cascade labels
  void processFindableCascades(aod::CascDatas const& casctable, aod::FindableV0sLinked const&, aod::V0Datas const& /*v0table*/, aod::McTrackLabels const&, aod::McParticles const& mcParticles)
  {
    generateCascadeMCinfo(casctable, mcParticles);
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // build kf cascade labels
  void processKFCascades(aod::KFCascDatas const& casctable, aod::V0s const&, aod::McTrackLabels const&, aod::McParticles const&)
  {
    for (auto& casc : casctable) {
      int lLabel = -1;

      // Acquire all three daughter tracks, please
      auto lBachTrack = casc.template bachelor_as<aod::McTrackLabels>();
      auto lNegTrack = casc.template negTrack_as<aod::McTrackLabels>();
      auto lPosTrack = casc.template posTrack_as<aod::McTrackLabels>();

      // Association check
      // There might be smarter ways of doing this in the future
      if (lNegTrack.has_mcParticle() && lPosTrack.has_mcParticle() && lBachTrack.has_mcParticle()) {
        auto lMCBachTrack = lBachTrack.template mcParticle_as<aod::McParticles>();
        auto lMCNegTrack = lNegTrack.template mcParticle_as<aod::McParticles>();
        auto lMCPosTrack = lPosTrack.template mcParticle_as<aod::McParticles>();
        // Step 1: check if the mother is the same, go up a level
        if (lMCNegTrack.has_mothers() && lMCPosTrack.has_mothers()) {
          for (auto& lNegMother : lMCNegTrack.template mothers_as<aod::McParticles>()) {
            for (auto& lPosMother : lMCPosTrack.template mothers_as<aod::McParticles>()) {
              if (lNegMother == lPosMother) {
                // if we got to this level, it means the mother particle exists and is the same
                // now we have to go one level up and compare to the bachelor mother too
                for (auto& lV0Mother : lNegMother.template mothers_as<aod::McParticles>()) {
                  for (auto& lBachMother : lMCBachTrack.template mothers_as<aod::McParticles>()) {
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
  void processTrackedCascades(aod::TraCascDatas const& casctable, aod::V0sLinked const&, aod::V0Datas const& /*v0table*/, aod::McTrackLabels const&, aod::McParticles const&)
  {
    for (auto& casc : casctable) {
      int lLabel = -1;

      // Acquire all three daughter tracks, please
      auto lBachTrack = casc.bachelor_as<aod::McTrackLabels>();
      auto lNegTrack = casc.negTrack_as<aod::McTrackLabels>();
      auto lPosTrack = casc.posTrack_as<aod::McTrackLabels>();

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
  void processBBTags(aod::CascDatas const& casctable, aod::V0sLinked const&, aod::V0Datas const& /*v0table*/, aod::McTrackLabels const&, aod::McParticles const&)
  {
    for (auto& casc : casctable) {
      bool bbTag = false; // bachelor-baryon correlation tag to pass

      // Acquire all three daughter tracks, please
      auto lBachTrack = casc.bachelor_as<aod::McTrackLabels>();
      auto lNegTrack = casc.negTrack_as<aod::McTrackLabels>();
      auto lPosTrack = casc.posTrack_as<aod::McTrackLabels>();

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
  PROCESS_SWITCH(cascademcbuilder, processFindableCascades, "Produce findable cascade label tables", false);
  PROCESS_SWITCH(cascademcbuilder, processKFCascades, "Produce KF cascade label tables", false);
  PROCESS_SWITCH(cascademcbuilder, processTrackedCascades, "Produce tracked cascade label tables", false);
  PROCESS_SWITCH(cascademcbuilder, processBBTags, "Produce cascade bach-baryon correlation tags", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascademcbuilder>(cfgc)};
}
