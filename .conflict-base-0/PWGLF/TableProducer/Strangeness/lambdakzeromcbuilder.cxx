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

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <map>
#include <utility>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

//*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
struct lambdakzeromcbuilder {
  Produces<aod::McV0Labels> v0labels;           // MC labels for V0s
  Produces<aod::V0MCCores> v0mccores;           // optionally aggregate information from MC side for posterior analysis (derived data)
  Produces<aod::V0CoreMCLabels> v0CoreMCLabels; // interlink V0Cores -> V0MCCores in asymmetric mode
  Produces<aod::V0MCCollRefs> v0mccollref;      // references collisions from V0MCCores

  Configurable<bool> populateV0MCCoresSymmetric{"populateV0MCCoresSymmetric", false, "populate V0MCCores table for derived data analysis, keep V0MCCores joinable with V0Cores"};
  Configurable<bool> populateV0MCCoresAsymmetric{"populateV0MCCoresAsymmetric", false, "populate V0MCCores table for derived data analysis, create V0Cores -> V0MCCores interlink. Saves only labeled V0s."};

  Configurable<bool> addGeneratedK0Short{"addGeneratedK0Short", false, "add V0MCCore entry for generated, not-recoed K0Short"};
  Configurable<bool> addGeneratedLambda{"addGeneratedLambda", false, "add V0MCCore entry for generated, not-recoed Lambda"};
  Configurable<bool> addGeneratedAntiLambda{"addGeneratedAntiLambda", false, "add V0MCCore entry for generated, not-recoed AntiLambda"};
  Configurable<bool> addGeneratedGamma{"addGeneratedGamma", false, "add V0MCCore entry for generated, not-recoed Gamma"};

  Configurable<bool> treatPiToMuDecays{"treatPiToMuDecays", true, "if true, will correctly capture pi -> mu and V0 label will still point to originating V0 decay in those cases. Nota bene: prong info will still be for the muon!"};

  Configurable<float> rapidityWindow{"rapidityWindow", 0.5, "rapidity window to save non-recoed candidates"};

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

    auto hK0s = histos.add<TH1>("hStatisticsK0s", "hBuildingStatisticsK0s", kTH1F, {{3, -0.5, 2.5f}});
    hK0s->GetXaxis()->SetBinLabel(1, "MC associated to reco.");
    hK0s->GetXaxis()->SetBinLabel(2, "Not originating from decay");
    hK0s->GetXaxis()->SetBinLabel(3, "MC with un-recoed V0");

    auto hLambda = histos.add<TH1>("hStatisticsLambda", "hBuildingStatisticsLambda", kTH1F, {{3, -0.5, 2.5f}});
    hLambda->GetXaxis()->SetBinLabel(1, "MC associated to reco.");
    hLambda->GetXaxis()->SetBinLabel(2, "Not originating from decay");
    hLambda->GetXaxis()->SetBinLabel(3, "MC with un-recoed V0");

    auto hAlambda = histos.add<TH1>("hStatisticsAlambda", "hBuildingStatisticsAlambda", kTH1F, {{3, -0.5, 2.5f}});
    hAlambda->GetXaxis()->SetBinLabel(1, "MC associated to reco.");
    hAlambda->GetXaxis()->SetBinLabel(2, "Not originating from decay");
    hAlambda->GetXaxis()->SetBinLabel(3, "MC with un-recoed V0");
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // Helper struct to contain V0MCCore information prior to filling
  struct mcV0info {
    int label = -1;
    int motherLabel = -1;
    int pdgCode = 0;
    int pdgCodeMother = 0;
    int pdgCodePositive = 0;
    int pdgCodeNegative = 0;
    int mcCollision = -1;
    bool isPhysicalPrimary = false;
    int processPositive = -1;
    int processNegative = -1;
    std::array<float, 3> xyz;
    std::array<float, 3> posP;
    std::array<float, 3> negP;
    std::array<float, 3> momentum;
  };
  mcV0info thisInfo;
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

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // build V0 labels
  void process(aod::V0Datas const& v0table, aod::McTrackLabels const&, aod::McParticles const& mcParticles)
  {
    // to be used if using the populateV0MCCoresAsymmetric mode, kept empty otherwise
    std::vector<mcV0info> mcV0infos;                               // V0MCCore information
    std::vector<bool> mcParticleIsReco(mcParticles.size(), false); // mc Particle not recoed by V0s

    for (auto& v0 : v0table) {
      thisInfo.label = -1;
      thisInfo.motherLabel = -1;
      thisInfo.pdgCode = 0;
      thisInfo.pdgCodeMother = 0;
      thisInfo.pdgCodePositive = 0;
      thisInfo.pdgCodeNegative = 0;
      thisInfo.mcCollision = -1;
      thisInfo.xyz[0] = thisInfo.xyz[1] = thisInfo.xyz[2] = 0.0f;
      thisInfo.posP[0] = thisInfo.posP[1] = thisInfo.posP[2] = 0.0f;
      thisInfo.negP[0] = thisInfo.negP[1] = thisInfo.negP[2] = 0.0f;
      thisInfo.momentum[0] = thisInfo.momentum[1] = thisInfo.momentum[2] = 0.0f;
      auto lNegTrack = v0.negTrack_as<aod::McTrackLabels>();
      auto lPosTrack = v0.posTrack_as<aod::McTrackLabels>();

      // Association check
      // There might be smarter ways of doing this in the future
      if (lNegTrack.has_mcParticle() && lPosTrack.has_mcParticle()) {
        auto lMCNegTrack = lNegTrack.mcParticle_as<aod::McParticles>();
        auto lMCPosTrack = lPosTrack.mcParticle_as<aod::McParticles>();

        thisInfo.pdgCodePositive = lMCPosTrack.pdgCode();
        thisInfo.pdgCodeNegative = lMCNegTrack.pdgCode();
        thisInfo.processPositive = lMCPosTrack.getProcess();
        thisInfo.processNegative = lMCNegTrack.getProcess();
        thisInfo.posP[0] = lMCPosTrack.px();
        thisInfo.posP[1] = lMCPosTrack.py();
        thisInfo.posP[2] = lMCPosTrack.pz();
        thisInfo.negP[0] = lMCNegTrack.px();
        thisInfo.negP[1] = lMCNegTrack.py();
        thisInfo.negP[2] = lMCNegTrack.pz();

        // check for pi -> mu + antineutrino decay
        // if present, de-reference original V0 correctly and provide label to original object
        // NOTA BENE: the prong info will still correspond to a muon, treat carefully!
        int negOriginating = -1, posOriginating = -1, particleForDecayPositionIdx = -1;
        negOriginating = getOriginatingParticle(lMCNegTrack, particleForDecayPositionIdx);
        posOriginating = getOriginatingParticle(lMCPosTrack, particleForDecayPositionIdx);

        if (negOriginating > -1 && negOriginating == posOriginating) {
          auto originatingV0 = mcParticles.rawIteratorAt(negOriginating);
          auto particleForDecayPosition = mcParticles.rawIteratorAt(particleForDecayPositionIdx);

          thisInfo.label = originatingV0.globalIndex();
          thisInfo.xyz[0] = particleForDecayPosition.vx();
          thisInfo.xyz[1] = particleForDecayPosition.vy();
          thisInfo.xyz[2] = particleForDecayPosition.vz();

          // MC pos. and neg. daughters are the same! Looking for replacement...
          // if (lMCPosTrack.globalIndex() == lMCNegTrack.globalIndex()) {
          //   auto const& daughters = lNegMother.daughters_as<aod::McParticles>();
          //   for (auto& ldau : daughters) {
          //     // check if the candidate originates from a decay
          //     // if not, this is not a suitable candidate for one of the decay daughters
          //     if (ldau.getProcess() != 4) // see TMCProcess.h
          //       continue;

          //     if (lMCPosTrack.pdgCode() < 0 && ldau.pdgCode() > 0) { // the positive track needs to be changed
          //       thisInfo.pdgCodePositive = ldau.pdgCode();
          //       thisInfo.processPositive = ldau.getProcess();
          //       thisInfo.posP[0] = ldau.px();
          //       thisInfo.posP[1] = ldau.py();
          //       thisInfo.posP[2] = ldau.pz();
          //       thisInfo.xyz[0] = ldau.vx();
          //       thisInfo.xyz[1] = ldau.vy();
          //       thisInfo.xyz[2] = ldau.vz();
          //     }
          //     if (lMCNegTrack.pdgCode() > 0 && ldau.pdgCode() < 0) { // the negative track needs to be changed
          //       thisInfo.pdgCodeNegative = ldau.pdgCode();
          //       thisInfo.processNegative = ldau.getProcess();
          //       thisInfo.negP[0] = ldau.px();
          //       thisInfo.negP[1] = ldau.py();
          //       thisInfo.negP[2] = ldau.pz();
          //     }
          //   }
          // }

          if (originatingV0.has_mcCollision()) {
            thisInfo.mcCollision = originatingV0.mcCollisionId(); // save this reference, please
          }

          // acquire information
          thisInfo.pdgCode = originatingV0.pdgCode();
          thisInfo.isPhysicalPrimary = originatingV0.isPhysicalPrimary();
          thisInfo.momentum[0] = originatingV0.px();
          thisInfo.momentum[1] = originatingV0.py();
          thisInfo.momentum[2] = originatingV0.pz();

          if (originatingV0.has_mothers()) {
            for (auto& lV0Mother : originatingV0.mothers_as<aod::McParticles>()) {
              thisInfo.pdgCodeMother = lV0Mother.pdgCode();
              thisInfo.motherLabel = lV0Mother.globalIndex();
            }
          }
        }

      } // end association check
      // Construct label table (note: this will be joinable with V0Datas!)
      v0labels(
        thisInfo.label, thisInfo.motherLabel);

      // Mark mcParticle as recoed (no searching necessary afterwards)
      if (thisInfo.label > -1) {
        mcParticleIsReco[thisInfo.label] = true;
      }

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
          thisInfo.negP[0], thisInfo.negP[1], thisInfo.negP[2],
          thisInfo.momentum[0], thisInfo.momentum[1], thisInfo.momentum[2]);
        v0mccollref(thisInfo.mcCollision);

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
          if (thisInfo.label == mcV0infos[ii].label && mcV0infos[ii].label > -1) {
            thisV0MCCoreIndex = ii;
            histos.fill(HIST("hBuildingStatistics"), 2.0f); // found
            break;                                          // this exists already in list
          }
        }
        if (thisV0MCCoreIndex < 0 && thisInfo.label > -1) {
          // this V0MCCore does not exist yet. Create it and reference it
          histos.fill(HIST("hBuildingStatistics"), 3.0f); // new
          thisV0MCCoreIndex = mcV0infos.size();
          mcV0infos.push_back(thisInfo);

          // For bookkeeping
          if (thisInfo.label > -1 && thisInfo.isPhysicalPrimary) {
            float ymc = 1e3;
            if (thisInfo.pdgCode == 310)
              ymc = RecoDecay::y(std::array{thisInfo.posP[0] + thisInfo.negP[0], thisInfo.posP[1] + thisInfo.negP[1], thisInfo.posP[2] + thisInfo.negP[2]}, o2::constants::physics::MassKaonNeutral);
            else if (TMath::Abs(thisInfo.pdgCode) == 3122)
              ymc = RecoDecay::y(std::array{thisInfo.posP[0] + thisInfo.negP[0], thisInfo.posP[1] + thisInfo.negP[1], thisInfo.posP[2] + thisInfo.negP[2]}, o2::constants::physics::MassLambda);

            if (thisInfo.pdgCode == 310 && TMath::Abs(ymc) < rapidityWindow) {
              histos.fill(HIST("hStatisticsK0s"), 0.0f); // found
              if (thisInfo.processPositive != 4 || thisInfo.processNegative != 4) {
                histos.fill(HIST("hStatisticsK0s"), 1.0f); // Not originating from decay
              }
            }
            if (thisInfo.pdgCode == 3122 && TMath::Abs(ymc) < rapidityWindow) {
              histos.fill(HIST("hStatisticsLambda"), 0.0f); // found
              if (thisInfo.processPositive != 4 || thisInfo.processNegative != 4) {
                histos.fill(HIST("hStatisticsLambda"), 1.0f); // Not originating from decay
              }
            }
            if (thisInfo.pdgCode == -3122 && TMath::Abs(ymc) < rapidityWindow) {
              histos.fill(HIST("hStatisticsAlambda"), 0.0f); // found
              if (thisInfo.processPositive != 4 || thisInfo.processNegative != 4) {
                histos.fill(HIST("hStatisticsAlambda"), 1.0f); // Not originating from decay
              }
            }
          }
        }
        v0CoreMCLabels(thisV0MCCoreIndex); // interlink index
      }
    }

    // now populate V0MCCores if in asymmetric mode
    if (populateV0MCCoresAsymmetric) {
      // first step: add any un-recoed v0mmcores that were requested
      for (auto& mcParticle : mcParticles) {
        thisInfo.label = -1;
        thisInfo.motherLabel = -1;
        thisInfo.pdgCode = 0;
        thisInfo.pdgCodeMother = -1;
        thisInfo.pdgCodePositive = -1;
        thisInfo.pdgCodeNegative = -1;
        thisInfo.mcCollision = -1;
        thisInfo.xyz[0] = thisInfo.xyz[1] = thisInfo.xyz[2] = 0.0f;
        thisInfo.posP[0] = thisInfo.posP[1] = thisInfo.posP[2] = 0.0f;
        thisInfo.negP[0] = thisInfo.negP[1] = thisInfo.negP[2] = 0.0f;
        thisInfo.momentum[0] = thisInfo.momentum[1] = thisInfo.momentum[2] = 0.0f;

        if (mcParticleIsReco[mcParticle.globalIndex()] == true)
          continue; // skip if already created in list

        if (TMath::Abs(mcParticle.y()) > rapidityWindow)
          continue; // skip outside midrapidity

        if (
          (addGeneratedK0Short && mcParticle.pdgCode() == 310) ||
          (addGeneratedLambda && mcParticle.pdgCode() == 3122) ||
          (addGeneratedAntiLambda && mcParticle.pdgCode() == -3122) ||
          (addGeneratedGamma && mcParticle.pdgCode() == 22)) {
          thisInfo.pdgCode = mcParticle.pdgCode();
          thisInfo.isPhysicalPrimary = mcParticle.isPhysicalPrimary();
          thisInfo.label = mcParticle.globalIndex();

          if (mcParticle.has_mcCollision()) {
            thisInfo.mcCollision = mcParticle.mcCollisionId(); // save this reference, please
          }

          //
          thisInfo.momentum[0] = mcParticle.px();
          thisInfo.momentum[1] = mcParticle.py();
          thisInfo.momentum[2] = mcParticle.pz();

          if (mcParticle.has_mothers()) {
            auto const& mother = mcParticle.mothers_first_as<aod::McParticles>();
            thisInfo.pdgCodeMother = mother.pdgCode();
            thisInfo.motherLabel = mother.globalIndex();
          }
          if (mcParticle.has_daughters()) {
            auto const& daughters = mcParticle.daughters_as<aod::McParticles>();

            for (auto& dau : daughters) {
              if (dau.getProcess() != 4)
                continue;

              if (dau.pdgCode() > 0) {
                thisInfo.pdgCodePositive = dau.pdgCode();
                thisInfo.processPositive = dau.getProcess();
                thisInfo.posP[0] = dau.px();
                thisInfo.posP[1] = dau.py();
                thisInfo.posP[2] = dau.pz();
                thisInfo.xyz[0] = dau.vx();
                thisInfo.xyz[1] = dau.vy();
                thisInfo.xyz[2] = dau.vz();
              }
              if (dau.pdgCode() < 0) {
                thisInfo.pdgCodeNegative = dau.pdgCode();
                thisInfo.processNegative = dau.getProcess();
                thisInfo.negP[0] = dau.px();
                thisInfo.negP[1] = dau.py();
                thisInfo.negP[2] = dau.pz();
              }
            }
          }

          // For bookkeeping
          float ymc = 1e3;
          if (mcParticle.pdgCode() == 310)
            ymc = RecoDecay::y(std::array{thisInfo.posP[0] + thisInfo.negP[0], thisInfo.posP[1] + thisInfo.negP[0], thisInfo.posP[2] + thisInfo.negP[2]}, o2::constants::physics::MassKaonNeutral);
          else if (TMath::Abs(mcParticle.pdgCode()) == 3122)
            ymc = RecoDecay::y(std::array{thisInfo.posP[0] + thisInfo.negP[0], thisInfo.posP[1] + thisInfo.negP[0], thisInfo.posP[2] + thisInfo.negP[2]}, o2::constants::physics::MassLambda);

          if (mcParticle.pdgCode() == 310 && mcParticle.isPhysicalPrimary() && TMath::Abs(ymc) < rapidityWindow) {
            histos.fill(HIST("hStatisticsK0s"), 2.0f); // found
          }
          if (mcParticle.pdgCode() == 3122 && mcParticle.isPhysicalPrimary() && TMath::Abs(ymc) < rapidityWindow) {
            histos.fill(HIST("hStatisticsLambda"), 2.0f); // found
          }
          if (mcParticle.pdgCode() == -3122 && mcParticle.isPhysicalPrimary() && TMath::Abs(ymc) < rapidityWindow) {
            histos.fill(HIST("hStatisticsAlambda"), 2.0f); // found
          }

          // if I got here, it means this MC particle was not recoed and is of interest. Add it please
          mcV0infos.push_back(thisInfo);
        }
      }

      for (auto info : mcV0infos) {
        v0mccores(
          info.label, info.pdgCode,
          info.pdgCodeMother, info.pdgCodePositive, info.pdgCodeNegative,
          info.isPhysicalPrimary, info.xyz[0], info.xyz[1], info.xyz[2],
          info.posP[0], info.posP[1], info.posP[2],
          info.negP[0], info.negP[1], info.negP[2],
          info.momentum[0], info.momentum[1], info.momentum[2]);
        v0mccollref(info.mcCollision);
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
