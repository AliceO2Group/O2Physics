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

/// \file phiphimctrigger.cxx
/// \brief Code for triggered doublephi resonance MC
/// \author Sawan

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h>
#include <TPDGCode.h>

#include <cmath>   // std::abs
#include <cstddef> // size_t / std::size_t
#include <vector>

using namespace o2;
using namespace o2::framework;

struct PhiCandidateTrigger {
  ROOT::Math::PxPyPzMVector p4;
  int t1Index = -1;
  int t2Index = -1;
  int mcPhiIndex = -1;

  // MC Truth attributes
  bool isValidMcPhi = false;
};

struct PhiCandidateGen {
  ROOT::Math::PxPyPzMVector p4;
};

struct PhiPhiMCTrigger {

  // Configurables
  struct : ConfigurableGroup {
    Configurable<bool> isSel8{"isSel8", true, "Event selection sel8"};
    Configurable<float> cutDCAxy{"cutDCAxy", 0.1f, "DCAxy max cut"};
    Configurable<float> cutDCAz{"cutDCAz", 0.1f, "DCAz max cut"};
    Configurable<float> cutNsigmaTPC{"cutNsigmaTPC", 3.0f, "TPC Nsigma cut for Kaon"};
    Configurable<bool> isApplyPIDSelection{"isApplyPIDSelection", true, "Apply PID selection for Kaons"};
    Configurable<float> trackEtaCut{"trackEtaCut", 0.8f, "Track Eta Cut"};
    Configurable<float> phiMassCutLow{"phiMassCutLow", 0.9f, "Lower Cut for Phi Mass"};
    Configurable<float> phiMassCutHigh{"phiMassCutHigh", 1.1f, "Upper Cut for Phi Mass"};
  } config;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  ConfigurableAxis phiMassAxis{"phiMassAxis", {200, 0.98, 1.22}, "Phi Mass Axis"};
  ConfigurableAxis phiPtAxis{"phiPtAxis", {200, 0.0, 100.0}, "Phi Pt Axis"};
  ConfigurableAxis phiPhiMassAxis{"phiPhiMassAxis", {200, 2.0, 4.0}, "Phi-Phi Mass Axis"};
  ConfigurableAxis phiPhiPtAxis{"phiPhiPtAxis", {200, 0.0, 50.0}, "Phi-Phi Pt Axis"};

  enum RecoSelectionStage {
    kLoop1Track = 0,
    kLoop1DCA,
    kLoop1PID,
    kLoop1MCLabel,
    kLoop1Eta,
    kLoop2Track,
    kLoop2Index,
    kLoop2Eta,
    kLoop2DCA,
    kLoop2PID,
    kLoop2MCLabel,
    kLoop2Charge,
    kLoop2TrackSelected,
    kLoop2KaonDaughters,
    kCommonMotherLoop,
    kCommonMother,
    kPhiMother,
    kPhiMassCut,
    kLoop2PhiPairs,
    kUniqueTracks,
    kDifferentPhiMothers,
    kValidPhiPhiPair
  };

  void init(InitContext const&)
  {

    // Event histogram
    registry.add("h1_TrackCount", "Track Count;Track Count;Counts", kTH1F, {{30, -0.5, 29.5}});
    registry.add("h1_EventCount", "Event Count;Event Count;Counts", kTH1F, {{15, -0.5, 14.5}});
    // Generated MC Histograms
    registry.add("Generated/h1_GenEventAndParticleCount", "Gen Level Summary;Category;Counts", kTH1F, {{4, -0.5, 3.5}});
    registry.add("Generated/h2_InvMassPhiPhi", "Gen Phi-Phi: Mass vs pT;Mass (GeV/#it{c}^{2});p_{T} (GeV/#it{c})", kTHnSparseF, {{phiPhiMassAxis}, {phiPhiPtAxis}});
    registry.add("Generated/h1_PhiPhi_Phi", "Gen Phi-Phi: #phi Angle;#phi;Counts", kTH1F, {{100, -o2::constants::math::PI, o2::constants::math::PI}});
    registry.add("Generated/h1_PhiPhi_Eta", "Gen Phi-Phi: #eta;#eta;Counts", kTH1F, {{100, -1.0, 1.0}});
    registry.add("Generated/h1_PhiPhi_Y", "Gen Phi-Phi: Rapidity y;y;Counts", kTH1F, {{100, -1.0, 1.0}});
    registry.add("Generated/h1_InvMassPhi", "Gen Phi: Mass;Mass (GeV/#it{c}^{2});Counts", kTH1F, {{phiMassAxis}});
    registry.add("Generated/h1_PhiPt", "Gen Phi: pT;p_{T} (GeV/#it{c});Counts", kTH1F, {{phiPtAxis}});
    registry.add("Generated/h1_PhiEta", "Gen Phi: #eta;#eta;Counts", kTH1F, {{100, -1.0, 1.0}});
    registry.add("Generated/h1_PhiY", "Gen Phi: Rapidity y;y;Counts", kTH1F, {{100, -1.0, 1.0}});

    // Reconstructed MC Histograms - 2D Mass vs pT of Mother Phi-Phi
    registry.add("Reconstructed/h2_PhiPhi_InvMass", "Rec Phi-Phi All Pairs: Mass vs pT;Mass (GeV/#it{c}^{2});p_{T} (GeV/#it{c})", kTHnSparseF, {{phiPhiMassAxis}, {phiPhiPtAxis}});
    registry.add("Reconstructed/h1_Phi1_Mass", "Rec Phi Mass: Mass;Mass (GeV/#it{c}^{2});Counts", kTH1F, {{phiMassAxis}});
    registry.add("Reconstructed/h1_Phi1_Pt", "Rec Phi pT: pT;p_{T} (GeV/#it{c});Counts", kTH1F, {{phiPtAxis}});
    registry.add("Reconstructed/h1_Phi1_Eta", "Rec Phi Eta: Eta;#eta;Counts", kTH1F, {{100, -1.0, 1.0}});
    registry.add("Reconstructed/h1_Phi1_Y", "Rec Phi Rapidity y;y;Counts", kTH1F, {{100, -1.0, 1.0}});
    registry.add("Reconstructed/h1_PhiPhi_Y", "Rec Phi Rapidity y;y;Counts", kTH1F, {{100, -1.0, 1.0}});
  }

  // PID Selection function
  bool checkPID(float ptcand, float nsigmaTPC, float nsigmaTOF, bool TOFHit)
  {
    constexpr float Radius2Max = 2.5f * 2.5f;
    const bool hasTOF = (TOFHit == 1);

    if (!hasTOF) {
      if (ptcand < 0.5f) {
        return std::abs(nsigmaTPC) < config.cutNsigmaTPC;
      } else {
        return nsigmaTPC > -2.0f && nsigmaTPC < config.cutNsigmaTPC;
      }
    } else {
      const float radius2 = nsigmaTPC * nsigmaTPC + nsigmaTOF * nsigmaTOF;
      if (radius2 >= Radius2Max) {
        return false;
      }
      if (ptcand < 2.0f) {
        return true;
      }

      float slope = 0.0f;
      float intercept = 0.0f;

      if (ptcand < 2.5f) {
        slope = 1.29f;
        intercept = 3.34f;
      } else if (ptcand < 3.0f) {
        slope = 0.94f;
        intercept = 2.45f;
      } else if (ptcand < 4.0f) {
        slope = 0.62f;
        intercept = 1.55f;
      } else if (ptcand < 5.0f) {
        slope = 0.32f;
        intercept = 1.45f;
      } else if (ptcand < 7.0f) {
        slope = 0.22f;
        intercept = 1.35f;
      } else {
        slope = 0.10f;
        intercept = 1.25f;
      }

      return nsigmaTPC < slope * nsigmaTOF + intercept;
    }
  }

  // Define data table bindings
  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels>;
  using TrackCandidatesMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::McTrackLabels, aod::pidTOFbeta, aod::TrackSelectionExtension>;

  ROOT::Math::PxPyPzEVector vecGen, vecDau1, vecDau2;
  ROOT::Math::PxPyPzMVector daughter1Gen, daughter2Gen, motherGen, daughter1Rec, daughter2Rec, motherRec;

  // Process Generated MC
  void processGen(aod::McCollision const&, aod::McParticles const& mcParticles)
  {
    // Fill total generated events counter
    registry.fill(HIST("Generated/h1_GenEventAndParticleCount"), 0.0);

    std::vector<PhiCandidateGen> phiCandidatesGen;
    for (const auto& mcParticle : mcParticles) {

      if (mcParticle.pdgCode() == o2::constants::physics::kPhi) {

        // Fill total phi produced counter
        registry.fill(HIST("Generated/h1_GenEventAndParticleCount"), 1.0);

        auto phiDaughters = mcParticle.daughters_as<aod::McParticles>();
        if (phiDaughters.size() != 2) {
          continue;
        }

        bool passPosKaon = false;
        bool passNegKaon = false;
        for (const auto& kCurrentDaughter : phiDaughters) {

          if (kCurrentDaughter.pdgCode() == PDG_t::kKPlus) {
            passPosKaon = true;
            daughter1Gen = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), o2::constants::physics::MassKPlus);

            // Fill total phi decay to K+
            registry.fill(HIST("Generated/h1_GenEventAndParticleCount"), 2.0);

          } else if (kCurrentDaughter.pdgCode() == PDG_t::kKMinus) {
            passNegKaon = true;
            daughter2Gen = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), o2::constants::physics::MassKPlus);

            // Fill total phi decay to K-
            registry.fill(HIST("Generated/h1_GenEventAndParticleCount"), 3.0);
          }
        }
        if (passPosKaon && passNegKaon) {
          motherGen = daughter1Gen + daughter2Gen;
          PhiCandidateGen phiCand;
          phiCand.p4 = motherGen;
          phiCandidatesGen.push_back(phiCand);

          registry.fill(HIST("Generated/h1_InvMassPhi"), motherGen.M());
          registry.fill(HIST("Generated/h1_PhiPt"), motherGen.Pt());
          registry.fill(HIST("Generated/h1_PhiEta"), motherGen.Eta());
          registry.fill(HIST("Generated/h1_PhiY"), motherGen.Rapidity());
        }
      }
    }

    for (size_t i = 0; i < phiCandidatesGen.size(); ++i) {
      for (size_t j = i + 1; j < phiCandidatesGen.size(); ++j) {
        const auto& phi1 = phiCandidatesGen[i];
        const auto& phi2 = phiCandidatesGen[j];

        ROOT::Math::PxPyPzMVector motherPhiPhi = phi1.p4 + phi2.p4;

        registry.fill(HIST("Generated/h2_InvMassPhiPhi"), motherPhiPhi.M(), motherPhiPhi.Pt());
        registry.fill(HIST("Generated/h1_PhiPhi_Phi"), motherPhiPhi.Phi());
        registry.fill(HIST("Generated/h1_PhiPhi_Eta"), motherPhiPhi.Eta());
        registry.fill(HIST("Generated/h1_PhiPhi_Y"), motherPhiPhi.Rapidity());
      }
    }
  }
  PROCESS_SWITCH(PhiPhiMCTrigger, processGen, "Process Generated MC", true);

  // Process Reconstructed MC
  void processRec(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&)
  {

    registry.fill(HIST("h1_EventCount"), 0.0);
    if (collision.sel8()) {
      registry.fill(HIST("h1_EventCount"), 1.0);
    }

    if (config.isSel8 && !collision.sel8()) {
      return;
    }

    std::vector<PhiCandidateTrigger> phiCandidates;

    // STEP 1: Reconstruct individual Phi candidates (Kaon track pairs)
    for (const auto& t1 : tracks) {
      registry.fill(HIST("h1_TrackCount"), kLoop1Track);

      if (std::abs(t1.dcaXY()) > config.cutDCAxy || std::abs(t1.dcaZ()) > config.cutDCAz) {
        continue;
      }
      registry.fill(HIST("h1_TrackCount"), kLoop1DCA);

      if (config.isApplyPIDSelection && !checkPID(t1.pt(), t1.tpcNSigmaKa(), t1.tofNSigmaKa(), t1.hasTOF())) {
        continue;
      }
      registry.fill(HIST("h1_TrackCount"), kLoop1PID);

      if (!t1.has_mcParticle()) {
        continue;
      }
      registry.fill(HIST("h1_TrackCount"), kLoop1MCLabel);

      if (std::abs(t1.eta()) > config.trackEtaCut) {
        continue;
      }
      registry.fill(HIST("h1_TrackCount"), kLoop1Eta);

      for (const auto& t2 : tracks) {
        registry.fill(HIST("h1_TrackCount"), kLoop2Track);

        if (t2.index() <= t1.index()) {
          continue;
        }
        registry.fill(HIST("h1_TrackCount"), kLoop2Index);

        if (std::abs(t2.eta()) > config.trackEtaCut) {
          continue;
        }
        registry.fill(HIST("h1_TrackCount"), kLoop2Eta);

        if (std::abs(t2.dcaXY()) > config.cutDCAxy || std::abs(t2.dcaZ()) > config.cutDCAz) {
          continue;
        }
        registry.fill(HIST("h1_TrackCount"), kLoop2DCA);

        if (config.isApplyPIDSelection && !checkPID(t2.pt(), t2.tpcNSigmaKa(), t2.tofNSigmaKa(), t2.hasTOF())) {
          continue;
        }
        registry.fill(HIST("h1_TrackCount"), kLoop2PID);

        if (!t2.has_mcParticle()) {
          continue;
        }
        registry.fill(HIST("h1_TrackCount"), kLoop2MCLabel);

        if (t1.sign() * t2.sign() > 0) {
          continue;
        }
        registry.fill(HIST("h1_TrackCount"), kLoop2Charge);

        registry.fill(HIST("h1_TrackCount"), kLoop2TrackSelected);

        auto mcKaTrk1 = t1.mcParticle();
        auto mcKaTrk2 = t2.mcParticle();

        if (std::abs(mcKaTrk1.pdgCode()) != PDG_t::kKPlus ||
            std::abs(mcKaTrk2.pdgCode()) != PDG_t::kKPlus) {
          continue;
        }
        registry.fill(HIST("h1_TrackCount"), kLoop2KaonDaughters);

        bool foundMother = false;

        // Search for common mother
        for (const auto& KaMotherPhi1 : mcKaTrk1.template mothers_as<aod::McParticles>()) {

          for (const auto& KaMotherPhi2 : mcKaTrk2.template mothers_as<aod::McParticles>()) {

            registry.fill(HIST("h1_TrackCount"), kCommonMotherLoop);

            if (KaMotherPhi1.globalIndex() != KaMotherPhi2.globalIndex())
              continue;

            registry.fill(HIST("h1_TrackCount"), kCommonMother);

            int motherPdg = KaMotherPhi1.pdgCode();

            if (motherPdg != o2::constants::physics::kPhi)
              continue;

            registry.fill(HIST("h1_TrackCount"), kPhiMother);

            daughter1Rec = ROOT::Math::PxPyPzMVector(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassKPlus);
            daughter2Rec = ROOT::Math::PxPyPzMVector(t2.px(), t2.py(), t2.pz(), o2::constants::physics::MassKPlus);
            motherRec = daughter1Rec + daughter2Rec;

            if (motherRec.M() < config.phiMassCutLow || motherRec.M() > config.phiMassCutHigh) {
              continue;
            }
            registry.fill(HIST("h1_TrackCount"), kPhiMassCut);

            PhiCandidateTrigger phiCand;

            phiCand.p4 = motherRec;
            phiCand.t1Index = t1.globalIndex();
            phiCand.t2Index = t2.globalIndex();
            phiCand.isValidMcPhi = true;
            phiCand.mcPhiIndex = KaMotherPhi1.globalIndex();

            phiCandidates.push_back(phiCand);
            foundMother = true;
            break;
          }
          if (foundMother) {
            break;
          }
        }
      }
    }

    // STEP 2: Fill Single Reconstructed Phi Histograms (Loop ONCE per candidate)
    for (const auto& phi : phiCandidates) {
      // All Reconstructed Phi Candidates
      if (phi.isValidMcPhi) {

        registry.fill(HIST("Reconstructed/h1_Phi1_Mass"), phi.p4.M());
        registry.fill(HIST("Reconstructed/h1_Phi1_Pt"), phi.p4.Pt());
        registry.fill(HIST("Reconstructed/h1_Phi1_Eta"), phi.p4.Eta());
        registry.fill(HIST("Reconstructed/h1_Phi1_Y"), phi.p4.Rapidity());
      }
    }

    // STEP 3: Reconstruct Phi-Phi Mother Candidates (Double Loop for Pairs Only)
    for (size_t i = 0; i < phiCandidates.size(); ++i) {
      const auto& phi1 = phiCandidates[i];

      for (size_t j = i + 1; j < phiCandidates.size(); ++j) {
        const auto& phi2 = phiCandidates[j];

        registry.fill(HIST("h1_TrackCount"), kLoop2PhiPairs);

        // Ensure all 4 kaon tracks are unique
        if (phi1.t1Index == phi2.t1Index || phi1.t1Index == phi2.t2Index ||
            phi1.t2Index == phi2.t1Index || phi1.t2Index == phi2.t2Index) {
          continue;
        }
        registry.fill(HIST("h1_TrackCount"), kUniqueTracks);

        if (phi1.mcPhiIndex == phi2.mcPhiIndex)
          continue;
        registry.fill(HIST("h1_TrackCount"), kDifferentPhiMothers);

        if (!phi1.isValidMcPhi || !phi2.isValidMcPhi) {
          continue;
        }
        registry.fill(HIST("h1_TrackCount"), kValidPhiPhiPair);

        // Combined 4-momentum of the Mother Phi-Phi System
        ROOT::Math::PxPyPzMVector phiPhiMother = phi1.p4 + phi2.p4;

        registry.fill(HIST("Reconstructed/h2_PhiPhi_InvMass"), phiPhiMother.M(), phiPhiMother.Pt());
        registry.fill(HIST("Reconstructed/h1_PhiPhi_Y"), phiPhiMother.Rapidity());
      }
    }
  }
  PROCESS_SWITCH(PhiPhiMCTrigger, processRec, "Process Reconstructed MC", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PhiPhiMCTrigger>(cfgc)};
}
