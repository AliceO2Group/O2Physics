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

/// \file HFCandidateCreatorBs.cxx
/// \brief Reconstruction of Bs candidates
/// \note Adapted from HFCandidateCreatorXicc
///
/// \author Panos Christakoglou <panos.christakoglou@cern.ch>, Nikhef

#include "Framework/AnalysisTask.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::aod::hf_cand_prong3;
using namespace o2::aod::hf_cand_bs;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, true, {"Perform MC matching."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// Reconstruction of Λb candidates
struct HFCandidateCreatorBs {
  Produces<aod::HfCandBsBase> rowCandidateBase;

  Configurable<double> magneticField{"magneticField", 20., "magnetic field"};
  Configurable<bool> b_propdca{"b_propdca", true, "create tracks version propagated to PCA"};
  Configurable<double> d_maxr{"d_maxr", 200., "reject PCA's above this radius"};
  Configurable<double> d_maxdzini{"d_maxdzini", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> d_minparamchange{"d_minparamchange", 1.e-3, "stop iterations if largest change of any Bs is smaller than this"};
  Configurable<double> d_minrelchi2change{"d_minrelchi2change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<double> ptPionMin{"ptPionMin", 0.5, "minimum pion pT threshold (GeV/c)"};

  OutputObj<TH1F> hMassDsToKKPi{TH1F("hMassDsToKKPi", "D_{s}^{#plus} candidates;inv. mass (K^{#plus}K^{#minus} #pi^{#plus}) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtDs{TH1F("hPtDs", "D_{s}^{#plus} candidates;D_{s}^{#plus} candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hPtPion{TH1F("hPtPion", "#pi^{#minus} candidates;#pi^{#minus} candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hCPADs{TH1F("hCPADs", "D_{s}^{#plus} candidates;D_{s}^{#plus} cosine of pointing angle;entries", 110, -1.1, 1.1)};
  OutputObj<TH1F> hMassBsToDsPi{TH1F("hMassBsToDsPi", "2-prong candidates;inv. mass (B_{s}^{0} #rightarrow D_{s}^{#plus}#pi^{#minus} #rightarrow K^{#plus}K^{#minus}#pi^{#plus}#pi^{#minus}) (GeV/#it{c}^{2});entries", 500, 3., 8.)};
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx position (cm^{2});entries", 100, 0., 0.2)};

  double massPi = RecoDecay::getMassPDG(kPiMinus);
  double massDs = RecoDecay::getMassPDG(pdg::Code::kDs);
  double massDsPi = 0.;

  Configurable<int> d_selectionFlagDs{"d_selectionFlagDs", 1, "Selection Flag for Ds"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Filter filterSelectCandidates = (aod::hf_selcandidate_ds::isSelDsKKpi >= d_selectionFlagDs || aod::hf_selcandidate_ds::isSelDspiKK >= d_selectionFlagDs);

  void process(aod::Collision const& collision,
               soa::Filtered<soa::Join<
                 aod::HfCandProng3,
                 aod::HFSelDsCandidate>> const& dsCands,
               aod::BigTracks const& tracks)
  {
    // 2-prong vertex fitter
    o2::vertexing::DCAFitterN<2> df2;
    df2.setBz(magneticField);
    df2.setPropagateToPCA(b_propdca);
    df2.setMaxR(d_maxr);
    df2.setMaxDZIni(d_maxdzini);
    df2.setMinParamChange(d_minparamchange);
    df2.setMinRelChi2Change(d_minrelchi2change);
    df2.setUseAbsDCA(true);

    // 3-prong vertex fitter (to rebuild Ds vertex)
    o2::vertexing::DCAFitterN<3> df3;
    df3.setBz(magneticField);
    df3.setPropagateToPCA(b_propdca);
    df3.setMaxR(d_maxr);
    df3.setMaxDZIni(d_maxdzini);
    df3.setMinParamChange(d_minparamchange);
    df3.setMinRelChi2Change(d_minrelchi2change);
    df3.setUseAbsDCA(true);

    // loop over Ds candidates
    for (auto& dsCand : dsCands) {
      if (!(dsCand.hfflag() & 1 << o2::aod::hf_cand_prong3::DecayType::DsToPiKK)) {
        continue;
      }
      if (dsCand.isSelDsKKpi() >= d_selectionFlagDs) {
        hMassDsToKKPi->Fill(InvMassDsKKpi(dsCand), dsCand.pt());
      }
      if (dsCand.isSelDspiKK() >= d_selectionFlagDs) {
        hMassDsToKKPi->Fill(InvMassDspiKK(dsCand), dsCand.pt());
      }
      hPtDs->Fill(dsCand.pt());
      hCPADs->Fill(dsCand.cpa());

      auto track0 = dsCand.index0_as<aod::BigTracks>();
      auto track1 = dsCand.index1_as<aod::BigTracks>();
      auto track2 = dsCand.index2_as<aod::BigTracks>();
      auto trackParVar0 = getTrackParCov(track0);
      auto trackParVar1 = getTrackParCov(track1);
      auto trackParVar2 = getTrackParCov(track2);
      auto collision = track0.collision();

      // reconstruct the 3-prong secondary vertex
      if (df3.process(trackParVar0, trackParVar1, trackParVar2) == 0) {
        continue;
      }
      const auto& secondaryVertex = df3.getPCACandidate();
      trackParVar0.propagateTo(secondaryVertex[0], magneticField);
      trackParVar1.propagateTo(secondaryVertex[0], magneticField);
      trackParVar2.propagateTo(secondaryVertex[0], magneticField);

      array<float, 3> pvecpK = {track0.px() + track1.px(), track0.py() + track1.py(), track0.pz() + track1.pz()};
      array<float, 3> pvecDs = {pvecpK[0] + track2.px(), pvecpK[1] + track2.py(), pvecpK[2] + track2.pz()};
      auto trackpK = o2::dataformats::V0(df3.getPCACandidatePos(), pvecpK, df3.calcPCACovMatrixFlat(),
                                         trackParVar0, trackParVar1, {0, 0}, {0, 0});
      auto trackDs = o2::dataformats::V0(df3.getPCACandidatePos(), pvecDs, df3.calcPCACovMatrixFlat(),
                                         trackpK, trackParVar2, {0, 0}, {0, 0});

      int index0Ds = track0.globalIndex();
      int index1Ds = track1.globalIndex();
      int index2Ds = track2.globalIndex();
      // int charge = track0.sign() + track1.sign() + track2.sign();

      for (auto& trackPion : tracks) {
        if (trackPion.pt() < ptPionMin) {
          continue;
        }
        if (trackPion.sign() > 0) {
          continue;
        }
        if (trackPion.globalIndex() == index0Ds || trackPion.globalIndex() == index1Ds || trackPion.globalIndex() == index2Ds) {
          continue;
        }
        hPtPion->Fill(trackPion.pt());
        array<float, 3> pvecPion;
        auto trackParVarPi = getTrackParCov(trackPion);

        // reconstruct the 3-prong Ds vertex
        if (df2.process(trackDs, trackParVarPi) == 0) {
          continue;
        }

        // calculate relevant properties
        const auto& secondaryVertexBs = df2.getPCACandidate();
        auto chi2PCA = df2.getChi2AtPCACandidate();
        auto covMatrixPCA = df2.calcPCACovMatrixFlat();

        df2.propagateTracksToVertex();
        df2.getTrack(0).getPxPyPzGlo(pvecDs);
        df2.getTrack(1).getPxPyPzGlo(pvecPion);

        auto primaryVertex = getPrimaryVertex(collision);
        auto covMatrixPV = primaryVertex.getCov();
        o2::dataformats::DCA impactParameter0;
        o2::dataformats::DCA impactParameter1;
        trackDs.propagateToDCA(primaryVertex, magneticField, &impactParameter0);
        trackParVarPi.propagateToDCA(primaryVertex, magneticField, &impactParameter1);

        hCovSVXX->Fill(covMatrixPCA[0]);
        hCovPVXX->Fill(covMatrixPV[0]);

        // get uncertainty of the decay length
        double phi, theta;
        getPointDirection(array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertexBs, phi, theta);
        auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
        auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

        int hfFlag = 1 << hf_cand_bs::DecayType::BsToDsPi;

        // fill the candidate table for the Bs here:
        rowCandidateBase(collision.globalIndex(),
                         collision.posX(), collision.posY(), collision.posZ(),
                         secondaryVertexBs[0], secondaryVertexBs[1], secondaryVertexBs[2],
                         errorDecayLength, errorDecayLengthXY,
                         chi2PCA,
                         pvecDs[0], pvecDs[1], pvecDs[2],
                         pvecPion[0], pvecPion[1], pvecPion[2],
                         impactParameter0.getY(), impactParameter1.getY(),
                         std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()),
                         dsCand.globalIndex(), trackPion.globalIndex(),
                         hfFlag);

        // calculate invariant mass
        auto arrayMomenta = array{pvecDs, pvecPion};
        massDsPi = RecoDecay::M(std::move(arrayMomenta), array{massDs, massPi});
        if (dsCand.isSelDsKKpi() > 0) {
          hMassBsToDsPi->Fill(massDsPi);
        }
        if (dsCand.isSelDspiKK() > 0) {
          hMassBsToDsPi->Fill(massDsPi);
        }
      } // pi- loop
    }   // Ds loop
  }     // process
};      // struct

/// Extends the base table with expression columns.
struct HFCandidateCreatorBsExpressions {
  Spawns<aod::HfCandBsExt> rowCandidateBs;
  void init(InitContext const&) {}
};

/// Performs MC matching.
struct HFCandidateCreatorBsMC {
  Produces<aod::HfCandBsMCRec> rowMCMatchRec;
  Produces<aod::HfCandBsMCGen> rowMCMatchGen;

  void process(aod::HfCandBs const& candidates,
               aod::HfCandProng3,
               aod::BigTracksMC const& tracks,
               aod::McParticles const& particlesMC)
  {
    int indexRec = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t debug = 0;

    // Match reconstructed candidates.
    for (auto& candidate : candidates) {
      // Printf("New rec. candidate");
      flag = 0;
      origin = 0;
      debug = 0;
      auto dsCand = candidate.index0();
      auto arrayDaughters = array{dsCand.index0_as<aod::BigTracksMC>(),
                                  dsCand.index1_as<aod::BigTracksMC>(),
                                  dsCand.index2_as<aod::BigTracksMC>(),
                                  candidate.index1_as<aod::BigTracksMC>()};
      auto arrayDaughtersDs = array{dsCand.index0_as<aod::BigTracksMC>(),
                                    dsCand.index1_as<aod::BigTracksMC>(),
                                    dsCand.index2_as<aod::BigTracksMC>()};
      // Bs → Ds+ π-
      // Printf("Checking Bs → Ds+ π-");
      indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kBs, array{+kKPlus, -kKPlus, +kPiPlus, -kPiPlus}, true, &sign, 2);
      if (indexRec > -1) {
        // Bs → Ds+ π-
        // Printf("Checking Bs → Ds+ π-");
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughtersDs, pdg::Code::kDs, array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 1);
        if (indexRec > -1) {
          flag = 1 << hf_cand_bs::DecayType::BsToDsPi;
        } else {
          debug = 1;
          LOGF(info, "WARNING: Λb in decays in the expected final state but the condition on the intermediate state is not fulfilled");
        }
      }
      rowMCMatchRec(flag, origin, debug);
    }

    // Match generated particles.
    for (auto& particle : particlesMC) {
      // Printf("New gen. candidate");
      flag = 0;
      origin = 0;
      // Bs → Ds+ π-
      if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdg::Code::kBs, array{int(pdg::Code::kDs), -kPiPlus}, true)) {
        // Match Ds+ -> φπ -> K+K-π
        auto DsCandMC = particlesMC.iteratorAt(particle.daughter0Id());

        // Printf("Checking Ds+ -> φπ -> K+K-π");
        if (RecoDecay::isMatchedMCGen(particlesMC, DsCandMC, int(pdg::Code::kDs), array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign)) {
          flag = sign * (1 << hf_cand_bs::DecayType::BsToDsPi);
        }
      }
      rowMCMatchGen(flag, origin);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HFCandidateCreatorBs>(cfgc, TaskName{"hf-cand-creator-bs"}),
    adaptAnalysisTask<HFCandidateCreatorBsExpressions>(cfgc, TaskName{"hf-cand-creator-bs-expressions"})};
  const bool doMC = cfgc.options().get<bool>("doMC");
  if (doMC) {
    workflow.push_back(adaptAnalysisTask<HFCandidateCreatorBsMC>(cfgc, TaskName{"hf-cand-creator-bs-mc"}));
  }
  return workflow;
}
