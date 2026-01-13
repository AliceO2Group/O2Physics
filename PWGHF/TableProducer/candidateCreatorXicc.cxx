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

/// \file candidateCreatorXicc.cxx
/// \brief Reconstruction of Xiccplusplus candidates
/// \note Extended from candidateCreator2Prong.cxx, candidateCreator3Prong.cxx, candidateCreatorX.cxx
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Luigi Dello Stritto <luigi.dello.stritto@cern.ch >, SALERNO
/// \author Mattia Faggin <mattia.faggin@cern.ch>, University and INFN PADOVA

#include "PWGHF/ALICE3/Core/DecayChannelsLegacy.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"

#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/Variant.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/V0.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec const optionDoMC{"doMC", VariantType::Bool, true, {"Perform MC matching."}};
  workflowOptions.push_back(optionDoMC);
}

#include <Framework/runDataProcessing.h>

/// Reconstruction of xicc candidates
struct HfCandidateCreatorXicc {
  Produces<aod::HfCandXiccBase> rowCandidateBase;

  // vertexing
  Configurable<double> bz{"bz", 5., "magnetic field"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  // selection
  Configurable<int> selectionFlagXic{"selectionFlagXic", 1, "Selection Flag for Xic"};
  Configurable<double> cutPtPionMin{"cutPtPionMin", 1., "min. pt pion track"};

  o2::vertexing::DCAFitterN<3> df3; // 3-prong vertex fitter to rebuild the Xic vertex
  o2::vertexing::DCAFitterN<2> df2; // 2-prong vertex fitter to build the Xicc vertex

  Filter filterSelectCandidates = (aod::hf_sel_candidate_xic::isSelXicToPKPi >= selectionFlagXic || aod::hf_sel_candidate_xic::isSelXicToPiKP >= selectionFlagXic);

  OutputObj<TH1F> hMassXic{TH1F("hMassXic", "xic candidates;inv. mass (#pi K #pi) (GeV/#it{c}^{2});entries", 500, 1.6, 2.6)};
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "3-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "3-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};

  void init(InitContext const&)
  {
    df3.setBz(bz);
    df3.setPropagateToPCA(propagateToPCA);
    df3.setMaxR(maxR);
    df3.setMaxDZIni(maxDZIni);
    df3.setMinParamChange(minParamChange);
    df3.setMinRelChi2Change(minRelChi2Change);
    df3.setUseAbsDCA(useAbsDCA);
    df3.setWeightedFinalPCA(useWeightedFinalPCA);

    df2.setBz(bz);
    df2.setPropagateToPCA(propagateToPCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(useAbsDCA);
    df2.setWeightedFinalPCA(useWeightedFinalPCA);
  }

  void process(aod::Collision const&,
               soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelXicToPKPi>> const& xicCands,
               aod::TracksWCov const& tracks)
  {
    for (const auto& xicCand : xicCands) {
      if ((xicCand.hfflag() & 1 << o2::aod::hf_cand_3prong::DecayType::XicToPKPi) == 0) {
        continue;
      }
      if (xicCand.isSelXicToPKPi() >= selectionFlagXic) {
        hMassXic->Fill(HfHelper::invMassXicToPKPi(xicCand), xicCand.pt());
      }
      if (xicCand.isSelXicToPiKP() >= selectionFlagXic) {
        hMassXic->Fill(HfHelper::invMassXicToPiKP(xicCand), xicCand.pt());
      }
      auto track0 = xicCand.prong0_as<aod::TracksWCov>();
      auto track1 = xicCand.prong1_as<aod::TracksWCov>();
      auto track2 = xicCand.prong2_as<aod::TracksWCov>();
      auto trackParVar0 = getTrackParCov(track0);
      auto trackParVar1 = getTrackParCov(track1);
      auto trackParVar2 = getTrackParCov(track2);
      auto collision = track0.collision();

      // reconstruct the 3-prong secondary vertex
      if (df3.process(trackParVar0, trackParVar1, trackParVar2) == 0) {
        continue;
      }
      const auto& secondaryVertex = df3.getPCACandidate();
      trackParVar0.propagateTo(secondaryVertex[0], bz);
      trackParVar1.propagateTo(secondaryVertex[0], bz);
      trackParVar2.propagateTo(secondaryVertex[0], bz);

      std::array<float, 3> const pvecpK = RecoDecay::pVec(track0.pVector(), track1.pVector());
      std::array<float, 3> pvecxic = RecoDecay::pVec(pvecpK, track2.pVector());
      auto trackpK = o2::dataformats::V0(df3.getPCACandidatePos(), pvecpK, df3.calcPCACovMatrixFlat(), trackParVar0, trackParVar1);
      auto trackxic = o2::dataformats::V0(df3.getPCACandidatePos(), pvecxic, df3.calcPCACovMatrixFlat(), trackpK, trackParVar2);

      int const index0Xic = track0.globalIndex();
      int const index1Xic = track1.globalIndex();
      int const index2Xic = track2.globalIndex();
      int const charge = track0.sign() + track1.sign() + track2.sign();

      for (const auto& trackpion : tracks) {
        if (trackpion.pt() < cutPtPionMin) {
          continue;
        }
        if (trackpion.sign() * charge < 0) {
          continue;
        }
        if (trackpion.globalIndex() == index0Xic || trackpion.globalIndex() == index1Xic || trackpion.globalIndex() == index2Xic) {
          continue;
        }
        std::array<float, 3> pvecpion{};
        auto trackParVarPi = getTrackParCov(trackpion);

        // reconstruct the 3-prong X vertex
        if (df2.process(trackxic, trackParVarPi) == 0) {
          continue;
        }

        // calculate relevant properties
        const auto& secondaryVertexXicc = df2.getPCACandidate();
        auto chi2PCA = df2.getChi2AtPCACandidate();
        auto covMatrixPCA = df2.calcPCACovMatrixFlat();

        df2.propagateTracksToVertex();
        df2.getTrack(0).getPxPyPzGlo(pvecxic);
        df2.getTrack(1).getPxPyPzGlo(pvecpion);

        auto primaryVertex = getPrimaryVertex(collision);
        auto covMatrixPV = primaryVertex.getCov();
        o2::dataformats::DCA impactParameter0;
        o2::dataformats::DCA impactParameter1;
        trackxic.propagateToDCA(primaryVertex, bz, &impactParameter0);
        trackParVarPi.propagateToDCA(primaryVertex, bz, &impactParameter1);

        // get uncertainty of the decay length
        double phi{}, theta{};
        getPointDirection(std::array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertexXicc, phi, theta);
        auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
        auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

        int hfFlag = 1 << aod::hf_cand_xicc::DecayType::XiccToXicPi;

        rowCandidateBase(collision.globalIndex(),
                         collision.posX(), collision.posY(), collision.posZ(),
                         secondaryVertexXicc[0], secondaryVertexXicc[1], secondaryVertexXicc[2],
                         errorDecayLength, errorDecayLengthXY,
                         chi2PCA,
                         pvecxic[0], pvecxic[1], pvecxic[2],
                         pvecpion[0], pvecpion[1], pvecpion[2],
                         impactParameter0.getY(), impactParameter1.getY(),
                         std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()),
                         xicCand.globalIndex(), trackpion.globalIndex(),
                         hfFlag);
      } // if on selected Xicc
    } // loop over candidates
  } // end of process
}; // end of struct

/// Extends the base table with expression columns.
struct HfCandidateCreatorXiccExpressions {
  Spawns<aod::HfCandXiccExt> rowCandidateXicc;

  void init(InitContext const&) {}
};

/// Performs MC matching.
struct HfCandidateCreatorXiccMc {
  Produces<aod::HfCandXiccMcRec> rowMcMatchRec;
  Produces<aod::HfCandXiccMcGen> rowMcMatchGen;

  void process(aod::HfCandXicc const& candidates,
               aod::HfCand3Prong const&,
               aod::TracksWMc const&,
               aod::McParticles const& mcParticles)
  {
    int8_t sign{};
    int8_t flag{};
    int8_t origin{};

    // Match reconstructed candidates.
    for (const auto& candidate : candidates) {
      int8_t debug = 0;
      flag = 0;
      origin = 0;
      auto xicCand = candidate.prong0();
      auto arrayDaughters = std::array{xicCand.prong0_as<aod::TracksWMc>(),
                                       xicCand.prong1_as<aod::TracksWMc>(),
                                       xicCand.prong2_as<aod::TracksWMc>(),
                                       candidate.prong1_as<aod::TracksWMc>()};
      auto arrayDaughtersXic = std::array{xicCand.prong0_as<aod::TracksWMc>(),
                                          xicCand.prong1_as<aod::TracksWMc>(),
                                          xicCand.prong2_as<aod::TracksWMc>()};
      // Ξcc±± → p± K∓ π± π±
      auto indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kXiCCPlusPlus, std::array{+kProton, -kKPlus, +kPiPlus, +kPiPlus}, true, &sign, 2);
      if (indexRec > -1) {
        // Ξc± → p± K∓ π±
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersXic, Pdg::kXiCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 1);
        if (indexRec > -1) {
          flag = 1 << aod::hf_cand_xicc::DecayType::XiccToXicPi;
        } else {
          debug = 1;
          LOGF(info, "WARNING: Ξc±± in decays in the expected final state but the condition on the intermediate state is not fulfilled");
        }
      }
      rowMcMatchRec(flag, origin, debug);
    }

    // Match generated particles.
    for (const auto& particle : mcParticles) {
      flag = 0;
      origin = 0;
      // Ξcc±± → Ξc± + π±
      if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kXiCCPlusPlus, std::array{static_cast<int>(Pdg::kXiCPlus), +kPiPlus}, true)) {
        // Ξc± → p± K∓ π±
        auto candXicMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
        if (RecoDecay::isMatchedMCGen(mcParticles, candXicMC, static_cast<int>(Pdg::kXiCPlus), std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign)) {
          flag = sign * (1 << aod::hf_cand_xicc::DecayType::XiccToXicPi);
        }
      }
      rowMcMatchGen(flag, origin);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HfCandidateCreatorXicc>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorXiccExpressions>(cfgc)};
  const bool doMC = cfgc.options().get<bool>("doMC");
  if (doMC) {
    workflow.push_back(adaptAnalysisTask<HfCandidateCreatorXiccMc>(cfgc));
  }
  return workflow;
}
