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

/// \file candidateCreatorLb.cxx
/// \brief Reconstruction of Lb candidates
/// \note Adapted from candidateCreatorXicc.cxx
///
/// \author Panos Christakoglou <panos.christakoglou@cern.ch>, Nikhef

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"

#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/V0.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_trkcandsel;
using namespace o2::hf_decay::hf_cand_beauty;

/// Reconstruction of Λb candidates
struct HfCandidateCreatorLb {
  Produces<aod::HfCandLbBase> rowCandidateBase;
  Produces<aod::HfCandLbProngs> rowCandidateProngs;
  // vertexing
  Configurable<float> bz{"bz", 20., "magnetic field"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<float> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<float> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<float> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any Lb is smaller than this"};
  Configurable<float> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  // selection
  Configurable<float> ptPionMin{"ptPionMin", 0.5, "minimum pion pT threshold (GeV/c)"};
  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<float> yCandMax{"yCandMax", -1., "max. cand. rapidity"};

  o2::vertexing::DCAFitterN<2> df2; // 2-prong vertex fitter
  o2::vertexing::DCAFitterN<3> df3; // 3-prong vertex fitter (to rebuild Lc vertex)

  double massLcPi{0.};

  Filter filterSelectCandidates = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLc || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLc);

  std::shared_ptr<TH1> hCandidatesLc, hCandidatesLb;
  HistogramRegistry registry{"registry"};

  OutputObj<TH1F> hMassLcToPKPi{TH1F("hMassLcToPKPi", "#Lambda_{c}^{#plus} candidates;inv. mass (pK^{#minus} #pi^{#plus}) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtLc{TH1F("hPtLc", "#Lambda_{c}^{#plus} candidates;#Lambda_{c}^{#plus} candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hPtPion{TH1F("hPtPion", "#pi^{#minus} candidates;#pi^{#minus} candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hCPALc{TH1F("hCPALc", "#Lambda_{c}^{#plus} candidates;#Lambda_{c}^{#plus} cosine of pointing angle;entries", 110, -1.1, 1.1)};
  OutputObj<TH1F> hMassLbToLcPi{TH1F("hMassLbToLcPi", "2-prong candidates;inv. mass (#Lambda_{b}^{0} #rightarrow #Lambda_{c}^{#plus}#pi^{#minus} #rightarrow pK^{#minus}#pi^{#plus}#pi^{#minus}) (GeV/#it{c}^{2});entries", 500, 3., 8.)};
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};

  void init(InitContext const&)
  {
    df2.setBz(bz);
    df2.setPropagateToPCA(propagateToPCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(useAbsDCA);
    df2.setWeightedFinalPCA(useWeightedFinalPCA);

    df3.setBz(bz);
    df3.setPropagateToPCA(propagateToPCA);
    df3.setMaxR(maxR);
    df3.setMaxDZIni(maxDZIni);
    df3.setMinParamChange(minParamChange);
    df3.setMinRelChi2Change(minRelChi2Change);
    df3.setUseAbsDCA(useAbsDCA);
    df3.setWeightedFinalPCA(useWeightedFinalPCA);

    /// candidate monitoring
    hCandidatesLc = registry.add<TH1>("hCandidatesLc", "Lc candidate counter", {HistType::kTH1D, {axisCands}});
    hCandidatesLb = registry.add<TH1>("hCandidatesLb", "B candidate counter", {HistType::kTH1D, {axisCands}});
    setLabelHistoCands(hCandidatesLc);
    setLabelHistoCands(hCandidatesLb);
  }

  void process(aod::Collision const&,
               soa::Filtered<soa::Join<
                 aod::HfCand3Prong,
                 aod::HfSelLc>> const& lcCands,
               aod::TracksWCov const& tracks)
  {
    // loop over Lc candidates
    for (const auto& lcCand : lcCands) {
      if ((lcCand.hfflag() & 1 << o2::aod::hf_cand_3prong::DecayType::LcToPKPi) == 0) {
        continue;
      }
      if (lcCand.isSelLcToPKPi() >= selectionFlagLc) {
        hMassLcToPKPi->Fill(HfHelper::invMassLcToPKPi(lcCand), lcCand.pt());
      }
      if (lcCand.isSelLcToPiKP() >= selectionFlagLc) {
        hMassLcToPKPi->Fill(HfHelper::invMassLcToPiKP(lcCand), lcCand.pt());
      }
      hPtLc->Fill(lcCand.pt());
      hCPALc->Fill(lcCand.cpa());

      auto track0 = lcCand.prong0_as<aod::TracksWCov>();
      auto track1 = lcCand.prong1_as<aod::TracksWCov>();
      auto track2 = lcCand.prong2_as<aod::TracksWCov>();
      auto trackParVar0 = getTrackParCov(track0);
      auto trackParVar1 = getTrackParCov(track1);
      auto trackParVar2 = getTrackParCov(track2);
      auto collision = lcCand.collision();

      // reconstruct the 3-prong secondary vertex
      hCandidatesLc->Fill(SVFitting::BeforeFit);
      try {
        if (df3.process(trackParVar0, trackParVar1, trackParVar2) == 0) {
          continue;
        }
      } catch (const std::runtime_error& error) {
        LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
        hCandidatesLc->Fill(SVFitting::Fail);
        continue;
      }
      hCandidatesLc->Fill(SVFitting::FitOk);

      const auto& secondaryVertex = df3.getPCACandidate();
      trackParVar0.propagateTo(secondaryVertex[0], bz);
      trackParVar1.propagateTo(secondaryVertex[0], bz);
      trackParVar2.propagateTo(secondaryVertex[0], bz);

      std::array<float, 3> const pvecpK = RecoDecay::pVec(track0.pVector(), track1.pVector());
      std::array<float, 3> pvecLc = RecoDecay::pVec(pvecpK, track2.pVector());
      auto trackpK = o2::dataformats::V0(df3.getPCACandidatePos(), pvecpK, df3.calcPCACovMatrixFlat(), trackParVar0, trackParVar1);
      auto trackLc = o2::dataformats::V0(df3.getPCACandidatePos(), pvecLc, df3.calcPCACovMatrixFlat(), trackpK, trackParVar2);

      int const index0Lc = track0.globalIndex();
      int const index1Lc = track1.globalIndex();
      int const index2Lc = track2.globalIndex();
      // int charge = track0.sign() + track1.sign() + track2.sign();

      for (const auto& trackPion : tracks) {
        if (trackPion.pt() < ptPionMin) {
          continue;
        }
        if (trackPion.sign() > 0) {
          continue;
        }
        if (trackPion.globalIndex() == index0Lc || trackPion.globalIndex() == index1Lc || trackPion.globalIndex() == index2Lc) {
          continue;
        }
        hPtPion->Fill(trackPion.pt());
        std::array<float, 3> pvecPion{};
        auto trackParVarPi = getTrackParCov(trackPion);

        // reconstruct the 3-prong Lc vertex
        hCandidatesLb->Fill(SVFitting::BeforeFit);
        try {
          if (df2.process(trackLc, trackParVarPi) == 0) {
            continue;
          }
        } catch (const std::runtime_error& error) {
          LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
          hCandidatesLb->Fill(SVFitting::Fail);
          continue;
        }
        hCandidatesLb->Fill(SVFitting::FitOk);

        // calculate relevant properties
        const auto& secondaryVertexLb = df2.getPCACandidate();
        auto chi2PCA = df2.getChi2AtPCACandidate();
        auto covMatrixPCA = df2.calcPCACovMatrixFlat();

        df2.propagateTracksToVertex();
        df2.getTrack(0).getPxPyPzGlo(pvecLc);
        df2.getTrack(1).getPxPyPzGlo(pvecPion);

        auto primaryVertex = getPrimaryVertex(collision);
        auto covMatrixPV = primaryVertex.getCov();
        o2::dataformats::DCA impactParameter0;
        o2::dataformats::DCA impactParameter1;
        trackLc.propagateToDCA(primaryVertex, bz, &impactParameter0);
        trackParVarPi.propagateToDCA(primaryVertex, bz, &impactParameter1);

        hCovSVXX->Fill(covMatrixPCA[0]);
        hCovPVXX->Fill(covMatrixPV[0]);

        // get uncertainty of the decay length
        double phi, theta;
        getPointDirection(std::array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertexLb, phi, theta);
        auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
        auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

        // fill the candidate table for the Lb here:
        rowCandidateBase(collision.globalIndex(),
                         collision.posX(), collision.posY(), collision.posZ(),
                         secondaryVertexLb[0], secondaryVertexLb[1], secondaryVertexLb[2],
                         errorDecayLength, errorDecayLengthXY,
                         chi2PCA,
                         pvecLc[0], pvecLc[1], pvecLc[2],
                         pvecPion[0], pvecPion[1], pvecPion[2],
                         impactParameter0.getY(), impactParameter1.getY(),
                         std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()));
        rowCandidateProngs(lcCand.globalIndex(), trackPion.globalIndex());
        // calculate invariant mass
        auto arrayMomenta = std::array{pvecLc, pvecPion};
        massLcPi = RecoDecay::m(arrayMomenta, std::array{MassLambdaCPlus, MassPiMinus});
        if (lcCand.isSelLcToPKPi() > 0) {
          hMassLbToLcPi->Fill(massLcPi);
        }
        if (lcCand.isSelLcToPiKP() > 0) {
          hMassLbToLcPi->Fill(massLcPi);
        }
      } // pi- loop
    } // Lc loop
  } // process
}; // struct

/// Extends the base table with expression columns.
struct HfCandidateCreatorLbExpressions {
  Spawns<aod::HfCandLbExt> rowCandidateLb;
  Produces<aod::HfCandLbMcRec> rowMcMatchRec;
  Produces<aod::HfCandLbMcGen> rowMcMatchGen;

  void init(InitContext const&) {}

  /// @brief dummy process function, to be run on data
  void process(aod::Tracks const&) {}

  void processMc(aod::HfCand3Prong const&,
                 aod::TracksWMc const&,
                 aod::McParticles const& mcParticles,
                 aod::HfCandLbProngs const& candsLb)
  {
    int indexRec = -1;
    int8_t sign = 0;
    int8_t flagChannelMain = 0;
    int8_t flagChannelReso = 0;
    int8_t origin = 0;
    int8_t debug = 0;

    // Match reconstructed candidates.
    for (const auto& candidate : candsLb) {
      flagChannelMain = 0;
      flagChannelReso = 0;
      origin = 0;
      debug = 0;
      auto lcCand = candidate.prong0();
      auto arrayDaughters = std::array{lcCand.prong0_as<aod::TracksWMc>(),
                                       lcCand.prong1_as<aod::TracksWMc>(),
                                       lcCand.prong2_as<aod::TracksWMc>(),
                                       candidate.prong1_as<aod::TracksWMc>()};
      auto arrayDaughtersLc = std::array{lcCand.prong0_as<aod::TracksWMc>(),
                                         lcCand.prong1_as<aod::TracksWMc>(),
                                         lcCand.prong2_as<aod::TracksWMc>()};
      // Λb → Λc+ π-
      indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kLambdaB0, std::array{+kProton, -kKPlus, +kPiPlus, -kPiPlus}, true, &sign, 2);
      if (indexRec > -1) {
        // Λb → Λc+ π-
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersLc, Pdg::kLambdaCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 1);
        if (indexRec > -1) {
          flagChannelMain = sign * DecayChannelMain::LbToLcPi;
        } else {
          debug = 1;
          LOGF(info, "WARNING: Λb in decays in the expected final state but the condition on the intermediate state is not fulfilled");
        }
      }
      rowMcMatchRec(flagChannelMain, flagChannelReso, origin, debug);
    }

    // Match generated particles.
    for (const auto& particle : mcParticles) {
      flagChannelMain = 0;
      flagChannelReso = 0;
      origin = 0;
      // Λb → Λc+ π-
      if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kLambdaB0, std::array{static_cast<int>(Pdg::kLambdaCPlus), -kPiPlus}, true)) {
        // Λc+ → p K- π+
        auto candLcMc = mcParticles.rawIteratorAt(particle.daughtersIds().front());
        if (RecoDecay::isMatchedMCGen(mcParticles, candLcMc, static_cast<int>(Pdg::kLambdaCPlus), std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign)) {
          flagChannelMain = sign * DecayChannelMain::LbToLcPi;
        }
      }
      rowMcMatchGen(flagChannelMain, flagChannelReso, origin);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorLbExpressions, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HfCandidateCreatorLb>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorLbExpressions>(cfgc)};
  return workflow;
}
