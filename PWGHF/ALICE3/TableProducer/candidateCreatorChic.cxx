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

/// \file candidateCreatorChic.cxx
/// \brief Reconstruction of chi_c candidates
/// \note Adapted from candidateCreatorX.cxx
///
/// \author Alessandro De Falco <alessandro.de.falco@ca.infn.it>, Cagliari University

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "ALICE3/DataModel/ECAL.h"
#include "Common/Core/trackUtilities.h"

#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"

#include <utility>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, true, {"Perform MC matching."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// Reconstruction of Chic candidates
struct HfCandidateCreatorChic {
  Produces<aod::HfCandChicBase> rowCandidateBase;

  // vertexing
  Configurable<double> bz{"bz", 20., "magnetic field"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  // selection
  Configurable<double> energyGammaMin{"energyGammaMin", 0.4, "minimum gamma energy threshold (GeV)"};
  Configurable<double> etaGammaMin{"etaGammaMin", -1.00, "minimum gamma pseudorapidity"};
  Configurable<double> etaGammaMax{"etaGammaMax", 1.00, "maximum gamma pseudorapidity"};
  Configurable<int> selectionFlagJpsi{"selectionFlagJpsi", 1, "Selection Flag for Jpsi"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};

  o2::vertexing::DCAFitterN<2> df2; // 2-prong vertex fitter (to rebuild Jpsi vertex)

  double massJpsi{0.};
  double massJpsiGamma{0.};

  Filter filterSelectCandidates = (aod::hf_sel_candidate_jpsi::isSelJpsiToEE >= selectionFlagJpsi || aod::hf_sel_candidate_jpsi::isSelJpsiToMuMu >= selectionFlagJpsi);

  OutputObj<TH1F> hMassJpsiToEE{TH1F("hMassJpsiToEE", "J/#psi candidates;inv. mass (e^{#plus} e^{#minus}) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hMassJpsiToMuMu{TH1F("hMassJpsiToMuMu", "J/#psi candidates;inv. mass (#mu^{#plus} #mu^{#minus}) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtJpsi{TH1F("hPtJpsi", "J/#psi candidates;candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hCPAJpsi{TH1F("hCPAJpsi", "J/#psi candidates;cosine of pointing angle;entries", 110, -1.1, 1.1)};
  OutputObj<TH1F> hMassChicToJpsiToEEGamma{TH1F("hMassChicToJpsiToEEGamma", "2-prong candidates;inv. mass (J/#psi (#rightarrow e+ e-) #gamma) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hMassChicToJpsiToMuMuGamma{TH1F("hMassChicToJpsiToMuMuGamma", "2-prong candidates;inv. mass (J/#psi (#rightarrow #mu+ #mu-) #gamma) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};

  void init(InitContext const&)
  {
    massJpsi = MassJPsi;

    df2.setBz(bz);
    df2.setPropagateToPCA(propagateToPCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(useAbsDCA);
    df2.setWeightedFinalPCA(useWeightedFinalPCA);
  }

  void process(aod::Collision const& collision,
               soa::Filtered<soa::Join<
                 aod::HfCand2Prong,
                 aod::HfSelJpsi>> const& jpsiCands,
               aod::TracksWCov const&,
               aod::ECALs const& ecals)
  {
    // loop over Jpsi candidates
    for (const auto& jpsiCand : jpsiCands) {
      if (!(jpsiCand.hfflag() & 1 << hf_cand_2prong::DecayType::JpsiToEE) && !(jpsiCand.hfflag() & 1 << hf_cand_2prong::DecayType::JpsiToMuMu)) {
        continue;
      }
      if (yCandMax >= 0. && std::abs(HfHelper::yJpsi(jpsiCand)) > yCandMax) {
        continue;
      }
      if (jpsiCand.isSelJpsiToEE() > 0) {
        hMassJpsiToEE->Fill(HfHelper::invMassJpsiToEE(jpsiCand));
      }
      if (jpsiCand.isSelJpsiToMuMu() > 0) {
        hMassJpsiToMuMu->Fill(HfHelper::invMassJpsiToMuMu(jpsiCand));
      }
      hPtJpsi->Fill(jpsiCand.pt());
      hCPAJpsi->Fill(jpsiCand.cpa());
      // create Jpsi track to pass to DCA fitter; use cand table + rebuild vertex
      const std::array<float, 3> vertexJpsi = {jpsiCand.xSecondaryVertex(), jpsiCand.ySecondaryVertex(), jpsiCand.zSecondaryVertex()};
      std::array<float, 3> pvecJpsi = jpsiCand.pVector();
      auto prong0 = jpsiCand.prong0_as<aod::TracksWCov>();
      auto prong1 = jpsiCand.prong1_as<aod::TracksWCov>();
      auto prong0TrackParCov = getTrackParCov(prong0);
      auto prong1TrackParCov = getTrackParCov(prong1);

      if (df2.process(prong0TrackParCov, prong1TrackParCov) == 0) {
        continue;
      }

      // propagate prong tracks to Jpsi vertex
      prong0TrackParCov.propagateTo(jpsiCand.xSecondaryVertex(), bz);
      prong1TrackParCov.propagateTo(jpsiCand.xSecondaryVertex(), bz);
      const std::array<float, 6> covJpsi = df2.calcPCACovMatrixFlat();
      // define the Jpsi track
      auto trackJpsi = o2::dataformats::V0(vertexJpsi, pvecJpsi, covJpsi, prong0TrackParCov, prong1TrackParCov); // FIXME: also needs covxyz???

      // -----------------------------------------------------------------
      // loop over gamma candidates

      for (const auto& ecal : ecals) {

        if (ecal.e() < energyGammaMin) {
          continue;
        }
        auto etagamma = RecoDecay::eta(std::array{ecal.px(), ecal.py(), ecal.pz()});
        if (etagamma < etaGammaMin || etagamma > etaGammaMax) { // calcolare la pseudorapidità da posz
          continue;
        }

        std::array<float, 3> pvecGamma{static_cast<float>(ecal.px()), static_cast<float>(ecal.py()), static_cast<float>(ecal.pz())};

        // get track impact parameters
        // This modifies track momenta!

        auto primaryVertex = getPrimaryVertex(collision);
        auto covMatrixPV = primaryVertex.getCov();
        hCovPVXX->Fill(covMatrixPV[0]);
        o2::dataformats::DCA impactParameter0;
        o2::dataformats::DCA impactParameter1;
        trackJpsi.propagateToDCA(primaryVertex, bz, &impactParameter0);

        // get uncertainty of the decay length
        // double phi, theta;
        // getPointDirection(std::array{collision.posX(), collision.posY(), collision.posZ()}, ChicsecondaryVertex, phi, theta);
        // auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
        // auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

        int hfFlag = 0;
        if (TESTBIT(jpsiCand.hfflag(), hf_cand_2prong::DecayType::JpsiToMuMu)) {
          SETBIT(hfFlag, hf_cand_chic::DecayType::ChicToJpsiToMuMuGamma); // dimuon channel
        }
        if (TESTBIT(jpsiCand.hfflag(), hf_cand_2prong::DecayType::JpsiToEE)) {
          SETBIT(hfFlag, hf_cand_chic::DecayType::ChicToJpsiToEEGamma); // dielectron channel
        }

        // fill the candidate table for the chi_c here:
        rowCandidateBase(collision.globalIndex(),
                         collision.posX(), collision.posY(), collision.posZ(),
                         0.f, 0.f, 0.f,               //    ChicsecondaryVertex[0], ChicsecondaryVertex[1], ChicsecondaryVertex[2],
                         0.f, 0.f,                    // errorDecayLength, errorDecayLengthXY,
                         df2.getChi2AtPCACandidate(), // chi2PCA of Jpsi
                         pvecJpsi[0], pvecJpsi[1], pvecJpsi[2],
                         pvecGamma[0], pvecGamma[1], pvecGamma[2],
                         impactParameter0.getY(), 0.f,                  // impactParameter1.getY(),
                         std::sqrt(impactParameter0.getSigmaY2()), 0.f, // std::sqrt(impactParameter1.getSigmaY2()),
                         jpsiCand.globalIndex(), ecal.globalIndex(),
                         hfFlag, HfHelper::invMassJpsiToMuMu(jpsiCand));

        // calculate invariant mass
        auto arrayMomenta = std::array{pvecJpsi, pvecGamma};
        massJpsiGamma = RecoDecay::m(std::move(arrayMomenta), std::array{massJpsi, 0.});
        if (jpsiCand.isSelJpsiToEE() > 0) {
          hMassChicToJpsiToEEGamma->Fill(massJpsiGamma);
        }
        if (jpsiCand.isSelJpsiToMuMu() > 0) {
          hMassChicToJpsiToMuMuGamma->Fill(massJpsiGamma);
        }
      } // ecal loop
    } // Jpsi loop
  } // process
}; // struct

/// Extends the base table with expression columns.
struct HfCandidateCreatorChicExpressions {
  Spawns<aod::HfCandChicExt> rowCandidateChic;

  void init(InitContext const&) {}
};

/// Performs MC matching.
struct HfCandidateCreatorChicMc {
  Produces<aod::HfCandChicMcRec> rowMcMatchRec;
  Produces<aod::HfCandChicMcGen> rowMcMatchGen;

  OutputObj<TH1F> hMassJpsiToMuMuMatched{TH1F("hMassChicToJpsiToMuMuMatched", "2-prong candidates;inv. mass (J/#psi (#rightarrow #mu+ #mu-)) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hMassEMatched{TH1F("hMassEMatched", "2-prong candidates;inv. mass (J/#psi (#rightarrow #mu+ #mu-)) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hEphotonMatched{TH1F("hEphotonMatched", "2-prong candidates;inv. mass (J/#psi (#rightarrow #mu+ #mu-)) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hMassChicToJpsiToMuMuGammaMatched{TH1F("hMassChicToJpsiToMuMuGammaMatched", "2-prong candidates;inv. mass (J/#psi (#rightarrow #mu+ #mu-) #gamma) (GeV/#it{c}^{2});entries", 500, 0., 5.)};

  void process(aod::HfCandChic const& candidates,
               aod::HfCand2Prong const&,
               aod::TracksWMc const&,
               aod::McParticles const& mcParticles,
               aod::ECALs const&)
  {
    int indexRec = -1;
    // int8_t sign = 0;
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t channel = 0;

    // Match reconstructed candidates.
    for (const auto& candidate : candidates) {
      flag = 0;
      origin = 0;
      channel = 0;
      auto jpsiTrack = candidate.prong0();
      auto daughterPosJpsi = jpsiTrack.prong0_as<aod::TracksWMc>();
      auto daughterNegJpsi = jpsiTrack.prong1_as<aod::TracksWMc>();
      auto arrayJpsiDaughters = std::array{daughterPosJpsi, daughterNegJpsi};

      // chi_c → J/ψ gamma
      indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayJpsiDaughters, Pdg::kJPsi, std::array{+kMuonPlus, -kMuonPlus}, true);
      if (indexRec > -1) {
        hMassJpsiToMuMuMatched->Fill(HfHelper::invMassJpsiToMuMu(candidate.prong0()));

        int indexMother = RecoDecay::getMother(mcParticles, mcParticles.rawIteratorAt(indexRec), Pdg::kChiC1);
        int indexMotherGamma = RecoDecay::getMother(mcParticles, mcParticles.rawIteratorAt(candidate.prong1().mcparticleId()), Pdg::kChiC1);
        if (indexMother > -1 && indexMotherGamma == indexMother && candidate.prong1().mcparticle().pdgCode() == kGamma) {
          auto particleMother = mcParticles.rawIteratorAt(indexMother);
          hEphotonMatched->Fill(candidate.prong1().e());
          hMassEMatched->Fill(RecoDecay::p(candidate.pVectorProng1()));
          if (particleMother.has_daughters()) {
            std::vector<int> arrAllDaughtersIndex;
            RecoDecay::getDaughters(particleMother, &arrAllDaughtersIndex, std::array{static_cast<int>(kGamma), static_cast<int>(Pdg::kJPsi)}, 1);
            if (arrAllDaughtersIndex.size() == 2) {
              flag = 1 << hf_cand_chic::DecayType::ChicToJpsiToMuMuGamma;
              hMassChicToJpsiToMuMuGammaMatched->Fill(HfHelper::invMassChicToJpsiGamma(candidate));
            }
          }
        }
      }
      if (flag != 0) {
        auto particle = mcParticles.rawIteratorAt(indexRec);
        origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle);
      }
      rowMcMatchRec(flag, origin, channel);
    }

    // Match generated particles.
    for (const auto& particle : mcParticles) {
      flag = 0;
      origin = 0;
      channel = 0;

      // chi_c → J/ψ gamma
      if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kChiC1, std::array{static_cast<int>(Pdg::kJPsi), static_cast<int>(kGamma)}, true)) {
        // Match J/psi --> e+e-
        std::vector<int> arrDaughter;
        RecoDecay::getDaughters(particle, &arrDaughter, std::array{static_cast<int>(Pdg::kJPsi)}, 1);
        auto jpsiCandMC = mcParticles.rawIteratorAt(arrDaughter[0]);
        if (RecoDecay::isMatchedMCGen(mcParticles, jpsiCandMC, Pdg::kJPsi, std::array{+kElectron, -kElectron}, true)) {
          flag = 1 << hf_cand_chic::DecayType::ChicToJpsiToEEGamma;
        }

        if (flag == 0) {
          if (RecoDecay::isMatchedMCGen(mcParticles, jpsiCandMC, Pdg::kJPsi, std::array{+kMuonPlus, -kMuonPlus}, true)) {
            flag = 1 << hf_cand_chic::DecayType::ChicToJpsiToMuMuGamma;
          }
        }
      }

      rowMcMatchGen(flag, origin, channel);
    } // candidate loop
  } // process
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HfCandidateCreatorChic>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorChicExpressions>(cfgc)};
  const bool doMC = cfgc.options().get<bool>("doMC");
  if (doMC) {
    workflow.push_back(adaptAnalysisTask<HfCandidateCreatorChicMc>(cfgc));
  }
  return workflow;
}
