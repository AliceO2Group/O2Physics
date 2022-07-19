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

/// \file HFCandidateCreatorChic.cxx
/// \brief Reconstruction of chi_c candidates
/// \note Adapted from HFCandidateCreatorX
///
/// \author Alessandro De Falco <alessandro.de.falco@ca.infn.it>, Cagliari University

#include "Framework/AnalysisTask.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "ALICE3/DataModel/ECAL.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand;
//using namespace o2::aod::alice3ecal;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::aod::hf_cand_chic;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, true, {"Perform MC matching."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// Reconstruction of chic candidates
struct HFCandidateCreatorChic {
  Produces<aod::HfCandChicBase> rowCandidateBase;

  Configurable<double> magneticField{"magneticField", 20., "magnetic field"};
  Configurable<bool> b_propdca{"b_propdca", true, "create tracks version propagated to PCA"};
  Configurable<double> d_maxr{"d_maxr", 200., "reject PCA's above this radius"};
  Configurable<double> d_maxdzini{"d_maxdzini", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> d_minparamchange{"d_minparamchange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> d_minrelchi2change{"d_minrelchi2change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<double> eneGammaMin{"eneGammaMin", 0.4, "minimum gamma energy threshold (GeV)"};
  Configurable<double> etaGammaMin{"etaGammaMin", -1.00, "minimum gamma pseudorapidity"};
  Configurable<double> etaGammaMax{"etaGammaMax", 1.00, "maximum gamma pseudorapidity"};

  OutputObj<TH1F> hMassJpsiToEE{TH1F("hMassJpsiToEE", "J/#psi candidates;inv. mass (e^{#plus} e^{#minus}) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hMassJpsiToMuMu{TH1F("hMassJpsiToMuMu", "J/#psi candidates;inv. mass (#mu^{#plus} #mu^{#minus}) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtJpsi{TH1F("hPtJpsi", "J/#psi candidates;candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hCPAJpsi{TH1F("hCPAJpsi", "J/#psi candidates;cosine of pointing angle;entries", 110, -1.1, 1.1)};
  OutputObj<TH1F> hMassChicToJpsiToEEGamma{TH1F("hMassChicToJpsiToEEGamma", "2-prong candidates;inv. mass (J/#psi (#rightarrow e+ e-) #gamma) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hMassChicToJpsiToMuMuGamma{TH1F("hMassChicToJpsiToMuMuGamma", "2-prong candidates;inv. mass (J/#psi (#rightarrow #mu+ #mu-) #gamma) (GeV/#it{c}^{2});entries", 500, 0., 5.)};

  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx position (cm^{2});entries", 100, 0., 0.2)};

  double massJpsi = RecoDecay::getMassPDG(pdg::Code::kJpsi);
  double massJpsiGamma = 0;

  Configurable<int> d_selectionFlagJpsi{"d_selectionFlagJpsi", 1, "Selection Flag for Jpsi"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Filter filterSelectCandidates = (aod::hf_selcandidate_jpsi::isSelJpsiToEE >= d_selectionFlagJpsi || aod::hf_selcandidate_jpsi::isSelJpsiToMuMu >= d_selectionFlagJpsi);

  void process(aod::Collision const& collision,
               soa::Filtered<soa::Join<
                 aod::HfCandProng2,
                 aod::HFSelJpsiCandidate>> const& jpsiCands,
               aod::BigTracks const& tracks,
               aod::ECALs const& ecals)
  {
    // 2-prong vertex fitter (to rebuild Jpsi vertex)
    o2::vertexing::DCAFitterN<2> df2;
    df2.setBz(magneticField);
    df2.setPropagateToPCA(b_propdca);
    df2.setMaxR(d_maxr);
    df2.setMaxDZIni(d_maxdzini);
    df2.setMinParamChange(d_minparamchange);
    df2.setMinRelChi2Change(d_minrelchi2change);
    df2.setUseAbsDCA(true);

    // loop over Jpsi candidates
    for (auto& jpsiCand : jpsiCands) {
      if (!(jpsiCand.hfflag() & 1 << hf_cand_prong2::DecayType::JpsiToEE) && !(jpsiCand.hfflag() & 1 << hf_cand_prong2::DecayType::JpsiToMuMu)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YJpsi(jpsiCand)) > cutYCandMax) {
        continue;
      }
      if (jpsiCand.isSelJpsiToEE() > 0) {
        hMassJpsiToEE->Fill(InvMassJpsiToEE(jpsiCand));
      }
      if (jpsiCand.isSelJpsiToMuMu() > 0) {
        hMassJpsiToMuMu->Fill(InvMassJpsiToMuMu(jpsiCand));
      }
      hPtJpsi->Fill(jpsiCand.pt());
      hCPAJpsi->Fill(jpsiCand.cpa());
      // create Jpsi track to pass to DCA fitter; use cand table + rebuild vertex
      const std::array<float, 3> vertexJpsi = {jpsiCand.xSecondaryVertex(), jpsiCand.ySecondaryVertex(), jpsiCand.zSecondaryVertex()};
      array<float, 3> pvecJpsi = {jpsiCand.px(), jpsiCand.py(), jpsiCand.pz()};
      auto prong0 = jpsiCand.index0_as<aod::BigTracks>();
      auto prong1 = jpsiCand.index1_as<aod::BigTracks>();
      auto prong0TrackParCov = getTrackParCov(prong0);
      auto prong1TrackParCov = getTrackParCov(prong1);

      if (df2.process(prong0TrackParCov, prong1TrackParCov) == 0) {
        continue;
      }

      // propogate prong tracks to Jpsi vertex
      prong0TrackParCov.propagateTo(jpsiCand.xSecondaryVertex(), magneticField);
      prong1TrackParCov.propagateTo(jpsiCand.xSecondaryVertex(), magneticField);
      const std::array<float, 6> covJpsi = df2.calcPCACovMatrixFlat();
      // define the Jpsi track
      auto trackJpsi = o2::dataformats::V0(vertexJpsi, pvecJpsi, covJpsi, prong0TrackParCov, prong1TrackParCov, {0, 0}, {0, 0}); //FIXME: also needs covxyz???

      // -----------------------------------------------------------------
      // loop over gamma candidates

      for (auto& ecal : ecals) {

        if (ecal.e() < eneGammaMin) {
          continue;
        }
        auto etagamma = RecoDecay::eta(array{ecal.px(), ecal.py(), ecal.pz()});
        if (etagamma < etaGammaMin || etagamma > etaGammaMax) { // calcolare la pseudorapidità da posz
          continue;
        }

        array<float, 3> pvecGamma{(float)ecal.px(), (float)ecal.py(), (float)ecal.pz()};

        // get track impact parameters
        // This modifies track momenta!

        auto primaryVertex = getPrimaryVertex(collision);
        auto covMatrixPV = primaryVertex.getCov();
        hCovPVXX->Fill(covMatrixPV[0]);
        o2::dataformats::DCA impactParameter0;
        o2::dataformats::DCA impactParameter1;
        trackJpsi.propagateToDCA(primaryVertex, magneticField, &impactParameter0);

        // get uncertainty of the decay length
        //double phi, theta;
        // getPointDirection(array{collision.posX(), collision.posY(), collision.posZ()}, ChicsecondaryVertex, phi, theta);
        //auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
        //auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

        int hfFlag = 0;
        if (TESTBIT(jpsiCand.hfflag(), hf_cand_prong2::DecayType::JpsiToMuMu)) {
          SETBIT(hfFlag, hf_cand_chic::DecayType::ChicToJpsiToMuMuGamma); // dimuon channel
        }
        if (TESTBIT(jpsiCand.hfflag(), hf_cand_prong2::DecayType::JpsiToEE)) {
          SETBIT(hfFlag, hf_cand_chic::DecayType::ChicToJpsiToEEGamma); // dielectron channel
        }

        // fill the candidate table for the chi_c here:
        rowCandidateBase(collision.globalIndex(),
                         collision.posX(), collision.posY(), collision.posZ(),
                         0.f, 0.f, 0.f,               //    ChicsecondaryVertex[0], ChicsecondaryVertex[1], ChicsecondaryVertex[2],
                         0.f, 0.f,                    // errorDecayLength, errorDecayLengthXY,
                         df2.getChi2AtPCACandidate(), //chi2PCA of Jpsi
                         pvecJpsi[0], pvecJpsi[1], pvecJpsi[2],
                         pvecGamma[0], pvecGamma[1], pvecGamma[2],
                         impactParameter0.getY(), 0.f,                  // impactParameter1.getY(),
                         std::sqrt(impactParameter0.getSigmaY2()), 0.f, // std::sqrt(impactParameter1.getSigmaY2()),
                         jpsiCand.globalIndex(), ecal.globalIndex(),
                         hfFlag, InvMassJpsiToMuMu(jpsiCand));

        // calculate invariant mass
        auto arrayMomenta = array{pvecJpsi, pvecGamma};
        massJpsiGamma = RecoDecay::m(std::move(arrayMomenta), array{massJpsi, 0.});
        if (jpsiCand.isSelJpsiToEE() > 0) {
          hMassChicToJpsiToEEGamma->Fill(massJpsiGamma);
        }
        if (jpsiCand.isSelJpsiToMuMu() > 0) {
          hMassChicToJpsiToMuMuGamma->Fill(massJpsiGamma);
        }
      } // ecal loop
    }   // Jpsi loop
  }     // process
};      // struct

/// Extends the base table with expression columns.
struct HFCandidateCreatorChicExpressions {
  Spawns<aod::HfCandChicExt> rowCandidateChic;
  void init(InitContext const&) {}
};

/// Performs MC matching.
struct HFCandidateCreatorChicMC {
  Produces<aod::HfCandChicMCRec> rowMCMatchRec;
  Produces<aod::HfCandChicMCGen> rowMCMatchGen;
  OutputObj<TH1F> hMassJpsiToMuMuMatched{TH1F("hMassChicToJpsiToMuMuMatched", "2-prong candidates;inv. mass (J/#psi (#rightarrow #mu+ #mu-)) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hMassEMatched{TH1F("hMassEMatched", "2-prong candidates;inv. mass (J/#psi (#rightarrow #mu+ #mu-)) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hEphotonMatched{TH1F("hEphotonMatched", "2-prong candidates;inv. mass (J/#psi (#rightarrow #mu+ #mu-)) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hMassChicToJpsiToMuMuGammaMatched{TH1F("hMassChicToJpsiToMuMuGammaMatched", "2-prong candidates;inv. mass (J/#psi (#rightarrow #mu+ #mu-) #gamma) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  void process(aod::HfCandChic const& candidates,
               aod::HfCandProng2 const&,
               aod::BigTracksMC const& tracks,
               aod::McParticles const& particlesMC,
               aod::ECALs const& ecals)
  {
    int indexRec = -1;
    //int8_t sign = 0;
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t channel = 0;

    // Match reconstructed candidates.
    for (auto& candidate : candidates) {
      flag = 0;
      origin = 0;
      channel = 0;
      auto jpsiTrack = candidate.index0();
      auto daughterPosJpsi = jpsiTrack.index0_as<aod::BigTracksMC>();
      auto daughterNegJpsi = jpsiTrack.index1_as<aod::BigTracksMC>();
      auto arrayJpsiDaughters = array{daughterPosJpsi, daughterNegJpsi};

      // chi_c → J/ψ gamma
      indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayJpsiDaughters, pdg::Code::kJpsi, array{+kMuonPlus, -kMuonPlus}, true);
      if (indexRec > -1) {
        hMassJpsiToMuMuMatched->Fill(InvMassJpsiToMuMu(candidate.index0()));

        int indexMother = RecoDecay::getMother(particlesMC, particlesMC.rawIteratorAt(indexRec), pdg::Code::kChic1);
        int indexMotherGamma = RecoDecay::getMother(particlesMC, particlesMC.rawIteratorAt(candidate.index1().mcparticleId()), pdg::Code::kChic1);
        if (indexMother > -1 && indexMotherGamma == indexMother && candidate.index1().mcparticle().pdgCode() == kGamma) {
          auto particleMother = particlesMC.rawIteratorAt(indexMother);
          hEphotonMatched->Fill(candidate.index1().e());
          hMassEMatched->Fill(sqrt(candidate.index1().px() * candidate.index1().px() + candidate.index1().py() * candidate.index1().py() + candidate.index1().pz() * candidate.index1().pz()));
          if (particleMother.has_daughters()) {
            std::vector<int> arrAllDaughtersIndex;
            RecoDecay::getDaughters(particleMother, &arrAllDaughtersIndex, array{(int)(kGamma), (int)(pdg::Code::kJpsi)}, 1);
            if (arrAllDaughtersIndex.size() == 2) {
              flag = 1 << hf_cand_chic::DecayType::ChicToJpsiToMuMuGamma;
              hMassChicToJpsiToMuMuGammaMatched->Fill(InvMassChicToJpsiGamma(candidate));
            }
          }
        }
      }
      if (flag != 0) {
        auto particle = particlesMC.rawIteratorAt(indexRec);
        origin = RecoDecay::getCharmHadronOrigin(particlesMC, particle);
      }
      rowMCMatchRec(flag, origin, channel);
    }

    // Match generated particles.
    for (auto& particle : particlesMC) {
      //Printf("New gen. candidate");
      flag = 0;
      origin = 0;
      channel = 0;

      // chi_c → J/ψ gamma
      if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdg::Code::kChic1, array{(int)(pdg::Code::kJpsi), (int)(kGamma)}, true)) {
        // Match J/psi --> e+e-
        std::vector<int> arrDaughter;
        RecoDecay::getDaughters(particle, &arrDaughter, array{(int)(pdg::Code::kJpsi)}, 1);
        auto jpsiCandMC = particlesMC.rawIteratorAt(arrDaughter[0]);
        if (RecoDecay::isMatchedMCGen(particlesMC, jpsiCandMC, pdg::Code::kJpsi, array{+kElectron, -kElectron}, true)) {
          flag = 1 << hf_cand_chic::DecayType::ChicToJpsiToEEGamma;
        }

        if (flag == 0) {
          if (RecoDecay::isMatchedMCGen(particlesMC, jpsiCandMC, pdg::Code::kJpsi, array{+kMuonPlus, -kMuonPlus}, true)) {
            flag = 1 << hf_cand_chic::DecayType::ChicToJpsiToMuMuGamma;
          }
        }
      }

      rowMCMatchGen(flag, origin, channel);
    } // candidate loop
  }   // process
};    // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HFCandidateCreatorChic>(cfgc, TaskName{"hf-cand-creator-chic"}),
    adaptAnalysisTask<HFCandidateCreatorChicExpressions>(cfgc, TaskName{"hf-cand-creator-chic-expressions"})};
  const bool doMC = cfgc.options().get<bool>("doMC");
  if (doMC) {
    workflow.push_back(adaptAnalysisTask<HFCandidateCreatorChicMC>(cfgc, TaskName{"hf-cand-creator-chic-mc"}));
  }
  return workflow;
}
