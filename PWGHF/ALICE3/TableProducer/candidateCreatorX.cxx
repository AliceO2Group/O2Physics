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

/// \file candidateCreatorX.cxx
/// \brief Reconstruction of X(3872) candidates
/// \note Adapted from candidateCreator3Prong.cxx
///
/// \author Rik Spijkers <r.spijkers@students.uu.nl>, Utrecht University
/// \author Luca Micheletti <luca.micheletti@to.infn.it>, INFN

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

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

/// Reconstruction of X candidates
struct HfCandidateCreatorX {
  Produces<aod::HfCandXBase> rowCandidateBase;

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
  Configurable<double> ptPionMin{"ptPionMin", 1., "minimum pion pT threshold (GeV/c)"};
  Configurable<int> selectionFlagJpsi{"selectionFlagJpsi", 1, "Selection Flag for Jpsi"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<double> diffMassJpsiMax{"diffMassJpsiMax", 0.07, "max. diff. between Jpsi rec. and PDG mass"};

  o2::vertexing::DCAFitterN<2> df2; // 2-prong vertex fitter (to rebuild Jpsi vertex)
  o2::vertexing::DCAFitterN<3> df3; // 3-prong vertex fitter

  double massPi{0.};
  double massJpsi{0.};
  double massJpsiPiPi{0.};

  Filter filterSelectCandidates = (aod::hf_sel_candidate_jpsi::isSelJpsiToEE >= selectionFlagJpsi || aod::hf_sel_candidate_jpsi::isSelJpsiToMuMu >= selectionFlagJpsi);

  OutputObj<TH1F> hMassJpsiToEE{TH1F("hMassJpsiToEE", "J/#psi candidates;inv. mass (e^{#plus} e^{#minus}) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hMassJpsiToMuMu{TH1F("hMassJpsiToMuMu", "J/#psi candidates;inv. mass (#mu^{#plus} #mu^{#minus}) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtJpsi{TH1F("hPtJpsi", "J/#psi candidates;candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hCPAJpsi{TH1F("hCPAJpsi", "J/#psi candidates;cosine of pointing angle;entries", 110, -1.1, 1.1)};
  OutputObj<TH1F> hMassXToJpsiToEEPiPi{TH1F("hMassXToJpsiToEEPiPi", "3-prong candidates;inv. mass (J/#psi (#rightarrow e+ e-) #pi+ #pi-) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hMassXToJpsiToMuMuPiPi{TH1F("hMassXToJpsiToMuMuPiPi", "3-prong candidates;inv. mass (J/#psi (#rightarrow #mu+ #mu-) #pi+ #pi-) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "3-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "3-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};

  void init(InitContext const&)
  {
    massPi = MassPiPlus;
    massJpsi = MassJPsi;

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
  }

  void process(aod::Collision const& collision,
               soa::Filtered<soa::Join<
                 aod::HfCand2Prong,
                 aod::HfSelJpsi>> const& jpsiCands,
               aod::TracksWCov const& tracks)
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
        if (std::abs(HfHelper::invMassJpsiToEE(jpsiCand) - massJpsi) > diffMassJpsiMax) {
          continue;
        }
        hMassJpsiToEE->Fill(HfHelper::invMassJpsiToEE(jpsiCand));
      }
      if (jpsiCand.isSelJpsiToMuMu() > 0) {
        if (std::abs(HfHelper::invMassJpsiToMuMu(jpsiCand) - massJpsi) > diffMassJpsiMax) {
          continue;
        }
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

      // used to check that prongs used for Jpsi and X reco are not the same prongs
      int index0Jpsi = jpsiCand.prong0Id();
      int index1Jpsi = jpsiCand.prong1Id();

      // loop over pi+ candidates
      for (const auto& trackPos : tracks) {
        if (trackPos.pt() < ptPionMin) {
          continue;
        }
        if (trackPos.sign() < 0) { // select only positive tracks - use partitions?
          continue;
        }
        if (trackPos.globalIndex() == index0Jpsi) {
          continue;
        }

        // loop over pi- candidates
        for (const auto& trackNeg : tracks) {
          if (trackNeg.pt() < ptPionMin) {
            continue;
          }
          if (trackNeg.sign() > 0) { // select only negative tracks - use partitions?
            continue;
          }
          if (trackNeg.globalIndex() == index1Jpsi) {
            continue;
          }

          auto trackParVarPos = getTrackParCov(trackPos);
          auto trackParVarNeg = getTrackParCov(trackNeg);
          std::array<float, 3> pvecPos;
          std::array<float, 3> pvecNeg;

          // reconstruct the 3-prong X vertex
          if (df3.process(trackJpsi, trackParVarPos, trackParVarNeg) == 0) {
            continue;
          }

          // calculate relevant properties
          const auto& XsecondaryVertex = df3.getPCACandidate();
          auto chi2PCA = df3.getChi2AtPCACandidate();
          auto covMatrixPCA = df3.calcPCACovMatrixFlat();
          hCovSVXX->Fill(covMatrixPCA[0]); // FIXME: Calculation of errorDecayLength(XY) gives wrong values without this line.

          df3.propagateTracksToVertex();          // propagate the pions and Jpsi to the X vertex
          df3.getTrack(0).getPxPyPzGlo(pvecJpsi); // update momentum of Jpsi at the X vertex
          df3.getTrack(1).getPxPyPzGlo(pvecPos);  // momentum of pi+ at the X vertex
          df3.getTrack(2).getPxPyPzGlo(pvecNeg);  // momentum of pi- at the X vertex

          // get track impact parameters
          // This modifies track momenta!
          auto primaryVertex = getPrimaryVertex(collision);
          auto covMatrixPV = primaryVertex.getCov();
          hCovPVXX->Fill(covMatrixPV[0]);
          o2::dataformats::DCA impactParameter0;
          o2::dataformats::DCA impactParameter1;
          o2::dataformats::DCA impactParameter2;
          trackJpsi.propagateToDCA(primaryVertex, bz, &impactParameter0);
          trackParVarPos.propagateToDCA(primaryVertex, bz, &impactParameter1);
          trackParVarNeg.propagateToDCA(primaryVertex, bz, &impactParameter2);

          // get uncertainty of the decay length
          double phi, theta;
          getPointDirection(std::array{collision.posX(), collision.posY(), collision.posZ()}, XsecondaryVertex, phi, theta);
          auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
          auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

          int hfFlag = 0;
          if (TESTBIT(jpsiCand.hfflag(), hf_cand_2prong::DecayType::JpsiToMuMu)) {
            SETBIT(hfFlag, hf_cand_x::DecayType::XToJpsiToMuMuPiPi); // dimuon channel
          }
          if (TESTBIT(jpsiCand.hfflag(), hf_cand_2prong::DecayType::JpsiToEE)) {
            SETBIT(hfFlag, hf_cand_x::DecayType::XToJpsiToEEPiPi); // dielectron channel
          }

          // fill the candidate table for the X here:
          rowCandidateBase(collision.globalIndex(),
                           collision.posX(), collision.posY(), collision.posZ(),
                           XsecondaryVertex[0], XsecondaryVertex[1], XsecondaryVertex[2],
                           errorDecayLength, errorDecayLengthXY,
                           chi2PCA,
                           pvecJpsi[0], pvecJpsi[1], pvecJpsi[2],
                           pvecPos[0], pvecPos[1], pvecPos[2],
                           pvecNeg[0], pvecNeg[1], pvecNeg[2],
                           impactParameter0.getY(), impactParameter1.getY(), impactParameter2.getY(),
                           std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()), std::sqrt(impactParameter2.getSigmaY2()),
                           jpsiCand.globalIndex(), trackPos.globalIndex(), trackNeg.globalIndex(),
                           hfFlag);

          // calculate invariant mass
          auto arrayMomenta = std::array{pvecJpsi, pvecPos, pvecNeg};
          massJpsiPiPi = RecoDecay::m(std::move(arrayMomenta), std::array{massJpsi, massPi, massPi});
          if (jpsiCand.isSelJpsiToEE() > 0) {
            hMassXToJpsiToEEPiPi->Fill(massJpsiPiPi);
          }
          if (jpsiCand.isSelJpsiToMuMu() > 0) {
            hMassXToJpsiToMuMuPiPi->Fill(massJpsiPiPi);
          }
        } // pi- loop
      } // pi+ loop
    } // Jpsi loop
  } // process
}; // struct

/// Extends the base table with expression columns.
struct HfCandidateCreatorXExpressions {
  Spawns<aod::HfCandXExt> rowCandidateX;

  void init(InitContext const&) {}
};

/// Performs MC matching.
struct HfCandidateCreatorXMc {
  Produces<aod::HfCandXMcRec> rowMcMatchRec;
  Produces<aod::HfCandXMcGen> rowMcMatchGen;

  void process(aod::HfCandX const& candidates,
               aod::HfCand2Prong const&,
               aod::TracksWMc const&,
               aod::McParticles const& mcParticles)
  {
    int indexRec = -1;
    int pdgCodeX = Pdg::kX3872;
    int pdgCodeJpsi = Pdg::kJPsi;
    int8_t sign = 0;
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
      auto arrayDaughters = std::array{candidate.prong1_as<aod::TracksWMc>(),
                                       candidate.prong2_as<aod::TracksWMc>(),
                                       daughterPosJpsi,
                                       daughterNegJpsi};

      // X → J/ψ π+ π−

      // J/ψ → e+ e−
      indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayJpsiDaughters, Pdg::kJPsi, std::array{+kElectron, -kElectron}, true);
      // X → π+ π− e+ e−
      if (indexRec > -1) {
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, pdgCodeX, std::array{+kPiPlus, -kPiPlus, +kElectron, -kElectron}, true, &sign, 2);
        if (indexRec > -1) {
          flag = 1 << hf_cand_x::DecayType::XToJpsiToEEPiPi;
        }
      }

      // J/ψ → μ+ μ−
      if (flag == 0) {
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayJpsiDaughters, Pdg::kJPsi, std::array{+kMuonPlus, -kMuonPlus}, true);
        // X → π+ π− μ+ μ−
        if (indexRec > -1) {
          indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, pdgCodeX, std::array{+kPiPlus, -kPiPlus, +kMuonPlus, -kMuonPlus}, true, &sign, 2);
          if (indexRec > -1) {
            flag = 1 << hf_cand_x::DecayType::XToJpsiToMuMuPiPi;
          }
        }
      }

      // Check whether the particle is non-prompt (from a b quark).
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

      // X → J/ψ π+ π−
      if (RecoDecay::isMatchedMCGen(mcParticles, particle, pdgCodeX, std::array{pdgCodeJpsi, +kPiPlus, -kPiPlus}, true)) {
        std::vector<int> arrDaughter;
        RecoDecay::getDaughters(particle, &arrDaughter, std::array{pdgCodeJpsi}, 1);
        auto jpsiCandMC = mcParticles.rawIteratorAt(arrDaughter[0]);
        // J/ψ → e+ e−
        if (RecoDecay::isMatchedMCGen(mcParticles, jpsiCandMC, pdgCodeJpsi, std::array{+kElectron, -kElectron}, true)) {
          flag = 1 << hf_cand_x::DecayType::XToJpsiToEEPiPi;
        }

        // J/ψ → μ+ μ−
        if (flag == 0) {
          if (RecoDecay::isMatchedMCGen(mcParticles, jpsiCandMC, pdgCodeJpsi, std::array{+kMuonPlus, -kMuonPlus}, true)) {
            flag = 1 << hf_cand_x::DecayType::XToJpsiToMuMuPiPi;
          }
        }
      }

      // Check whether the particle is non-prompt (from a b quark).
      if (flag != 0) {
        origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle);
      }

      rowMcMatchGen(flag, origin, channel);
    } // candidate loop
  } // process
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HfCandidateCreatorX>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorXExpressions>(cfgc)};
  const bool doMC = cfgc.options().get<bool>("doMC");
  if (doMC) {
    workflow.push_back(adaptAnalysisTask<HfCandidateCreatorXMc>(cfgc));
  }
  return workflow;
}
