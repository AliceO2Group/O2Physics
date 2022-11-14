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

/// \file candidateCreatorB0.cxx
/// \brief Reconstruction of B0 candidates
/// \note Adapted from candidateCreatorXicc.cxx
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::aod::hf_cand_prong3;
using namespace o2::aod::hf_cand_b0; // from CandidateReconstructionTables.h
using namespace o2::framework::expressions;

/// Reconstruction of B0 candidates
struct HfCandidateCreatorB0 {
  Produces<aod::HfCandB0Base> rowCandidateBase; // table defined in CandidateReconstructionTables.h

  // vertexing
  Configurable<double> bz{"bz", 20., "magnetic field"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any B0 is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  // selection
  Configurable<double> ptPionMin{"ptPionMin", 0.5, "minimum pion pT threshold (GeV/c)"};
  Configurable<int> selectionFlagD{"selectionFlagD", 1, "Selection Flag for D"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};

  double massPi = RecoDecay::getMassPDG(kPiPlus);
  double massD = RecoDecay::getMassPDG(pdg::Code::kDMinus);
  double massDPi = 0.;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagD); // FIXME

  OutputObj<TH1F> hMassDToPiKPi{TH1F("hMassB0ToPiKPi", "D^{#minus} candidates;inv. mass (p^{#minus} K^{#plus} #pi^{#minus}) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtD{TH1F("hPtD", "D^{#minus} candidates;D^{#minus} candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hPtPion{TH1F("hPtPion", "#pi^{#plus} candidates;#pi^{#plus} candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hCPAD{TH1F("hCPAD", "D^{#minus} candidates;D^{#minus} cosine of pointing angle;entries", 110, -1.1, 1.1)};
  OutputObj<TH1F> hMassB0ToDPi{TH1F("hMassB0ToDPi", "2-prong candidates;inv. mass (B^{0} #rightarrow D^{#minus}#pi^{#plus} #rightarrow #pi^{#minus}K^{#plus}#pi^{#minus}#pi^{#plus}) (GeV/#it{c}^{2});entries", 500, 3., 8.)};
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};

  void process(aod::Collision const& collision,
               soa::Filtered<soa::Join<
                 aod::HfCandProng3,
                 aod::HfSelDplusToPiKPi>> const& dCands,
               aod::BigTracks const& tracks)
  {
    // Initialise fitter for B vertex (2-prong vertex filter)
    o2::vertexing::DCAFitterN<2> df2;
    df2.setBz(bz);
    df2.setPropagateToPCA(propagateToPCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(true);

    // Initial fitter to redo D-vertex to get extrapolated daughter tracks (3-prong vertex filter)
    o2::vertexing::DCAFitterN<3> df3;
    df3.setBz(bz);
    df3.setPropagateToPCA(propagateToPCA);
    df3.setMaxR(maxR);
    df3.setMaxDZIni(maxDZIni);
    df3.setMinParamChange(minParamChange);
    df3.setMinRelChi2Change(minRelChi2Change);
    df3.setUseAbsDCA(true);

    // loop over D candidates
    for (auto const& dCand : dCands) {
      if (!TESTBIT(dCand.hfflag(), hf_cand_prong3::DecayType::DplusToPiKPi)) {
        continue;
      }
      if (dCand.isSelDplusToPiKPi() >= selectionFlagD) {
        hMassDToPiKPi->Fill(InvMassDPlus(dCand), dCand.pt());
      }
      hPtD->Fill(dCand.pt());
      hCPAD->Fill(dCand.cpa());

      // track0 <-> pi, track1 <-> K, track2 <-> pi
      auto track0 = dCand.prong0_as<aod::BigTracks>();
      auto track1 = dCand.prong1_as<aod::BigTracks>();
      auto track2 = dCand.prong2_as<aod::BigTracks>();
      auto trackParVar0 = getTrackParCov(track0);
      auto trackParVar1 = getTrackParCov(track1);
      auto trackParVar2 = getTrackParCov(track2);
      auto collision = track0.collision();

      // reconstruct 3-prong secondary vertex (D±)
      if (df3.process(trackParVar0, trackParVar1, trackParVar2) == 0) {
        continue;
      }

      const auto& secondaryVertex = df3.getPCACandidate();
      trackParVar0.propagateTo(secondaryVertex[0], bz);
      trackParVar1.propagateTo(secondaryVertex[0], bz);
      trackParVar2.propagateTo(secondaryVertex[0], bz);

      // D∓ → π∓ K± π∓
      array<float, 3> pVecpiK = {track0.px() + track1.px(), track0.py() + track1.py(), track0.pz() + track1.pz()};
      array<float, 3> pVecD = {pVecpiK[0] + track2.px(), pVecpiK[1] + track2.py(), pVecpiK[2] + track2.pz()};
      auto trackParVarPiK = o2::dataformats::V0(df3.getPCACandidatePos(), pVecpiK, df3.calcPCACovMatrixFlat(),
                                                trackParVar0, trackParVar1, {0, 0}, {0, 0});
      auto trackParVarD = o2::dataformats::V0(df3.getPCACandidatePos(), pVecD, df3.calcPCACovMatrixFlat(),
                                              trackParVarPiK, trackParVar2, {0, 0}, {0, 0});

      int index0D = track0.globalIndex();
      int index1D = track1.globalIndex();
      int index2D = track2.globalIndex();
      // int charge = track0.sign() + track1.sign() + track2.sign();

      // loop on D-
      // D- → π- K+ π-
      // we don't have direct access to D sign so we use the sign of the daughters (the pion track0 here)
      if (track0.sign() < 0) {
        // loop over pions
        for (auto const& trackPion : tracks) {
          // minimum pT selection
          if (trackPion.pt() < ptPionMin) {
            continue;
          }
          // we reject pions that are D daughters
          if (trackPion.globalIndex() == index0D || trackPion.globalIndex() == index1D || trackPion.globalIndex() == index2D) {
            continue;
          }
          // we only keep pi+ to combine them with Dminus and reconstruct B0
          if (trackPion.sign() < 0) {
            continue;
          }

          hPtPion->Fill(trackPion.pt());
          array<float, 3> pVecPion;
          auto trackParVarPi = getTrackParCov(trackPion);

          // ---------------------------------
          // reconstruct the 2-prong B0 vertex
          if (df2.process(trackParVarD, trackParVarPi) == 0) {
            continue;
          }

          // calculate relevant properties
          const auto& secondaryVertexB0 = df2.getPCACandidate();
          auto chi2PCA = df2.getChi2AtPCACandidate();
          auto covMatrixPCA = df2.calcPCACovMatrixFlat();

          df2.propagateTracksToVertex();
          df2.getTrack(0).getPxPyPzGlo(pVecD);
          df2.getTrack(1).getPxPyPzGlo(pVecPion);

          auto primaryVertex = getPrimaryVertex(collision);
          auto covMatrixPV = primaryVertex.getCov();
          o2::dataformats::DCA impactParameter0;
          o2::dataformats::DCA impactParameter1;
          trackParVarD.propagateToDCA(primaryVertex, bz, &impactParameter0);
          trackParVarPi.propagateToDCA(primaryVertex, bz, &impactParameter1);

          hCovSVXX->Fill(covMatrixPCA[0]);
          hCovPVXX->Fill(covMatrixPV[0]);

          // get uncertainty of the decay length
          double phi, theta;
          getPointDirection(array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertexB0, phi, theta);
          auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
          auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

          int hfFlag = BIT(hf_cand_b0::DecayType::B0ToDPi);

          // fill the candidate table for the B0 here:
          rowCandidateBase(collision.globalIndex(),
                           collision.posX(), collision.posY(), collision.posZ(),
                           secondaryVertexB0[0], secondaryVertexB0[1], secondaryVertexB0[2],
                           errorDecayLength, errorDecayLengthXY,
                           chi2PCA,
                           pVecD[0], pVecD[1], pVecD[2],
                           pVecPion[0], pVecPion[1], pVecPion[2],
                           impactParameter0.getY(), impactParameter1.getY(),
                           std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()),
                           dCand.globalIndex(), trackPion.globalIndex(),
                           hfFlag);

          // calculate invariant mass
          auto arrayMomenta = array{pVecD, pVecPion};
          massDPi = RecoDecay::m(std::move(arrayMomenta), array{massD, massPi});
          if (dCand.isSelDplusToPiKPi() > 0) {
            hMassB0ToDPi->Fill(massDPi);
          }
        } // pi+ loop
      }   // if D-

      // loop on D+
      // D+ → π+ K- π+
      // we now loop on the D+
      if (track0.sign() > 0) {
        // loop over pions
        for (auto const& trackPion : tracks) {
          // minimum pT selection
          if (trackPion.pt() < ptPionMin) {
            continue;
          }
          // we reject pions that are D daughters
          if (trackPion.globalIndex() == index0D || trackPion.globalIndex() == index1D || trackPion.globalIndex() == index2D) {
            continue;
          }
          // we onnly keep pi- to combine them with D+ and reconstruct B0bar
          if (trackPion.sign() > 0) {
            continue;
          }

          hPtPion->Fill(trackPion.pt());
          array<float, 3> pVecPion;
          auto trackParVarPi = getTrackParCov(trackPion);

          // ---------------------------------
          // reconstruct the 2-prong B0bar vertex
          if (df2.process(trackParVarD, trackParVarPi) == 0) {
            continue;
          }

          // calculate relevant properties
          const auto& secondaryVertexB0 = df2.getPCACandidate();
          auto chi2PCA = df2.getChi2AtPCACandidate();
          auto covMatrixPCA = df2.calcPCACovMatrixFlat();

          df2.propagateTracksToVertex();
          df2.getTrack(0).getPxPyPzGlo(pVecD);
          df2.getTrack(1).getPxPyPzGlo(pVecPion);

          auto primaryVertex = getPrimaryVertex(collision);
          auto covMatrixPV = primaryVertex.getCov();
          o2::dataformats::DCA impactParameter0;
          o2::dataformats::DCA impactParameter1;
          trackParVarD.propagateToDCA(primaryVertex, bz, &impactParameter0);
          trackParVarPi.propagateToDCA(primaryVertex, bz, &impactParameter1);

          hCovSVXX->Fill(covMatrixPCA[0]);
          hCovPVXX->Fill(covMatrixPV[0]);

          // get uncertainty of the decay length
          double phi, theta;
          getPointDirection(array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertexB0, phi, theta);
          auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
          auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

          int hfFlag = BIT(hf_cand_b0::DecayType::B0ToDPi);

          // fill the candidate table for the B0 here:
          rowCandidateBase(collision.globalIndex(),
                           collision.posX(), collision.posY(), collision.posZ(),
                           secondaryVertexB0[0], secondaryVertexB0[1], secondaryVertexB0[2],
                           errorDecayLength, errorDecayLengthXY,
                           chi2PCA,
                           pVecD[0], pVecD[1], pVecD[2],
                           pVecPion[0], pVecPion[1], pVecPion[2],
                           impactParameter0.getY(), impactParameter1.getY(),
                           std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()),
                           dCand.globalIndex(), trackPion.globalIndex(),
                           hfFlag);

          // calculate invariant mass
          auto arrayMomenta = array{pVecD, pVecPion};
          massDPi = RecoDecay::m(std::move(arrayMomenta), array{massD, massPi});
          if (dCand.isSelDplusToPiKPi() > 0) {
            hMassB0ToDPi->Fill(massDPi);
          }
        } // pi- loop
      }   // if D+
    }     // D loop
  }       // process
};        // struct

/// Extends the base table with expression columns.
struct HfCandidateCreatorB0Expressions {
  Spawns<aod::HfCandB0Ext> rowCandidateB0;

  void init(InitContext const&) {}
};

/// Performs MC matching.
struct HfCandidateCreatorB0Mc {
  Produces<aod::HfCandB0MCRec> rowMCMatchRec; // table defined in CandidateReconstructionTables.h
  Produces<aod::HfCandB0MCGen> rowMCMatchGen; // table defined in CandidateReconstructionTables.h

  void processMc(aod::HfCandB0 const& candidates,
                 aod::HfCandProng3 const&,
                 aod::BigTracksMC const& tracks,
                 aod::McParticles const& particlesMC)
  {
    int indexRec = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t debug = 0;

    // Match reconstructed candidates.
    for (auto const& candidate : candidates) {
      // Printf("New rec. candidate");
      flag = 0;
      origin = 0;
      debug = 0;
      auto candD = candidate.prong0();
      auto arrayDaughters = array{candD.prong0_as<aod::BigTracksMC>(),
                                  candD.prong1_as<aod::BigTracksMC>(),
                                  candD.prong2_as<aod::BigTracksMC>(),
                                  candidate.prong1_as<aod::BigTracksMC>()};
      auto arrayDaughtersD = array{candD.prong0_as<aod::BigTracksMC>(),
                                   candD.prong1_as<aod::BigTracksMC>(),
                                   candD.prong2_as<aod::BigTracksMC>()};
      // B0 → D- π+ → (π- K+ π-) π+
      // Printf("Checking B0 → D- π+");
      indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kB0, array{-kPiPlus, +kKPlus, -kPiPlus, +kPiPlus}, true, &sign, 2);
      if (indexRec > -1) {
        // D- → π- K+ π-
        // Printf("Checking D- → π- K+ π-");
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughtersD, pdg::Code::kDMinus, array{-kPiPlus, +kKPlus, -kPiPlus}, true, &sign, 1);
        if (indexRec > -1) {
          flag = sign * BIT(hf_cand_b0::DecayType::B0ToDPi);
        } else {
          debug = 1;
          LOGF(info, "WARNING: B0 in decays in the expected final state but the condition on the intermediate state is not fulfilled");
        }
      }
      rowMCMatchRec(flag, origin, debug);
    }

    // Match generated particles.
    for (auto const& particle : particlesMC) {
      // Printf("New gen. candidate");
      flag = 0;
      origin = 0;
      // B0 → D- π+
      if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdg::Code::kB0, array{-int(pdg::Code::kDPlus), +kPiPlus}, true)) {
        // Match D- -> π- K+ π-
        auto candDMC = particlesMC.rawIteratorAt(particle.daughtersIds().front());
        // Printf("Checking D- -> π- K+ π-");
        if (RecoDecay::isMatchedMCGen(particlesMC, candDMC, -int(pdg::Code::kDPlus), array{-kPiPlus, +kKPlus, -kPiPlus}, true, &sign)) {
          flag = sign * BIT(hf_cand_b0::DecayType::B0ToDPi);
        }
      }
      rowMCMatchGen(flag, origin);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorB0Mc, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HfCandidateCreatorB0>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorB0Expressions>(cfgc)};
  workflow.push_back(adaptAnalysisTask<HfCandidateCreatorB0Mc>(cfgc));
  return workflow;
}